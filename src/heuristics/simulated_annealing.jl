export simulated_annealing

global opt_fit = 1
global sort_stations = false

function simulated_annealing(initial_solution::Solution, τ⁰::Float64=300.0, τˢ::Float64=10., α::Float64=0.98, Ι::Int64=35, β::Float64=0.5)
    #keep track of the starting time
    sa_start_time = time()

    current_solution = deepcopy(initial_solution)
    current_cost = ECS_objective_function(current_solution)
    best_solution = deepcopy(current_solution)
    best_cost = current_cost
    global stations_requests = get_station_requests_from_solution(sol)
    τ = τ⁰
    while τ > τˢ
        for _ in 1:Ι  # Number of iterations at each temperature
           
            neighbor_solution = sample_neighbor(current_solution, β)
            
            #@info "sa = neighbor_solution = $neighbor_solution"

            #clean_up_cars_number!(neighbor_solution)
            neighbor_cost = ECS_objective_function(neighbor_solution)
            #= if !is_feasible_solution(neighbor_solution)
                @warn "infeasible solution"
                τ = 0
                break
            end =#

            if neighbor_cost <= current_cost
                current_solution = neighbor_solution
                current_cost = neighbor_cost
                if neighbor_cost < best_cost
                    best_solution = deepcopy(neighbor_solution)
                    best_cost = neighbor_cost
                end
            else
                #total_tried += 1
                acceptance_probability = exp((current_cost - neighbor_cost) / τ)
                #@info "Δ = $(current_cost - neighbor_cost), prob = $acceptance_probability"
                if rand(rng) < acceptance_probability
                    current_solution = neighbor_solution
                    current_cost = neighbor_cost
                    #total_accepted += 1
                end

            end

        end
        τ *= α
        #@info "current cost: $current_cost, best cost: $best_cost, temperature: $τ"
    end

    @info "best_cost = $best_cost, gap = $(round((best_cost - opt_fit )/ opt_fit * 100, digits=2))% time = $( time() - sa_start_time)"
    return best_solution, best_cost, (time() - sa_start_time)
end

##### Neighborhood functions #####
function sample_neighbor(sol::Solution, β::Float64=0.5)::Solution
    return ruin_recreate(sol)
    #=  if rand(rng) < β
        return open_station_neighborhood(sol)
    else
        return close_station_neighborhood(sol)
    end =#

end

function open_station_neighborhood(sol::Solution)
   
    ## flip state of one random station
    global rng
    global online_request_serving
    global scenario_list
    global stations_capacity

    #sol = load_sol("sol.jls")

    # step1: copy the solution
    neigh_sol = deepcopy(sol)

    closed_stations = findall(!, neigh_sol.open_stations_state)
    if length(closed_stations) == 0
        return neigh_sol
    end
    station_id = closed_stations[rand(rng, 1:length(closed_stations))] 
   

    # step2: select randomly the station and the scenario to be used
    #station_id = rand(rng, 1:length(neigh_sol.open_stations_state))
    scenario_id = rand(rng, 1:length(scenario_list))

    #make sure that the station is open
    neigh_sol.open_stations_state[station_id] = true

    if online_request_serving
        # step4: update the station initial number of cars
        neigh_sol.initial_cars_number[station_to_open] = floor(stations_capacity[station_to_open] / 2)
    else
        #in offline mode we need to serve new requests
        #serve_new_requests!(neigh_sol, scenario_list, collect(1:length(scenario_list)), station_id)
        serve_new_requests!(neigh_sol, scenario_list, [scenario_id], station_id)
    end

    return neigh_sol
end

function close_station_neighborhood(sol::Solution)

    ## flip state of one random station
    global rng
    global online_request_serving
    global scenario_list
    global stations_capacity

    # step 1: copy the solution
    neigh_sol = deepcopy(sol)

    # step 2: get list of opened stations and select randomly the station to close
    opened_stations = findall(neigh_sol.open_stations_state)
    if length(opened_stations) == 0
        return neigh_sol
    end
    stations_to_close = [opened_stations[rand(rng, 1:length(opened_stations))]]

    neigh_sol.open_stations_state[stations_to_close] .= false
    neigh_sol.initial_cars_number[stations_to_close] .= 0

    if !online_request_serving
        # step 3: get the requests List to lost
        stations_to_check, lost_requests = clean_up_trips!(neigh_sol, scenario_list, stations_to_close)

        # step 4: reassign the lost requests if we can
        # assigne_requests!(neigh_sol, scenario_list, lost_requests)
        for st in stations_to_check
            serve_new_requests!(neigh_sol, scenario_list, collect(1:length(scenario_list)), st)
        end
    end

    return neigh_sol
end
##### Utils functions #####

function clean_up_trips!(sol::Solution, scenario_list::Array{Scenario,1}, stations_ids::Vector{Int64})
    #sol = load_sol("sol.jls")
    #stations_ids = [7]

    unseleceted_requests = Vector{Vector{Int64}}()
    rechecked_stations = Int64[]
    stations_nodes_ids = get_potential_locations()[stations_ids]

    for scenario in scenario_list
        # scenario = scenario_list[2];
        # get requests to lost and unselect their trips
        all_stations_trips_ids = vcat(station_trips_ids[scenario.scenario_id][stations_ids]...)
        curr_scenario_requests_to_lost = Int64[]
        stations_to_recheck = Int64[]
        for trip_id in all_stations_trips_ids
            if sol.selected_paths[scenario.scenario_id][trip_id]
                trip = scenario.feasible_paths[trip_id, :]
                push!(curr_scenario_requests_to_lost, trip.req)
                sol.selected_paths[scenario.scenario_id][trip_id] = false
                origin_station_id, destination_station_id = locations_dict[trip.origin_station], locations_dict[trip.destination_station]
                origin_station_id ∉ stations_nodes_ids && push!(stations_to_recheck, origin_station_id)
                destination_station_id ∉ stations_nodes_ids && push!(stations_to_recheck, destination_station_id)
            end
        end

        unique!(stations_to_recheck)

        for station_id in stations_to_recheck
            # station_id = stations_to_recheck[4]
            additional_station_to_recheck, requests_to_lost_from_station = recheck_station!(sol, scenario, station_id)
            length(additional_station_to_recheck) > 0 && push!(stations_to_recheck, additional_station_to_recheck...)
            length(requests_to_lost_from_station) > 0 && push!(curr_scenario_requests_to_lost, requests_to_lost_from_station...)
        end

        push!(unseleceted_requests, curr_scenario_requests_to_lost)
        push!(rechecked_stations, stations_to_recheck ...)
    end
    
    return rechecked_stations, unseleceted_requests
end

function recheck_station!(sol::Solution, scenario::Scenario, station_id::Int64)
    global stations_capacity

    stations_to_recheck, requests_to_unserve = Int64[], Int64[]
    
    # Step 1: check the new bound of the station and creect it if the solution is not feasible
    sc_bounds = get_station_cars_bounds(scenario, sol, station_id)

    if sc_bounds[1] > sc_bounds[2]
        stations_to_recheck, requests_to_unserve = correct_station_trips(scenario, station_id, sol, sol.initial_cars_number[station_id])
        sc_bounds = get_station_cars_bounds(scenario, sol, station_id)
       
    end

    #check the other scenarios
    if length(scenario_list) == 1
        sol.initial_cars_number[station_id] = sc_bounds[1]
        return (unique!(stations_to_recheck), unique!(requests_to_unserve))
    end

    bounds = []
    for sc_id in eachindex(scenario_list)
        #sc_id = 2
        if sc_id == scenario.scenario_id
            continue
        end
        if isempty(bounds)
            bounds = get_station_cars_bounds(scenario_list[sc_id], sol, station_id)
        else
            bounds = hcat(bounds, get_station_cars_bounds(scenario_list[sc_id], sol, station_id))
        end
    end

    overall_bound = [maximum(bounds[1, :]), minimum(bounds[2, :])]
    bounds_intersection = intersect(overall_bound[1]:overall_bound[2], sc_bounds[1]:sc_bounds[2])
    if isempty(bounds_intersection)
        #the station is not feasible compared to the other scenarios
        sol.initial_cars_number[station_id] = overall_bound[1]
        stations_to_recheck1, requests_to_unserve1 = correct_station_trips(scenario, station_id, sol, sol.initial_cars_number[station_id])
        push!(stations_to_recheck, stations_to_recheck1...)
        push!(requests_to_unserve, requests_to_unserve1...)
    else
        #it is important to set the cars number to lower value of intersection 
        sol.initial_cars_number[station_id] = bounds_intersection[1] 
    end

    unique!(stations_to_recheck)
    unique!(requests_to_unserve)

    #return the stations to recheck and the requests to unserved as well as the new cars number of the station
    return stations_to_recheck, requests_to_unserve
end

function assigne_requests!(sol::Solution, scenario_list::Array{Scenario,1}, requests_list::Vector{Vector{Int64}})

    #step 1: order all the requests by their revenue
    requests_list_df = vcat([DataFrame(reqId=requests_list[sc_id],
        Rev=scenario_list[sc_id].request_list.Rev[requests_list[sc_id]],
        scenario_id=ones(Int, length(requests_list[sc_id])) .* sc_id)
                             for sc_id in eachindex(scenario_list)]...)

    sort!(requests_list_df, [:scenario_id, :Rev], rev=[false, true])

    #loop over the requests
    for req in eachrow(requests_list_df)
        #req = requests_list_df[3, :]
        #get the feasible trips for the current request
        #curr_req_feasible_trips = filter(x -> x.req == req.reqId, scenario_list[req.scenario_id].feasible_paths)
        curr_req_feasible_trips = scenario_list[req.scenario_id].feasible_paths[request_feasible_trips_ids[req.scenario_id][req.reqId], :]

        origin_stations_ids = unique([locations_dict[x] for x in curr_req_feasible_trips.origin_station])
        destination_stations_ids = unique([locations_dict[x] for x in curr_req_feasible_trips.destination_station])
        # sort the station to less restrictive first
        #origin_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in origin_stations_ids]
        #destination_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in destination_stations_ids]

        #= origin_station_order = sortperm([origin_stations_bounds[i][2] - #= upper =#
                                         origin_stations_bounds[i][1] #= upper =#
                                         for i in eachindex(origin_stations_bounds)], rev=true) =#

        #= destination_station_order = sortperm([destination_stations_bounds[i][2] - #= upper =#
                                              destination_stations_bounds[i][1] #= upper =#
                                              for i in eachindex(destination_stations_bounds)], rev=true) =#

        # no sorting strategy for now
        origin_station_order = collect(1:length(origin_stations_ids))
        destination_station_order = collect(1:length(destination_stations_ids))
        origin_can_serve, des_can_serve = false, false

        for o_id in origin_stations_ids[origin_station_order]
            #o_id = origin_stations_ids[3]
            origin_trips = filter(x -> x.origin_station == get_potential_locations()[o_id], curr_req_feasible_trips)
            origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario_list, req.scenario_id, sol, o_id, origin_trips[1, :])

            if !origin_can_serve
                continue
            end

            for d_id in destination_stations_ids[destination_station_order]
                trips = filter(x -> x.destination_station == get_potential_locations()[d_id], origin_trips)
                if isempty(trips)
                    continue
                end
                trip = trips[1, :]
                des_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario_list, req.scenario_id, sol, d_id, trip)

                if !des_can_serve
                    continue
                end

                sol.initial_cars_number[o_id] = origin_new_cars
                sol.initial_cars_number[d_id] = destination_new_cars

                sol.selected_paths[req.scenario_id][trip.fp_id] = true # we can serve the trip
                break # we skip all the other trips of the request 
            end

            #here we served the request so no need to check the other origin stations
            origin_can_serve && des_can_serve && break
        end
    end
end

function correct_station_trips(scenario::Scenario, station_id::Int64, sol::Solution, initial_number_of_cars)
    global stations_capacity
    global station_trips_ids

   
    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]

    trips = scenario.feasible_paths[station_trips_ids[scenario.scenario_id][station_id], :]

    
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : trip.arriving_time + 0.1 for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time)

    # the results of the correction
    stations_to_recheck, requests_to_unserve = Int64[], Int64[]

    # we set the cars number to lower bound of all the other scenarios
    cars_number = initial_number_of_cars

    for curr_trip_id in trips_order
        #curr_trip_id = trips_order[1]
        if !sol.selected_paths[scenario.scenario_id][trips.fp_id[curr_trip_id]] 
            continue
        end

        cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1

        if cars_number < 0 || cars_number > station_capacity
            #we can not serve the requests anymore
            sol.selected_paths[scenario.scenario_id][trips.fp_id[curr_trip_id]] = false
            #@info "new method : unserve $(trips.req[curr_trip_id])"
            push!(requests_to_unserve, trips.req[curr_trip_id])
            station_to_recheck = (trips.origin_station[curr_trip_id] == station_node_id) ? trips.destination_station[curr_trip_id] : trips.origin_station[curr_trip_id]
            push!(stations_to_recheck, locations_dict[station_to_recheck])

            cars_number -= (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        end
    end

    return stations_to_recheck, requests_to_unserve
end

function get_station_cars_bounds(scenario::Scenario, sol::Solution, station_id::Int64)
    
    #the lower bound of station is the minimum number of cars needed to serve the requests (all the cars are used)
    #the upper bound is the maximum number of cars needed to serve the requests (ther may be unused cars)
    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]

    #get the list of trips of the scenario where station intervenes 
    trips_ids = station_trips_ids[scenario.scenario_id][station_id][findall(sol.selected_paths[scenario.scenario_id][station_trips_ids[scenario.scenario_id][station_id]])]
    trips = scenario.feasible_paths[trips_ids, :]
    #trips = scenario.feasible_paths[station_trips_ids[scenario.scenario_id][station_id], :]
    #trips = filter(x-> sol.selected_paths[1][x.fp_id] && (x.origin_station == station_node_id || x.destination_station == station_node_id),scenario.feasible_paths)
    
    # for the time order we privilige the trips where the station is the origin station (for that we added .1 to the arriving time)
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : (trip.arriving_time + 0.1) for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time, alg=InsertionSort)

    # check if the nbr of cars is always positive and less than the station capacity 
    curr_cars_number = sol.initial_cars_number[station_id]
    cars_number_list = [curr_cars_number]

    for curr_trip_id in trips_order
        #curr_trip_id = trips_order[3]
        if !sol.selected_paths[scenario.scenario_id][trips.fp_id[curr_trip_id]]
            continue
        end
            
        curr_cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        push!(cars_number_list, curr_cars_number)
    end

    lower_bound = sol.initial_cars_number[station_id] - minimum(cars_number_list)
    upper_bound = sol.initial_cars_number[station_id] + station_capacity - maximum(cars_number_list)


    return [lower_bound, upper_bound]
end

function can_serve_and_get_cars_number(scenario_list::Vector{Scenario}, scenario_id::Int64, sol::Solution, station_id::Int64, trip::DataFrameRow)::Tuple{Bool,Int64}
    # scenario_id, station_id = 1, 4
    #first check if the station is open
    !sol.open_stations_state[station_id] && return (false, -1)
    
    #let set the trip as selected
    sol.selected_paths[scenario_id][trip.fp_id] = true
    
    # check the feasibility of serving the trip in the conserned scenario
    sc_bounds = get_station_cars_bounds(scenario_list[scenario_id], sol, station_id)

    if sc_bounds[1] > sc_bounds[2]
        #we can not serve the trip
        sol.selected_paths[scenario_id][trip.fp_id] = false
        return (false, -1)
    end

    if length(scenario_list) == 1
        sol.selected_paths[scenario_id][trip.fp_id] = false
        return (true, sc_bounds[1])
    end

    # here we can serve the trip in this station se we check the other scenarios
    bounds = hcat([get_station_cars_bounds(scenario, sol, station_id) for scenario in scenario_list]...)

    # the feasibility of serving the trip in the other scenarios is checked if all the bounds intersects
    overall_bound = maximum(bounds[1,:]), minimum(bounds[2,:])
    if overall_bound[1] > overall_bound[2]
        #we can not serve the trip
        sol.selected_paths[scenario_id][trip.fp_id] = false
        return (false, -1)
    end

    # here we can serve the trip without fuck the other scenarios
    sol.selected_paths[scenario_id][trip.fp_id] = false
    return (true, overall_bound[1])
end


function get_all_stations_bounds(sol)
    stations_bounds = []
    for st in eachindex(sol.open_stations_state)
        bounds = []
        for i in eachindex(scenario_list)
            if isempty(bounds)
                bounds = e.get_station_cars_bounds(e.scenario_list[i], sol, st)
            else
                bounds = hcat(bounds, e.get_station_cars_bounds(e.scenario_list[i], sol, st))
            end

        end
        push!(stations_bounds, [maximum(bounds[1, :]), minimum(bounds[2, :])])
    end
    return stations_bounds
end

function station_all_scenario_bounds(sol::Solution, scenario_list::Array{Scenario,1}, station_id::Int64)
    #sol, station_id = generate_random_solution(), 15 
    bounds = Matrix{Int}(undef, 2, length(scenario_list))
    for i in eachindex(scenario_list)
        #i = 2
        bounds[:, i] = get_station_cars_bounds(scenario_list[i], sol, station_id)
    end
    return maximum(bounds[1,:]), minimum(bounds[2,:])
end

function serve_new_requests!(sol::Solution, scenario_list::Array{Scenario,1}, selected_scenarios::Vector{Int64}, station_id::Int64)
    #save_sol(sol, "sol.jls")
    #sol, station_id = load_sol("sol.jls"), 42; 
    global stations_capacity
    global request_feasible_trips_ids
    global locations_dict

    station_node_id = get_potential_locations()[station_id]

    # Step 1: get all the trips for all the scenario where station intervene
    trips = [filter(x -> x.origin_station == station_node_id || x.destination_station == station_node_id,
        scenario.feasible_paths) for scenario in scenario_list[selected_scenarios]]

    for i in eachindex(trips)
        sc_id = selected_scenarios[i]
        trips[i].scenario_id = ones(Int, nrow(trips[i])) .* sc_id
    end
    trips = vcat(trips...)

    requests_list = unique(trips, [:req, :scenario_id])
   
    sort!(requests_list, [:scenario_id, :Rev], rev=[false, true])
   
    for req in eachrow(requests_list)

        # check if the requests is not served yet
        is_served = !isnothing(findfirst(sol.selected_paths[req.scenario_id][request_feasible_trips_ids[req.scenario_id][req.req]]))
        if is_served
            #the trip is already served
            continue
        end
        #the request is not served yet

        #check if the request can be served from the station id
        station_can_serve, station_new_cars = can_serve_and_get_cars_number(scenario_list, req.scenario_id, sol, station_id, req)

        if !station_can_serve
            continue
        end

        #get the trips of the request
        req_trips = filter(x -> x.origin_station == station_node_id || x.destination_station == station_node_id,
            scenario_list[req.scenario_id].feasible_paths[request_feasible_trips_ids[req.scenario_id][req.req], :])

        other_stations = locations_dict[req_trips.origin_station[1]] == station_id ?
                         [locations_dict[x] for x in req_trips.destination_station] :
                         [locations_dict[x] for x in req_trips.origin_station]

        #get_other stations_bounds to privilege the less ristricted stations
        if sort_stations
            other_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, other_stations[i]) for i in eachindex(other_stations)]
            other_station_order = sortperm([other_stations_bounds[i][2] - #= upper =#
                                            other_stations_bounds[i][1] #= upper =#
                                            for i in eachindex(other_stations)], rev=true)
        else
            other_station_order = collect(1:length(other_stations))
        end
        #other_station_order = collect(1:length(other_stations))
        #check if we can serve the trip
        for i in other_station_order
            # i = other_station_order[1] 
            trip = req_trips[i, :]
            other_station_id = other_stations[i]
            other_station_can_serve, other_station_new_cars = can_serve_and_get_cars_number(scenario_list, req.scenario_id, sol, other_station_id, trip)

            if other_station_can_serve
                #here we are sure that we can serve the trip
                sol.initial_cars_number[station_id] = station_new_cars
                sol.initial_cars_number[other_station_id] = other_station_new_cars
                sol.selected_paths[req.scenario_id][trip.fp_id] = true # we can serve the trip
                break
            end
        end
    end

end

function correct_station_trips_new(scenario::Scenario, station_id::Int64, sol::Solution, initial_number_of_cars)
    #@info "############# corection station $station_id ###############"
    #sol, station_id = load_sol("old_sol.jls"), 34;
    #initial_number_of_cars = sol.initial_cars_number[station_id]
     
    global stations_capacity
    global station_trips_ids
    # the results of the correction
    stations_to_recheck, requests_to_unserve = Int64[], Int64[]


    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]
 
    #get the trips of the station 
    trips_ids = station_trips_ids[scenario.scenario_id][station_id][findall(sol.selected_paths[scenario.scenario_id][station_trips_ids[scenario.scenario_id][station_id]])]
    trips = scenario.feasible_paths[trips_ids, :]
    
    #sort the trips according to their interesting time for this station
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : trip.arriving_time + 0.1 for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time)
    trips = trips[trips_order, :]
    
    not_correct = true
    while not_correct
        not_correct = false
        cars_number = initial_number_of_cars
        pickup_trips = []
        dropoff_trips = []
        for i in 1:nrow(trips)
            #i = 15
            trip = trips[i, :]
            cars_number += (trip.origin_station == station_node_id) ? -1 : 1
            (trip.origin_station == station_node_id) ? push!(pickup_trips, i) : push!(dropoff_trips, i)
            
            if cars_number < 0 

                not_correct = true

                # delete the less revenue trip
                #selected_trip_id is the index of the selected trip in the trips dataframe
                selected_trip_id = pickup_trips[argmin(trips.Rev[pickup_trips])]
                sol.selected_paths[scenario.scenario_id][trips.fp_id[selected_trip_id]] = false
                push!(requests_to_unserve, trips.req[selected_trip_id])
                station_to_recheck = locations_dict[trips.destination_station[selected_trip_id]]
                push!(stations_to_recheck, station_to_recheck)
                
                deleteat!(trips, selected_trip_id)
                
                break
            elseif cars_number > station_capacity
                not_correct = true

                # delete the less revenue trip
                #selected_trip_id is the index of the selected trip in the trips dataframe
                selected_trip_id = dropoff_trips[argmin(trips.Rev[dropoff_trips])]
                sol.selected_paths[scenario.scenario_id][trips.fp_id[selected_trip_id]] = false
                push!(requests_to_unserve, trips.req[selected_trip_id])
                station_to_recheck = locations_dict[trips.origin_station[selected_trip_id]]
                push!(stations_to_recheck, station_to_recheck)
                
                deleteat!(trips, selected_trip_id)
                
                break
            end
        end
    end
   #=  bounds = get_station_cars_bounds(scenario, sol, station_id)
    if bounds[1] > bounds[2]
        @warn "infeasible solution after correction"
    end =#
    return stations_to_recheck, requests_to_unserve
end

