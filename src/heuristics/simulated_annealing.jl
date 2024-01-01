export simulated_annealing

global request_feasible_trips_ids = [] #to get feasible trips for req i on scenario s : request_feasible_trips_ids[s][i]

function set_trips_to_requets_var()
    global request_feasible_trips_ids = []

    for scenario in scenario_list
        fp = [Int[] for _ in eachindex(scenario.request_list.reqId)]
        for i in eachindex(scenario.feasible_paths.req)
            push!(fp[scenario.feasible_paths.req[i]], i)
        end
        push!(request_feasible_trips_ids, fp)
    end
end

function simulated_annealing(initial_solution::Solution, τ⁰::Float64=329., τˢ::Float64=0.1, α::Float64=0.9998, Ι::Int64=1775, β::Float64=0.5)
    
    set_trips_to_requets_var()
    nbr_try = 0
    accepted_worse = 0
    current_solution = deepcopy(initial_solution)
    current_cost = E_carsharing_sim(current_solution)
    best_solution = deepcopy(current_solution)
    best_cost = current_cost
    τ = τ⁰
    while τ > τˢ
        for _ in 1:Ι  # Number of iterations at each temperature
            neighbor_solution = sample_neighbor(current_solution, β)
            #clean_up_cars_number!(neighbor_solution)
            neighbor_cost = ECS_objective_function(neighbor_solution)
            #= fit = E_carsharing_sim(neighbor_solution)
            if fit != neighbor_cost
                println("fit: $fit, neighbor_cost: $neighbor_cost")
                break
            end =#

            if neighbor_cost < current_cost
                current_solution = neighbor_solution
                current_cost = neighbor_cost
                if neighbor_cost < best_cost
                    best_solution = deepcopy(neighbor_solution)
                    best_cost = neighbor_cost
                end
            else
                acceptance_probability = exp((current_cost - neighbor_cost) / τ)
                nbr_try += 1
                if rand(rng) < acceptance_probability
                    accepted_worse += 1
                    current_solution = neighbor_solution
                    current_cost = neighbor_cost
                end
            end
        end
        τ *= α
        #@info "current cost: $current_cost, best cost: $best_cost, temperature: $τ"
    end
    #@info "nbr of try: $nbr_try, nbr of accepted worse: $accepted_worse"
    return best_solution, best_cost
end

##### Neighborhood functions #####
function sample_neighbor(sol::Solution, β::Float64=0.5)::Solution
    
    if rand(rng) < β
        return open_station_neighborhood(sol)
    else
        return close_station_neighborhood(sol)
    end

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
    
    # step2: get list of closed stations and select randomly the station to be opened
    closed_stations = findall(!, neigh_sol.open_stations_state)
    if length(closed_stations) == 0
        return neigh_sol
    end
    station_to_open = closed_stations[rand(rng, 1:length(closed_stations))]
    
    neigh_sol.open_stations_state[station_to_open] = true

    if online_request_serving
        # step4: update the station initial number of cars
        neigh_sol.initial_cars_number[station_to_open] = floor(stations_capacity[station_to_open] / 2)
    else
        #in offline mode we need to serve new requests
        serve_new_requests!(neigh_sol, scenario_list, station_to_open)
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
    station_to_close = opened_stations[rand(rng, 1:length(opened_stations))]
    
    neigh_sol.open_stations_state[station_to_close] = false
    neigh_sol.initial_cars_number[station_to_close] = 0

    if !online_request_serving
        # step 3: get the requests List to lost
        lost_requests = clean_up_trips!(neigh_sol, scenario_list, station_to_close)
        
        # step 4: reassign the lost requests if we can
        assigne_requests!(neigh_sol, scenario_list, lost_requests)

    end
    
    return neigh_sol
end

##### Utils functions #####
function serve_new_requests!(sol::Solution, scenario_list::Array{Scenario,1}, station_id::Int64)
    #sol, station_id = load_sol("sol_error.jls"), 43;
    #sol.open_stations_state[station_id] = true
    global stations_capacity
    global request_feasible_trips_ids
    global locations_dict
    station_node_id = get_potential_locations()[station_id]

    for scenario in scenario_list
        # scenario = scenario_list[1]
        #step 1 : get the list of trips where station intervenes as a pickup station or a dropoff station
        trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id, scenario.feasible_paths)

        #step 2 : order the trips by their revenue
        sort!(trips, [:Rev], rev = true)
 
        for trip in eachrow(trips)
            # trip = trips[6,:]
            # check if the requests is not served yet
            is_served = !isnothing(findfirst(sol.selected_paths[scenario.scenario_id][request_feasible_trips_ids[scenario.scenario_id][trip.req]]))
            if is_served
                #the trip is already served
                continue
            end
            
            #the request is not served yet
            destination_station_id = locations_dict[trip.destination_station]
            origin_station_id = locations_dict[trip.origin_station]
            #check if we can serve the trip
            #origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario, sol, origin_station_id, trip)
            origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, origin_station_id, trip)
            if origin_can_serve
                des_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, destination_station_id, trip)
                
                if des_can_serve                     
                    #here we are sure that we can serve the trip
                    sol.initial_cars_number[origin_station_id] = origin_new_cars
                    sol.initial_cars_number[destination_station_id] = destination_new_cars
                        
                    sol.selected_paths[scenario.scenario_id][trip.fp_id] = true # we can serve the trip
                    #@info "we can serve the trip $(trip.fp_id)"
                end
            end
            
        end
        
    end
   
end

function clean_up_trips!(sol::Solution, scenario_list::Array{Scenario,1}, station_id::Int64)
    
    unseleceted_requests = Vector{Vector{Int64}}()
    station_node_id = get_potential_locations()[station_id]

    for scenario in scenario_list
        #scenario = scenario_list[1]
        # get requests to lost and unselect their trips
        trips = filter(x-> x.origin_station == station_node_id ||
                           x.destination_station == station_node_id, scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])

        curr_scenario_requests_to_lost = unique(trips.req)
        sol.selected_paths[scenario.scenario_id][trips.fp_id] .= false
        
        # Second check: if a car which start from station is used to serve a request in other station
        stations_to_recheck = map(x-> locations_dict[x], filter(x-> x != station_node_id, union(trips.origin_station, trips.destination_station)))

        for station_id in stations_to_recheck
            additional_station_to_recheck, requests_to_lost_from_station, sol.initial_cars_number[station_id] = recheck_station!(sol, scenario, station_id)
            length(additional_station_to_recheck) > 0 && push!(stations_to_recheck, additional_station_to_recheck...)
            length(requests_to_lost_from_station) > 0 && push!(curr_scenario_requests_to_lost, requests_to_lost_from_station...)
        end
        
        push!(unseleceted_requests, curr_scenario_requests_to_lost)
    end
    return unseleceted_requests
end

function recheck_station!(sol::Solution, scenario::Scenario, station_id::Int64)
    # sol = load_sol("sol_error.jls")
    # scenario = scenario_list[1]
    # station_id = findall(sol.open_stations_state)[1]
    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]
    
    stations_to_recheck, requests_to_unserve = Int64[], Int64[]

    # Step 1: check the new bound of the station
    sc_bounds = get_station_cars_bounds(scenario, sol, station_id)
    bounds = hcat([get_station_cars_bounds(scenario, sol, station_id) for scenario in scenario_list]...)
    overall_bound = [maximum(bounds[1,:]), minimum(bounds[2,:])]
    #@info "recheck station $station_id :  starting correction"
    #save_sol(sol, "sol_error.jls")
    while overall_bound[1] > overall_bound[2]
        # sure the current scenario that causes the problem
        # correction is needed
        trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id,
                 scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])
   
        trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : trip.arriving_time for trip in eachrow(trips)]
        trips_order = sortperm(trips_interesting_time)
        
        # check if the nbr of cars is always positive and less than the station capacity
        cars_number = sc_bounds[2] 
        cars_number_list = [cars_number]

        for curr_trip_id in trips_order
            #curr_trip_id = trips_order[1]
            cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
            if cars_number < 0 || cars_number > station_capacity
                #we can not serve the requests anymore
                sol.selected_paths[scenario.scenario_id][trips.fp_id[curr_trip_id]] = false
                #@info "new method : unserve $(trips.req[curr_trip_id])"
                push!(requests_to_unserve, trips.req[curr_trip_id])
                station_to_recheck = (trips.origin_station[curr_trip_id] == station_node_id) ? trips.destination_station[curr_trip_id] : trips.origin_station[curr_trip_id]
                push!(stations_to_recheck, locations_dict[station_to_recheck])
                
                cars_number -= (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
            else
                push!(cars_number_list, cars_number)
            end
        end

        sc_bounds[1] = sc_bounds[1] - minimum(cars_number_list)
        sc_bounds[2] = sc_bounds[1] + station_capacity - maximum(cars_number_list)

        bounds[1, scenario.scenario_id] = sc_bounds[1]
        bounds[2, scenario.scenario_id] = sc_bounds[2]
        overall_bound = [maximum(bounds[1,:]), minimum(bounds[2,:])]
    end
    #@info "recheck station $station_id : END"
    #return the stations to recheck and the requests to unserved as well as the new cars number of the station
    return stations_to_recheck, requests_to_unserve, overall_bound[1]
end

function assigne_requests!(sol::Solution, scenario_list::Array{Scenario,1}, requests_list::Vector{Vector{Int64}})
    #sol, requests_list = load_sol("sol.jls"), [[271, 421, 772]]
    for sc_id in eachindex(scenario_list)
        scenario = scenario_list[sc_id]
        #loop over the requests
        for req_id in requests_list[sc_id]
            # req_id = requests_list[sc_id][2]
            #get the feasible trips for the current request
            curr_req_feasible_trips = filter(x-> x.req == req_id, scenario.feasible_paths)
            
            for trip in eachrow(curr_req_feasible_trips)
                #trip = curr_req_feasible_trips[2,:]
                origin_station_id = locations_dict[trip.origin_station]
                destination_station_id = locations_dict[trip.destination_station]
                
                origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, origin_station_id, trip)
            
                if origin_can_serve
                    des_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, destination_station_id, trip)
                    
                    if des_can_serve 
                        #here we are sure that we can serve the trip
                        sol.initial_cars_number[origin_station_id] = origin_new_cars
                        sol.initial_cars_number[destination_station_id] = destination_new_cars
        
                        sol.selected_paths[scenario.scenario_id][trip.fp_id] = true # we can serve the trip
                        break # we skip all the other trips of the request 
                    end
                end
                
            end
        end

    end
end

function get_station_cars_bounds(scenario::Scenario, sol::Solution, station_id::Int64)
    
    #the lower bound of station is the minimum number of cars needed to serve the requests (all the cars are used)
    #the upper bound is the maximum number of cars needed to serve the requests (ther may be unused cars)
    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]
    
    #get the list of trips of the scenario where station intervenes as a dropoff station
    trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id,
                 scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])
    # for the time order we privilige the trips where the station is the origin station (for that we added .1 to the arriving time)
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : (trip.arriving_time + .1) for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time)
     
    # check if the nbr of cars is always positive and less than the station capacity 
    curr_cars_number = sol.initial_cars_number[station_id]
    cars_number_list = [curr_cars_number]
    
    for curr_trip_id in trips_order
        #curr_trip_id = trips_order[2]
        curr_cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        push!(cars_number_list, curr_cars_number)
    end
    lower_bound = sol.initial_cars_number[station_id] - minimum(cars_number_list)
    upper_bound = sol.initial_cars_number[station_id] + station_capacity - maximum(cars_number_list) 
    
    return [lower_bound, upper_bound]
end

function can_serve_and_get_cars_number(scenario_list::Vector{Scenario},scenario_id::Int64, sol::Solution, station_id::Int64, trip::DataFrameRow)::Tuple{Bool, Int64}
    
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

    # here we can serve the trip in this station se we check the other scenarios
    bounds = hcat([get_station_cars_bounds(scenario, sol, station_id) for scenario in scenario_list]...)

    # the feasibility of serving the trip in the other scenarios is checked if all the bounds intersects
    overall_bound = [maximum(bounds[1,:]), minimum(bounds[2,:])]
    if overall_bound[1] > overall_bound[2]
        #we can not serve the trip
        sol.selected_paths[scenario_id][trip.fp_id] = false
        return (false, -1)
    end

    # here we can serve the trip without fuck the other scenarios
    sol.selected_paths[scenario_id][trip.fp_id] = false
    return (true , overall_bound[1])
end

