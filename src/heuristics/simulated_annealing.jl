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
    #sol = generate_random_solution();
    global stations_capacity
    global request_feasible_trips_ids
    global locations_dict
    station_node_id = get_potential_locations()[station_id]

    cars_number_bounds = hcat([get_station_cars_bounds(scenario, sol, station_id) for scenario in scenario_list]...)
    for scenario in scenario_list
        # scenario = scenario_list[1]
        #step 1 : get the list of trips where station intervenes as a pickup station or a dropoff station
        trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id, scenario.feasible_paths)

        #step 2 : order the trips by their revenue
        sort!(trips, [:Rev], rev = true)
 
        for trip in eachrow(trips)
            # trip = trips[1,:]
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
            origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario, sol, origin_station_id, trip)
            
            if origin_can_serve
                des_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario, sol, destination_station_id, trip)
                
                if des_can_serve 

                    #check the other scenarios
                    if (origin_new_cars != sol.initial_cars_number[origin_station_id]) || (destination_new_cars != sol.initial_cars_number[destination_station_id])
                        origin_old_cars, destination_old_cars = sol.initial_cars_number[origin_station_id], sol.initial_cars_number[destination_station_id]
                        sol.initial_cars_number[origin_station_id] = origin_new_cars
                        sol.initial_cars_number[destination_station_id] = destination_new_cars
                        if !recheck_scenarios()
                            sol.initial_cars_number[origin_station_id] = origin_old_cars
                            sol.initial_cars_number[destination_station_id] = destination_old_cars
                            continue
                        end
                    end
                    
                    #here we are sure that we can serve the trip
                    sol.selected_paths[scenario.scenario_id][trip.fp_id] = true # we can serve the trip
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
            additional_station_to_recheck, requests_to_lost_from_station = recheck_station!(sol, scenario, station_id)
            length(additional_station_to_recheck) > 0 && push!(stations_to_recheck, additional_station_to_recheck...)
            length(requests_to_lost_from_station) > 0 && push!(curr_scenario_requests_to_lost, requests_to_lost_from_station...)
        end
        
        push!(unseleceted_requests, curr_scenario_requests_to_lost)
    end
    return unseleceted_requests
end

function recheck_station!(sol::Solution, scenario::Scenario, station_id::Int64)
    
    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]
    
    #get the list of trips of the scenario where station intervenes as a dropoff station
    trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id,
                 scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])
                 
   
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : trip.arriving_time for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time)
     
    # check if the nbr of cars is always positive and less than the station capacity
    cars_number = sol.initial_cars_number[station_id]
    requests_to_unserve = Int64[]
    stations_to_recheck = Int64[]
    for curr_trip_id in trips_order
        #curr_trip_id = trips_order[1]
        cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        if cars_number < 0 || cars_number > station_capacity
            #we can not serve the requests anymore
            sol.selected_paths[scenario.scenario_id][trips.fp_id[curr_trip_id]] = false
            
            push!(requests_to_unserve, trips.req[curr_trip_id])
            station_to_recheck = (trips.origin_station[curr_trip_id] == station_node_id) ? trips.destination_station[curr_trip_id] : trips.origin_station[curr_trip_id]
            push!(stations_to_recheck, locations_dict[station_to_recheck])
            
            cars_number -= (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        end
    end
    return stations_to_recheck, requests_to_unserve
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
                
                
                origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario, sol, origin_station_id, trip)
            
                if origin_can_serve
                    des_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario, sol, destination_station_id, trip)
                    
                    if des_can_serve 
    
                        #check the other scenarios
                        if (origin_new_cars != sol.initial_cars_number[origin_station_id]) || (destination_new_cars != sol.initial_cars_number[destination_station_id])
                            origin_old_cars, destination_old_cars = sol.initial_cars_number[origin_station_id], sol.initial_cars_number[destination_station_id]
                            sol.initial_cars_number[origin_station_id] = origin_new_cars
                            sol.initial_cars_number[destination_station_id] = destination_new_cars
                            if !recheck_scenarios()
                                sol.initial_cars_number[origin_station_id] = origin_old_cars
                                sol.initial_cars_number[destination_station_id] = destination_old_cars
                                continue
                            end
                        end
                        
                        #here we are sure that we can serve the trip
                        sol.selected_paths[scenario.scenario_id][trip.fp_id] = true # we can serve the trip
                        break # we skip all the other trips of the request 
                    end
                end
                
            end
        end

    end
end

function can_serve_and_get_cars_number(scenario::Scenario, sol::Solution, station_id::Int64, trip::DataFrameRow)::Tuple{Bool, Int64}
    
    global stations_capacity
    
    #first check if the station is open
    !sol.open_stations_state[station_id] && return (false, -1)

    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]
    
    #get the list of trips of the scenario where station intervenes as a dropoff station
    trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id,
                 scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])
                 
    # add the trip to check the feasibilty
    push!(trips, trip)
    
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : trip.arriving_time for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time)
     
    # check if the nbr of cars is always positive and less than the station capacity
    cars_number = sol.initial_cars_number[station_id]
    
    curr_cars_number = sol.initial_cars_number[station_id]
    cars_number_list = [cars_number]
    
    for curr_trip_id in trips_order
        #curr_trip_id = trips_order[2]
        curr_cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        push!(cars_number_list, curr_cars_number)
    end
    
    # substruct the minimum number of cars so the minimum number of cars is 0
    cars_number_list .-= minimum(cars_number_list)  

    # first check if the capacity constraints is respected
    if isnothing(findfirst(cars_number_list .> station_capacity))
        # get the new number of cars
        new_cars_number = cars_number_list[1] 
        return (true, new_cars_number) # we can serve the trip without incrementing the initial number of cars
    else
        return (false, -1) # we can not serve the trip
    end
end


function recheck_scenarios(sol::Solution, scenario_list::Array{Scenario,1}, station_id::Int64, bounds::Array{Int64,2})
    return true
   #=  # just we make sure that the number of cars in the other scenarios are within the upper and lower bounds
    #sol = generate_random_solution();

    # station_id = findall(sol.open_stations_state)[2]
    bounds = hcat([get_station_cars_bounds(scenario, sol, station_id) for scenario in scenario_list]...)
    
    #check that the car is within the bounds
    bounds[:][1] =#
end

function get_station_cars_bounds(scenario::Scenario, sol::Solution, station_id::Int64)
    #the lower bound of station is the minimum number of cars needed to serve the requests (all the cars are used)
    #the upper bound is the maximum number of cars needed to serve the requests (ther may be unused cars)
    station_node_id = get_potential_locations()[station_id]
    station_capacity = stations_capacity[station_id]
    
    #get the list of trips of the scenario where station intervenes as a dropoff station
    trips = filter(x-> x.origin_station == station_node_id || x.destination_station == station_node_id,
                 scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])
    
    trips_interesting_time = [(trip.origin_station == station_node_id) ? trip.start_driving_time : trip.arriving_time for trip in eachrow(trips)]
    trips_order = sortperm(trips_interesting_time)
     
    # check if the nbr of cars is always positive and less than the station capacity
    cars_number = sol.initial_cars_number[station_id]
    
    curr_cars_number = sol.initial_cars_number[station_id]
    cars_number_list = [cars_number]
    
    for curr_trip_id in trips_order
        #curr_trip_id = trips_order[2]
        curr_cars_number += (trips.origin_station[curr_trip_id] == station_node_id) ? -1 : 1
        push!(cars_number_list, curr_cars_number)
    end
    lower_bound = sol.initial_cars_number[station_id] - minimum(cars_number_list)
    upper_bound = sol.initial_cars_number[station_id] + station_capacity - maximum(cars_number_list) 
    
    return [lower_bound, upper_bound]
end

#convert 

