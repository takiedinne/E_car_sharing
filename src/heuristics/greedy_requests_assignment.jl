export greedy_assign_requests
global adjacent_stations
global ruin_depth = 0.012 # the percentage of stations to be closed
global γ = 0.1 # the blink probability

global use_adjacent_selection = true

function fill_adjacent_stations()
    # for each station: get a list od the station sorted by their distance
    global adjacent_stations = Matrix{Int64}(undef, length(get_potential_locations()), length(get_potential_locations()))
    for i in 1:length(get_potential_locations())
        st = get_potential_locations()[i]
        distances =[get_walking_distance(st, st2) for st2 in get_potential_locations()]
        adjacent_stations[i, :] = sortperm(distances)
    end
    #delete the first column
    adjacent_stations = adjacent_stations[:, 2:end]
    if !use_adjacent_selection
        for i in 1:length(get_potential_locations())
            adjacent_stations[i, :] = shuffle(adjacent_stations[i, :])
        end
    end

end

function greedy_assign_requests()
    sa_start_time = time()
    global scenario_list

    scenario = scenario_list[1]
    sol = Solution()

    #get list of feasible requests
    feasible_requests = scenario.request_list[unique(scenario_list[1].feasible_paths.req), [:reqId, :Rev]]
    sort!(feasible_requests, [:Rev], rev=true)
   #=  i = 1  =#
    for req in eachrow(feasible_requests)
        #= @info "i = $i"
        i +=1 =#
        #req = feasible_requests[116, :]
        #get origin stations for the request
        curr_req_trips_ids = request_feasible_trips_ids[scenario.scenario_id][req.reqId]
        curr_req_trips = scenario.feasible_paths[curr_req_trips_ids, :]

        #get list of origin stations
        origin_stations_ids = unique([locations_dict[x] for x in curr_req_trips.origin_station])
        destination_stations_ids = unique([locations_dict[x] for x in curr_req_trips.destination_station])

        # try to serve the request from the origin stations

        # step 1: get the stations that are open
        origin_open_stations = origin_stations_ids[findall(sol.open_stations_state[origin_stations_ids])]
        origin_can_serve, origin_new_cars, origin_station = false, -1, -1

        if !isempty(origin_open_stations)
            #there is(are) station(s) open
            origin_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in origin_open_stations]
            origin_station_order = sortperm([origin_stations_bounds[i][2] - origin_stations_bounds[i][1] for i in eachindex(origin_stations_bounds)], rev=true)

            for st in origin_open_stations[origin_station_order]
                trip_id = findfirst(x -> x.origin_station == get_potential_locations()[st], eachrow(curr_req_trips))
                trip = curr_req_trips[trip_id, :]
                origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)
                if origin_can_serve
                    origin_station = st
                    break
                end
            end
        end

        if !origin_can_serve
            #try to open new station
            station_to_open = origin_stations_ids[findall(.!sol.open_stations_state[origin_stations_ids])]
            
            if !isempty(station_to_open)
                #there is(are) station(s) open
                stations_cost = [get_station_total_cost(st) for st in scenario.stations[station_to_open]]

                origin_station_order = sortperm(stations_cost)

                for st in station_to_open[origin_station_order]
                    #st = 53
                    trip_id = findfirst(x -> x.origin_station == get_potential_locations()[st], eachrow(curr_req_trips))
                    trip = curr_req_trips[trip_id, :]
                    sol.open_stations_state[st] = true
                    origin_can_serve, origin_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)
                    sol.open_stations_state[st] = false
                    if origin_can_serve
                        origin_station = st
                        break
                    end
                end
            end
        end

        if !origin_can_serve
            #we can not serve the request from the origin stations
            continue
        end

        destination_open_stations = destination_stations_ids[findall(sol.open_stations_state[destination_stations_ids])]
        destination_can_serve, destination_new_cars, destination_station = false, -1, -1

        if !isempty(destination_open_stations)
            #there is(are) station(s) open
            destination_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in destination_open_stations]
            destination_station_order = sortperm([destination_stations_bounds[i][2] - destination_stations_bounds[i][1] for i in eachindex(destination_stations_bounds)], rev=true)

            for st in destination_open_stations[destination_station_order]
                trip_id = findfirst(x -> x.destination_station == get_potential_locations()[st] && 
                                        x.origin_station == get_potential_locations()[origin_station], eachrow(curr_req_trips))
                if isnothing(trip_id)
                    #there is no feasible path origin and dst stations
                    continue
                end
                trip = curr_req_trips[trip_id, :]
                
                destination_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)
                if destination_can_serve
                    destination_station = st
                    break
                end
            end
        end

        if !destination_can_serve
            #try to open new station
            station_to_open = destination_stations_ids[findall(.!sol.open_stations_state[destination_stations_ids])]

            if !isempty(station_to_open)
                #there is(are) station(s) open
                stations_cost = [get_station_total_cost(st) for st in scenario.stations[station_to_open]]

                destination_station_order = sortperm(stations_cost)

                for st in station_to_open[destination_station_order]
                    #st = 23
                    trip_id = findfirst(x -> x.destination_station == get_potential_locations()[st] && 
                                            x.origin_station == get_potential_locations()[origin_station], eachrow(curr_req_trips))
                    if isnothing(trip_id)
                        #there is no feasible path origin and dst stations
                        continue
                    end
                    trip = curr_req_trips[trip_id, :]
                    sol.open_stations_state[st] = true
                    destination_can_serve, destination_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)
                    sol.open_stations_state[st] = false
                    if destination_can_serve
                        destination_station = st
                        break
                    end
                end
            end
            
        end
        if !destination_can_serve
            #we can not serve the request from the origin stations
            continue
        end

        #here we can serve the request
        sol.open_stations_state[origin_station] = true
        sol.initial_cars_number[origin_station] = origin_new_cars
        sol.open_stations_state[destination_station] = true
        sol.initial_cars_number[destination_station] = destination_new_cars
        
        #println("$(req.reqId), $origin_station, $destination_station")
        
        trip_id = findfirst(x -> x.origin_station == get_potential_locations()[origin_station] && x.destination_station == get_potential_locations()[destination_station], eachrow(curr_req_trips))
        sol.selected_paths[1][curr_req_trips.fp_id[trip_id]] = true

    end
    total_cpu = (time() - sa_start_time)
    global total_time
    
    return sol, E_carsharing_sim(sol), total_cpu  
end

function get_station_total_cost(station::Station)
    #=  max_number_of_charging_points
    max_power
    charging_station_base_cost
    charging_point_cost_slow
    max_charging_rate_per_charging_point_slow
    charging_point_cost_fast
    max_charging_rate_per_charging_point_fast =#
    total_cost = station.max_number_of_charging_points * station.charging_point_cost_fast + station.charging_station_base_cost
    return total_cost
end

############ Ruin procedure ##########################

function adjacent_ruin!(sol::Solution)
    # global variables
    global rng
    global adjacent_stations
    global scenario_list

    open_stations = findall(sol.open_stations_state)
    
    if length(open_stations) == 0
        # all the stations are closed so we can not close any station
        return
    end

    # calculate the number of stations to be closed
    n_stations_to_close = calculate_station_to_close(sol)
    
    #select the station that is used as seed for the ruin
    st_seed = rand(rng, 1:length(open_stations))
    
    #get the adjacent stations
    closed_stations_counter = 0 #the index of stations to close
    stations_to_close = Int64[]
    for adj_st in adjacent_stations[st_seed, :]
        
        if sol.open_stations_state[adj_st]
            sol.open_stations_state[adj_st] = false
            sol.initial_cars_number[adj_st] = 0
            push!(stations_to_close, adj_st)
            closed_stations_counter += 1
        end
        
        #check if we heve already closed engough stations
        if closed_stations_counter == n_stations_to_close
            break
        end
    end
    
    #clean the requests assignments
    clean_up_trips!(sol, scenario_list, stations_to_close)
    
end

function calculate_station_to_close(sol)
    global ruin_depth
    max_length = round(Int64, length(findall(sol.open_stations_state)) * ruin_depth) + 1
    
    return rand(rng, 1:max_length)
end

############ Recreate procedure #######################

function greedy_recreate!(sol)
    #sol = load_sol("sol.jls")
    #save_sol(sol, "sol.jls")
    global scenario_list
    global request_feasible_trips_ids
    global γ = 0. # the blink probability
    
    #get unserved requests
    unserved_requests = get_unserved_requests(sol, scenario_list)
    sort!(unserved_requests, [:Rev], rev=true, alg=InsertionSort)
    
    for req in eachrow(unserved_requests) 
        #req = unserved_requests[91, :]
        #get origin stations for the request
        request_trips_ids = request_feasible_trips_ids[req.scenario_id][req.reqId] 
        curr_req_trips = scenario_list[req.scenario_id].feasible_paths[request_trips_ids, :]
        scenario = scenario_list[req.scenario_id]

        #get list of origin stations
        origin_stations_ids = unique([locations_dict[x] for x in curr_req_trips.origin_station])
        destination_stations_ids = unique([locations_dict[x] for x in curr_req_trips.destination_station])

        # step 1: select an origin station
        origin_open_stations = origin_stations_ids[findall(sol.open_stations_state[origin_stations_ids])]
        origin_can_serve, origin_new_cars, origin_station = false, -1, -1
        
        #sort the stations by their ristricted cars bounds
        origin_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in origin_open_stations]
        origin_station_order = sortperm([origin_stations_bounds[i][2] - origin_stations_bounds[i][1] for i in eachindex(origin_stations_bounds)], rev=true)
        origin_open_stations = origin_open_stations[origin_station_order]
        
        for ost in origin_open_stations
            #the blink mechanism
            if rand(rng) < γ
                continue
            end
            trip_id = findfirst(x -> x.origin_station == get_potential_locations()[ost], eachrow(curr_req_trips))
            trip = curr_req_trips[trip_id, :]
            
            this_st_can_serve, this_st_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, ost, trip)
            
            if !origin_can_serve && this_st_can_serve
                origin_can_serve = true
                origin_station = ost
                origin_new_cars = this_st_new_cars

                #if this station can serve and we gained in terms of cars number so this is the best alternative
                if this_st_new_cars == sol.initial_cars_number[ost]
                    #accept directly this station
                    break
                end
            elseif this_st_can_serve 
                #here we found befor a feasible station, however we have to increament number of cars 
                # so we need to check if selecting this stations is better "accepted when no need to increament the cars number"
                if this_st_new_cars == sol.initial_cars_number[ost]
                    origin_station = ost
                    origin_new_cars = this_st_new_cars
                    break
                end
            end
        end
        
        if !origin_can_serve
            #we can not serve the request from an open a station
            #try to open new station
            origin_closed_stations = origin_stations_ids[findall(.!sol.open_stations_state[origin_stations_ids])]

            if !isempty(origin_closed_stations)
                #there is(are) station(s) open
                stations_cost = [get_station_total_cost(st) for st in scenario.stations[origin_closed_stations]]
                origin_closed_stations = origin_closed_stations[sortperm(stations_cost)]

                for st in origin_closed_stations
                    
                    # the blink mechanism
                    if rand(rng) < γ
                        continue
                    end

                    #no need to check if we can serve the trip as we are opening a new station
                    origin_can_serve, origin_station, origin_new_cars = true, st, 1
                    break
                end
            end
        end

        if !origin_can_serve
            #we can not serve the request from the origin stations either by already open stations or by opening a new station
            continue
        end
        
        #step 2: select a destination station
        destination_open_stations = destination_stations_ids[findall(sol.open_stations_state[destination_stations_ids])]
        destination_can_serve, destination_new_cars, destination_station = false, -1, -1
        
        #sort the stations by their ristricted cars bounds
        destination_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in destination_open_stations]
        destination_station_order = sortperm([destination_stations_bounds[i][2] - destination_stations_bounds[i][1] for i in eachindex(destination_stations_bounds)], rev=true)
        destination_open_stations = destination_open_stations[destination_station_order]
        
        for st in destination_open_stations
            #the blink mechanism
            if rand(rng) < γ
                continue
            end
          
            trip_id = findfirst(x -> x.origin_station == get_potential_locations()[origin_station] &&
                                    x.destination_station == get_potential_locations()[st], eachrow(curr_req_trips))
            
            if isnothing(trip_id)
                #there is no feasible path origin and dst stations
                continue
            end

            trip = curr_req_trips[trip_id, :]

            this_st_can_serve, this_st_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)
           
            if !destination_can_serve && this_st_can_serve
                #temporally accept this station
                destination_can_serve = true
                destination_station = st
                destination_new_cars = this_st_new_cars

                if this_st_new_cars < sol.initial_cars_number[st]
                    #accept directly current station
                    break
                end
            elseif this_st_can_serve 
                # here we found befor a feasible station, however we did not gain cars
                # so we need to check if selecting this stations is better "accepted when no need to decrease the cars number"
                if this_st_new_cars == sol.initial_cars_number[st]
                    destination_station = st
                    destination_new_cars = this_st_new_cars
                    break
                end
            end
        end

        if !destination_can_serve
            # try to open new station
            destination_stations_to_open = destination_stations_ids[findall(.!sol.open_stations_state[destination_stations_ids])]

            if !isempty(destination_stations_to_open)
                #there is(are) station(s) open
                stations_cost = [get_station_total_cost(st) for st in scenario.stations[destination_stations_to_open]]

                destination_stations_to_open = destination_stations_to_open[sortperm(stations_cost)]
                
                for st in destination_stations_to_open
                    #the blink mechanism
                    if rand(rng) < γ
                        continue
                    end
                    
                    # we need to check if a trip origin -> st is a valid trip
                    trip_id = findfirst(x -> x.origin_station == get_potential_locations()[origin_station] &&
                                    x.destination_station == get_potential_locations()[st], eachrow(curr_req_trips))
                    if isnothing(trip_id)
                        #there is no feasible path origin and dst stations
                        #@info "no feasible path origin and dst stations"
                        continue
                    end

                    destination_can_serve, destination_station, destination_new_cars = true, st, 0
                    break
                end
            end
        end

        if !destination_can_serve
            #we can not serve the request from the destination stations either by already open stations or by opening a new station
            continue
        end
        
        
        # here we can serve the request

        sol.open_stations_state[origin_station] = true
        sol.initial_cars_number[origin_station] = origin_new_cars
        sol.open_stations_state[destination_station] = true
        sol.initial_cars_number[destination_station] = destination_new_cars
        
        trip_id = findfirst(x -> x.origin_station == get_potential_locations()[origin_station] && x.destination_station == get_potential_locations()[destination_station], eachrow(curr_req_trips))
        #println("$(req.reqId), $origin_station, $destination_station")
        sol.selected_paths[req.scenario_id][curr_req_trips.fp_id[trip_id]] = true
    end
    return sol
end

function get_unserved_requests(sol::Solution, scenario_list::Vector{Scenario})
    global feasible_requests_masks
    
    unserved_requests = []

    for scenario in scenario_list
        #scenario = scenario_list[2];
        curr_sc_served_requests_mask = falses(nrow(scenario.request_list))
        curr_sc_served_requests_mask[scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :].req] .= true
    
        curr_sc_unserved_requests_ids = findall(xor.(feasible_requests_masks[scenario.scenario_id], curr_sc_served_requests_mask))
        curr_sc_unserved_requests = scenario.request_list[curr_sc_unserved_requests_ids, :]
        curr_sc_unserved_requests.scenario_id .= scenario.scenario_id
        push!(unserved_requests, curr_sc_unserved_requests)
    end

    return vcat(unserved_requests...)
   
end

function ruin_recreate(sol)
    new_sol = deepcopy(sol)
    adjacent_ruin!(new_sol)
    greedy_recreate!(new_sol)
    return new_sol
end
