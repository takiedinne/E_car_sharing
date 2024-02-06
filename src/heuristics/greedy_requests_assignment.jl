export greedy_assign_requests


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
    @info "adjacent_stations is filled"
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
            #=  origin_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in origin_open_stations]
            origin_station_order = sortperm([origin_stations_bounds[i][2] - origin_stations_bounds[i][1] for i in eachindex(origin_stations_bounds)], rev=true)
            =#         
            origin_station_order = shuffle(1:length(origin_open_stations))
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
            #= destination_stations_bounds = [station_all_scenario_bounds(sol, scenario_list, st) for st in destination_open_stations]
            destination_station_order = sortperm([destination_stations_bounds[i][2] - destination_stations_bounds[i][1] for i in eachindex(destination_stations_bounds)], rev=true)
             =#
            destination_station_order = shuffle(1:length(destination_open_stations))
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
