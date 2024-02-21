export ruin_recreate

global adjacent_stations
global ruin_depth = 0.15 # the percentage of stations to be closed
global γ = 0.05 # the blink probability
global use_adjacent_selection = true

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
    global rng
    global scenario_list
    global request_feasible_trips_ids
    global γ
    #sol = load_sol("sol.jls")
    #get unserved requests

    unserved_requests = get_unserved_requests(sol, scenario_list)
    # sort according to revenue and number of feasible paths and random
    sorting_probabilities = [1.0, 0.0, 0.0]
    p = rand(rng)
    if p < sorting_probabilities[1]
        sort!(unserved_requests, [:Rev], rev=[true])
    elseif p < sorting_probabilities[1] + sorting_probabilities[2]
        unserved_requests.nbr_fp = map(length, unserved_requests.fp)
        sort!(unserved_requests, [:nbr_fp, :Rev], rev=[false, true])
    else
        shuffle!(unserved_requests)
    end

    for req in eachrow(unserved_requests)
       
        #req = unserved_requests[1, :]
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
        for ost in origin_open_stations
            #ost = origin_open_stations[1]
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
        for st in destination_open_stations
            #st = destination_open_stations[1]
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
        #scenario = scenario_list[1];
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
    #ruin!(new_sol)

    greedy_recreate!(new_sol)
    #greedy_cars_use_recreate!(new_sol)
    return new_sol
end
#############################################@

function greedy_cars_use_recreate!(sol)
    global rng
    global scenario_list
    global request_feasible_trips_ids
    global γ

    #sol = load_sol("sol.jls")
    #get unserved requests

    unserved_requests = get_unserved_requests(sol, scenario_list)
    # sort according to revenue and number of feasible paths and random
    sort!(unserved_requests, [:Rev], rev=[true])

    request_ids = unserved_requests.reqId
    for req_id in request_ids

        req = filter(x -> x.reqId == req_id, unserved_requests)[1, :]
        #@info "greedy requests : $(req.reqId)"
        #req = unserved_requests[2, :]
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

        for ost in origin_open_stations
            #ost = origin_open_stations[1]
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
                break
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

        for st in destination_open_stations
            #st = destination_open_stations[1]
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
                break
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

        sol.selected_paths[req.scenario_id][curr_req_trips.fp_id[trip_id]] = true

        filter!(x -> x != curr_req_trips.req[trip_id], request_ids)

        #=  if E_carsharing_sim(sol) > 100000
            @info "1"
            break
        end =#
        continue_from_station!(scenario_list, scenario.scenario_id, sol, request_ids, destination_station)
        #=  if E_carsharing_sim(sol) > 100000
             @info "2"
             save_sol(sol, "sol12.jls")
             break
         end =#
    end

    return sol
end

function continue_from_station!(scenario_list, scenario_id, sol, unserved_requests, station_id)
    #save_sol(sol, "sol.jls")
    #sol , scenario_id, station_id = load_sol("sol.jls"), 1, 83

    scenario = scenario_list[scenario_id]
    # println("##########################")
    while true
        #@info "Continue from station : $(station_id)"
        station_node_id = get_potential_locations()[station_id]

        #@info "Continue from station : $(station_id)"
        # get the unserved request where the origin station is the current station
        curr_station_trips_ids = station_trips_ids[scenario.scenario_id][station_id]

        # normally I should change the filter to take into account the arriving time of the car from last request 
        curr_station_unserved_trips = filter(x -> x.req in unserved_requests && x.origin_station == station_node_id, scenario.feasible_paths[curr_station_trips_ids, :])

        # sort the trips according to their revenues
        sort!(curr_station_unserved_trips, [:Rev], rev=[true])

        curr_station_unserved_requests = scenario.request_list[unique(curr_station_unserved_trips.req), :]


        next_station = nothing
        for req in eachrow(curr_station_unserved_requests)
            #req = curr_station_unserved_requests[1, :]

            curr_req_trips = curr_station_unserved_trips[curr_station_unserved_trips.req.==req.reqId, :]

            curr_st_can_serve, curr_st_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, station_id, curr_req_trips[1, :])

            if !curr_st_can_serve
                # see the next trip 
                continue
            end

            # step 2: check the destination stations
            destination_stations_ids = [locations_dict[x] for x in curr_req_trips.destination_station]

            destination_open_stations = destination_stations_ids[findall(sol.open_stations_state[destination_stations_ids])]
            destination_can_serve, destination_new_cars, destination_station = false, -1, -1

            for st in destination_open_stations
                #st = destination_open_stations[1]
                #the blink mechanism
                if rand(rng) < γ
                    continue
                end

                trip_id = findfirst(curr_req_trips.destination_station .== get_potential_locations()[st])
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
                        trip_id = findfirst(x -> x.destination_station == get_potential_locations()[st], eachrow(curr_req_trips))
                        trip = curr_req_trips[trip_id, :]


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
            sol.initial_cars_number[station_id] = curr_st_new_cars
            sol.open_stations_state[destination_station] = true
            sol.initial_cars_number[destination_station] = destination_new_cars

            trip_id = findfirst(x -> x.origin_station == station_node_id &&
                    x.destination_station == get_potential_locations()[destination_station],
                eachrow(curr_req_trips))
            #println("$(req.reqId), $origin_station, $destination_station")
            sol.selected_paths[scenario.scenario_id][curr_req_trips.fp_id[trip_id]] = true

            #= if E_carsharing_sim(sol) > 100000
                @info "2"
                break
            end =#
            #delete the served request
            filter!(x -> x != curr_req_trips.req[trip_id], unserved_requests)

            next_station = destination_station
            break
        end

        if isnothing(next_station)
            break
        else
            station_id = next_station
        end

    end
end


function ruin!(sol)
    
    # global variables
    global rng
    global scenario_list
    
    served_trips = [scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :] for scenario in scenario_list]
    for i in eachindex(served_trips)
        served_trips[i].scenario_id .= i
    end
    served_trips = vcat(served_trips...)
    nbr_of_trips_to_delete = calculate_nbr_of_trips_to_delete(served_trips)
    #@info "nbr_of_trips_to_delete : $nbr_of_trips_to_delete"
    nbr_deleted_trips = 0
    
    trips_to_delete =  rand(rng, 1:nrow(served_trips), nbr_of_trips_to_delete)
    for trip_to_delete_id in trips_to_delete
        
        trip_to_delete_id = rand(rng, 1:nrow(served_trips))
        trip_to_delete = served_trips[trip_to_delete_id, :]

        # delete the trip
        sol.selected_paths[trip_to_delete.scenario_id][trip_to_delete.fp_id] = false
       
        stations_to_recheck = [locations_dict[trip_to_delete.origin_station], locations_dict[trip_to_delete.destination_station]]
        deleted_requests = [trip_to_delete.req]
        for station_id in stations_to_recheck
            # station_id = stations_to_recheck[2]
            additional_station_to_recheck, requests_to_lost_from_station = recheck_station!(sol, scenario_list[trip_to_delete.scenario_id], station_id)
            length(additional_station_to_recheck) > 0 && push!(stations_to_recheck, additional_station_to_recheck...)
            length(deleted_requests) > 0 && push!(deleted_requests, requests_to_lost_from_station...)
        end

        filter!(x -> !(x.req in deleted_requests), served_trips)
        nbr_deleted_trips += length(deleted_requests)
    end
    
    #check if there is a station that does not serve any request so we close it
    used_stations = [locations_dict[x] for x in unique(union(served_trips.origin_station, served_trips.destination_station)) ]
    used_st_mask = falses(length(sol.open_stations_state))
    used_st_mask[used_stations] .= true

    open_stations = findall(sol.open_stations_state)
    open_st_mask = falses(length(sol.open_stations_state))
    open_st_mask[open_stations] .= true

    sol.open_stations_state = used_st_mask .& open_st_mask


    #@info "Deleted trips : $deleted_requests"
    return sol

end

function calculate_nbr_of_trips_to_delete(trips::DataFrame)
    global ruin_depth
    max_trips_nbr = round(Int64, nrow(trips) * 0.3) + 1

    return rand(rng, 1:max_trips_nbr)
end

global adjacent_requests= []

function fill_adjacent_requests()
    # for each station: get a list od the station sorted by their distance
    scenario = scenario_list[1];
    global adjacent_requests = Matrix{Int64}(undef, 1000, 1000)
    for i in 1:1000
        
        rq = scenario.request_list.ON[i][1]
        distances =[get_walking_distance(rq, req1[1]) for req1 in scenario.request_list.ON]
        adjacent_requests[i, :] = sortperm(distances)
    end
    #delete the first column
    adjacent_requests = adjacent_requests[:, 2:end]
end