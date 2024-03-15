export ruin_recreate

global prob_sorting = [0.4, 0.2, 0.3, 0.1]
global adjacent_stations
global ruin_depth = 0.05 # the percentage of stations to be closed
global γ = 0.00 # the blink probability
global use_adjacent_selection = true

global number_of_lost_requests = DataFrame(calculated_nbr = Int64[], real_nbr = Int64[])
global stations_requests = []
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
    st_seed = open_stations[rand(rng, 1:length(open_stations))]

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
    return sol
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
    global prob_sorting
    #get unserved requests
    #= rng = MersenneTwister(1584)
    sol = greedy_assign_requests()[1]
    adjacent_ruin!(sol) 
    =#

    unserved_requests = get_unserved_requests(sol, scenario_list)
    # sort according to revenue and number of feasible paths and random

    p = rand(rng)
    if p < prob_sorting[1]
        sort!(unserved_requests, [:Rev], rev=[true])
    elseif p < prob_sorting[1] + prob_sorting[2]
        sort!(unserved_requests, [:ST], rev=[false])
    elseif p < prob_sorting[1] + prob_sorting[2] + prob_sorting[3]
        #sort according to the number of feasible paths
        unserved_requests.nbr_fps = [length(request_feasible_trips_ids[req.scenario_id][req.reqId]) for req in eachrow(unserved_requests)]
        sort!(unserved_requests, [:nbr_fps, :Rev], rev=[false, true])
    else
        shuffle!(rng, unserved_requests)
    end

    for req in eachrow(unserved_requests)

        #req = unserved_requests[4, :]
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
        # scenario = scenario_list[1];
        
        curr_sc_served_requests_mask = falses(nrow(scenario.request_list))
        curr_sc_served_requests_mask[scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :].req] .= true
        
        curr_sc_unserved_requests_ids = findall(xor.(feasible_requests_masks[scenario.scenario_id], curr_sc_served_requests_mask))
    
        curr_sc_unserved_requests = scenario.request_list[curr_sc_unserved_requests_ids, :]
    
        curr_sc_unserved_requests.scenario_id .= scenario.scenario_id
    
        push!(unserved_requests, curr_sc_unserved_requests)
        
       
    end

    return vcat(unserved_requests...)
    
end

function get_served_requests(sol::Solution, scenario_list::Vector{Scenario})
    global feasible_requests_masks
    served_requests = []

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
############### Ruin Recreate main function ################
function ruin_recreate(sol)
    
    new_sol = deepcopy(sol)
    #adjacent_ruin!(new_sol)
    adjacent_requests_ruin!(new_sol)
    greedy_recreate!(new_sol)
    
    return new_sol
end
#############################################
"""
    Greedy recreate function that continue from the stations of the distination of the selected requests
    it uses the Function continue_from_station()
"""
function recreate_with_continuing!(sol)
    global rng
    global scenario_list
    global request_feasible_trips_ids
    global γ
    #sol = load_sol("sol.jls")
    #get unserved requests

    unserved_requests = get_unserved_requests(sol, scenario_list)
    # sort according to revenue and number of feasible paths and random
    p = rand(rng)
    if p < prob_sorting[1]
        sort!(unserved_requests, [:Rev], rev=[true])
    elseif p < prob_sorting[1] + prob_sorting[2]
        sort!(unserved_requests, [:ST], rev=[false])
    elseif p < prob_sorting[1] + prob_sorting[2] + prob_sorting[3]
        #sort according to the number of feasible paths
        unserved_requests.nbr_fps = [length(request_feasible_trips_ids[req.scenario_id][req.reqId]) for req in eachrow(unserved_requests)]
        sort!(unserved_requests, [:nbr_fps, :Rev], rev=[false, true])
    else
        shuffle!(rng, unserved_requests)
    end

    while !isempty(unserved_requests)
        req = unserved_requests[1, :]
        #@info " iteration $i processing request $(req.reqId)"

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
            deleteat!(unserved_requests, 1)
            continue
        end

        #step 2: select a destination station
        destination_can_serve, destination_new_cars, destination_station = false, -1, -1
        best_cost = Inf64
        best_sol = nothing
        best_requests_ids = nothing
        selected_trip = nothing
        #sort the stations by their ristricted cars bounds
        for st in destination_stations_ids
            #st = destination_stations_ids[1]
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
            old_station_state = sol.open_stations_state[st]
            sol.open_stations_state[st] = true
            this_st_can_serve, this_st_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)
            sol.open_stations_state[st] = old_station_state

            if this_st_can_serve
                #calculate the cost 
                new_sol = deepcopy(sol)

                # the destination new configuration
                new_sol.open_stations_state[st] = true
                new_sol.initial_cars_number[st] = this_st_new_cars

                new_sol.open_stations_state[origin_station] = true
                new_sol.initial_cars_number[origin_station] = origin_new_cars

                new_sol.selected_paths[req.scenario_id][trip.fp_id] = true

                requests_ids = filter(x -> x.scenario_id == req.scenario_id, unserved_requests).reqId
                filter!(x -> x != req.reqId, requests_ids)

                continue_from_station!(scenario_list, req.scenario_id, new_sol, requests_ids, st, trip.arriving_time)

                new_cost = ECS_objective_function(new_sol)

                if (!destination_can_serve || new_cost < best_cost)

                    destination_can_serve = true
                    destination_station = st
                    destination_new_cars = this_st_new_cars
                    selected_trip = trip.fp_id
                    best_cost = new_cost
                    best_sol = deepcopy(new_sol)
                    best_requests_ids = copy(requests_ids)

                end

            end
        end

        if !destination_can_serve
            #we can not serve the request from the destination stations either by already open stations or by opening a new station
            deleteat!(unserved_requests, 1)
            continue
        end


        #@info "[greedy recreate] assigned request $(req.reqId)"
        # here we can serve the request
        sol.open_stations_state = best_sol.open_stations_state
        sol.initial_cars_number = best_sol.initial_cars_number
        sol.selected_paths = best_sol.selected_paths


        #delete the served requests
        filter!(x -> x.reqId ∈ best_requests_ids, unserved_requests)
    end

    return sol
end

function continue_from_station!(scenario_list, scenario_id, sol, unserved_requests_ids, station_id, curr_time)
    global station_trips_ids
    #= sol = load_sol("sol.jls")
    scenario_id, station_id, curr_time = 2, 4, 113 
    unserved_requests_ids = get_unserved_requests(sol, scenario_list).reqId=#
    scenario = scenario_list[scenario_id]
    #@info "continue from station $station_id at scenario $scenario_id at time $time"
    #is_feasible_solution(sol)

    while true

        #step 1: get all the unserved requests where their origin station is the current station
        station_node_id = get_potential_locations()[station_id]
        curr_station_all_trips_ids = station_trips_ids[scenario.scenario_id][station_id]

        curr_station_unserved_trips = filter(x -> x.req in unserved_requests_ids &&
                                                      x.origin_station == station_node_id &&
                                                      x.start_driving_time > curr_time,
            scenario.feasible_paths[curr_station_all_trips_ids, :])

        # sort the trips according to their revenues
        sort!(curr_station_unserved_trips, [:Rev], rev=[true])
        candidate_requests_ids = unique(curr_station_unserved_trips.req)

        next_station = nothing

        for req in candidate_requests_ids

            curr_req_trips = curr_station_unserved_trips[curr_station_unserved_trips.req.==req, :]
            curr_st_can_serve, curr_st_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, station_id, curr_req_trips[1, :])

            if !curr_st_can_serve
                # see the next trip 
                continue
            end

            # step 2: check the destination stations
            destination_stations_ids = [locations_dict[x] for x in curr_req_trips.destination_station]

            destination_open_stations = destination_stations_ids[findall(sol.open_stations_state[destination_stations_ids])]
            destination_can_serve, destination_new_cars, destination_station = false, -1, -1
            sol_trip_id = nothing

            # sort the destination stations by the nbr of potential trips that start after the arriving time of the current trip
            #TO DO
            for st in destination_open_stations
                # st = destination_open_stations[1]
                trip_id = findfirst(curr_req_trips.destination_station .== get_potential_locations()[st])
                trip = curr_req_trips[trip_id, :]

                this_st_can_serve, this_st_new_cars = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, st, trip)

                if this_st_can_serve
                    destination_can_serve, destination_new_cars, destination_station = true, this_st_new_cars, st
                    sol_trip_id = trip.fp_id
                    break
                end
            end

            # we don't try to open a new station
            if !destination_can_serve
                #we can not serve the request from the destination stations either by already open stations or by opening a new station
                continue
            end

            # here we can serve the request
            sol.initial_cars_number[station_id] = curr_st_new_cars
            sol.initial_cars_number[destination_station] = destination_new_cars
            sol.open_stations_state[destination_station] = true #eventhough it is already open
            sol.selected_paths[scenario.scenario_id][sol_trip_id] = true

            #delete the served request
            filter!(x -> x != req, unserved_requests_ids)
            #@info "[continue] we are assigning request $req"
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

#############################################
"""
    Adjacent Rui thatn not close the stations however it delete some requests from the stations
    and if the all the requests of the stations are deleted then the station is closed
"""
function adjacent_requests_ruin!(sol::Solution)
    #sol = generate_random_solution()
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
    n_requests = calculate_requests_to_delete(sol)

    #select the station that is used as seed for the ruin
    st_seed = open_stations[rand(rng, 1:length(open_stations))]

    #get the adjacent stations
    deleted_requests_counter = 0 #the index of stations to close
    stations_to_recheck = Int64[]
    adjacent_station_counter = 1
    st = st_seed
    lost_requests = [[] for _ in eachindex(scenario_list)]
    while deleted_requests_counter < n_requests
        trips = []
        for scenario in scenario_list
            #scenario = scenario_list[1]
            curr_scenario_trips_ids = station_trips_ids[scenario.scenario_id][st][findall(sol.selected_paths[scenario.scenario_id][station_trips_ids[scenario.scenario_id][st]])]
            curr_trips = scenario.feasible_paths[curr_scenario_trips_ids, :]
            curr_trips.scenario_id .= scenario.scenario_id
            push!(trips, curr_trips)
        end
        trips = vcat(trips...)

        if nrow(trips) < n_requests - deleted_requests_counter
            #here delete all the trips that start from this station
            for trip in eachrow(trips)
                #trip = trips[1, :]
                sol.selected_paths[trip.scenario_id][trip.fp_id] = false
                trip.origin_station == get_potential_locations()[st] ? push!(stations_to_recheck, locations_dict[trip.destination_station]) : push!(stations_to_recheck, locations_dict[trip.origin_station])
                deleted_requests_counter += 1
                push!(lost_requests[trip.scenario_id], trip.req)
            end
            sol.open_stations_state[st] = false
            sol.initial_cars_number[st] = 0
        else
            trips_ids = sample(rng, 1:nrow(trips), n_requests - deleted_requests_counter, replace=false)
            for trip in eachrow(trips[trips_ids, :])
                #trip = trips[1, :]
                sol.selected_paths[trip.scenario_id][trip.fp_id] = false
                trip.origin_station == get_potential_locations()[st] ? push!(stations_to_recheck, locations_dict[trip.destination_station]) : push!(stations_to_recheck, locations_dict[trip.origin_station])
                deleted_requests_counter += 1
                push!(lost_requests[trip.scenario_id], trip.req)
            end
        end

        st = adjacent_stations[st_seed, adjacent_station_counter]
        adjacent_station_counter += 1
    end
    #E_carsharing_sim(sol)
    #recheck the stations that are affected by the ruin
    unique!(stations_to_recheck)
    for station_id in stations_to_recheck
        for scenario in scenario_list
            # station_id = stations_to_recheck[4]
            additional_station_to_recheck, requests_to_lost_from_station = recheck_station!(sol, scenario, station_id)
            length(additional_station_to_recheck) > 0 && push!(stations_to_recheck, additional_station_to_recheck...)
            length(requests_to_lost_from_station) > 0 && push!(lost_requests[scenario.scenario_id], requests_to_lost_from_station...)
        end
    end

    #global number_of_lost_requests
    push!(number_of_lost_requests, (n_requests, sum([length(lost_requests[i]) for i in eachindex(lost_requests)])))
end

function calculate_requests_to_delete(sol)
    global ruin_depth
    max_request_nbr_to_delete = sum([nrow(scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :]) for scenario in scenario_list]) * ruin_depth
    max_request_nbr_to_delete = round(Int64, max_request_nbr_to_delete) + 1

    return rand(rng, 1:max_request_nbr_to_delete)
end


#############################################
"""
    Implement teh idea of selecting the destination stations accoording to the number of requests that can have in the future
"""
function greedy_recreate_with_future_requests!(sol)
    global rng
    global scenario_list
    global request_feasible_trips_ids
    global γ
    global prob_sorting
    #get unserved requests
    #= 
    rng = MersenneTwister(1584)
    sol = greedy_assign_requests()[1]
    adjacent_ruin!(sol) 
    =#

    unserved_requests = get_unserved_requests(sol, scenario_list)
    # sort according to revenue and number of feasible paths and random

    p = rand(rng)
    if p < prob_sorting[1]
        sort!(unserved_requests, [:Rev], rev=[true])
    elseif p < prob_sorting[1] + prob_sorting[2]
        sort!(unserved_requests, [:ST], rev=[false])
    elseif p < prob_sorting[1] + prob_sorting[2] + prob_sorting[3]
        #sort according to the number of feasible paths
        unserved_requests.nbr_fps = [length(request_feasible_trips_ids[req.scenario_id][req.reqId]) for req in eachrow(unserved_requests)]
        sort!(unserved_requests, [:nbr_fps, :Rev], rev=[false, true])
    else
        shuffle!(rng, unserved_requests)
    end

    for req in eachrow(unserved_requests)

        #req = unserved_requests[4, :]
        #get origin stations for the request
        request_trips_ids = request_feasible_trips_ids[req.scenario_id][req.reqId]
        curr_req_trips = scenario_list[req.scenario_id].feasible_paths[request_trips_ids, :]
        scenario = scenario_list[req.scenario_id];
        
        is_accepted = Dict{Int64, Tuple{Bool, Int64}}()
        accepted_trips = []
        for trip in eachrow(curr_req_trips)
            #trip = curr_req_trips[1, :]
            #blinking mechanism
            if rand(rng) < γ
                continue
            end
            
            ori_st = locations_dict[trip.origin_station]
            des_st = locations_dict[trip.destination_station]
            
            # step 1: check is the origin can serve
            if !haskey(is_accepted, ori_st)
                # we need to check
                if !sol.open_stations_state[ori_st]
                    is_accepted[ori_st] = (true, 1)
                else
                    is_accepted[ori_st] = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, ori_st, trip)
                end
            end
            # Here we are 100% sure that we checked this station
            if !is_accepted[ori_st][1]
                continue
            end

            # step 2: check is the destination can serve
            if !haskey(is_accepted, des_st)
                # we need to check
                if !sol.open_stations_state[des_st]
                    is_accepted[des_st] = (true, 0)
                else
                    is_accepted[des_st] = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, sol, des_st, trip)
                end
            end

            if !is_accepted[des_st][1]
                continue
            end

            # add the trip to the list of accepted trips
            push!(accepted_trips, trip)

        end
        #!isempty(accepted_trips) && @info " req $(req.reqId) could be served"
        
        if isempty(accepted_trips)
            continue
        end

        # select which requests
        # First calculate 
        heur_values = []
        prop_dict = Dict{Int64, Float64}()
        for trip in accepted_trips
            ori_st = locations_dict[trip.origin_station]
            des_st = locations_dict[trip.destination_station]
            #!haskey(prop_dict, ori_st) && (prop_dict[ori_st] = get_station_trips_proportion(sol, scenario, ori_st, req.reqId))
            !haskey(prop_dict, des_st) && (prop_dict[des_st] = get_station_trips_proportion(sol, scenario, des_st, req.reqId))
            
            hv = sol.open_stations_state[ori_st] ? 2 : 0
            hv += sol.open_stations_state[des_st] ? 2 : 0
            hv += sol.initial_cars_number[ori_st] == is_accepted[ori_st][2] ? 1 : 0
            hv += sol.initial_cars_number[des_st] < is_accepted[des_st][2] ? 1 : 0
            hv += #= (1 - prop_dict[ori_st]) + =# prop_dict[des_st]

            push!(heur_values, hv)
        end

        #selecte the trip with the maximum heuristic value
        selected_trip = accepted_trips[argmax(heur_values)]
        origin_station = locations_dict[selected_trip.origin_station]
        destination_station = locations_dict[selected_trip.destination_station]

        sol.open_stations_state[origin_station] = true
        sol.initial_cars_number[origin_station] = is_accepted[origin_station][2]
        sol.open_stations_state[destination_station] = true
        sol.initial_cars_number[destination_station] = is_accepted[destination_station][2]

        sol.selected_paths[req.scenario_id][selected_trip.fp_id] = true
    end

    return sol
end

function get_nbr_future_requests_from(sol, scenario, station_id, req_id)
    global station_trips_ids

    #get the trips ids of the stations
    trips_ids = station_trips_ids[scenario.scenario_id][station_id]
    trips_ids = trips_ids[.!sol.selected_paths[scenario.scenario_id][trips_ids]]
    
    #calculate the time point
    station_node_id = get_potential_locations()[station_id]
    interesting_trip_id = findfirst(scenario.feasible_paths[trips_ids, :req] .== req_id #= .&&
                                    scenario.feasible_paths[trips_ids, :origin_station] .== station_node_id =#)
    interesting_trip = scenario.feasible_paths[trips_ids[interesting_trip_id], :]

    is_origin = (interesting_trip.origin_station == station_node_id)
    time_point = is_origin ? interesting_trip.start_driving_time : interesting_trip.arriving_time
   
    trips_ids = trips_ids[scenario.feasible_paths[trips_ids, :start_driving_time] .> time_point .&&
                          scenario.feasible_paths[trips_ids, :origin_station] .== station_node_id]

    nbr_requests = length(unique(scenario.feasible_paths[trips_ids, :].req))
    return nbr_requests
end

#STUDENT-22EU20kqt18903
function get_station_trips_proportion(sol, scenario, station_id, req_id)
    global station_trips_ids

    #get the unserved trips ids of the stations
    trips_ids = station_trips_ids[scenario.scenario_id][station_id]
    trips_ids = trips_ids[.!sol.selected_paths[scenario.scenario_id][trips_ids]]
    
    #calculate the time point from where we should calculate the proportion
    station_node_id = get_potential_locations()[station_id]
    interesting_trip_id = findfirst(scenario.feasible_paths[trips_ids, :req] .== req_id)
    interesting_trip = scenario.feasible_paths[trips_ids[interesting_trip_id], :]
    is_origin = (interesting_trip.origin_station == station_node_id)
    time_point = is_origin ? interesting_trip.start_driving_time : interesting_trip.arriving_time
   
     
    trips_from_ids = trips_ids[ scenario.feasible_paths[trips_ids, :req] .!= req_id .&& 
                            scenario.feasible_paths[trips_ids, :start_driving_time] .>= time_point .&&
                            scenario.feasible_paths[trips_ids, :origin_station] .== station_node_id]
    trips_to_ids = trips_ids[ scenario.feasible_paths[trips_ids, :req] .!= req_id .&& 
                            scenario.feasible_paths[trips_ids, :arriving_time] .>= time_point .&&
                            scenario.feasible_paths[trips_ids, :destination_station] .== station_node_id]

    # propotion represents the percentage of departure trips from the station aftre time_point
    #if no further trips then the propotion is 0 to encourage taking this request frim this station
    isempty(trips_from_ids) && isempty(trips_to_ids) && return -1.0 

    proportion = length(trips_from_ids) / (length(trips_from_ids) + length(trips_to_ids))
    return proportion
end

function create_tree!(sol, scenario, station_id, st_time=0)

    global stations_trips_ids
    served_requets = scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :req]
    #scenario  = scenario_list[1]
    g = MetaDiGraph()
    props_dict = Dict(:station_id => station_id,
    :time => st_time, 
    :solution => sol, 
    :served_requests => served_requets)
    
    add_vertex!(g, props_dict)
    
    #create a graph in BFS manner
    queue = [1]
    leaves_nodes = []
    while !isempty(queue)

        cur_node = popfirst!(queue)
        node_props = props(g, cur_node)
        cur_st, starting_time, cur_sol, served_requets = node_props[:station_id], node_props[:time], node_props[:solution], node_props[:served_requests]
                  
        cur_st_node_id = get_potential_locations()[cur_st]
        #get the unserved trips that start from the current station
        trips_ids = station_trips_ids[scenario.scenario_id][cur_st]
        trips_ids = trips_ids[.!cur_sol.selected_paths[scenario.scenario_id][trips_ids]]
        trips = scenario.feasible_paths[trips_ids, :]
        cur_st_unserved_trips = trips[trips.origin_station .== cur_st_node_id .&&
                                    trips.start_driving_time .>= starting_time, :] 
        
        # filter the trips of already served requests
        cur_st_unserved_trips = cur_st_unserved_trips[[req ∉ served_requets for req in cur_st_unserved_trips.req] , :]
        
        is_accepted = Dict{Int64, Tuple{Bool, Int64}}()
        nbr_acepted_trips = 0
        for trip in eachrow(cur_st_unserved_trips)
            #trip = cur_st_unserved_trips[1, :]
           
            if !haskey(is_accepted, trip.req)
                is_accepted[trip.req] = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, cur_sol, cur_st, trip)
            end
            ori_can_serve, ori_new_cars = is_accepted[trip.req]
            if !ori_can_serve
                continue
            end
            #check on the destination station
            des_st = locations_dict[trip.destination_station]
            des_can_serve, des_cars_nbr = can_serve_and_get_cars_number(scenario_list, scenario.scenario_id, cur_sol, des_st, trip)
            
            if !des_can_serve
                continue
            end

            #add the trip to the graph
            new_sol = deepcopy(cur_sol)
            new_sol.open_stations_state[cur_st] = true
            new_sol.initial_cars_number[cur_st] = ori_new_cars
            new_sol.open_stations_state[des_st] = true
            new_sol.initial_cars_number[des_st] = des_cars_nbr
            new_sol.selected_paths[scenario.scenario_id][trip.fp_id] = true
            
            node_props_dict = Dict(:station_id => locations_dict[trip.destination_station],
                              :time => trip.arriving_time + 1, #1 is the charging time
                              :solution => new_sol,
                              :served_requests => [served_requets; trip.req])
            
            add_vertex!(g, node_props_dict)
            new_node = nv(g)
            
            #cost of opening new stations
            stations_costs = 0 
            cur_sol.open_stations_state[cur_st] && (stations_costs += get_station_total_cost(scenario.stations[cur_st]))
            cur_sol.open_stations_state[des_st] && (stations_costs += get_station_total_cost(scenario.stations[des_st]))
            stations_costs /= cost_factor
            
            edge_props_dict = Dict(:trip_id => trip.fp_id, :weight => (-1 * trip.Rev + stations_costs))
            add_edge!(g, cur_node, new_node, edge_props_dict)

            push!(queue, new_node)
            nbr_acepted_trips += 1
        end

        if nbr_acepted_trips == 0
            push!(leaves_nodes, cur_node)
            continue
        end

    end
    
    #calculate the shortest path
    shortest_paths = bellman_ford_shortest_paths(g, 1) 
    distances = [shortest_paths.dists[i] for i in leaves_nodes]
    leaf_node_with_longest_rev = leaves_nodes[argmin(distances)]
    path = enumerate_paths(shortest_paths, leaf_node_with_longest_rev)
    trips_ids = []
    for i in eachindex(path)
        if i == 1
            continue
        end
        edge = Edge(path[i-1], path[i])
        trip_id = get_prop(g, edge, :trip_id)
        push!(trips_ids, trip_id)
    end

    new_sol = get_prop(g, leaf_node_with_longest_rev, :solution)
    
    #update the input solution
    sol.open_stations_state = new_sol.open_stations_state
    sol.initial_cars_number = new_sol.initial_cars_number
    sol.selected_paths = new_sol.selected_paths

    return new_sol, trips_ids
end

function recreate_with_best_tree_path!(sol)
    global rng
    global scenario_list
    global request_feasible_trips_ids
    global γ
   
    unserved_requests = get_unserved_requests(sol, scenario_list)
    # sort according to revenue and number of feasible paths and random

    p = rand(rng)
    if p < prob_sorting[1]
        sort!(unserved_requests, [:Rev], rev=[true])
    elseif p < prob_sorting[1] + prob_sorting[2]
        sort!(unserved_requests, [:ST], rev=[false])
    elseif p < prob_sorting[1] + prob_sorting[2] + prob_sorting[3]
        #sort according to the number of feasible paths
        unserved_requests.nbr_fps = [length(request_feasible_trips_ids[req.scenario_id][req.reqId]) for req in eachrow(unserved_requests)]
        sort!(unserved_requests, [:nbr_fps, :Rev], rev=[false, true])
    else
        shuffle!(rng, unserved_requests)
    end
    requests_list = unserved_requests.reqId
    
    for req in requests_list

    end
    return sol
end
