export applyHeuristic

global mutationHeuristics = [1, 2]
global localSearchHeuristics = [3, 4, 5]
global ruinRecreateHeuristics = [6, 7]
global crossoverHeuristics = [8]
""" 
    this file implments the heuristics for the carsharing problem
    all the heuristic need to have the sam schema of input and output
    function myheuristic(sol::Solution, kwargs...)
        #do something without altering the solution in input
        return new_sol, new_fit
    end
"""

"""
    Mutation heuristics:
    take a solution and number of mutation parameter, perform a random aleration to
    the soltion and return the new fit and solution
"""

"""
    flipStationState(sol::Solution; mutation_number::Int64=1)

    Select randomly a number of stations and flip their state (open -> close or close -> open)

    @sol: the solution to be mutated
    @mutation_number: the number of stations to flip

    return the new fit and solution
"""
function flipStationState(sol::Solution; mutation_number::Int64=1)
    
    global rng
    global online_request_serving
    
    # step1: copy the solution
    curr_sol = deepcopy(sol)
    
    # step2: get the list of station to mutation (open -> close or close -> open)
    station_to_mutate = sample(rng, 1:length(curr_sol.open_stations_state), mutation_number, replace=false)

    # step3: flip the state of the station
    closed_stations = Int64[]
    opened_stations = Int64[]
    for st in station_to_mutate
        curr_sol.open_stations_state[st] = !curr_sol.open_stations_state[st]

        # keep track of the stations that were closed and opened
        if curr_sol.open_stations_state[st]
            push!(opened_stations, st)
            
            # set the initial cars number to zero
            init_cars_nbr = floor(stations[st].max_number_of_charging_points / 2)
            curr_sol.initial_cars_number[st] = init_cars_nbr
        else
            push!(closed_stations, st)
            # set the initial cars number to zero
            curr_sol.initial_cars_number[st] = 0
        end
    end
    
    # treat the selected_paths if we are in offline mode
    
    if !online_request_serving
        # step4: clean up the solution
        if length(closed_stations) > 0 
            fitness, _ = clean_up_selected_paths!(curr_sol)
        end
        # step5: serve the requests
        if length(opened_stations) > 0
            fitness, _ = serve_requests_after_opening_station(curr_sol, opened_stations)
        end
    else
        fitness = E_carsharing_sim(curr_sol)
    end
    #return sol
    fitness, curr_sol
end

"""
    interchangeStations(sol::Solution; mutation_number::Int64=1)

    Select randomly a number of stations pairs, each pair consists of a closed station and an open station
    Then we close the open station and open the closed station

    @sol: the solution to be mutated
    @mutation_number: the number of stations to flip

    return the new fit and solution
"""
function interchangeStations(sol::Solution; mutation_number::Int64=1)
    # global variables
    global rng
    global online_request_serving

    # step1: copy the solution
    curr_sol = deepcopy(sol)
    
    opened_stations = findall(curr_sol.open_stations_state)
    closed_stations = findall(.!curr_sol.open_stations_state)
    
    # step2: make sure that mutation number is feasible
    mutation_number = min(mutation_number, length(opened_stations), length(closed_stations))

    # step3: get the list of station pairs to interchange
    station_to_open = sample(rng, 1:length(closed_stations), mutation_number, replace=false)
    station_to_close = sample(rng, 1:length(opened_stations), mutation_number, replace=false)
    
    # step4: interchange the stations states 
    curr_sol.open_stations_state[station_to_close] .= false
    curr_sol.initial_cars_number[station_to_close] .= 0
    curr_sol.open_stations_state[station_to_open] .= true
    curr_sol.initial_cars_number[station_to_open] .= floor.(
        [stations[st].max_number_of_charging_points for st in station_to_open] 
        ./ 2)

    # step5: treat the selected_paths if we are in offline mode
    if !online_request_serving
        # step6: clean up the solution
        fitness, _ = clean_up_selected_paths!(curr_sol)
        # step7: serve the requests
        fitness, _ = serve_requests_after_opening_station(curr_sol, station_to_open)
    else
        fitness = E_carsharing_sim(curr_sol)
    end
    
    return fitness, curr_sol
end

""" 
    Local search methods
    they take a solution and try to move to a better one
"""

"""
    greedyAdding(sol::Solution; local_search_depth::Int64=1)

    first improvement variable neighboorhoud local search. it tries to open a station in a greedy manner
    if no station can be opened, it tries to open two stations at a time

    @sol: the solution that we start the search from
    @local_search_depth: the number of iteration for the local search
"""
function greedyAdding(sol::Solution; local_search_depth::Int64=4)
    global online_request_serving
    
    #step1: copy the solution
    curr_sol = deepcopy(sol)
    curr_fit = E_carsharing_sim(curr_sol)
    
    #get the list of closed stations
    closed_stations = findall(.!curr_sol.open_stations_state)
    
    #sort the closed station by their number of requests
    requests_to_serve_per_station = [sum(
                                        [length(unique(get_trips_station(st, scenario_list[i].feasible_paths).req)) 
                                         for i in eachindex(scenario_list)])
                                    for st in closed_stations]
    closed_stations = closed_stations[sortperm(requests_to_serve_per_station, rev=true)]
    
    #make sure that the depth of search is less or equal to the number of closed stations
    local_search_depth = min(local_search_depth, length(closed_stations))
    opened_stations = Int64[] # keep track of the opened stations
    for _ in 1:local_search_depth
        find_better_solution = false
        for closed_st in closed_stations
            if closed_st in opened_stations
                continue #avoid the stations already decided to be open
            end
            # save the current solution (useful to reset)
            new_sol = deepcopy(curr_sol)
            # open the station
            new_sol.open_stations_state[closed_st] = true
            new_sol.initial_cars_number[closed_st] = floor(stations[closed_st].max_number_of_charging_points / 2)
            
            #evalute the new solution
            if !online_request_serving
                new_fit, _ = serve_requests_after_opening_station(new_sol, [closed_st]) 
            else
                new_fit = E_carsharing_sim(new_sol)
            end

            if new_fit < curr_fit
                curr_sol = new_sol
                curr_fit = new_fit
                push!(opened_stations, closed_st)
                
                find_better_solution = true
                break
            end
        end

        if !find_better_solution
            #we couldn't enhance the solution by opening a station so we stop the local search
            # to continue Add the second depth serearch
            for (closed_st1, closed_st2) in combinations(closed_stations, 2)
                
                if closed_st1 in opened_stations || closed_st2 in opened_stations
                    continue #avoid the stations already decided to be open
                end

                # save the current solution (useful to reset)
                new_sol = deepcopy(curr_sol)
                # open the station
                new_sol.open_stations_state[closed_st1] = true
                new_sol.open_stations_state[closed_st2] = true

                new_sol.initial_cars_number[closed_st1] = floor(stations[closed_st1].max_number_of_charging_points / 2)
                new_sol.initial_cars_number[closed_st2] = floor(stations[closed_st2].max_number_of_charging_points / 2)
                #evalute the new solution
                if !online_request_serving
                    new_fit, _ = serve_requests_after_opening_station(new_sol, [closed_st1, closed_st2]) 
                else
                    new_fit = E_carsharing_sim(new_sol)
                end

                if new_fit < curr_fit
                    curr_sol = new_sol
                    curr_fit = new_fit
                    push!(opened_stations, closed_st1)
                    push!(opened_stations, closed_st2)

                    find_better_solution = true
                    break
                end
            end
            # if we couldn't find a better solution by opening two stations, we stop the local search
            if !find_better_solution
                break
            end
        end
    end
    #delete the no used cars in the solution
    if !online_request_serving
        curr_fit = clean_up_cars_number!(curr_sol)
    end
    return curr_fit, curr_sol
end


"""
    greedyClosing(sol::Solution; local_search_depth::Int64=1)

    first improvement variable neighboorhoud local search. it tries to close a station in a greedy manner
    if no station can be closed, it tries to close two stations at a time

    @sol: the solution that we start the search from
    @local_search_depth: the number of iteration for the local search

    return the new fit and solution
"""
function greedyClosing(sol::Solution; local_search_depth::Int64=1)
    global online_request_serving
    
    #step1: copy the solution
    curr_sol = deepcopy(sol)
    curr_fit = E_carsharing_sim(curr_sol)
    #get the list of opened stations
    opened_stations = findall(curr_sol.open_stations_state)

    #make sure that the depth of search is less or equal to the number of opened stations
    local_search_depth = min(local_search_depth, length(opened_stations))

    #sort the stations by their number of requests
    requests_to_loose_per_station = [sum(
                                        [nrow(get_trips_station(st,
                                         scenario_list[i].feasible_paths[curr_sol.selected_paths[i], :])) 
                                         for i in eachindex(scenario_list)])
                                    for st in opened_stations]
    opened_stations = opened_stations[sortperm(requests_to_loose_per_station)]
    
    #step2: start the local search
    closed_stations = Int64[] # keep track of the closed stations
    for _ in 1:local_search_depth
        find_better_solution = false
        for opened_st in opened_stations
            
            if opened_st in closed_stations
                continue
            end

            # save the current solution (useful to reset)
            new_sol = deepcopy(curr_sol)
            
            # close the station
            new_sol.open_stations_state[opened_st] = false
            new_sol.initial_cars_number[opened_st] = 0

            #evalute the new solution
            if !online_request_serving 
                # get the requests that were served by the station
                requests_served_by_station = [get_trips_station(opened_st, 
                                                scenario_list[i].feasible_paths[curr_sol.selected_paths[i], :]).req 
                                                for i in 1:length(scenario_list)]
                #clean the selected paths
                new_fit, _ = clean_up_selected_paths!(new_sol)

                #try to serve the requests that were served by the station
                if sum(x->length(x), requests_served_by_station) > 0
                    new_fit, _, _ = serve_requests!(new_sol, requests_served_by_station)
                end
            else
                new_fit = E_carsharing_sim(new_sol)
            end

            if new_fit < curr_fit
                curr_sol = new_sol
                curr_fit = new_fit
                push!(closed_stations, opened_st)
                find_better_solution = true

                break
            end
        end
        if !find_better_solution
            #we couldn't enhance the solution by closing a station so we stop the local search
            break
        end
    end  
    
    return curr_fit, curr_sol  
end

"""
    addCarsLS(sol::Solution; local_search_depth::Int64=1)

    first improvement variable neighboorhoud local search. it tries to add a car to a station in a greedy manner

    @sol: the solution that we start the search from
    @local_search_depth: the number of iteration for the local search

    return the new fit and solution
"""
function addCarsLS(sol::Solution; local_search_depth::Int64=1)
    #global information
    global online_request_serving

    #step1: copy the solution
    curr_sol = deepcopy(sol)
    curr_fit = E_carsharing_sim(curr_sol)

    #get the list of opened stations
    opened_stations = findall(curr_sol.open_stations_state)
    #make sure that the depth of search is less or equal to the number of opened stations
    local_search_depth = min(local_search_depth, length(opened_stations))

    # see how can we enhance this heuristic by sorting the stations

    for _ in 1:local_search_depth
        find_better_solution = false
        for opened_st in opened_stations
            # save the current solution (useful to reset)
            new_sol = deepcopy(curr_sol)
            
            station_capacity = stations[opened_st].max_number_of_charging_points
            
            if new_sol.initial_cars_number[opened_st] >= station_capacity
                continue # we cannot add more cars to this station so we skip it
            end

            # add a car to the station
            new_sol.initial_cars_number[opened_st] += 1

            # evalute the new solution
            if !online_request_serving
                new_fit, _ = serve_requests_after_opening_station(new_sol, [opened_st])
            else
                new_fit = E_carsharing_sim(new_sol)
            end

            if new_fit < curr_fit
                curr_sol = new_sol
                curr_fit = new_fit
                find_better_solution = true
                break
            end
        end
        
        if !find_better_solution
            #we couldn't enhance the solution by adding a car so we stop the local search
            break
        end
    end

    return curr_fit, curr_sol
end

"""
    Ruin and recreate heuristics
    they perform on two phases: ruin where they dustruct the solution on hand 
    and the second phase where yhey recreate the solution accordding to some strategy
"""

"""
    FIFSSelectedPaths(sol::Solution; ruin_recreate_depth::Int64=1)

    it destroys the solution by closing a number of stations, Then it handles all the requests in a FIFO manner

    @sol: the solution that we start the search from
    @ruin_recreate_depth: the number of stations to close (Ruine phase) and open (Recreate)

    return the new fit and solution
"""
function FIFSSelectedPaths(sol::Solution; ruin_recreate_depth::Int64=1)
    # global variables
    global online_request_serving
    global rng
    global online_selected_paths # the selected paths that were used to serve the requests in the online mode

    #step1: copy the solution
    curr_sol = deepcopy(sol)

    opened_stations = findall(curr_sol.open_stations_state)
    
    #make sure that ruin_recreate_depth is less than the number of opened stations
    ruin_recreate_depth = min(ruin_recreate_depth, length(opened_stations))

    #step2: ruin phase
    stations_to_close = sample(rng, opened_stations, ruin_recreate_depth, replace=false)
    curr_sol.open_stations_state[stations_to_close] .= false
    curr_sol.initial_cars_number[stations_to_close] .= 0
    
    #step3: recreate phase
    closed_stations = findall(.!curr_sol.open_stations_state)
    stations_to_open = sample(rng, closed_stations, ruin_recreate_depth, replace=false)
    curr_sol.open_stations_state[stations_to_open] .= true
    curr_sol.initial_cars_number[stations_to_open] .= floor.(
        [stations[st].max_number_of_charging_points for st in stations_to_open] 
        ./ 2)

    #step4: recreate phase (handle the selected_paths)
    if !online_request_serving
       online_request_serving = true
       new_fit = E_carsharing_sim(curr_sol)
       curr_sol.selected_paths = deepcopy(online_selected_paths)
       online_request_serving = false
    else
        new_fit = E_carsharing_sim(curr_sol)
    end

    return new_fit, curr_sol
end

"""
    greedyRuinRecreate(sol::Solution; ruin_recreate_depth::Int64=1)

    it destroys the solution by closing a number of stations. Then it open the same number of stations
    in a greedy manner
"""
function greedyRuinRecreate(sol::Solution; ruin_recreate_depth::Int64=1)
    # global variables
    global online_request_serving
    global rng

    #step1: copy the solution
    curr_sol = deepcopy(sol)

    opened_stations = findall(curr_sol.open_stations_state)
    #make sure that ruin_recreate_depth is less than the number of opened stations
    ruin_recreate_depth = min(ruin_recreate_depth, length(opened_stations))
    
    #step2: ruin phase (randomly close the stations)
    stations_to_close = sample(rng, opened_stations, ruin_recreate_depth, replace=false)
    curr_sol.open_stations_state[stations_to_close] .= false
    curr_sol.initial_cars_number[stations_to_close] .= 0

    # clean up the selected paths if we are in offline mode
    if !online_request_serving
        new_fit, _ = clean_up_selected_paths!(curr_sol)
    end

    #step3: recreate phase (open the stations in a greedy manner)
    closed_stations = findall(.!curr_sol.open_stations_state)
    #sort the closed station by their number of requests
    requests_to_serve_per_station = [sum(
                                        [length(unique(get_trips_station(st, scenario_list[i].feasible_paths).req)) 
                                        for i in eachindex(scenario_list)])
                                    for st in closed_stations]
    closed_stations = closed_stations[sortperm(requests_to_serve_per_station, rev=true)]
    # we simply take the first stations as they are ordred according to their number of requests
    stations_to_open = closed_stations[1:ruin_recreate_depth]
    curr_sol.open_stations_state[stations_to_open] .= true
    curr_sol.initial_cars_number[stations_to_open] .= floor.(
                            [stations[st].max_number_of_charging_points for st in stations_to_open] 
                            ./ 2)
    #handle the selected paths
    if !online_request_serving
        new_fit, _ = serve_requests_after_opening_station(curr_sol, stations_to_open)
        new_fit = clean_up_cars_number!(curr_sol)
    else
        new_fit = E_carsharing_sim(curr_sol)
    end

    return new_fit, curr_sol
end

"""
    Crossover methods 

    they take two solutions and try to combine them to get a better solution
"""

"""
        oneXCrossover(sol1::Solution, sol2::Solution)

        it takes two solutions and combine them to get a new solution according one poin strategy
"""
function oneXCrossover(sol1::Solution, sol2::Solution)
    # global variables
    global online_request_serving
    global rng
    
    #step1: decide whether the order of parents
    if rand(rng) < 0.5
        parent1 = sol1
        parent2 = sol2
    else
        parent1 = sol2
        parent2 = sol1
    end

    xPoint = rand(rng, 1:length(parent1.open_stations_state))
    #create the child solution
    child_sol = deepcopy(parent1)
    child_sol.open_stations_state[1:xPoint] .= parent1.open_stations_state[1:xPoint]
    child_sol.open_stations_state[xPoint+1:end] .= parent2.open_stations_state[xPoint+1:end]
    child_sol.initial_cars_number[1:xPoint] .= parent1.initial_cars_number[1:xPoint]
    child_sol.initial_cars_number[xPoint+1:end] .= parent2.initial_cars_number[xPoint+1:end]

    if !online_request_serving
        new_fit, child_sol = FIFSSelectedPaths(child_sol, ruin_recreate_depth=0)
    else
        new_fit = E_carsharing_sim(child_sol)
    end
    return new_fit, child_sol
end

function applyHeuristic(heurId::Int64, sol::Solution, mutation_num::Int64, depth_of_search::Int64)
   
    if heurId == 1
        fit, new_sol = flipStationState(sol, mutation_number=mutation_num)
    elseif heurId == 2
        fit, new_sol = interchangeStations(sol, mutation_number=mutation_num)
    elseif heurId == 3
        fit, new_sol = greedyAdding(sol, local_search_depth=depth_of_search)
    elseif heurId == 4
        fit, new_sol = greedyClosing(sol, local_search_depth=depth_of_search)
    elseif heurId == 5
        fit, new_sol = addCarsLS(sol, local_search_depth=depth_of_search)
    elseif heurId == 6
        fit, new_sol = FIFSSelectedPaths(sol, ruin_recreate_depth=mutation_num)
    elseif heurId == 7
        fit, new_sol = greedyRuinRecreate(sol, ruin_recreate_depth=mutation_num)
    else
        error("heuristic id must be in $(mutationHeuristics) ∪ $(localSearchHeuristics) ∪ $(ruinRecreateHeuristics)")
    end
    return fit, new_sol
end

function applyHeuristic(heurId::Int64, sol1::Solution, sol2::Solution)
    
    if heurId == 8
        return oneXCrossover(sol1, sol2)
    else
        error("heuristic id must be in $(crossoverHeuristics)")
    end
end

###############################################################################################################
##################################### Utils functions #########################################################
###############################################################################################################

"""
    serve_requests_after_opening_station(sol::Solution, stations_idx::Array{Int64,1})
    
    this function tries to serve new requests after opening new station for each scenario
    inputs:
        @sol: the solution
        @station_id: the id of the station that we recently oponed
    outputs:
        basicaly it return the new seletced paths ad assigne it as well to the onling_selected_paths so we 
        could get it from there as well
"""
function serve_requests_after_opening_station(sol::Solution, stations_idx::Array{Int64,1})
    #@info "*********** serve_requests_after_opening_station ***************"
    #= @show stations_idx =#
    @assert !online_request_serving && all(x -> x, sol.open_stations_state[stations_idx]) "we are in wrong mode or the station is closed"
    
    # declare some global variables
    global online_request_serving
    global failed

    #get the nodes ids for the stations that we recently opened
    station_nodes_idx = get_potential_locations()[stations_idx]
    sol_stations_nodes_idx = get_potential_locations()[sol.open_stations_state]
    
    stations_to_increment_cars = Int64[]

    #loop over each scenario and try to serve new requests
    @threads for sc_id in eachindex(scenario_list)
        #sc_id = 1
        scenario = scenario_list[sc_id] # handle one scenario a time

        # get the already served requests 
        served_requests = scenario.feasible_paths.req[sol.selected_paths[sc_id]]

        potential_feasible_paths = filter(scenario.feasible_paths) do fp
            (fp.origin_station in station_nodes_idx && fp.destination_station in sol_stations_nodes_idx) ||
                (fp.destination_station in station_nodes_idx && fp.origin_station in sol_stations_nodes_idx)
        end

        # keep only the requests that are not served yet
        filter!(potential_feasible_paths) do fp
            fp.req ∉ served_requests
        end

        #sort the requests according to their revenue
        sort!(potential_feasible_paths, [:Rev, :req], rev=true)

        new_served_requests = Int64[]

        for curr_fp in eachrow(potential_feasible_paths)
            #curr_fp = potential_feasible_paths[2, :]
            #check if we already served the request
            curr_fp.req in new_served_requests && continue
            
            #check if we can serve the request without trying to serve it
            can_serve, can_serve_if_increment_cars = can_use_trip(sc_id, curr_fp, sol)
            if can_serve
                sol.selected_paths[sc_id][curr_fp.fp_id] = true
                push!(new_served_requests, curr_fp.req)
            end
            
            if can_serve_if_increment_cars
                station_id = findfirst(x -> x == curr_fp.origin_station, get_potential_locations())
                sol.initial_cars_number[station_id] += 1
                push!(stations_to_increment_cars, station_id)
            end

        end
    end

    global online_selected_paths = sol.selected_paths
    #@info "*********** end of serve_requests_after_opening_station ***************"
    E_carsharing_sim(sol), stations_to_increment_cars
end

"""
    can_use_trip(sc_id::Int64, trip::DataFrameRow, sol::Solution)

    check if we can serve the trip in the solution
    inputs:
        @sc_id: the id of the scenario
        @trip: the trip to be served
        @sol: the solution

    outputs:
        true if we can serve the trip, false otherwise
"""
function can_use_trip(sc_id::Int64, trip::DataFrameRow, sol::Solution)
    #list of global variables
    stations = scenario_list[sc_id].stations
    
    #list of served trips 
    sol_feasile_trips = scenario_list[sc_id].feasible_paths[sol.selected_paths[sc_id], :]

    #general information about the new trip
    trip_origin_station = trip.origin_station
    trip_destination_station = trip.destination_station

    initial_car_origin_station = sol.initial_cars_number[findfirst(get_potential_locations() .== trip_origin_station)]
    origin_station_capacity = stations[findfirst(get_potential_locations() .== trip_origin_station)].max_number_of_charging_points
    initial_car_destination_station = sol.initial_cars_number[findfirst(get_potential_locations() .== trip_destination_station)]
    destination_station_capacity = stations[findfirst(get_potential_locations() .== trip_destination_station)].max_number_of_charging_points

    #get all the trips involving the origin station and add the new trips and sort them
    origin_station_trips = filter(row -> row.destination_station == trip_origin_station ||
            row.origin_station == trip_origin_station, sol_feasile_trips)
    #add a column to the dataframe to keep track of the taking or parking time to make the sorting easy
    origin_station_trips.taking_or_parking_time = [row.origin_station == trip_origin_station ?
                                                   row.start_driving_time : row.arriving_time for row in eachrow(origin_station_trips)]

    #we add the trip and see if always the cars at this station are >= 0 
    #(if negative we stop and we can not serve this request)
    trip_as_df = DataFrame(trip)
    trip_as_df.taking_or_parking_time = [trip.start_driving_time]
    origin_station_trips = vcat(origin_station_trips, trip_as_df)
    sort!(origin_station_trips, :taking_or_parking_time)

    cars_at_station = initial_car_origin_station
    can_serve = true
    car_arriving_time = zeros(initial_car_origin_station) #keep track when cars are arrived to the station

    for row in eachrow(origin_station_trips)
        
        if row.origin_station == trip_origin_station
            cars_at_station -= 1
            if cars_at_station < 0
                can_serve = false
                break
            end
            start_charging = pop!(car_arriving_time)
            if start_charging == row.start_driving_time
                can_serve = false
                break
            end
        else
            cars_at_station += 1
            pushfirst!(car_arriving_time, row.arriving_time)
        end
    end
    can_serve_if_we_increment_cars = false
    if !can_serve
        #=  # we can not handle the current request so we try to increment the number of cars in the station
        #vars to investigate if we add a car to origin station we can serve the request
        can_serve = true
        cars_at_station = initial_car_origin_station + 1
        car_arriving_time = zeros(initial_car_origin_station + 1)
        #check the capacity
        if cars_at_station > origin_station_capacity
            #we cannot serve in both cases
            return (false, false) 
        end
        for row in eachrow(origin_station_trips)
        
            if row.origin_station == trip_origin_station
                cars_at_station -= 1
                if cars_at_station < 0
                    can_serve = false
                    break
                end
                start_charging = pop!(car_arriving_time)
                if start_charging == row.start_driving_time
                    can_serve = false
                    break
                end
            else
                cars_at_station += 1
                pushfirst!(car_arriving_time, row.arriving_time)
                if cars_at_station > origin_station_capacity
                    can_serve = false
                    break
                end
            end
        end
        if !can_serve
            return (false, false)
        else
            can_serve_if_we_increment_cars = true
            
        end =#
        return (false, false)
    end

    destination_station_trips = filter(row -> row.destination_station == trip_destination_station ||
            row.origin_station == trip_destination_station,
        sol_feasile_trips)
    #add a column to the dataframe to keep track of the taking or parking time to make the sorting easy
    destination_station_trips.taking_or_parking_time = [row.origin_station == trip_destination_station ?
                                                        row.start_driving_time : row.arriving_time for row in eachrow(destination_station_trips)]

    #we add the trip and see if always the free parking spots this station are >= 0 
    #(if negative we stop and we can not serve this request)
    trip_as_df.taking_or_parking_time = [trip.arriving_time]
    destination_station_trips = vcat(destination_station_trips, trip_as_df)
    sort!(destination_station_trips, :taking_or_parking_time)

    free_spots = destination_station_capacity - initial_car_destination_station

    for row in eachrow(destination_station_trips)
        row.origin_station == trip_destination_station ? free_spots += 1 : free_spots -= 1
        if free_spots < 0
            can_serve = false
            break
        end
    end

    return can_serve, can_serve_if_we_increment_cars
end


"""
    clean_up_cars_number!(sol::Solution)
    
    This function clean up the solution to serve the requests using only 
    minimal number of cars it alters the initial number of cars in the solution
    inputs:
        @sol: the solution
    outputs:
        the new objective function
"""
function clean_up_cars_number!(sol::Solution)
    global used_cars
    
    # the way that the cars are created are basically when we instantiate the stations
    E_carsharing_sim(sol)
    new_initial_number_of_cars = Array{Int64}(undef, length(sol.initial_cars_number))
    car_id = 1

    for station_id in eachindex(sol.initial_cars_number)
        old_init_car_numbers = car_id:car_id+sol.initial_cars_number[station_id]-1
        car_id += sol.initial_cars_number[station_id]

        new_initial_number_of_cars[station_id]  = count([x in used_cars for x in old_init_car_numbers])
    end

    sol.initial_cars_number = new_initial_number_of_cars

    #return the new objective function
    E_carsharing_sim(sol)

end

"""
    clean Up the selecte paths after closing (a) station(s).
    inputs:
        @sol: the solution
    outputs:
        the new objective function
"""
function clean_up_selected_paths_MT_old!(sol::Solution)
    #println("************** cleaning up selected paths **************")
    global failed
    global current_scenario_id
    global trips_to_unselect
    
    sol_selected_path_lock = ReentrantLock() 
    #first the trivial case: unselect all the trips that contains a closed station
    @threads for sc in eachindex(sol.selected_paths)
        scenario = scenario_list[sc]
        for fp in eachindex(sol.selected_paths[sc])
            
            if sol.selected_paths[sc][fp]
                #get the trip information
                origin_station_id = findfirst(get_potential_locations() .== scenario.feasible_paths.origin_station[fp])
                destination_station_id = findfirst(get_potential_locations() .== scenario.feasible_paths.destination_station[fp])

                if !sol.open_stations_state[origin_station_id] || !sol.open_stations_state[destination_station_id]
                    Threads.lock(sol_selected_path_lock)
                    try
                        sol.selected_paths[sc][fp] = false
                    finally
                        Threads.unlock(sol_selected_path_lock)
                    end
                end
            end
        end
    end
    
    # second step is to run the simulation and see if everything is okay.
    # if there still trips that cause infeasible solution we unselect them
    failed = Threads.Atomic{Bool}(true)
    fit_value = 10^16
    set_online_mode(false)
    while failed[]
        findall(x->length(x)>0, trips_to_unselect) 
        trips_to_unselect = [Int64[] for _ in eachindex(scenario_list)]
        fit_value = E_carsharing_sim(sol)
        
        for i in eachindex(scenario_list)
            sol.selected_paths[i][trips_to_unselect[i]] .= false
        end
    end


    fit_value, sol.selected_paths
end

"""
    clean_up_selected_paths!(sol::Solution)
    clean Up the selecte paths after closing (a) station(s).
    inputs:
        @sol: the solution
    outputs:
        the new objective function
        the new list of the sletced paths
"""
function clean_up_selected_paths!(sol::Solution)
    
    #println("************** cleaning up selected paths **************")
    global failed
    global current_scenario_id
    global trips_to_unselect
        
    sol_selected_path_lock = ReentrantLock() 
    
    @threads for sc in eachindex(sol.selected_paths)
        scenario = scenario_list[sc]
        stations_to_recheck = Int64[]
        for fp in eachindex(sol.selected_paths[sc])
            #first Step: unselect all the trips that contains a closed station and mark the closed stations
            if sol.selected_paths[sc][fp]
                #get the trip information
                origin_station_id = findfirst(get_potential_locations() .== scenario.feasible_paths.origin_station[fp])
                destination_station_id = findfirst(get_potential_locations() .== scenario.feasible_paths.destination_station[fp])

                #check if one of its stations is closed
                if !sol.open_stations_state[origin_station_id] || !sol.open_stations_state[destination_station_id]
                    #unselect the trip
                    Threads.lock(sol_selected_path_lock)
                    try
                        sol.selected_paths[sc][fp] = false
                    finally
                        Threads.unlock(sol_selected_path_lock)
                    end
                    #keep track of the opened station to recheck them
                    if sol.open_stations_state[origin_station_id]
                        push!(stations_to_recheck, origin_station_id)
                    elseif sol.open_stations_state[destination_station_id]
                        push!(stations_to_recheck, destination_station_id)
                    end
                end
            end

        end
        #save_sol(sol, "sol1.jls")
        #second step: recheck the opened stations
        unique!(stations_to_recheck)
        for station_id in stations_to_recheck
            additional_station_to_recheck = recheck_station!(sol, scenario, station_id, sol_selected_path_lock)
            length(additional_station_to_recheck) > 0 && push!(stations_to_recheck, additional_station_to_recheck...)
        end
    end
    #re-evaluate the solution
    fit_value = E_carsharing_sim(sol)
    fit_value, sol.selected_paths
end

"""
    recheck_station!(sol::Solution, scenario::Scenario, station_id::Int64, sol_lock::ReentrantLock)

    when we close a station, some other stations are serving requests by cars came from the closed station
    or after freeing a parking place after sending a car to the closed station. This function check if a station
    still serve the requests after closing a station. If not, it unselect the trip.

    inputs:
        @sol: the solution
        @scenario: the scenario
        @station_id: the id of the station to recheck
        @sol_lock: the lock to protect the solution

    outputs:
        the list of stations to recheck again
"""
function recheck_station!(sol::Solution, scenario::Scenario, station_id::Int64, sol_lock::ReentrantLock)
    # sol, scenario, station_id = load_sol("sol1.jls"), scenario_list[4], 79;
    #list of served trips 
    selected_trips = scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :]
    
    initial_car_nbr = sol.initial_cars_number[station_id]
    station_capacity = scenario.stations[station_id].max_number_of_charging_points
    station_node_id = get_potential_locations()[station_id]
    
    #get all the trips involving the origin station and add the new trips and sort them
    station_trips = filter(row -> row.destination_station == station_node_id ||
            row.origin_station == station_node_id, selected_trips)
    
    #add a column to the dataframe to keep track of the taking or parking time to make the sorting easy
    station_trips.taking_or_parking_time = [row.origin_station == station_node_id ?
                                                   row.start_driving_time : row.arriving_time for row in eachrow(station_trips)]

    #sort the trips according to their taking or parking time
    sort!(station_trips, :taking_or_parking_time)

    cars_at_station = initial_car_nbr
    car_arriving_time = zeros(initial_car_nbr) #keep track when cars are arrived to the station

    add_stations_to_recheck = Int64[]

    for row in eachrow(station_trips)
        #row = station_trips[5, :]
        if row.origin_station == station_node_id #we are taking a car from the station
            cars_at_station -= 1
            if cars_at_station < 0
                #unselect the trip
                Threads.lock(sol_lock)
                try
                    sol.selected_paths[scenario.scenario_id][row.fp_id] = false
                    cars_at_station += 1
                    push!(add_stations_to_recheck, findfirst(get_potential_locations() .== row.destination_station))
                finally
                    Threads.unlock(sol_lock)
                    continue
                end
            end
            
            start_charging = pop!(car_arriving_time)
            if start_charging == row.start_driving_time
                Threads.lock(sol_lock)
                try
                    sol.selected_paths[scenario.scenario_id][row.fp_id] = false
                    cars_at_station += 1
                    push!(add_stations_to_recheck, findfirst(get_potential_locations() .== row.destination_station))
                    push!(car_arriving_time, start_charging)
                finally
                    Threads.unlock(sol_lock)
                end
            end
        else
            cars_at_station += 1
            if cars_at_station > station_capacity
                Threads.lock(sol_lock)
                try
                    sol.selected_paths[scenario.scenario_id][row.fp_id] = false
                    cars_at_station -= 1
                    push!(add_stations_to_recheck, findfirst(get_potential_locations() .== row.origin_station))
                finally
                    Threads.unlock(sol_lock)
                end
            end
            pushfirst!(car_arriving_time, row.arriving_time)
        end

    end

    return add_stations_to_recheck
end

"""
    serve_requests(sol::Solution, requests_list::Array{Int64, 1})
    given a list of requests try to serve them using theopened stations in the solution

    inputs:
        @sol: the solution
        @requests_list: the list of requests to serve for each scenario
    outputs:
        the new objective function
"""
function serve_requests!(sol::Solution, requests_list::Vector{Vector{Int64}})
    #@info "************** serving requests **************"
    
    opend_stations_node_ids = get_potential_locations()[sol.open_stations_state]
    stations_to_increment_cars = Int64[]

    @threads for sc_id in eachindex(requests_list)
        curr_scenario = scenario_list[sc_id]
        #loop over the requests
        for req_id in requests_list[sc_id]
            #get the feasible trips for the current request
            curr_req_feasible_trips_id = findall(x -> curr_scenario.feasible_paths.req[x] == req_id &&
                                                    curr_scenario.feasible_paths.origin_station[x] in opend_stations_node_ids &&
                                                    curr_scenario.feasible_paths.destination_station[x] in opend_stations_node_ids,
                1:nrow(curr_scenario.feasible_paths))

            for trip_id in curr_req_feasible_trips_id
                #check if we can serve the request without trying to serve it
                curr_fp = curr_scenario.feasible_paths[trip_id, :]
                can_serve, can_serve_if_increment_cars = can_use_trip(sc_id, curr_fp, sol)
                if can_serve
                    sol.selected_paths[sc_id][curr_fp.fp_id] = true
                end

                if can_serve_if_increment_cars
                    station_id = findfirst(x -> x == curr_fp.origin_station, get_potential_locations())
                    sol.initial_cars_number[station_id] += 1
                    push!(stations_to_increment_cars, station_id)
                end
            end
        end
    end
    #return the new objective function
    E_carsharing_sim(sol), sol.selected_paths, stations_to_increment_cars
end

"""
    get_trips_station(station_id::Int64, feasible_trips::DataFrame)

    given the list of feasible trips this function return the trips where station_id intervened
    inputs:
        @station_id: the id of the station
        @feasible_trips: the list of feasible trips

    outputs:
        the trips where station_id intervened
"""
function get_trips_station(station_id::Int64, feasible_trips::DataFrame)
    station_node_id = get_potential_locations()[station_id]
    trips = filter(x -> x.origin_station == station_node_id || x.destination_station == station_node_id, feasible_trips)
    trips
end

function serve_requests_old!(sol::Solution, requests_list::Vector{Vector{Int64}})
    #@info "************** serving requests **************"
   
    opend_stations_node_ids = get_potential_locations()[sol.open_stations_state]
    stations_to_increment_cars = Int64[]

    for sc_id in eachindex(requests_list)
        curr_scenario = scenario_list[sc_id]
        #loop over the requests
        for req_id in requests_list[sc_id]
            #get the feasible trips for the current request
            curr_req_feasible_trips_id = findall(x -> curr_scenario.feasible_paths.req[x] == req_id &&
                                                    scenario_list[sc_id].feasible_paths.origin_station[x] in opend_stations_node_ids &&
                                                    scenario_list[sc_id].feasible_paths.destination_station[x] in opend_stations_node_ids,
                1:nrow(scenario_list[sc_id].feasible_paths))

            for trip_id in curr_req_feasible_trips_id
                #check if we can serve the request without trying to serve it
                
                curr_fp = scenario_list[sc_id].feasible_paths[trip_id, :]
                can_serve, can_serve_if_increment_cars = can_use_trip(sc_id, curr_fp, sol)
                if can_serve
                    sol.selected_paths[sc_id][curr_fp.fp_id] = true
                end

                if can_serve_if_increment_cars
                    station_id = findfirst(x -> x == curr_fp.origin_station, get_potential_locations())
                    sol.initial_cars_number[station_id] += 1
                    push!(stations_to_increment_cars, station_id)
                end
            end
        end
    end
    #return the new objective function
    E_carsharing_sim(sol), sol.selected_paths, stations_to_increment_cars
end