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