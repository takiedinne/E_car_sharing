#= #################################################################### =#
#= ############# Global variables for all scenarios#################### =#

#create the graph
# manhaten_city_graph = create_graph_from_XML(Manhatten_network_details_file, save_file = "Data/manhatten_graph.mg")
# manhaten_city_graph = create_graph_from_XML("Tests/test_graph.xml", save_file = "Data/test_graph.mg")

#to load the graph use the following instruction ===> it is better to load the graph in terme of computational performance
global manhaten_city_graph = loadgraph(Manhatten_network_Metagraph_file, MGFormat())

global non_directed_manhaten_city_graph = MetaGraph(manhaten_city_graph) # for walking purposes we don't take into account the directed edges 

#read all requests details
global all_request_df = CSV.read(all_request_details_path, DataFrame)

#initialization
global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph


#= #################################################################### =#
#= ############# Global variables for each scenario #################### =#
#result of the pre-processing step
global all_feasible_paths_scenario = Array{DataFrame,1}()
global all_feasible_paths 

global failed = false #for the offline mode it is needed to stop the simulation

global revenues = 0 #

global stations = Array{Station, 1}() # list of stations

global scenario # the requests lists

@resumable function request_arrival_process_offline_mode(env::Environment, scenario::DataFrame, sol::Solution)
    #browse all the requests (the requests are sorted according to their arrival time)
    for req in eachrow(scenario)
        
        if failed
            break
        end

        # waiting until a new request is arrived 
        sleeping_time = req.ST - now(env)
        sleeping_time > 0 && @yield timeout(env, sleeping_time)
        print_simulation && println("Customer [", req.reqId, "]: request arrive at ", now(env))

        # in the offline mode we can check if the request is decided to be served according to the decision variables
        if isempty(req.feasible_paths)
            # the resuest is refused by the decision variables
            print_simulation && println("Customer [", req.reqId, "]: the requests is rejected according to the decision variables")
            continue
        else
            pickup_station_id = findfirst(sol.open_stations_ids .== req.feasible_paths.origin_station[1])
            drop_off_station_id = findfirst(sol.open_stations_ids .== req.feasible_paths.destination_station[1])
            selected_car_id, parking_place_id  =  -1, -1
            # check the availabilty of paths to serve the request
            if pickup_station_id != -1 && drop_off_station_id != -1
                # for the online serving we check all the variables otherwise we nee only to check the stations
                @process perform_the_trip_process(env, req, sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
            else
                # we couldn't serve the request there is no feasible path
                print_simulation && printstyled(stdout, "Customer [", req.reqId, "]: We can not serve the request there is no feasible path\n", color = :yellow)
                global failed = true
                break
            end
        end
    end
end

@resumable function request_arrival_process_online_mode(env::Environment, scenario::DataFrame, sol::Solution)
    #browse all the requests (the requests are sorted according to their arrival time)
    for req in eachrow(scenario)
       
        # in the offline mode we check if we failed to serve a request so we stop the simulation
       
        # waiting until a new request is arrived 
        sleeping_time = req.ST - now(env)
        sleeping_time > 0 && @yield timeout(env, sleeping_time)
        print_simulation && println("Customer [", req.reqId, "]: request arrive at ", now(env))

        # get the trip information (pickup station, drop off station, the selected car id, parking place)
        # if the request can not be served, this function will return -1 in one of the information variables
        
        # in the offline mode we can check if the request is decided to be served according to the decision variables
        if !online_request_serving
            if isempty(req.feasible_paths)
                # the resuest is refused by the decision variables
                print_simulation && println("Customer [", req.reqId, "]: the requests is rejected according to the decision variables")
                continue
            else
                pickup_station_id = findfirst(sol.open_stations_ids .== req.feasible_paths.origin_station[1])
                drop_off_station_id = findfirst(sol.open_stations_ids .== req.feasible_paths.destination_station[1])
                selected_car_id, parking_place_id  =  -1, -1
            end
        else
            pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id = get_trip_info_for_request(req, sol, now(env)) # one row Dataframe
        end
        
        # check the availabilty of paths to serve the request
        if pickup_station_id != -1 && drop_off_station_id != -1 && ( online_request_serving && selected_car_id != -1 && parking_place_id != -1 || !online_request_serving)
            # for the online serving we check all the variables otherwise we nee only to check the stations
            @process perform_the_trip_process(env, req, sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
        else
            #we couldn't serve the request there is no feasible path
            # the cause message
            if print_simulation
                if (pickup_station_id == -1 || drop_off_station_id == -1)
                    #there is no feasible path for the request
                    cause_message = "there is no feasible path"
                elseif print_simulation && selected_car_id == -1
                    cause_message = "there is no available car"
                elseif print_simulation && parking_place_id == -1
                    cause_message = "there is no parking place"
                end
            end
                 
            print_simulation && printstyled(stdout, "Customer [", req.reqId, "]: We can not serve the request there is no feasible path (",cause_message,")\n", color = :yellow)
        end
    end

end

@resumable function perform_the_trip_process_offline_mode(env::Environment, req::DataFrameRow, sol::Solution, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    
    #some information
    walking_duration = get_walking_time(req.ON, sol.open_stations_ids[pickup_station_id])
    walking_duration = work_with_time_slot ? ceil(Integer, walking_duration / time_slot_length) : walking_duration

    start_walking_time = now(env)

    driving_duration = get_trip_duration(sol.open_stations_ids[pickup_station_id], sol.open_stations_ids[drop_off_station_id])
    driving_duration = work_with_time_slot ? ceil(Integer, driving_duration / time_slot_length) : driving_duration
    
    print_simulation && println("Customer [", req.reqId, "]: the request is accepted") #
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    walking_duration > 0 && @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station ", sol.open_stations_ids[pickup_station_id])

    # check and select a car to performe the trip
    refrech_battery_levels(pickup_station_id, now(env)) #refrech the battery leveles 
    available_car_df = filter(row -> row.status == CAR_PARKED &&
            row.last_battery_level >= get_battery_level_needed((pickup_station_id, drop_off_station_id)),
        stations[pickup_station_id].cars)
    if isempty(available_car_df)
        #double check if there is another request (at this moment) that can help to serve the current request
        @yield timeout(env, 0.0, priority = -1) # put this process at the end of the heap for the current time
        available_car_df = filter(row -> row.status == CAR_PARKED &&
                row.last_battery_level >= get_battery_level_needed((pickup_station_id, drop_off_station_id)),
            stations[pickup_station_id].cars)
    end

    if !isempty(available_car_df)
        # simply we take the first found car
        selected_car_id = available_car_df.car_id[1]
        stations[pickup_station_id].cars.status[findfirst(x -> x == selected_car_id, stations[pickup_station_id].cars.car_id)] = CAR_RESERVED
        # add the car in the drop_off station 
        car = DataFrame(available_car_df[1, :]) # copy
        car.status[1] = CAR_ON_WAY
        car.start_charging_time[1] = NaN
        car.expected_arrival_time[1] = now(env) + driving_duration
        car.last_battery_level[1] -= get_trip_battery_consumption(sol.open_stations_ids[pickup_station_id],
            sol.open_stations_ids[drop_off_station_id],
            car.car_type[1])
        append!(stations[drop_off_station_id].cars, car)
    else
        print_simulation && printstyled(stdout, "Customer [$(req.reqId)]: (offline mode) there is no available car at the station \n"; color=:light_red)
        global failed = true
        return
    end
   
    print_simulation && println("Customer [", req.reqId, "]: start ridding the car number $selected_car_id at $(now(env))")

    # start event with a higher priority to always favour the start driving event over the drop event
    start_ev = Event(env)
    @callback start_driving(start_ev, pickup_station_id, selected_car_id) # this function put failed true if there is not a car with the consumption constraints
    succeed(start_ev, priority = 10)

    #simulate the riding time
    driving_duration > 0 && @yield timeout(env, driving_duration)
    
    # for offline mode we have to see if there is a parking place
    #get the available places
    parking_and_cars_df = leftjoin(stations[drop_off_station_id].parking_places, stations[drop_off_station_id].cars, on=:cars => :car_id, makeunique=true)
    available_places_df = filter(row -> row.status == P_FREE, stations[drop_off_station_id].parking_places)
    if isempty(available_places_df)
        #double check if there is another request that can help to serve the current request
        @yield timeout(env, 0.0, priority = -1) # put this process at the end of the heap for the current time
        available_places_df = filter(row -> row.status == P_FREE, stations[drop_off_station_id].parking_places)
    end
    if !isempty(available_places_df)
        # simply we select the first free place
        parking_place_id = available_places_df.p_id[1]
        # reserve the place 
        stations[drop_off_station_id].parking_places.status[parking_place_id] = P_RESERVED
    else
        print_simulation && printstyled(stdout, "Customer [$(req.reqId)]: (offline mode) there is no free parking place at station that has id  $drop_off_station_id \n"; color=:light_red)
        global failed = true
        return
    end
    
    # drop the car
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", sol.open_stations_ids[drop_off_station_id], " at ", now(env))
    #@process drop_car(env, drop_off_station_id, selected_car_id, parking_place_id) # this function put failed variable to true if we couldn't find a parking place 
    
    drop_ev = Event(env)
    @callback drop_car(drop_ev, drop_off_station_id, selected_car_id, parking_place_id, req.reqId)
    succeed(drop_ev, priority = 10)

    # make the payment
    global revenues += req.Rev
    
end

@resumable function perform_the_trip_process(env::Environment, req::DataFrameRow, sol::Solution, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    
    #some information
    walking_duration = get_walking_time(req.ON, sol.open_stations_ids[pickup_station_id])
    walking_duration = work_with_time_slot ? ceil(Integer, walking_duration / time_slot_length) : walking_duration

    start_walking_time = now(env)

    driving_duration = get_trip_duration(sol.open_stations_ids[pickup_station_id], sol.open_stations_ids[drop_off_station_id])
    driving_duration = work_with_time_slot ? ceil(Integer, driving_duration / time_slot_length) : driving_duration

    #book the trip (the car + parking palce ,etc) only for the 
    if online_request_serving
        book_trip(sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id, start_walking_time + walking_duration, start_walking_time + walking_duration + driving_duration)
    end
    print_simulation && println("Customer [", req.reqId, "]: the request is accepted") #
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    walking_duration > 0 && @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station ", sol.open_stations_ids[pickup_station_id])

    #for the offline mode we have to check if there is a car 
    if !online_request_serving
        # check and select a car to performe the trip
        refrech_battery_levels(pickup_station_id, now(env)) #refrech the battery leveles 
        available_car_df = filter(row -> row.status == CAR_PARKED &&
                row.last_battery_level >= get_battery_level_needed((pickup_station_id, drop_off_station_id)),
            stations[pickup_station_id].cars)
        if isempty(available_car_df)
            #double check if there is another request that can help to serve the current request
            @yield timeout(env, 0.0, priority = -1) # put this process at the end of the heap for the current time
            available_car_df = filter(row -> row.status == CAR_PARKED &&
                    row.last_battery_level >= get_battery_level_needed((pickup_station_id, drop_off_station_id)),
                stations[pickup_station_id].cars)
        end
        if !isempty(available_car_df)
            # simply we take the first found car
            selected_car_id = available_car_df.car_id[1]
            stations[pickup_station_id].cars.status[findfirst(x -> x == selected_car_id, stations[pickup_station_id].cars.car_id)] = CAR_RESERVED
            # add the car in the drop_off station 
            car = DataFrame(available_car_df[1, :]) # copy
            car.status[1] = CAR_ON_WAY
            car.start_charging_time[1] = NaN
            car.expected_arrival_time[1] = now(env) + driving_duration
            car.last_battery_level[1] -= get_trip_battery_consumption(sol.open_stations_ids[pickup_station_id],
                sol.open_stations_ids[drop_off_station_id],
                car.car_type[1])
            append!(stations[drop_off_station_id].cars, car)
        else
            print_simulation && printstyled(stdout, "Customer [$(req.reqId)]: (offline mode) there is no available car at the station \n"; color=:light_red)
            global failed = true
        end
    end

    if !failed # we could start the trip 
        print_simulation && println("Customer [", req.reqId, "]: start ridding the car number $selected_car_id at $(now(env))")

        # start event with a higher priority to always favour the start driving event over the drop event
       
        start_ev = Event(env)
        @callback start_driving(start_ev, pickup_station_id, selected_car_id) # this function put failed true if there is not a car with the consumption constraints
        succeed(start_ev, priority = 10)

        #simulate the riding time
        driving_duration > 0 && @yield timeout(env, driving_duration)
        # for offline mode we have to see if there is a parking place
        if !online_request_serving
            
            #get the available places
            parking_and_cars_df = leftjoin(stations[drop_off_station_id].parking_places, stations[drop_off_station_id].cars, on=:cars => :car_id, makeunique=true)

            #get the available places
            available_places_df = filter(row -> row.status == P_FREE, stations[drop_off_station_id].parking_places)
            if isempty(available_places_df)
                #double check if there is another request that can help to serve the current request
                @yield timeout(env, 0.0, priority = -1) # put this process at the end of the heap for the current time
                if req.reqId in [697, 714]
                    println("find parking place second tentative $(req.reqId)")
                end
                available_places_df = filter(row -> row.status == P_FREE, stations[drop_off_station_id].parking_places)
            end
            if !isempty(available_places_df)
                # simply we select the first free place
                parking_place_id = available_places_df.p_id[1]
                # reserve the place 
                stations[drop_off_station_id].parking_places.status[parking_place_id] = P_RESERVED
            else
                print_simulation && printstyled(stdout, "Customer [$(req.reqId)]: (offline mode) there is no free parking place at station that has id  $drop_off_station_id \n"; color=:light_red)
                global failed = true
            end
        end
        if !failed # there is no place to drop the car
            # drop the car
            print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", sol.open_stations_ids[drop_off_station_id], " at ", now(env))
            #@process drop_car(env, drop_off_station_id, selected_car_id, parking_place_id) # this function put failed variable to true if we couldn't find a parking place 
            
            drop_ev = Event(env)
            @callback drop_car(drop_ev, drop_off_station_id, selected_car_id, parking_place_id, req.reqId)
            succeed(drop_ev, priority = 10)

            # make the payment
            global revenues += req.Rev
        end
    end
end

function start_driving(env::AbstractEvent, station_id, car_id)
   
    car_indx = findfirst(stations[station_id].cars.car_id .== car_id)
    parking_place_indx = findfirst(stations[station_id].parking_places.cars .== car_id)

    if car_indx === nothing || parking_place_indx === nothing
        printstyled(stdout, "Error: the car or the parking place are not found \n"; color=:light_red)
        global failed = true
        return
    end

    # free the parking place
    if stations[station_id].parking_places.pending_reservation[parking_place_indx] > 0
        stations[station_id].parking_places.pending_reservation[parking_place_indx] -= 1
        stations[station_id].parking_places.status[parking_place_indx] = P_RESERVED
    else
        stations[station_id].parking_places.status[parking_place_indx] = P_FREE
    end
    stations[station_id].parking_places.cars[parking_place_indx] = -1 # there is no cars

    # delete the car from the station 
    filter!(row -> (row.car_id != car_id), stations[station_id].cars)
    #= if station_id == 30 
        println("pickup at $(now(env)) nbr cars: $(nrow(stations[station_id].cars))")
        #println(stations[station_id].cars)
    end =#
end

function drop_car(env::AbstractEvent, drop_off_station_id, car_id, parking_place_id, req)
    current_time = now(env.bev.env)
    # check if the parking place is free (basically it will be free but just we check if there is a problem)
    if stations[drop_off_station_id].parking_places.status[parking_place_id] == P_OCCUPIED
        printstyled(stdout, "Error: the parking place is occupied there is problem in the simulatiuon logic station $(stations[drop_off_station_id]) \n"; color=:light_red)
        global failed = true
        return
    end

    # occupy the parking place
    stations[drop_off_station_id].parking_places.status[parking_place_id] = P_OCCUPIED
    stations[drop_off_station_id].parking_places.cars[parking_place_id] = car_id

    # get the car index inside the data frame
    car_indx = findfirst(x -> x == car_id, stations[drop_off_station_id].cars.car_id)

    # change the status of the car
    stations[drop_off_station_id].cars.status[car_indx] = CAR_PARKED
    stations[drop_off_station_id].cars.start_charging_time[car_indx] = current_time

    #set the reseravtion and expected time 
    stations[drop_off_station_id].cars.start_reservation_time[car_indx] = NaN
    stations[drop_off_station_id].cars.expected_arrival_time[car_indx] = NaN

end

function initialize_sim(sol::Solution, scenario::DataFrame)
    global failed = false
    global revenues = 0
    
    global stations = Array{Station, 1}()
    #set the stations
    car_id = 1
    for i in 1:length(sol.open_stations_ids)
        initial_cars_number = sol.initial_cars_number[i]
        push!(stations, Station(sol.open_stations_ids[i], initial_cars_number, car_id))
        car_id += initial_cars_number
    end

    selected_paths = !online_request_serving ? all_feasible_paths[sol.selected_paths, :] : all_feasible_paths

    # put each feasible path beside its request
    scenario.feasible_paths = [filter( x-> x.req == i, selected_paths) for i in scenario.reqId]


end

function E_carsharing_sim(sol::Solution, scenario_id::Integer)
    # check the feasibilty of the solution
    scenario =  scenario_list[scenario_id]
    global all_feasible_paths = all_feasible_paths_scenario[scenario_id]

    if is_feasible_solution(sol)
        
        initialize_sim(sol, scenario)
        
        sim = Simulation()
        if online_request_serving
            @process request_arrival_process_online_mode(sim, scenario, sol)
        else
            @process request_arrival_process_offline_mode(sim, scenario, sol)
        end
        run(sim)
        if failed
            print_simulation && printstyled(stdout, "The simulation is stopped because a request couldn't be served\n", color=:light_red)
            return penality
        else
            # count the objective function
            print_simulation && println("counting the objective function")
            total_cars_cost = 0
            for i in 1:length(stations)
                if nrow(stations[i].cars) > 0
                    total_cars_cost += sum([vehicle_specific_values[stations[i].cars.car_type[j]][:car_cost] for j in 1:nrow(stations[i].cars)])
                end
            end
        
            #= total_station_cost = sum([station.charging_station_base_cost + 
                                        station.max_number_of_charging_points * station.charging_point_cost_fast
                                        for station in stations]) =#
            total_station_cost = sum([station.charging_station_base_cost for station in stations])
        
            return revenues - (total_cars_cost + total_station_cost) / cost_factor
        end
    else
        print_simulation && println("The solution is not feasible")
        return penality
    end
end
