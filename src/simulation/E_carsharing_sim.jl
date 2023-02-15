#= #################################################################### =#
#= ############# Global variables for all scenarios#################### =#

#create the graph
# manhaten_city_graph = create_graph_from_XML(Manhatten_network_details_file, save_file = "Data/manhatten_graph.mg")
# manhaten_city_graph = create_graph_from_XML(Manhatten_network_details_file, save_file = "Data/test_graph.mg")

#to load the graph use the following instruction ===> it is better to load the graph in terme of computational performance
global manhaten_city_graph = loadgraph(Manhatten_network_Metagraph_file, MGFormat())
global non_directed_manhaten_city_graph = MetaGraph(manhaten_city_graph) # for walking purposes we don't take into account the directed edges 
global potential_locations = [] # the index of all potential stations 
#read all requests details
global all_request_df = CSV.read(all_request_details_path, DataFrame)
global scenario_list  # contain all scenario instances as dataframe

#initialization
global shortest_car_paths = Dict{Int64,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
global shortest_walking_paths = Dict{Int64,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph


#= #################################################################### =#
#= ############# Global variables for each scenario #################### =#

global failed = false #for the offline mode it is needed to stop the simulation
global revenues = 0 #the revenue of serving customers

global stations = Array{Station,1}() # list of open_stations_state

global online_selected_paths = Array{Array{Bool,1}, 1}() #useful for the online mode to see which paths is chosen
global number_of_served_requests = 0 # a couter for the numlber of requests that the simulation served
################################ Online mode ####################################
@resumable function request_arrival_process_online_mode(env::Environment, scenario::Scenario, sol::Solution)
    #browse all the requests (the requests are sorted according to their arrival time)
    for req in eachrow(scenario.request_list)
        if failed
            break
        end

        # waiting until a new request is arrived 
        sleeping_time = req.ST - now(env)
        sleeping_time > 0 && @yield timeout(env, sleeping_time)
        print_simulation && println("Customer [", req.reqId, "]: request arrive at ", now(env))

        # get the trip information (pickup station, drop off station, the selected car id, parking place)
        # if the request can not be served, this function will return -1 in one of the information variables
        pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id = get_trip_info_for_request(req, sol, scenario, now(env)) # one row Dataframe
        
        # check the availabilty of paths to serve the request
        if pickup_station_id != -1 && drop_off_station_id != -1 && selected_car_id != -1 && parking_place_id != -1 
            #book the trip (the car + parking palce ,etc) only for the 
            walking_duration = get_walking_time(req.ON, potential_locations[pickup_station_id])
            start_walking_time = now(env)
            driving_duration = get_trip_duration(potential_locations[pickup_station_id], potential_locations[drop_off_station_id])

            book_trip(pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id, start_walking_time + walking_duration, start_walking_time + walking_duration + driving_duration)

            # for the online serving we check all the variables otherwise we have only to check the stations
            @process perform_the_trip_process_online_mode(env, req, sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
        else
            #we couldn't serve the request there is no feasible path
            # the cause message
            if print_simulation
                if (pickup_station_id == -1 || drop_off_station_id == -1)
                    #there is no feasible path for the request
                    cause_message = "there is no feasible path"
                elseif selected_car_id == -1
                    cause_message = "there is no available car"
                elseif parking_place_id == -1
                    cause_message = "there is no parking place"
                end
            end

            print_simulation && printstyled(stdout, "Customer [", req.reqId, "]: We can not serve the request there is no feasible path (", cause_message, ")\n", color=:yellow)
        end
    end

end

@resumable function perform_the_trip_process_online_mode(env::Environment, req::DataFrameRow, sol::Solution, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    #get the different duration for the trip.
    walking_duration = get_walking_time(req.ON, potential_locations[pickup_station_id])
    start_walking_time = now(env)
    driving_duration = get_trip_duration(potential_locations[pickup_station_id], potential_locations[drop_off_station_id])

    print_simulation && println("Customer [", req.reqId, "]: the request is accepted") #
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    walking_duration > 0 && @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station N°", potential_locations[pickup_station_id], " at ", now(env))

    print_simulation && println("Customer [", req.reqId, "]: start ridding the car number $selected_car_id at $(now(env))")

    # take the car and free the parking space
    start_driving(pickup_station_id, selected_car_id) # this function put failed true if there is not a car with the consumption constraints
    
    failed && print_simulation && printstyled(stdout, "Customer [", req.reqId, "]: Problem while starting the drive\n", color=:light_red)
    

    #simulate the driving time
    driving_duration > 0 && @yield timeout(env, driving_duration)
    
    # drop the car
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", potential_locations[drop_off_station_id], " at ", now(env))
    drop_car(drop_off_station_id, selected_car_id, parking_place_id, now(env))

    # make the payment
    global revenues += req.Rev
    
end

################################ Offline mode######################################
@resumable function request_arrival_process_offline_mode(env::Environment, scenario::Scenario, sol::Solution)
    #browse all the requests (the requests are sorted according to their arrival time)
    for req in eachrow(scenario.request_list)
        if failed
            break
        end

        # waiting until a new request is arrived 
        sleeping_time = req.ST - now(env)
        sleeping_time > 0 && @yield timeout(env, sleeping_time)
        print_simulation && println("Customer [", req.reqId, "]: request arrive at ", now(env))

        # in the offline mode we can check if the request is decided to be served according to the decision variables
        if isempty(req.fp)
            # the resuest is refused by the decision variables
            print_simulation && println("Customer [", req.reqId, "]: the requests is rejected according to the decision variables")
            continue
        else
            path = scenario.feasible_paths[req.fp[1], :]
            pickup_station_id = findfirst(get_potential_locations() .== path.origin_station[1])
            drop_off_station_id = findfirst(get_potential_locations() .== path.destination_station[1])
            selected_car_id, parking_place_id = -1, -1
            # check the availabilty of paths to serve the request
            if isempty(pickup_station_id) || isempty(drop_off_station_id)
                # we couldn't serve the request there is no feasible path
                print_simulation && printstyled(stdout, "Customer [", req.reqId, "]: We can not serve the request there is no feasible path\n", color=:light_red)
                global failed = true
                break
            else
                @process perform_the_trip_process_offline_mode(env, req, sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
            end
        end

    end
end

@resumable function perform_the_trip_process_offline_mode(env::Environment, req::DataFrameRow, sol::Solution, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    global potential_locations #tp get the node id of station
    
    walking_duration = get_walking_time(req.ON, potential_locations[pickup_station_id])
    driving_duration = get_trip_duration(potential_locations[pickup_station_id], potential_locations[drop_off_station_id])

    print_simulation && println("Customer [", req.reqId, "]: the request is accepted")
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    walking_duration > 0 && @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station N° ", potential_locations[pickup_station_id], " at ", now(env))

    # check and select a car to performe the trip
    refrech_battery_levels(pickup_station_id, now(env)) #refrech the battery leveles 
    available_car_df = filter(row -> row.status == CAR_PARKED &&
            row.last_battery_level >= get_battery_level_needed((pickup_station_id, drop_off_station_id)),
        stations[pickup_station_id].cars)

    if isempty(available_car_df)
        #double check if there is another request (at this moment) that can help to serve the current request
        # to do so only we need to handle this request again after we done with all the requests at this moment 
        @yield timeout(env, 0.0, priority=-1) # put this process at the end of the heap for the current time
        #check again
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
        car.last_battery_level[1] -= get_trip_battery_consumption(potential_locations[pickup_station_id],
            potential_locations[drop_off_station_id],
            car.car_type[1])
        append!(stations[drop_off_station_id].cars, car)
    else
        print_simulation && printstyled(stdout, "Customer [$(req.reqId)]: (offline mode) there is no available car at the station \n"; color=:light_red)
        global failed = true
        return
    end

    print_simulation && println("Customer [", req.reqId, "]: start ridding the car number $selected_car_id at $(now(env))")
    start_driving(pickup_station_id, selected_car_id) # this function put failed true if there is not a car with the consumption constraints
    
    #simulate the riding time
    driving_duration > 0 && @yield timeout(env, driving_duration)

    #get the available parking places
    available_places_df = filter(row -> row.status == P_FREE, stations[drop_off_station_id].parking_places)
    if isempty(available_places_df)
        #double check if there is another request that can help to serve the current request
        @yield timeout(env, 0.0, priority=-1) # put this process at the end of the heap for the current time
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
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", potential_locations[drop_off_station_id], " at ", now(env))
    drop_car(drop_off_station_id, selected_car_id, parking_place_id, now(env))
    # make the payment
    global revenues += req.Rev

end

############################### Genrale functions ###############################
function start_driving(station_id, car_id)

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
end

function drop_car(drop_off_station_id, car_id, parking_place_id, current_time)
   
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

    global number_of_served_requests +=1
end

function initialize_sim(sol::Solution, scenario_id::Int64)
    scenario = scenario_list[scenario_id]
    global failed = false
    global revenues = 0
    global stations = Array{Station,1}()

    #set the stations
    car_id = 1
    for i in 1:length(sol.open_stations_state)
        initial_cars_number = sol.initial_cars_number[i]
        push!(stations, Station(get_potential_locations()[i], initial_cars_number, car_id))
        car_id += initial_cars_number
    end

    if !online_request_serving
        #create a new column in the scenario.requets_list to get the selected path
        fp = [Int[] for i in eachindex(scenario.request_list.reqId)]
        for i in eachindex(sol.selected_paths[scenario_id])
            if sol.selected_paths[scenario_id][i]
                req = scenario.feasible_paths.req[i]
                fp[findfirst(==(req), scenario.request_list.reqId)] = [i]
            end
        end
        scenario.request_list.fp = fp
    end
end

function E_carsharing_sim(sol::Solution, scenario_id::Int64)
    scenario = scenario_list[scenario_id]
    
    initialize_sim(sol, scenario_id)

    sim = Simulation()
    if online_request_serving
        @process request_arrival_process_online_mode(sim, scenario, sol)
    else
        @process request_arrival_process_offline_mode(sim, scenario, sol)
    end
    run(sim)

    if failed
        print_simulation && printstyled(stdout, "Fatal Error: The simulation is stopped because there is a problem\n", color=:light_red)
        return Inf
    else
        # count the objective function
        print_simulation && println("counting the objective function")
        
        total_cars_cost = sum(Array{Float64}([vehicle_specific_values[i][:car_cost] for i in vcat([stations[i].cars.car_type for i in eachindex(stations)]...)]))

        #= total_station_cost = sum([station.charging_station_base_cost + 
                                    station.max_number_of_charging_points * station.charging_point_cost_fast
                                    for station in stations]) =#
        
        total_station_cost = sum(Array{Float64}([station.charging_station_base_cost for station in stations[sol.open_stations_state]]))

        return revenues - (total_cars_cost + total_station_cost) / cost_factor
    end
end

function E_carsharing_sim(sol::Solution)
    # set the counter to keep post on the number of requets that we served
    global number_of_served_requests = 0;

    #check the feasibilty of the solution
    if is_feasible_solution(sol)
        if online_request_serving
            # we have to reset the selected paths inside the solution only for the selected paths
            global online_selected_paths = [Array{Bool, 1}(falses(nrow(scenario_list[sc].feasible_paths))) for sc in eachindex(scenario_list)]
        end
        
        f_x  = 0
        for i in eachindex(scenario_list)
            f_x_sc = E_carsharing_sim(sol, i)
            if f_x_sc == Inf
                # penality expersion is as follow: 10^(16 - % of served request/10) so if 
                return 10 ^16 * 10 ^ (-1 *  (number_of_served_requests /sum(nrow(sc.request_list) for sc in scenario_list) * 100 / 10))
            end
            f_x += f_x_sc
        end
        return -1 * f_x/length(scenario_list)
    else
        print_simulation && println("The solution is not feasible")
        return (10 ^(16 - (number_of_served_requests /sum(nrow(sc.request_list) for sc in scenario_list) * 100 / 10)))
    end
end

function save_sol(sol::Solution)
    println("save solution")
    serialize("/Users/taki/Desktop/Preparation doctorat ERM/Projects/GIHH_V2.0/sol.jls", sol)
end

