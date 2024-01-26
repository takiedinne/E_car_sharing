using Base.Threads

#= #################################################################### =#
#= ############# Global variables for all scenarios#################### =#

# creating the different Graphs if they are not created before
!isfile(Manhatten_network_driving_graph_file) && global manhaten_city_driving_graph = create_graph_from_XML(Manhatten_network_details_file,
    save_file=Manhatten_network_driving_graph_file,
    GraphType=MetaDiGraph, weights_from="driving_time")
!isfile(Manhatten_network_length_graph_file) && global manhaten_city_length_graph = create_graph_from_XML(Manhatten_network_details_file,
    save_file=Manhatten_network_length_graph_file,
    GraphType=MetaDiGraph, weights_from="length")

# loading the graphs (100% they exists)
global manhaten_city_driving_graph = (multiple_driving_speeds) ? loadgraph(Manhatten_network_driving_graph_file, MGFormat()) : loadgraph(Manhatten_network_length_graph_file, MGFormat())
global manhaten_city_length_graph = loadgraph(Manhatten_network_length_graph_file, MGFormat())

global potential_locations = Int64[] # the index of all potential stations 
locations_dict = Dict{Int64,Int64}()
for i in eachindex(get_potential_locations())
    locations_dict[get_potential_locations()[i]] = i
end
#read all requests details
global all_request_df = CSV.read(all_request_details_path, DataFrame)

global shortest_car_paths = Dict{Int64,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
global shortest_car_paths_lock = ReentrantLock()

global shortest_walking_paths = Dict{Int64,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph
global shortest_walking_paths_lock = ReentrantLock()

global scenario_list = Array{Scenario,1}() # contain all scenario instances as dataframe

global online_selected_paths = Array{Array{Bool,1},1}() #useful for the online mode to see which paths is chosen

# global variable for trak the used cars
global used_cars = Array{Int64,1}()
global used_cars_lock = ReentrantLock()
#= #################################################################### =#
#= ############# Global variables for each scenario #################### =#

global failed = Threads.Atomic{Bool}(false) #for the offline mode it is needed to stop the simulation


global number_of_served_requests = Threads.Atomic{Int64}(0) # a couter for the numlber of requests that the simulation served

global current_scenario_id = 1
global trips_to_unselect = []
global trips_to_unselect_lock = ReentrantLock()

global stations_capacity = [get_prop(manhaten_city_driving_graph, potential_locations[i], :max_number_of_charging_points) for i in eachindex(get_potential_locations())]

global request_feasible_trips_ids = []
############################### Multi threading ##################################

function E_carsharing_sim(sol::Solution)
    
    #check if we did initialize the scenarios
    isempty(scenario_list) &&  @warn "you need to initialize the scenarios !"
    
    # set the counter to keep post on the number of requets that we served
    global number_of_served_requests = Threads.Atomic{Int64}(0)
    
    global used_cars = Int64[]
    global failed = Threads.Atomic{Bool}(false)
    #check the feasibilty of the solution
    
    if  is_feasible_solution(sol)
        
        if online_request_serving 
            # we have to reset the selected paths inside the solution only for the selected paths
            global online_selected_paths = [Array{Bool,1}(falses(nrow(scenario_list[sc].feasible_paths))) for sc in eachindex(scenario_list)]
        end
        
        sc_fit = Array{Float64,1}(undef, length(scenario_list))
        
        for i in eachindex(scenario_list)
            
            if failed[]
                break
            end
            
            f_x_sc = E_carsharing_sim(sol, i)
           
            sc_fit[i] = f_x_sc
        end
        if failed[]
            #println("The solution is not feasible")
            #penality expersion is as follow: 10^(16 - % of served request/10) so if 
            
            return 10^16 * 10^(-1 * (number_of_served_requests[] / 
                    sum(nrow(sc.request_list) for sc in scenario_list) * 100 / 10))
        else
            return -1 * sum(sc_fit) / length(scenario_list)
        end
    else
        print_simulation && println("The solution is not feasible")
        return (10^(16 - (number_of_served_requests[] / sum(nrow(sc.request_list) for sc in scenario_list) * 100 / 10)))
    end
end

function E_carsharing_sim(sol::Solution, scenario_id::Int64)
    
    scenario = scenario_list[scenario_id]

    initialize_sim(sol, scenario) 

    if !is_feasible_solution(sol)
        print_simulation && @warn "Error: The solution is not feasible"
        return Inf
    end

    sim = Simulation()
    if online_request_serving
        @process request_arrival_process_online_mode(sim, scenario, sol)
    else
        @process request_arrival_process_offline_mode(sim, scenario, sol)
    end
    run(sim)

    if failed[]
        print_simulation && @warn "Fatal Error: The simulation is stopped because there is a problem"
        return Inf
    else
        # count the objective function
        print_simulation && println("counting the objective function")
        
        total_cars_cost = sum(Float64[vehicle_specific_values[ct][:car_cost] for ct in 
                                vcat([scenario.stations[st].cars.car_type for st in findall(sol.open_stations_state)]...)])
        
        total_station_cost = sum(Float64[station.charging_station_base_cost +
                                  station.max_number_of_charging_points * station.charging_point_cost_fast
                                  for station in scenario.stations[sol.open_stations_state]])
                                    
       
        return scenario.revenue - (total_cars_cost + total_station_cost) / cost_factor
    end
end

function initialize_sim_old(sol::Solution, scenario::Scenario)
    
    scenario.stations = Array{Station,1}()
    scenario.revenue = 0.0
    
    #set the stations
   
    car_id = 1
    for i in eachindex(sol.open_stations_state)
        initial_cars_number = sol.initial_cars_number[i]
        push!(scenario.stations, Station(get_potential_locations()[i], initial_cars_number, car_id))
        car_id += initial_cars_number
    end

    # we set the feasible trips Ids in fp column
    fp = [Int[] for i in eachindex(scenario.request_list.reqId)]
    if online_request_serving
        # case 1: we set all feasible paths for each requests so we can take info from it later
        for i in eachindex(scenario.feasible_paths.req)
            push!(fp[scenario.feasible_paths.req[i]], i)
        end
    else
        for i in eachindex(scenario.feasible_paths.req)
            sol.selected_paths[scenario.scenario_id][i] && push!(fp[scenario.feasible_paths.req[i]], i)
        end
    end
    scenario.request_list.fp = fp
end

function initialize_sim(sol::Solution, scenario::Scenario)
    
    if isempty(scenario.stations)
        scenario.revenue = 0.0
        #set the stations
        car_id = 1
        for i in eachindex(sol.open_stations_state)
            initial_cars_number = sol.initial_cars_number[i]
            push!(scenario.stations, Station(get_potential_locations()[i], initial_cars_number, car_id))
            car_id += initial_cars_number
        end
    else
        
        car_id = 1
        @inbounds for i in eachindex(sol.open_stations_state)
            
            if !sol.open_stations_state[i]
                scenario.stations[i].cars = DataFrame()
                scenario.stations[i].parking_places = DataFrame()
            else

                initial_cars_number = sol.initial_cars_number[i]
                
                scenario.stations[i].cars = DataFrame(car_id = collect(car_id:(car_id + initial_cars_number -1)),
                            car_type = [Smart_ED for _ in 1:initial_cars_number],
                            status = [CAR_PARKED for _ in 1:initial_cars_number],
                            last_battery_level = [100.0 for _ in 1:initial_cars_number],
                            start_charging_time = [0.0 for _ in 1:initial_cars_number],
                            start_reservation_time = [NaN for _ in 1:initial_cars_number],
                            pending_reservation = [0 for _ in 1:initial_cars_number],
                            expected_arrival_time = [NaN for _ in 1:initial_cars_number])
                
                
                #create the parking places
                total_parking_places = scenario.stations[i].max_number_of_charging_points
                scenario.stations[i].parking_places = DataFrame(p_id = collect(1:total_parking_places), 
                                    status = vcat([P_OCCUPIED for _ in 1:initial_cars_number],
                                                    [P_FREE for _ in 1:(total_parking_places - initial_cars_number)]), 
                                    cars = vcat(collect(car_id:(car_id + initial_cars_number -1)),
                                                        [-1 for _ in 1:(total_parking_places - initial_cars_number)]),
                                    pending_reservation = [0 for _ in 1:total_parking_places] )

                car_id += initial_cars_number
            end  
        end
       
    end
    
    # we set the feasible trips Ids in fp column
    fp = [Int[] for _ in eachindex(scenario.request_list.reqId)]
    if online_request_serving
        # case 1: we set all feasible paths for each requests so we can take info from it later
        @inbounds for i in eachindex(scenario.feasible_paths.req)
            push!(fp[scenario.feasible_paths.req[i]], i)
        end
    else
        @inbounds for i in eachindex(scenario.feasible_paths.req)
            sol.selected_paths[scenario.scenario_id][i] && push!(fp[scenario.feasible_paths.req[i]], i)
        end
    end

    scenario.request_list.fp = fp
    scenario.revenue = 0.0
end

function initialize_scenarios(scenario_idx::Array{Int64,1}; nbr_requests_per_scenario::Union{Nothing,Int64}=nothing)
    #global variables
    global scenarios_paths
    global number_of_requests_per_scenario
    
    if !isnothing(nbr_requests_per_scenario)
        number_of_requests_per_scenario = nbr_requests_per_scenario
    end
    global scenario_list = [initialize_scenario(scenarios_paths[scenario_idx[i]], i) for i in eachindex(scenario_idx)]

    set_trips_to_requets_var()
end

function initialize_scenario(scenario_path::String, id::Int64=-1; check_file::Bool=true)
    global maximum_walking_time
    global number_of_requests_per_scenario
    # Replace "scenario_txt_files" with "scenarios_obj"
    serialized_file = replace(scenario_path, "scenario_txt_files" => "scenarios_objects")
    # Replace ".txt" with ".jls"
    serialized_file = replace(serialized_file, r"\.txt$" => "_$(number_of_requests_per_scenario)_requests_$(maximum_walking_time)_walking_time.jls")

    if check_file && isfile(serialized_file)
        sc = deserialize(serialized_file)
        sc.scenario_id = id
    else

        # construct the requests lists 
        requests_list = requests_as_dataframe(scenario_path)

        # preprossesing procedure
        afp = get_feasible_paths(requests_list, get_potential_locations(), maximum_walking_time)

        # id of each trips
        afp.fp_id = collect(1:nrow(afp))

        # link the feasible paths to their corresponding requests
        feasible_paths_ranges = Array{Array{Int64,1},1}()
        i = 1
        for r in 1:nrow(requests_list)
            start = nothing
            last = nothing
            while i <= nrow(afp) && afp.req[i] == requests_list.reqId[r]
                isnothing(start) && (start = i)
                i += 1
            end
            last = i - 1
            isnothing(start) ? push!(feasible_paths_ranges, Int[]) : push!(feasible_paths_ranges, collect(start:last))
        end
        requests_list.fp = feasible_paths_ranges
        
        stations = [Station(get_potential_locations()[i], 0, 0) for i in eachindex(get_potential_locations())]
        sc = Scenario(id, requests_list, afp, stations, 0)
        
        
        # save the file
        !isdir(dirname(serialized_file)) && mkpath(dirname(serialized_file))
        serialize(serialized_file, sc)
    end

    
    #return the scenario
    sc
end
############################### Offline functions ###############################
@resumable function request_arrival_process_offline_mode(env::Environment, scenario::Scenario, sol::Solution)
    #browse all the requests (the requests are sorted according to their arrival time)
    @inbounds for req in eachrow(scenario.request_list)
        
        if failed[]
            #@info "scenario : $(scenario.scenario_id) end Thread: $(Threads.threadid()) "
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

            pickup_station_id = locations_dict[path.origin_station[1]]
            drop_off_station_id = locations_dict[path.destination_station[1]]
        
            selected_car_id, parking_place_id = -1, -1
            @process perform_the_trip_process_offline_mode(env, req, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id, scenario)
        end

    end
    
end

@resumable function perform_the_trip_process_offline_mode(env::Environment,
    req::DataFrameRow, pickup_station_id, drop_off_station_id, selected_car_id,
    parking_place_id, scenario::Scenario)


    global potential_locations
    pickup_station, drop_off_station = scenario.stations[pickup_station_id], scenario.stations[drop_off_station_id]
    
    #simulate teh different steps of the trip
    walking_duration = get_walking_time(req.ON, potential_locations[pickup_station_id])
    work_with_time_slot && walking_duration != Inf && (walking_duration = ceil(Int64, walking_duration / time_slot_length))

    start_walking_time = now(env)
    start_driving_time = start_walking_time + walking_duration

    driving_duration = get_trip_duration(potential_locations[pickup_station_id], potential_locations[drop_off_station_id], start_driving_time)
    work_with_time_slot && driving_duration != Inf && (driving_duration = ceil(Int64, driving_duration / time_slot_length))

    print_simulation && println("Customer [", req.reqId, "]: the request is accepted")
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    walking_duration > 0 && @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station N° ", potential_locations[pickup_station_id], " at ", now(env))

    # check and select a car to performe the trip
    cars = pickup_station.cars
    refrech_battery_levels!(cars, now(env)) #refrech the battery leveles 
    
    battery_needed = get_battery_level_needed((pickup_station_id, drop_off_station_id))
    available_cars_ids = findall(cars.status .== CAR_PARKED .&& cars.last_battery_level .>= battery_needed)

    if isempty(available_cars_ids)
        #double check if there is another request (at this moment) that can help to serve the current request
        # to do so only we need to handle this request again after we done with all the requests at this moment 
        @yield timeout(env, 0.0, priority=-1) # put this process at the end of the heap for the current time
        #check again
        available_cars_ids = findall(cars.status .== CAR_PARKED .&& cars.last_battery_level .>= battery_needed)
    end

    if !isempty(available_cars_ids)
        # simply we take the last found car
        selected_car_id = cars.car_id[available_cars_ids[end]]
        # put the selected cars to used_cars to track the used them
        Threads.lock(used_cars_lock)
        try
            if !in(selected_car_id, used_cars)
                push!(used_cars, selected_car_id)
            end
        finally
            Threads.unlock(used_cars_lock)
        end

        # update the status of the car
        pickup_station.cars.status[available_cars_ids[end]] = CAR_RESERVED

        # add the car in the drop_off station 
        new_battery_level = cars.last_battery_level[available_cars_ids[end]] - get_trip_battery_consumption(potential_locations[pickup_station_id],
            potential_locations[drop_off_station_id], cars.car_type[available_cars_ids[end]])
        push!(drop_off_station.cars, (selected_car_id, cars.car_type[available_cars_ids[end]], CAR_ON_WAY, new_battery_level,
                NaN, NaN, 0, now(env) + driving_duration )) # copy
    else
        print_simulation && @warn "Customer [$(req.reqId)]: (offline mode) there is no available car at the station "
        failed[] = true

        Threads.lock(trips_to_unselect_lock)
        try
            push!(trips_to_unselect[scenario.scenario_id], req.fp[1])
        finally
            Threads.unlock(trips_to_unselect_lock)
        end

        return #stop the simulation
    end

    print_simulation && println("Customer [", req.reqId, "]: start ridding the car number $selected_car_id at $(now(env))")
    start_driving(pickup_station, selected_car_id) # this function put failed true if there is not a car with the consumption constraints
    
    #simulate the riding time
    driving_duration > 0 && @yield timeout(env, driving_duration)

    #get the available parking places
    available_places_ids = findall(drop_off_station.parking_places.status .== P_FREE)
    if isempty(available_places_ids)
        #double check if there is another request that can help to serve the current request
        @yield timeout(env, 0.0, priority=-1) # put this process at the end of the heap for the current time
        available_places_ids = findall(drop_off_station.parking_places.status .== P_FREE)
    end
    if !isempty(available_places_ids)
        # simply we select the first free place
        parking_place_id = drop_off_station.parking_places.p_id[available_places_ids[1]]
        # reserve the place 
        drop_off_station.parking_places.status[parking_place_id] = P_RESERVED
    else
        #show the stations vars
        print_simulation && @warn "Customer [$(req.reqId)]: (offline mode) there is no free parking place at station that has id  $drop_off_station_id "
        failed[] = true
        #@info "failed is set to true"
        Threads.lock(trips_to_unselect_lock)
        try
            push!(trips_to_unselect[scenario.scenario_id], req.fp[1])
        finally
            Threads.unlock(trips_to_unselect_lock)
        end
        return
    end

    # drop the car
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", potential_locations[drop_off_station_id], " at ", now(env))
    drop_car(drop_off_station, selected_car_id, parking_place_id, now(env))
    # make the payment
    scenario.revenue += req.Rev

end

############################### Online functions ###############################
@resumable function request_arrival_process_online_mode(env::Environment, scenario::Scenario, sol::Solution)
    
    #browse all the requests (the requests are sorted according to their arrival time)
    @inbounds for req in eachrow(scenario.request_list)
        
        if failed[]
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
        if  !isnothing(selected_car_id) && !isnothing(parking_place_id)
            #book the trip (the car + parking palce ,etc) only for the 
            walking_duration = get_walking_time(req.ON, potential_locations[pickup_station_id])
            work_with_time_slot && walking_duration != Inf && (walking_duration = ceil(Int64, walking_duration / time_slot_length))

            start_walking_time = now(env)
            start_driving_time = start_walking_time + walking_duration
            
            driving_duration = get_trip_duration(potential_locations[pickup_station_id], potential_locations[drop_off_station_id], start_driving_time)
            work_with_time_slot && driving_duration != Inf && (driving_duration = ceil(Int64, driving_duration / time_slot_length))

            book_trip(scenario.stations[pickup_station_id], scenario.stations[drop_off_station_id],
                 pickup_station_id, drop_off_station_id,
                 selected_car_id, parking_place_id, start_driving_time,
                  start_walking_time + walking_duration + driving_duration)

            # for the online serving we check all the variables otherwise we have only to check the stations
            @process perform_the_trip_process_online_mode(env, req, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id, scenario)
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

            print_simulation && println("Customer [$(req.reqId)]: We can not serve the request there is no feasible path ($cause_message)")
        end
    end

end

@resumable function perform_the_trip_process_online_mode(env::Environment, req::DataFrameRow, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id, scenario)
    
    pickup_station, drop_off_station = scenario.stations[pickup_station_id], scenario.stations[drop_off_station_id]
    
    #get the different duration for the trip.
    walking_duration = get_walking_time(req.ON, potential_locations[pickup_station_id])
    work_with_time_slot && walking_duration != Inf && (walking_duration = ceil(Int64, walking_duration / time_slot_length))

    start_walking_time = now(env)
    start_driving_time = start_walking_time + walking_duration

    driving_duration = get_trip_duration(potential_locations[pickup_station_id], potential_locations[drop_off_station_id], start_driving_time)
    work_with_time_slot && driving_duration != Inf && (driving_duration = ceil(Int64, driving_duration / time_slot_length))

    print_simulation && println("Customer [", req.reqId, "]: the request is accepted") #
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    walking_duration > 0 && @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station N°", potential_locations[pickup_station_id], " at ", now(env))

    print_simulation && println("Customer [", req.reqId, "]: start ridding the car number $selected_car_id at $(now(env))")

    Threads.lock(used_cars_lock)
    try
        if !in(selected_car_id, used_cars)
            push!(used_cars, selected_car_id)
        end
    finally
        Threads.unlock(used_cars_lock)
    end
    
    # take the car and free the parking space
    start_driving(pickup_station, selected_car_id) # this function put failed true if there is not a car with the consumption constraints

    failed[] && print_simulation && @warn "Customer [$(req.reqId)]: Problem while starting the drive"


    #simulate the driving time
    driving_duration > 0 && @yield timeout(env, driving_duration)

    # drop the car
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", potential_locations[drop_off_station_id], " at ", now(env))
    #@info "drop_off scenario: $(scenario.scenario_id), request: $(req.reqId)"
    drop_car(drop_off_station, selected_car_id, parking_place_id, now(env))

    # make the payment
    scenario.revenue += req.Rev

end

############################### Genrale functions ###############################
function start_driving(station::Station, car_id::Int64)
    car_indx = findfirst(station.cars.car_id .== car_id)
    parking_place_indx = findfirst(station.parking_places.cars .== car_id)
    if car_indx === nothing || parking_place_indx === nothing
        @warn "Error: the car or the parking place are not found "
        global failed[] = true
        return
    end

    # free the parking place
    if station.parking_places.pending_reservation[parking_place_indx] > 0
        station.parking_places.pending_reservation[parking_place_indx] -= 1
        station.parking_places.status[parking_place_indx] = P_RESERVED
    else
        station.parking_places.status[parking_place_indx] = P_FREE
    end
    
    station.parking_places.cars[parking_place_indx] = -1 # there is no cars

    # delete the car from the station
    filter!(row -> (row.car_id != car_id), station.cars)
end

function drop_car(drop_off_station, car_id, parking_place_id, current_time)
    # check if the parking place is free (basically it will be free but just we check if there is a problem)
    if drop_off_station.parking_places.status[parking_place_id] == P_OCCUPIED
        println("************************************")
        @show drop_off_station.cars
        @show drop_off_station.parking_places
        println("************************************")
        #@warn "Error: the parking place is occupied there is problem in the simulatiuon logic station $(drop_off_station) "
        global failed[] = true
        return
    end

    # occupy the parking place
    drop_off_station.parking_places.status[parking_place_id] = P_OCCUPIED
    drop_off_station.parking_places.cars[parking_place_id] = car_id

    # get the car index inside the data frame
    car_indx = findfirst(x -> x == car_id, drop_off_station.cars.car_id)

    # change the status of the car
    drop_off_station.cars.status[car_indx] = CAR_PARKED
    drop_off_station.cars.start_charging_time[car_indx] = current_time

    #set the reseravtion and expected time 
    drop_off_station.cars.start_reservation_time[car_indx] = NaN
    drop_off_station.cars.expected_arrival_time[car_indx] = NaN

    global number_of_served_requests[] += 1
end

function is_feasible_solution_old(sol::Solution)
    
    #check if the dimention of the solution fields are correctly defined
    if length(get_potential_locations()) != length(sol.open_stations_state) != length(sol.initial_cars_number)
        print_simulation && printstyled(stdout, "the solution fields are not correctly defined\n", color=:light_red)
        failed[] = true
        return false
    end
    #check the initial number of cars (constraint 7 and 8)
    for i in eachindex(sol.open_stations_state)
        if sol.initial_cars_number[i] > (sol.open_stations_state[i] ? get_prop(manhaten_city_driving_graph, potential_locations[i], :max_number_of_charging_points) : 0)
            print_simulation && println("the initial number of cars in the station ", potential_locations[i], " is greater the the total allowed number (or a stations contains cars despite it is closed")
            failed[] = true
            return false
        end
    end

    if !online_request_serving
        global scenario_list
        for s in eachindex(sol.selected_paths)
            selected_paths = scenario_list[s].feasible_paths[sol.selected_paths[s], :]

            # check if the customer is served by  at most one trip (constraint 2)
            if nrow(selected_paths) != length(unique(selected_paths.req))
                print_simulation && println("there is at least a customer served by more than one trip in scenario $s")
                failed[] = true
                return false
            end
            #check if each station in the selected trips is open (constraint 3)

            #get the oponed stations as dataframes (the number of node in the graph)
            open_stations_df = DataFrame(station_ids=get_potential_locations()[sol.open_stations_state]) # usuful for innerjoin

            #i used inner join to take select the opened stations used in the selected paths and then compared with the number of stations used
            if nrow(innerjoin(selected_paths, open_stations_df, on=:origin_station => :station_ids)) != nrow(selected_paths) ||
               nrow(innerjoin(selected_paths, open_stations_df, on=:destination_station => :station_ids)) != nrow(selected_paths)

                print_simulation && println("there is at least one closed station in the feasible paths")
                failed[] = true
                return false
            end
        end
        # constraint 4, 5 and 6  are to be checked in the simulation
    end

    return true
end

function is_feasible_solution(sol::Solution)
    #global variables
    global stations_capacity
    global scenario_list

    #check if the dimension of the solution fields are correctly defined
    if length(get_potential_locations()) != length(sol.open_stations_state) != length(sol.initial_cars_number)
        print_simulation && printstyled(stdout, "the solution fields are not correctly defined\n", color=:light_red)
        failed[] = true
        return false
    end

    #check the initial number of cars (constraint 7 and 8)
    capacity_mask = sol.open_stations_state .* stations_capacity
    if sol.initial_cars_number > capacity_mask 
        print_simulation && println("the initial number of cars in the station ", potential_locations[i], " is greater the the total allowed number (or a stations contains cars despite it is closed")
        failed[] = true
        return false
    end
    
    if !online_request_serving
        stations_nodes = potential_locations[sol.open_stations_state]
        for s in eachindex(sol.selected_paths)
            selected_paths = scenario_list[s].feasible_paths[sol.selected_paths[s], :]
            # check if the customer is served by  at most one trip (constraint 2)
            if !allunique(selected_paths.req)
                print_simulation && println("there is at least a customer served by more than one trip in scenario $s")
                failed[] = true
                return false
            end

            #check if each station in the selected trips is open (constraint 3)

            #get the oponed stations as dataframes (the number of node in the graph)
            
            if !all(in.(selected_paths.origin_station,  Ref(stations_nodes))) ||
                    !all(in.(selected_paths.destination_station,  Ref(stations_nodes)))
                
                print_simulation && println("there is at least one closed station in the feasible paths")
                failed[] = true
                return false
            end
        end
        # constraint 4, 5 and 6  are to be checked in the simulation
    end

    return true
end
"""
    description: GET ALL THE INFORMATION RELATIVE TO THE REQUEST
    inputs: 
        @req => the request to be served
        @sol => the solution to get the list of the opened stations
        @current_time => the surrent time
    output:
        pickup_station_id => the id of the pickup station 
        drop_off_station_id => the id the drop off station
        selected_car_id => the id of the car to be used to perform the trip
        parking_place_id => the id of the parking place in the drop off station
"""
function get_trip_info_for_request(req, sol::Solution, scenario::Scenario, current_time)
    # vars to return 
    pickup_station_id = nothing
    drop_off_station_id = nothing
    selected_car_id = nothing
    parking_place_id = nothing
    
    for path_id in req.fp
        path = scenario.feasible_paths[path_id, :]
        # reset the vars to be returned
        selected_car_id = nothing # the id of the car to be used to perform the trip
        parking_place_id = nothing # the id of the parking place in the drop off station
        
        pickup_station_id = locations_dict[path.origin_station]
        drop_off_station_id = locations_dict[path.destination_station]
        
        pickup_station, drop_off_station = scenario.stations[pickup_station_id], scenario.stations[drop_off_station_id]
        #check if the stations are opened
        if !sol.open_stations_state[pickup_station_id] || !sol.open_stations_state[drop_off_station_id]
            # at least one of the stations is close
            continue
        end

        battery_level_needed = get_battery_level_needed(path) # always 100%

        walking_duration = get_walking_time(req.ON, path.origin_station[1])
        work_with_time_slot && walking_duration != Inf && (walking_duration = ceil(Int64, walking_duration / time_slot_length))

        expected_start_riding_time = current_time + walking_duration
        
        #get the list (as DataFrame) of cars that are available for the customer (parked cars + expected to arrive befor the starting time)
       
        available_cars_ids = findall(pickup_station.cars.status .== CAR_PARKED)
        if isempty(available_cars_ids)
            # check on way cars
            available_cars_ids = findall(pickup_station.cars.status .== CAR_ON_WAY 
                                    .|| pickup_station.cars.expected_arrival_time .<= expected_start_riding_time) 
        end
    
        if !isempty(available_cars_ids)
            
            # first we count the actual battery levels
            refrech_battery_levels!(pickup_station.cars, current_time) # count the actual battery levels

            # second, count the expected battery level at the time of departure of the trip.
            expected_battery_levels = get_expected_battery_levels(pickup_station.cars, expected_start_riding_time)[available_cars_ids]
            
            #finaly, get the list of cars that meet the consumption constraint 
            car_index = findfirst(expected_battery_levels .>= battery_level_needed)

            if !isnothing(car_index)
                # here at least there is a car that meets the consumption constraint
                selected_car_id = pickup_station.cars.car_id[available_cars_ids[car_index]]
            end
        end

        #check if we could select a car from the pickup station
        if isnothing(selected_car_id)
            # the current path could not be used because we didn't find an available car for the trip
            continue # see the next path
        end

        #check the availability of the parking place (free places + expected to be free places)
        trip_duration = get_trip_duration(path.origin_station[1], path.destination_station[1], expected_start_riding_time)
        work_with_time_slot && trip_duration != Inf && (trip_duration = ceil(Int64, trip_duration / time_slot_length))

        expected_arrival_time = expected_start_riding_time + trip_duration
        
        #get an available parking slot
        parking_place_id = findfirst(drop_off_station.parking_places.status .== P_FREE)
        if isnothing(parking_place_id)
            car_to_leave_id = findfirst(drop_off_station.cars.status .== CAR_RESERVED 
                .&& drop_off_station.cars.start_reservation_time .<= expected_arrival_time
                #= .&& drop_off_station.cars.pending_reservation .== 0 =#)

            if !isnothing(car_to_leave_id)
                parking_place_id = findfirst(drop_off_station.parking_places.cars 
                                                .== drop_off_station.cars.car_id[car_to_leave_id] .&&
                                            drop_off_station.parking_places.pending_reservation .== 0)
            else
                #the path can not be used because there is no available place
                continue
            end
        end

        # we need to memorize the selected requests in the online mode.
        global online_selected_paths[scenario.scenario_id][path_id] = true
        break
    end

    return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)

end

"""
    in this function we don't do any check about the paramaters. the necessary check is done befor

        inputs:
        @pickup_station_id: the station where the car to pick up is parked
        @drop_off_station_id: the station where the car will be parked after performing the trip
        @car_id: the id of the car to be taken
    outputs:
            Nothing
"""
function book_trip(pickup_station, drop_off_station, pickup_station_id, drop_off_station_id, car_id,
    parking_place_id, start_trip_time, expected_arriving_time)
    # get the car index inside the data frame
    car_indx = findfirst(pickup_station.cars.car_id .== car_id)
    
    if isnothing(car_indx)
        printstyled(stdout, "Error: We can not book the trip the selected car is not parked in the station\n", color=:light_red)
        global failed[] = true
        return
    end

    #  reserve the car
    if pickup_station.cars.status[car_indx] == CAR_RESERVED
        printstyled(stdout, "Error: We can not book the trip the selected car is already reserved\n", color=:light_red)
        global failed[] = true
        return
    end

    if pickup_station.cars.status[car_indx] == CAR_PARKED
        pickup_station.cars.status[car_indx] = CAR_RESERVED
    else # CAR_ON_WAY
        #just we increment the pending reservation 
        pickup_station.cars.pending_reservation[car_indx] += 1
    end
    pickup_station.cars.start_reservation_time[car_indx] = start_trip_time

    #decrease the battery level
    pickup_station.cars.last_battery_level[car_indx] -= get_trip_battery_consumption(potential_locations[pickup_station_id], potential_locations[drop_off_station_id], pickup_station.cars.car_type[car_indx])

    # reserve the parking space
    #check if the parking place is occupied or resereved ( it will be free by the arriving time)
    if drop_off_station.parking_places.status[parking_place_id] in [P_OCCUPIED, P_RESERVED]
        drop_off_station.parking_places.pending_reservation[parking_place_id] += 1
    else
        # the place is free so we reserved it directely
        drop_off_station.parking_places.status[parking_place_id] = P_RESERVED
    end
    
    push!(drop_off_station.cars, (pickup_station.cars.car_id[car_indx],
            pickup_station.cars.car_type[car_indx],
            CAR_ON_WAY, 
            pickup_station.cars.last_battery_level[car_indx],
            NaN, 
            NaN,
            0,
            expected_arriving_time)) # copy
end


"""
    generate a random solution multi threaded manner:
        1- decide how much station to open if it is not precised in @open_stations_number
        2- open random stations 
        3- set random initial number of cars for each station
    inputs:
        @open_stations_number (optional): the number of station to open
    outputs:
        sol:: Solution 
"""
function generate_random_solution(; open_stations_number=-1)
    @assert length(scenario_list) > 0 "Error: you have the initialize the scenarios first ... "

    global online_request_serving
    global rng

    sol = Solution()
    potential_locations = get_potential_locations()
    #decide the number of stations to open if it is not precised by the user
    if open_stations_number == -1
        mean_value = 43  # (1+85)/2 = 43
        std_deviation = 10  # Adjust this value to control the spread of the distribution

        # Create a truncated normal distribution between 1 and 85
        truncated_dist = Truncated(Distributions.Normal(mean_value, std_deviation), 1, 85)

        # Generate a random integer from the truncated distribution
        open_stations_number = round(Int64, rand(rng, truncated_dist))
    end
    #randomly open stations
    sol.open_stations_state[sample(rng, 1:length(potential_locations), open_stations_number, replace=false)] .= true

    #set initial car number for each station
    for i in eachindex(sol.open_stations_state)
        max_number_of_charging_points = sol.open_stations_state[i] ? get_prop(manhaten_city_driving_graph, get_potential_locations()[i], :max_number_of_charging_points) : 0
        sol.initial_cars_number[i] = rand(rng, 0:max_number_of_charging_points)
    end

    # get the selected paths according to the FIFS policy
    old_online_serving_value = online_request_serving

    set_online_mode(true)
    E_carsharing_sim(sol)
    sol.selected_paths = deepcopy(online_selected_paths)

    set_online_mode(old_online_serving_value)
    #sol = load_sol("/Users/taki/Desktop/Preparation doctorat ERM/Projects/E_car_sharing/Data/other/GIHH_sol.jls")
    clean_up_cars_number!(sol)
    sol

end


function ECS_objective_function(sol::Solution)
    @assert length(scenario_list) > 0 "Error: you have the initialize the scenarios first ... "
    
    total_cars_cost = sum(sol.initial_cars_number) * vehicle_specific_values[Smart_ED][:car_cost]

    scenario = scenario_list[1]
    
    total_station_cost = sum(Float64[station.charging_station_base_cost +
        station.max_number_of_charging_points * station.charging_point_cost_fast
        for station in scenario.stations[sol.open_stations_state]])
    
    
    revenues_as_list = Vector{Vector{Float64}}([scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :Rev] for scenario in scenario_list])
    revenues = sum(vcat(revenues_as_list ...))
    
    return -1 * (revenues / length(scenario_list) - (total_cars_cost + total_station_cost) / cost_factor)
end

