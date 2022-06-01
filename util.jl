#=
    description: create the graph of the city to solve the path problems using djisktra as an example
    inputs:
        @xml_path => the file where the data are stored (e.g., longitude, latitude, etc)
        @save_file (optional) => the path where the MetaDiGraph object will be stored
    outputs:
       MetaDiGraph object => the created network
=#
function create_graph_from_XML(xml_path::String; save_file::String = "")
    doc = readxml(xml_path)

    #read the city nodes
    nodes = elements(elements(elements(doc.root)[2])[1])

    #initialize the MetaDiGraph
    graph = MetaDiGraph(length(nodes))

    #assigne the meta data to each node
    curr_index = 1
    for current_node in nodes
        #add_vertex!(graph)
        current_node_elements = elements(current_node)

        set_prop!(
            graph,
            curr_index,
            :latitude,
            parse(Float64, nodecontent(current_node_elements[1]))
        )
        set_prop!(
            graph,
            curr_index,
            :longitude,
            parse(Float64, nodecontent(current_node_elements[2]))
        )
        set_prop!(
            graph,
            curr_index,
            :type,
            parse(Int64,  current_node["type"])
        )
        if current_node["type"] == "2" || current_node["type"] == "3"
            custom_elements = elements(current_node_elements[3])

            set_prop!(
                graph,
                curr_index,
                :max_number_of_charging_points,
                parse(Int64, nodecontent(custom_elements[1]))
            )
            set_prop!(
                graph,
                curr_index,
                :max_power,
                parse(Int64, nodecontent(custom_elements[2]))
            )
            set_prop!(
                graph,
                curr_index,
                :charging_station_base_cost,
                parse(Int64, nodecontent(custom_elements[3]))
            )
            set_prop!(
                graph,
                curr_index,
                :charging_point_cost_slow,
                parse(Int64, nodecontent(custom_elements[4]))
            )
            set_prop!(
                graph,
                curr_index,
                :max_charging_rate_per_charging_point_slow,
                parse(Int64, nodecontent(custom_elements[5]))
            )
            set_prop!(
                graph,
                curr_index,
                :charging_point_cost_fast,
                parse(Int64, nodecontent(custom_elements[6]))
            )
            set_prop!(
                graph,
                curr_index,
                :max_charging_rate_per_charging_point_fast,
                parse(Int64, nodecontent(custom_elements[7]))
            )
        end

        curr_index += 1
    end

    #set the arcs and their informations
    links = elements(elements(elements(doc.root)[2])[2])

    for current_link in links
        tail = parse(Int64, current_link["tail"]) 
        head = parse(Int64, current_link["head"]) 
        length = parse(Float64, nodecontent(elements(current_link)[1]))
        frc = parse(Float64, nodecontent(elements(elements(current_link)[2])[1]))
        maxSpeed = parse(Float64, nodecontent(elements(elements(current_link)[2])[1]))
        # add the add_edge 
        Graphs.add_edge!(graph, tail + 1, head + 1)
        set_prop!(graph, tail + 1, head + 1, :weight, length)
        set_prop!(graph, tail + 1, head + 1, :frc, frc)
        set_prop!(graph, tail + 1, head + 1, :maxSpeed, maxSpeed)
    end
    if save_file != ""
        savegraph(save_file, graph)
    end
    graph
end

#=
    inputs: 
        @scenario_path => the path of the file wich contain the scenarion details
    output:
    scenario_df => list of requests with ther details as dataFrame
=#
function scenario_as_dataframe(scenario_path::String)
    # the scenario file contains only the id of the requests
    scenario_request_list_df = CSV.read(scenario_path, DataFrame, header = false)
    rename!(scenario_request_list_df, :Column1 => :reqId)

    #perform the join instruction to get all the details of the requests 
    scenario_df = innerjoin(all_request_df, scenario_request_list_df, on = :reqId)

    # sort the request according to their arriving time
    sort!(scenario_df, [:ST]) #sort the requests according to the starting time

    #as julia start indexing from 1 we have to add 1 to all origin and distination nodes
    scenario_df.ON .+= 1
    scenario_df.DN .+= 1
    return scenario_df
end


#=
    accept the request according to the decision variables
=#
function accept_request(req, sol::Solution)
    if !online_request_serving
        # check if the req is in one of the selected paths
        return req.reqId in all_feasible_paths[sol.selected_paths, :].req
    else 
        #online mode so there is no need to see the decision variables 
        return true
    end
end

function get_trip_info_for_request(req, sol::Solution, current_time)
    # information to return 
    pickup_station_id = -1 # the index of the station inside the solution
    drop_off_station_id = -1 # the index of the station inside the solution
    selected_car_id = -1 # the id of the car to be used to perform the trip
    parking_place_id =-1 # the id of the parking place in the drop off station
    
    if !online_request_serving
        #check if the request is selected to be served according to the decision variables
        path = filter(row -> row.req == req.reqId, all_feasible_paths[sol.selected_paths, :] )
        if isempty(path)
             print_simulation && println("Customer [", req.reqId, "]:the requests is rejected according to the decision variables")
            return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
        end
    else
        #get the path 
        paths = get_all_feasible_solution_for_request()
    end
         
    pickup_station_id = findfirst(in(path.origin_station), sol.open_stations_ids)
    drop_off_station_id = findfirst(in(path.destination_station), sol.open_stations_ids)
    battery_level_needed = get_battery_level_needed(path) # always 100%
    expected_starting_time = current_time + get_walking_time(req.reqId, path.origin_station[1])
    
    #get the list (as DataFrame) of cars that are available for the customer (parked cars + expected to arrive befor the starting time)
    available_car_df = filter(row -> row.status == CAR_PARKED || 
                                        (row.status == CAR_ON_WAY && row.expected_arrival_time <= expected_starting_time ),
                                        stations[pickup_station_id].cars)
    if !isempty(available_car_df)
        #check the battery consumption

        # first we count the actual battery levels
        refrech_battery_levels(sol, pickup_station_id, current_time) # count the actual battery levels
        
        # secondly, count the expected battery level at the time of departure of the trip.
        expected_battery_levels = get_expected_battery_levels(available_car_df, current_time + get_walking_time(req.reqId, path.origin_station[1]))
        
        # get the list of cars that meet the consumption constraint 
        car_index = findall(x -> x >= battery_level_needed, expected_battery_levels)
        
        if !isnothing(car_index)
            # here at least there is a car that meet the consumption constraint

            # first we privilege a parked car
            potential_selected_parked_cars_df = filter(row -> row.status == CAR_PARKED, available_car_df[car_index, :]  )
            if !isempty(potential_selected_parked_cars_df)
                # simply we select the first one --> or we can privilege  the car that has the maximum battery level
                selected_car_id = potential_selected_parked_cars_df[1, :].car_id
            else
                # or we select a comming car
                selected_car_id = available_car_ids[1, :].car_id
            end
        end
    end 

    #check if we could select a car from the pickup station
    if selected_car_id == -1
        print_simulation && println("Customer [", req.reqId, "]:we couldn't find any car in the pickup station so we can not serve the request")
        return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    end
    
    #check the availability of the parking place (free places + expected to be free places)
    expected_arrival_time = current_time + get_walking_time(req.reqId, path.origin_station[1]) + 
                                            get_trip_duration(path.origin_station[1], path.destination_station[1])
    # first we group all the information in one Dataframe
    parking_and_cars_df = leftjoin(stations[drop_off_station_id].parking_places, stations[drop_off_station_id].cars, on = :cars => :car_id, makeunique=true)
    
    #get the available places
    available_places_df = filter(row-> row.status == P_FREE || 
                                (!ismissing(row.status_1) && row.status_1 == CAR_RESERVED && row.start_reservation_time <= expected_arrival_time),
                            parking_and_cars_df )
    
    if !isempty(available_places_df)
        # here we are sure that there is a place
        potential_selected_parked_place_df = filter(row -> row.status == P_FREE, available_places_df)
        if !isempty(potential_selected_parked_place_df)
            # simply we select the first one --> or we can privilege  the car that has the maximum battery level
            parking_place_id = potential_selected_parked_place_df[1, :].p_id
        else
            # or we select a comming car
            parking_place_id = available_places_df[1, :].p_id
        end
    else
        print_simulation && println("Customer [", req.reqId, "]:there is no parking place")
        return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    end
    

    #everything went well so we return the pick up station the car id the dropoff station and the parking place id 
    return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)

end


#=
    preprocessing step
=#
function all_requests_feasible_paths(scenario::DataFrame, sol::Solution)
    paths = DataFrame( req = Integer[], origin_station = Integer[], destination_station = Integer[])
    
    for req in eachrow(scenario)
        
        print_preprocessing && println("we are with req " , req.reqId)
        #check all possible starting station
        for origin_station_id in sol.open_stations_ids
            
            walking_time_to_origin_station = get_walking_time(req.ON, origin_station_id)
            print_preprocessing && println("--> walking duration from ", req.ON ," to the station $origin_station_id is $walking_time_to_origin_station")
            #check accessibilty to the origin station (walking time)
            if walking_time_to_origin_station <= maximum_walking_time
                
                print_preprocessing && println("---->the origin station $origin_station_id is accesible")
                #check the ending station
                for destination_station_id in sol.open_stations_ids

                    walking_time_to_destination_node = get_walking_time(destination_station_id, req.DN)
                    print_preprocessing && println("---->walking duration from station $destination_station_id to distination node ", req.DN ," is $walking_time_to_destination_node")
            
                    if walking_time_to_destination_node <= maximum_walking_time

                        print_preprocessing && println("------>the destination node is accesible from the station $destination_station_id")
                        #check the total length
                        trip_duration = get_trip_duration(origin_station_id, destination_station_id) 
                        total_time_for_trip = walking_time_to_destination_node + walking_time_to_origin_station + trip_duration
                        print_preprocessing &&  println("-------->the total trip duration considering start station $origin_station_id and destination station $destination_station_id is $total_time_for_trip")
                        # check the total length of the trip constraints
                        if total_time_for_trip <= req.TrT
                            print_preprocessing && println("---------->the trip is accepted")
                            push!(paths, (req.reqId, origin_station_id, destination_station_id))
                        end
                    end
                end
            end 

        end
    end
    paths
end

function get_walking_time(id_node1, id_node2)
    #= 
    #the distances are in metres 
    langitude_dis = haversine((get_prop(manhaten_city_graph, id_node1, :latitude ), -74) , (get_prop(manhaten_city_graph, id_node2, :latitude ), -74 ) )
    longitude_distance = haversine((40, get_prop(manhaten_city_graph, id_node1, :longitude )) , (40, get_prop(manhaten_city_graph, id_node2, :longitude )) )

    walking_time = (langitude_dis + longitude_distance) * walking_speed / 60 #convert it to minutes
     =#
    if id_node1 ∉ keys(shortest_walking_paths)
        merge!(shortest_walking_paths, Dict(id_node1 => dijkstra_shortest_paths(non_directed_manhaten_city_graph, [id_node1]) ) )
    end
    walking_time = shortest_walking_paths[id_node1].dists[id_node2] / walking_speed / 60 #converte it to minutes
    
    walking_time

end

function get_trip_duration(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_graph, [id_node1]) ) )
    end
    driving_duration = shortest_car_paths[id_node1].dists[id_node2] / driving_speed / 1000 * 60 # convert it to minutes 
    driving_duration
end

function get_trip_distance(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_graph, [id_node1]) ) )
    end
    distance = shortest_car_paths[id_node1].dists[id_node2]
    distance
end

function get_trip_battery_consumption(id_node1, id_node2)
    return 10
end

function get_potential_locations()
    collect(filter_vertices(manhaten_city_graph, :type, 2))
end

function book_the_car(sol::Solution, req, current_time, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    
    #get the estimated starting time
    expected_starting_time = current_time + get_walking_time(req.ON, sol.open_stations_ids[pickup_station_id])

    #get the index of the cars
    car_indx = findfirst(in([selected_car_id]), stations[pickup_station_id].cars.car_id)
    
    #reserve the car and set the time of the reservation
    stations[pickup_station_id].cars.status[car_indx] = CAR_RESERVED
    stations[pickup_station_id].cars.start_reservation_time[car_indx] = expected_starting_time

    #reserve the parking place
    stations[drop_off_station_id].parking_places.status[parking_place_id] = P_RESERVED # reserved
end

function get_battery_level_needed(path)
    return 100
end

function get_expected_battery_levels(available_cars::DataFrame, time)
    expected_battery_level = copy(available_cars.last_battery_level)
    
    for i in 1:nrow(available_cars)
        car = available_cars[i,:]
        #count the expected battery level 
        if car.status in [CAR_PARKED, CAR_RESERVED]
            charging_amount =  100 * (time - car.start_charging_time) / 
                    ( battery_capacity / charging_rate / 60   #= to covenrt to minutes =#)
            expected_battery_level[i] = min(100, car.last_battery_level + charging_amount)
        end
    end
    expected_battery_level
end

#=
    this function recounts the battery level of the cars in the station 
=#
function refrech_battery_levels(sol::Solution, station_id::Integer, current_time)
    
    for car in eachrow(stations[station_id].cars)
        #count the new 
        if car.status in [CAR_PARKED, CAR_RESERVED]
            charging_amount =  100 * (current_time - car.start_charging_time) / 
                    ( battery_capacity / charging_rate / 60   #= to covenrt to minutes =#)
            car.last_battery_level = min(100, car.last_battery_level + charging_amount)
            car.start_charging_time = current_time
        end
    end
end

#=
    inputs:
        @pickup_station_id: the station where the car to pick up is parked
        @drop_off_station_id: the station where the car will be parked after performing the trip
        @car_id: the id of the car to be taken
    outputs:
            Nothing
=#

function book_trip(sol::Solution, pickup_station_id, drop_off_station_id, car_id,
                     parking_place_id, start_trip_time, expected_arriving_time )
    # get the car index inside the data frame
    car_indx = findfirst(x -> x == car_id, stations[pickup_station_id].cars.car_id)
    
    #  reserve the car
    stations[pickup_station_id].cars.status[car_indx] = CAR_RESERVED
    stations[pickup_station_id].cars.start_reservation_time[car_indx] = start_trip_time

    
    #decrease the battery level
    stations[pickup_station_id].cars.last_battery_level[car_indx] -= get_trip_battery_consumption(sol.open_stations_ids[pickup_station_id], sol.open_stations_ids[drop_off_station_id])
    
   
    # reserve the parking space
    #check if the parking place is occupied or resereved ( it will be free by the arriving time)
    if stations[drop_off_station_id].parking_places[parking_place_id, :].status in [P_OCCUPIED, P_RESERVED]
        stations[drop_off_station_id].parking_places.pending_reservation[parking_place_id] += 1
    else
        # the place is free so we reserved it directely
        stations[drop_off_station_id].parking_places.status[parking_place_id] = P_RESERVED
    end

    # change the status of the car and precise to the drop off station that it is in its way comming 
    car =DataFrame(stations[pickup_station_id].cars[car_indx, :])
    car.expected_arrival_time[1] = expected_arriving_time
    car.status[1] = CAR_ON_WAY
    append!(stations[drop_off_station_id].cars, car)
    
    return
end

function free_parking_place(parking_place)
    if parking_place[:pending_reservation] > 0
        pending_reservation =  parking_place[:pending_reservation] - 1
        status = P_RESERVED
    else
        status, pending_reservation  = P_FREE, 0
    end
    status, pending_reservation, -1# always there is no car
end

function start_riding(station_id, car_id)
    #free the parking place
    transform!(stations[station_id].parking_places, :, 
                AsTable(:) => ByRow(x -> (x.cars == car_id) ? free_parking_place(x) : 
                (x.status, x.pending_reservation, x.cars)) 
                => [:status, :pending_reservation, :cars])

    # delete the car from the station 
    filter!(row -> (row.car_id != car_id),stations[station_id].cars)
    
    return #Nothing
end

function drop_car(drop_off_station_id, car_id, parking_place_id, current_time)
    
    # check if the parking place is free (basically it will be free but just we check if there is a problem)
    if stations[drop_off_station_id].parking_places[parking_place_id, :].status == P_OCCUPIED
        print_simulation && printstyled(stdout, "the parking place is occupied there is problem in the simulatiuon logic\n"; color=:light_red )
    end
    stations[drop_off_station_id].parking_places[parking_place_id, :].status = P_OCCUPIED
    stations[drop_off_station_id].parking_places[parking_place_id, :].cars = car_id

    # get the car index inside the data frame
    car_indx = findfirst(x -> x == car_id, stations[drop_off_station_id].cars.car_id)
    
    # change the status of the car
    stations[drop_off_station_id].cars.status[car_indx] = CAR_PARKED
    stations[drop_off_station_id].cars.start_charging_time[car_indx] = current_time

    #set the reseravtion and expected time 
    stations[drop_off_station_id].cars.start_reservation_time[car_indx] = 0.0
    stations[drop_off_station_id].cars.expected_arrival_time[car_indx] = 0.0

    return # Nothing
end

