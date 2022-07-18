#=
    description: create the graph of the city to solve the path problems using djisktra as an example
    inputs:
        @xml_path => the file where the data are stored (e.g., longitude, latitude, etc)
        @save_file (optional) => the path where the MetaDiGraph object will be stored
    outputs:
       MetaDiGraph object => the created network
=#
function create_graph_from_XML(xml_path::String; save_file::String="")
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
            parse(Int64, current_node["type"])
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
    scenario_request_list_df = CSV.read(scenario_path, DataFrame, header=false)
    rename!(scenario_request_list_df, :Column1 => :reqId)

    #perform the join instruction to get all the details of the requests 
    scenario_df = innerjoin(all_request_df, scenario_request_list_df, on=:reqId)
    # sort the request according to their arriving time
    sort!(scenario_df, [:ST])
    #as julia start indexing from 1 we have to add 1 to all origin and distination nodes
    scenario_df.ON .+= 1
    scenario_df.DN .+= 1
    
    scenario_df.reqId = collect(1:nrow(scenario_df))
    
    return scenario_df
end

function get_trip_info_for_request(req, sol::Solution, current_time) # information to return 
    pickup_station_id = -1 # the id of the station (inside the solution object) where we can take a car
    drop_off_station_id = -1 # the id the station (inside the solution object) where we can drop off the car
    selected_car_id = -1 # the id of the car to be used to perform the trip
    parking_place_id = -1 # the id of the parking place in the drop off station

    #= if !online_request_serving
        #check if the request is selected to be served according to the decision variables
        paths = filter(row -> row.req == req.reqId, all_feasible_paths[sol.selected_paths, :])

    else
        #get the path 
        paths = paths = filter(row -> row.req == req.reqId, all_feasible_paths)
    end =#
    paths = req.feasible_paths

    for path in eachrow(paths)
        # if we are in the offline mode necessarilly paths contains only one path 

        # reset the information to be returned
        selected_car_id = -1 # the id of the car to be used to perform the trip
        parking_place_id = -1 # the id of the parking place in the drop off station
        pickup_station_id = findfirst(in(path.origin_station), sol.open_stations_ids)
        drop_off_station_id = findfirst(in(path.destination_station), sol.open_stations_ids)

        battery_level_needed = get_battery_level_needed(path) # always 100%
        expected_start_riding_time = work_with_time_slot ? ceil(Integer, current_time + get_walking_time(req.ON, path.origin_station[1])) : current_time + get_walking_time(req.ON, path.origin_station[1])
        

        #get the list (as DataFrame) of cars that are available for the customer (parked cars + expected to arrive befor the starting time)
        available_car_df = filter(row -> row.status == CAR_PARKED ||
                (row.status == CAR_ON_WAY && row.expected_arrival_time <= expected_start_riding_time),
            stations[pickup_station_id].cars)
        if !isempty(available_car_df)
            #check the battery consumption constraint
            # first we count the actual battery levels
            refrech_battery_levels(pickup_station_id, current_time) # count the actual battery levels

            # second, count the expected battery level at the time of departure of the trip.
            expected_battery_levels = get_expected_battery_levels(available_car_df, expected_start_riding_time)

            #finaly, get the list of cars that meet the consumption constraint 
            car_index = findall(x -> x >= battery_level_needed, expected_battery_levels)

            if !isempty(car_index)
                # here at least there is a car that meet the consumption constraint

                # first we privilege a parked car
                potential_selected_parked_cars_df = filter(row -> row.status == CAR_PARKED, available_car_df[car_index, :])
                if !isempty(potential_selected_parked_cars_df)
                    # simply we select the first one --> or we can privilege  the car that has the maximum battery level
                    selected_car_id = potential_selected_parked_cars_df[1, :].car_id
                else
                    # or we select a comming car
                    selected_car_id = available_car_df[1, :].car_id
                end
            end
        end

        #check if we could select a car from the pickup station
        if selected_car_id == -1
            # the current path could not be used because we didn't find an available car for the trip
            continue # see the next path
        end

        #check the availability of the parking place (free places + expected to be free places)
        expected_arrival_time = expected_start_riding_time + get_trip_duration(path.origin_station[1], path.destination_station[1])
        expected_arrival_time = work_with_time_slot ? ceil(Integer, expected_arrival_time) : expected_arrival_time
        # first we group all the information in one Dataframe
        parking_and_cars_df = leftjoin(stations[drop_off_station_id].parking_places, stations[drop_off_station_id].cars, on=:cars => :car_id, makeunique=true)

        #get the available places
        available_places_df = filter(row -> row.status == P_FREE ||
                (!ismissing(row.status_1) && row.status_1 == CAR_RESERVED && row.start_reservation_time <= expected_arrival_time
                 && row.pending_reservation == 0),
            parking_and_cars_df)

        if !isempty(available_places_df)
            # here we are sure that there is a place
            immediately_free_parking_place = filter(row -> row.status == P_FREE, available_places_df)
            if !isempty(immediately_free_parking_place)
                # simply we select the first free place
                parking_place_id = immediately_free_parking_place[1, :].p_id
            else
                # or we select any place
                parking_place_id = available_places_df[1, :].p_id
            end
        else
            #the path can not be used because there is no available place
            continue
        end

    end

    return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)

end

#=
    the preprocessing procedure
=#
function get_all_requests_feasible_paths(requests_list::DataFrame, stations::Vector{T}, maximum_walking_time) where T<: Integer
    paths = DataFrame(req=Integer[], origin_station=Integer[], destination_station=Integer[], start_driving_time = [], arriving_time = [], Rev = [])

    for req in eachrow(requests_list)

        print_preprocessing && println("we are with req ", req.reqId)
        #check all possible starting station
        for origin_station_id in stations

            walking_time_to_origin_station =  get_walking_time(req.ON, origin_station_id) 
            print_preprocessing && println("--> walking duration from ", req.ON, " to the station $origin_station_id is $walking_time_to_origin_station")
            #check accessibilty to the origin station (walking time)
            if walking_time_to_origin_station <= maximum_walking_time

                print_preprocessing && println("---->the origin station $origin_station_id is accesible")
                #check the ending station
                for destination_station_id in stations
                    if origin_station_id != destination_station_id
                        walking_time_to_destination_node = get_walking_time(destination_station_id, req.DN)
                        print_preprocessing && println("---->walking duration from station $destination_station_id to distination node ", req.DN, " is $walking_time_to_destination_node")

                        if walking_time_to_destination_node <= maximum_walking_time
                            print_preprocessing && println("------>the destination node is accesible from the station $destination_station_id")
                            #check the total length
                            trip_duration = get_trip_duration(origin_station_id, destination_station_id)
                            total_time_for_trip = walking_time_to_destination_node + walking_time_to_origin_station + trip_duration
                            print_preprocessing && println("-------->the total trip duration considering start station $origin_station_id and destination station $destination_station_id is $total_time_for_trip")
                            # check the total length of the trip constraints
                            if total_time_for_trip <= 1.1 * get_trip_duration(req.ON, req.DN)
                                print_preprocessing && println("---------->the trip is accepted")
                                start_time = walking_time_to_origin_station + req.ST * 5 # as the times in the request file representes the time slots
                                arriving_time = start_time + trip_duration

                                if work_with_time_slot
                                    start_time  = ceil(Int, start_time / time_slot_length)
                                    arriving_time = start_time + ceil(Int, trip_duration / time_slot_length)
                                end
                                push!(paths, (req.reqId, origin_station_id, destination_station_id, start_time, arriving_time, req.Rev))
                            end
                        end
                    end
                end
            end

        end
    end
    paths
end

function get_walking_time(id_node1, id_node2)
        
        if id_node1 ∉ keys(shortest_walking_paths)
            merge!(shortest_walking_paths, Dict(id_node1 => dijkstra_shortest_paths(non_directed_manhaten_city_graph, [id_node1])))
        end
        walking_time = shortest_walking_paths[id_node1].dists[id_node2] / walking_speed / 60 #converte it to minutes
        walking_time
end

function get_trip_duration(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_graph, [id_node1])))
    end
    driving_duration = shortest_car_paths[id_node1].dists[id_node2] / driving_speed / 1000 * 60 # convert it to minutes 
    driving_duration
end

function get_trip_distance(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_graph, [id_node1])))
    end
    distance = shortest_car_paths[id_node1].dists[id_node2]
    distance
end

function get_trip_battery_consumption(id_node1::Integer, id_node2::Integer, Tcar::car_type)
    α = vehicle_specific_values[Tcar][:α]
    β = vehicle_specific_values[Tcar][:β]
    γ = vehicle_specific_values[Tcar][:γ]
    vij = driving_speed * 1000 / 3600
    battery_capacity = vehicle_specific_values[Tcar][:battery_capacity]

    distance = get_trip_distance(id_node1, id_node2)
    percentage_energy_consumption = distance * (α + β * vij^2 + γ / vij) / battery_capacity * 100

    percentage_energy_consumption
    #total_distance = battery_capacity / (α + β * driving_speed^2 + γ / driving_speed)
end

function get_potential_locations()
    collect(filter_vertices(manhaten_city_graph, :type, 2))
end

function get_battery_level_needed(path)
    return 100
end

#=
    this function calculate the expected battery level at time @btime
=#
function get_expected_battery_levels(available_cars::DataFrame, time)
    expected_battery_level = copy(available_cars.last_battery_level)

    for i in 1:nrow(available_cars)
        current_car = available_cars[i, :]
        #count the expected battery level 
        if current_car.status in [CAR_PARKED, CAR_RESERVED]
            charging_time = (time - current_car.start_charging_time) * 60 # convert it to seconds
            if work_with_time_slot
                charging_time *= time_slot_length
            end
            battery_capacity = vehicle_specific_values[current_car.car_type][:battery_capacity]
            charging_rate = vehicle_specific_values[current_car.car_type][:fast_charging_rate]

            percentage_charging_amount = 100 * charging_time / (battery_capacity / charging_rate)

            expected_battery_level[i] = min(100, current_car.last_battery_level + percentage_charging_amount)
        end
    end
    expected_battery_level
end

#=
    this function recounts the battery level of the cars in the station 
=#
function refrech_battery_levels(station_id::Integer, current_time)

    for car in eachrow(stations[station_id].cars)

        if car.status in [CAR_PARKED, CAR_RESERVED] && car.last_battery_level < 100  # in both these cases the cars is supposed to be plugged-in the charger spot
            charging_time = (current_time - car.start_charging_time) * 60 # convert it to seconds
            if work_with_time_slot
                charging_time *= time_slot_length
            end
            battery_capacity = vehicle_specific_values[car.car_type][:battery_capacity]
            charging_rate = vehicle_specific_values[car.car_type][:fast_charging_rate]

            percentage_charging_amount = 100 * charging_time / (battery_capacity / charging_rate)

            car.last_battery_level = min(100, car.last_battery_level + percentage_charging_amount)
            car.start_charging_time = current_time
        end
    end
end

#=
    in this function we don't do any check about the paramaters. the necessary check is done befor

        inputs:
        @pickup_station_id: the station where the car to pick up is parked
        @drop_off_station_id: the station where the car will be parked after performing the trip
        @car_id: the id of the car to be taken
    outputs:
            Nothing
=#
function book_trip(sol::Solution, pickup_station_id, drop_off_station_id, car_id,
    parking_place_id, start_trip_time, expected_arriving_time)
    # get the car index inside the data frame
    car_indx = findfirst(x -> x == car_id, stations[pickup_station_id].cars.car_id)

    if car_indx === nothing
        printstyled(stdout, "Error: We can not book the trip the selected car is not parked in the station\n", color=:light_red)
        global failed = true
        return
    end

    #  reserve the car
    if stations[pickup_station_id].cars.status[car_indx] == CAR_RESERVED
        printstyled(stdout, "Error: We can not book the trip the selected car is already reserved\n", color=:light_red)
        global failed = true
        return
    end

    if stations[pickup_station_id].cars.status[car_indx] != CAR_ON_WAY
        stations[pickup_station_id].cars.status[car_indx] = CAR_RESERVED
    else
        stations[pickup_station_id].cars.pending_reservation[car_indx] += 1
    end

    stations[pickup_station_id].cars.start_reservation_time[car_indx] = start_trip_time


    #decrease the battery level
    stations[pickup_station_id].cars.last_battery_level[car_indx] -= get_trip_battery_consumption(sol.open_stations_ids[pickup_station_id], sol.open_stations_ids[drop_off_station_id], stations[pickup_station_id].cars.car_type[car_indx])


    # reserve the parking space
    #check if the parking place is occupied or resereved ( it will be free by the arriving time)
    if stations[drop_off_station_id].parking_places[parking_place_id, :].status in [P_OCCUPIED, P_RESERVED]
        stations[drop_off_station_id].parking_places.pending_reservation[parking_place_id] += 1
    else
        # the place is free so we reserved it directely
        stations[drop_off_station_id].parking_places.status[parking_place_id] = P_RESERVED
    end

    # change the status of the car and precise to the drop off station that it is in its way comming 
    car = DataFrame(stations[pickup_station_id].cars[car_indx, :]) # copy
    car.expected_arrival_time[1] = expected_arriving_time
    car.status[1] = CAR_ON_WAY
    car.start_charging_time = NaN
    append!(stations[drop_off_station_id].cars, car)

    return
end

function free_parking_place(parking_place)
    if parking_place[:pending_reservation] > 0
        pending_reservation = parking_place[:pending_reservation] - 1
        status = P_RESERVED
    else
        status, pending_reservation = P_FREE, 0
    end
    status, pending_reservation, -1# always there is no car
end

function draw_graph_and_scenario(manhaten_city_graph, scenario)
    cols = []
    node_size = Array{Real,1}()
    node_marker = [:circle for i in 0:nv(manhaten_city_graph)-1]

    for i in 1:nv(manhaten_city_graph)
        if get_prop(manhaten_city_graph, i, :type) == 2
            push!(node_size, 40)
            push!(cols, :blue)
            node_marker[i] = :rect
        else
            push!(node_size, 20)
            push!(cols, :black)
        end
    end
    req_colors = [:red, :yellow, :green, :gray]
    for req in eachrow(scenario)
        if node_marker[req.ON] == :circle
            node_marker[req.ON] = :star5
            node_size[req.ON] = 40
            cols[req.ON] = req_colors[req.reqId]
        end
        if node_marker[req.DN] == :circle
            node_marker[req.DN] = :diamond
            node_size[req.DN] = 40
            cols[req.DN] = req_colors[req.reqId]
        end
    end

    edge_labels = String[]
    for e in collect(edges(manhaten_city_graph))
        weight = get_prop(manhaten_city_graph, e.src, e.dst, :weight)
        push!(edge_labels, string(weight))
    end

    f, ax, p = graphplot(manhaten_city_graph, edge_width=[3 for i in 1:ne(manhaten_city_graph)],
        node_size=node_size,
        node_color=cols,
        nlabels=[string(i) for i in 1:nv(manhaten_city_graph)],
        elabels=edge_labels,
        node_marker=node_marker)

    deregister_interaction!(ax, :rectanglezoom)

    register_interaction!(ax, :nhover, NodeHoverHighlight(p))
    register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
    register_interaction!(ax, :ndrag, NodeDrag(p))

    f
end

function initilaize_scenario(scenario_path::String, sol::Solution)
    # construct the requests lists 
    scenario = scenario_as_dataframe(scenario_path)
    #= global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
    global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph
     =#
    # preprossesing procedure 
    global all_feasible_paths = get_all_requests_feasible_paths(scenario, get_potential_locations(), maximum_walking_time)
    
    selected_paths = !online_request_serving ? all_feasible_paths[sol.selected_paths, :] : all_feasible_paths

    # put each feasible path beside its request
    scenario.feasible_paths = [filter( x-> x.req == i, selected_paths) for i in scenario.reqId]

    scenario # return
end