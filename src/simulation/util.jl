
######################### General functions ############################
"""
    create_graph_from_XML(xml_path::String; save_file::String="", GraphType=MetaGraph, weights_from="length")
    
description: create the graph of the city to solve the path problems using djisktra as an example
inputs:
    @xml_path => the file where the data are stored (e.g., longitude, latitude, etc)
    @save_file (optional) => the path where the MetaDiGraph object will be stored
    @GraphType (optional) => the type of graph to be created (MetaGraph or MetaDiGraph)
    @weights_from (optional) => the type of weights to be used (length, walking_time, driving_time)
outputs:
       Meta(Di)Graph object => the created network
"""
function create_graph_from_XML(xml_path::String; save_file::String="", GraphType=MetaGraph, weights_from="length")
    @assert GraphType in [MetaGraph, MetaDiGraph] "the graph Type need to be Meta(Di)Graph "
    @assert weights_from in ["length", "walking_time", "driving_time"] "the weights_from need to be length, walking_time or driving_time"

    doc = readxml(xml_path)

    #read the city nodes
    nodes = findall("//node", doc)

    #initialize the MetaDiGraph
    graph = GraphType(length(nodes))

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
    links = findall("//link", doc)

    for current_link in links
        tail = parse(Int64, current_link["tail"]) + 1
        head = parse(Int64, current_link["head"]) + 1

        child_elems = elements(current_link)
        length = parse(Float64, nodecontent(child_elems[1]))
        frc = parse(Float64, nodecontent(elements(child_elems[2])[1]))
        maxSpeed = parse(Float64, nodecontent(elements(child_elems[2])[2]))
        # here we see which weight we use
        if weights_from == "length"
            weight = length
        elseif weights_from == "walking_time"
            weight = ceil(Int, length / 1.34) # is in second
        elseif weights_from == "driving_time"
            weight = ceil(Int, length / maxSpeed / 1000 * 3600) # convert to seconds
        end
        properties = Dict(:weight => weight, :frc => frc, :maxSpeed => maxSpeed)
        # add the add_edge 
        if !add_edge!(graph, tail, head, properties)
            if has_edge(graph, tail, head)
                #=if true  properties[:weight] < get_prop(graph, tail, head, :weight)  && properties[:weight] != 0 =#
                set_prop!(graph, tail, head, :weight, properties[:weight])
                #end
            end
        end
    end

    if save_file != ""
        savegraph(save_file, graph)
    end
    graph
end

function draw_graph_and_scenario(manhaten_city_driving_graph, scenario)
    cols = []
    node_size = Array{Real,1}()
    node_marker = [:circle for i in 0:nv(manhaten_city_driving_graph)-1]

    for i in 1:nv(manhaten_city_driving_graph)
        if get_prop(manhaten_city_driving_graph, i, :type) == 2
            push!(node_size, 40)
            push!(cols, :blue)
            node_marker[i] = :rect
        else
            push!(node_size, 20)
            push!(cols, :black)
        end
    end
    req_colors = [:red, :yellow, :green, :gray]
    for req in eachrow(scenario.request_list)
        if node_marker[req.ON] == :circle
            node_marker[req.ON] = :star5
            node_size[req.ON] = 40
            cols[req.ON] = req_colors[1]
        end
        if node_marker[req.DN] == :circle
            node_marker[req.DN] = :diamond
            node_size[req.DN] = 40
            cols[req.DN] = req_colors[2]
        end
    end

    edge_labels = String[]
    for e in collect(edges(manhaten_city_driving_graph))
        weight = get_prop(manhaten_city_driving_graph, e.src, e.dst, :weight)
        push!(edge_labels, string(weight))
    end

    f, ax, p = graphplot(manhaten_city_driving_graph, edge_width=[3 for i in 1:ne(manhaten_city_driving_graph)],
        node_size=node_size,
        node_color=cols,
        nlabels=[string(i) for i in 1:nv(manhaten_city_driving_graph)],
        elabels=edge_labels,
        node_marker=node_marker)

    deregister_interaction!(ax, :rectanglezoom)

    register_interaction!(ax, :nhover, NodeHoverHighlight(p))
    register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
    register_interaction!(ax, :ndrag, NodeDrag(p))

    f
end

"""
    load_travel_speed() charge the travel speeds from the Instance file

    it return Dict((frc, maxSpeed) => travel_speeds) where for each tuple of (frc, max_speed) => travel_speeds)
    
    frc: functional road class 
    max_speed: the maximum speed allowed on the road
    travel_speeds: the travel speed for each hour of the day 

    it is useful to calculate the driving time
"""
function load_travel_speed()
    # the path where the travel_speeds file is located (or will be stored)
    global travel_speed_file_path

    if isfile(travel_speed_file_path)
        travel_speeds = deserialize(travel_speed_file_path)
    else
        #read the instance  XML file
        doc = readxml(Manhatten_network_details_file)
        # point to the travel_speeds nodes
        travel_speeds_nodes = findall("//travel_speed", doc)
        # the dict that will store the traveling speeds for each hour in a day 
        travel_speeds = Dict{Tuple,Array{Float64,1}}()
        # iterate over the travel_speeds nodes
        for current_node in travel_speeds_nodes
            travel_speed_as_str = nodecontent(current_node)
            travel_speed = parse.(Float64, split(travel_speed_as_str))[1:24] # all the days are equal I check it !!
            frc = parse(Int64, current_node["frc"])
            max_speed = parse(Int64, current_node["maxSpeed"])
            # push the travleing time to the dictionnary
            travel_speeds[(frc, max_speed)] = travel_speed
        end
        serialize(travel_speed_file_path, travel_speeds)
    end
    travel_speeds
end
######################## Scenario functions ###########################
"""
    requests_as_dataframe(scenario_path::String)
    inputs: 
        @scenario_path => the path of the file wich contain the scenario details
    output:
        scenario_df => list of requests with ther details as dataFrame
"""
function requests_as_dataframe(scenario_path::String)

    global number_of_requests_per_scenario # this variables indicates how much customers to take 
    # the scenario file contains only the id of the requests
    scenario_requests_list_df = CSV.read(scenario_path, DataFrame, header=false)
    rename!(scenario_requests_list_df, :Column1 => :reqId)

    # take only number of requests 
    scenario_requests_list_df = scenario_requests_list_df[1:number_of_requests_per_scenario, :]
    #perform the join instruction to get all the details of the requests 
    scenario_df = innerjoin(all_request_df, scenario_requests_list_df, on=:reqId)
    # sort the request according to their arriving time
    sort!(scenario_df, [:ST])
    #as julia start indexing from 1 we have to add 1 to all origin and distination nodes
    scenario_df.ON .+= 1
    scenario_df.DN .+= 1

    scenario_df.reqId = collect(1:nrow(scenario_df))

    return scenario_df
end

"""
    the preprocessing procedure
"""

function get_feasible_paths(requests_list::DataFrame, stations::Vector{Int64}, β_w)
    global work_with_time_slot

    paths = DataFrame(req=Int64[], origin_station=Int64[], destination_station=Int64[], start_driving_time=[], arriving_time=[], Rev=[])

    for req in eachrow(requests_list)
        print_preprocessing && println("we are with req ", req.reqId)
        γ = ceil(Int, 1.1 * (req.TrT * 60)) # gamma for this request

        #check all possible starting station
        for origin_station_id in stations
            walking_time_to_origin_station = get_walking_time(req.ON, origin_station_id, unit="seconds")

            # the expected starting time 
            start_driving_time = (req.ST * 5 * 60) + walking_time_to_origin_station # seconds

            print_preprocessing && println("--> walking duration from ", req.ON, " to the station $origin_station_id is $walking_time_to_origin_station")

            #push!(jl_df, (req.ON, req.DN, origin_station_id, walking_time_to_origin_station))

            #check accessibilty to the origin station (walking time)
            if walking_time_to_origin_station <= (β_w * 60)

                print_preprocessing && println("---->the origin station $origin_station_id is accesible")

                #check the ending station
                for destination_station_id in stations
                    #skip the same station
                    if destination_station_id == origin_station_id
                        continue
                    end
                    walking_time_to_destination_node = get_walking_time(destination_station_id, req.DN, unit="seconds")

                    print_preprocessing && println("---->walking duration from station $destination_station_id to distination node ", req.DN, " is $walking_time_to_destination_node")

                    if walking_time_to_destination_node <= β_w * 60


                        print_preprocessing && println("------>the destination node is accesible from the station $destination_station_id")
                        #check the total length
                        driving_duration = get_trip_duration(origin_station_id, destination_station_id, start_driving_time, unit="seconds")

                        total_trip_duration = driving_duration #= walking_time_to_destination_node + walking_time_to_origin_station + =#
                        print_preprocessing && println("-------->the total trip duration considering start station $origin_station_id and destination station $destination_station_id is $total_trip_duration")
                        # check the total length of the trip constraints

                        if total_trip_duration <= γ
                            print_preprocessing && println("---------->the trip is accepted")

                            if work_with_time_slot
                                # special traitement to get time slots indexes
                                start_driving_slot_id = ceil(Int, start_driving_time / 60 / time_slot_length) + 1 # converted it to minutes and then to time slots index
                                end_driving_slot_id = ceil(Int, driving_duration / 60 / time_slot_length) + start_driving_slot_id
                                push!(paths, (req.reqId, origin_station_id, destination_station_id, start_driving_slot_id, end_driving_slot_id, req.Rev))
                            else
                                arriving_time = start_driving_time + driving_duration
                                push!(paths, (req.reqId, origin_station_id, destination_station_id, start_driving_time, arriving_time, req.Rev))
                            end
                        end
                    end
                end
            end

        end
    end
    paths
end

"""
    description: create the scenario object which gather the requests list and their feasible paths
    inputs:
        @scenario_path: the path to the file where the request od the scenario are stored
        @id: the id of the scenario
        @check_file: if true the function will check if the scenario is already serialized 
                      and if so it will return the serialized object
    outputs:
        scenario: the object of type scenario
"""
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
        sc = Scenario(id, requests_list, afp)

        # save the file
        !isdir(dirname(serialized_file)) && mkpath(dirname(serialized_file))
        serialize(serialized_file, sc)
    end

    #return the scenario
    sc
end

function initialize_scenarios(scenario_idx::Array{Int64,1}; nbr_requests_per_scenario::Union{Nothing,Int64}=nothing)
    global scenarios_paths
    global number_of_requests_per_scenario
    if !isnothing(nbr_requests_per_scenario)
        number_of_requests_per_scenario = nbr_requests_per_scenario
    end
    global scenario_list = [initialize_scenario(scenarios_paths[scenario_idx[i]], i) for i in eachindex(scenario_idx)]
end
########################### Simulation functions ###################################
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
    pickup_station_id = -1
    drop_off_station_id = -1
    selected_car_id = -1
    parking_place_id = -1

    for path_id in req.fp
        path = scenario.feasible_paths[path_id, :]
        # reset the vars to be returned
        selected_car_id = -1 # the id of the car to be used to perform the trip
        parking_place_id = -1 # the id of the parking place in the drop off station
        pickup_station_id = findfirst(in(path.origin_station), potential_locations)
        drop_off_station_id = findfirst(in(path.destination_station), potential_locations)

        #check if the stations are opened
        if !sol.open_stations_state[pickup_station_id] || !sol.open_stations_state[drop_off_station_id]
            # at least one of the stations is close
            continue
        end

        battery_level_needed = get_battery_level_needed(path) # always 100%

        walking_duration = get_walking_time(req.ON, path.origin_station[1])
        work_with_time_slot && walking_duration != Inf && (walking_duration = ceil(Int64, walking_duration / time_slot_length))

        expected_start_riding_time = current_time + get_walking_time(req.ON, path.origin_station[1]) #= path.start_driving_time =#
        #= if req.reqId == 99
            println("expected_start_riding_time = ", expected_start_riding_time)
            walking_duration = get_walking_time(req.ON, path.origin_station[1])
            work_with_time_slot && walking_duration != Inf && (walking_duration = ceil(Int64, walking_duration / time_slot_length))
            @show current_time + walking_duration
        end =#
        #get the list (as DataFrame) of cars that are available for the customer (parked cars + expected to arrive befor the starting time)
        available_car_df = filter(row -> row.status == CAR_PARKED ||
                (row.status == CAR_ON_WAY && row.expected_arrival_time <= expected_start_riding_time),
            stations[pickup_station_id].cars)
        if !isempty(available_car_df)

            # first we count the actual battery levels
            refrech_battery_levels(pickup_station_id, current_time) # count the actual battery levels

            # second, count the expected battery level at the time of departure of the trip.
            expected_battery_levels = get_expected_battery_levels(available_car_df, expected_start_riding_time)

            #finaly, get the list of cars that meet the consumption constraint 
            car_index = findall(x -> x >= battery_level_needed, expected_battery_levels)

            if !isempty(car_index)

                # here at least there is a car that meets the consumption constraint
                selected_car_id = available_car_df.car_id[minimum(car_index)]
                #= # first we privilege a parked car
                potential_selected_parked_cars_df = filter(row -> row.status == CAR_PARKED, available_car_df[car_index, :])
                if !isempty(potential_selected_parked_cars_df)
                    # simply we select the first one --> or we can privilege  the car that has the maximum battery level
                    selected_car_id = potential_selected_parked_cars_df[1, :].car_id
                else
                    # or we select a comming car
                    selected_car_id = available_car_df[1, :].car_id
                end =#
            end
        end

        #check if we could select a car from the pickup station
        if selected_car_id == -1
            # the current path could not be used because we didn't find an available car for the trip
            continue # see the next path
        end

        #check the availability of the parking place (free places + expected to be free places)
        expected_arrival_time = expected_start_riding_time + get_trip_duration(path.origin_station[1], path.destination_station[1], expected_start_riding_time)
        # first we group all the information in one Dataframe
        parking_and_cars_df = leftjoin(stations[drop_off_station_id].parking_places, stations[drop_off_station_id].cars, on=:cars => :car_id, makeunique=true)

        #get the available places
        available_places_df = filter(row -> row.status == P_FREE ||
                (!ismissing(row.status_1) && row.status_1 == CAR_RESERVED &&
                 row.start_reservation_time <= expected_arrival_time
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
        if pickup_station_id != -1 && drop_off_station_id != -1 && selected_car_id != -1 && parking_place_id != -1
            # we need to memorize the selected requests in the online mode.
            global online_selected_paths[scenario.scenario_id][path_id] = true
            break
        end
    end

    return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)

end
"""
    get_walking_time(id_node1, id_node2)
    returning the walking time between two nodes in minutes
    it uses the manhatten_city_walking_graph

    if multiple_driving_speeds is true, then the walking time is rounded up to the next integer
    else the walking time is computed using the walking speed and no rounding is applied
"""
function get_walking_time(id_node1, id_node2; unit::String="minutes")
    @assert unit in ["minutes", "seconds"] "unit must be either minutes or seconds"
    if id_node1 ∉ keys(shortest_walking_paths)
        merge!(shortest_walking_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_length_graph, [id_node1])))
    end
    conversion_factor = unit == "minutes" ? 60 : 1
    if shortest_walking_paths[id_node1].dists[id_node2] != Inf
        walking_time = (multiple_driving_speeds == true) ? ceil(Int, (shortest_walking_paths[id_node1].dists[id_node2] / walking_speed) / conversion_factor) :
                       (shortest_walking_paths[id_node1].dists[id_node2] / walking_speed / conversion_factor)

        walking_time
    else
        Inf
    end
end

"""
    get_trip_duration(id_node1, id_node2, strat_driving_time)

    get the driving time between two nodes in minutes
    it uses the manhatten_city_driving_graph
    if global use_dynamic_speed is true then the driving time is computed using the dynamic speed
    else the driving time is computed using the fixed_driving_speed

    if multiple_driving_speeds is true, then the driving time is rounded up to the next integer
    else the driving time is computed using the a fixed driving speed and no rounding is applied
"""
function get_trip_duration(id_node1, id_node2, strat_driving_time; unit::String="minutes")
    @assert unit in ["minutes", "seconds"] "unit must be either minutes or seconds"
    # check if we count the driving time using the MaxSpeed or dynamic speeds
    if use_dynamic_speeds
        return get_trip_duration_with_dynamic_speed(id_node1, id_node2, strat_driving_time, unit=unit)
    end
    # 
    if id_node1 ∉ keys(shortest_car_paths)
        shortest_car_paths[id_node1] = dijkstra_shortest_paths(manhaten_city_driving_graph, [id_node1])
    end
    conversion_factor = unit == "minutes" ? 60 : 1
    if shortest_car_paths[id_node1].dists[id_node2] != Inf
        driving_duration = (multiple_driving_speeds == true) ? ceil(shortest_car_paths[id_node1].dists[id_node2] / conversion_factor) :
                           shortest_car_paths[id_node1].dists[id_node2] / fixed_driving_speed / 1000 * conversion_factor

        driving_duration
    else
        Inf
    end
end

"""
    spath(dist, shortest_path_state, src)
get the shortest path between two nodes src -> dist

    @inputs:
        dist::Int64 : the id of the destination node
        shortest_path_state:: the results of Dijsktra algoritme after applying it to the graph and src node
        src::Int64 : the id of the origin node

"""
spath(dist, shortest_path_state, src) = dist == src ? dist : [spath(shortest_path_state.parents[dist], shortest_path_state, src) dist]

"""
    get_trip_duration_with_frc(id_node1::Int64, id_node2::Int64, start_time::Number; unit::String="minutes)
    
    get traveling time between two nodes by not using a fixed driving speed but a variant one depending on the time
    it uses the matrices from the xml file which give the speed for each road segment at each time.
    in this first version we take the shortest path in terms of length and then we compute the time using the speed matrix

    input:
        id_node1::Int64 : the id of the origin node
        id_node2::Int64 : the id of the distination node
        start_driving_time::Float64 : the time when the trip starts in seconds
        unit::
"""

function get_trip_duration_with_dynamic_speed(id_node1::Int64, id_node2::Int64, start_driving_time::Number; unit::String="minutes")
    if id_node1 ∉ keys(shortest_car_paths)
        shortest_car_paths[id_node1] = dijkstra_shortest_paths(manhaten_city_length_graph, [id_node1])
    end
    #check if the id_node2 can be achived from id_node1
    if shortest_car_paths[id_node1].dists[id_node2] == Inf
        return Inf
    else
        path = spath(id_node2, shortest_car_paths[id_node1], id_node1)
        #convert matrix to Array
        path = [path[1, i] for i in 1:size(path)[2]]

        current_hour = start_driving_time == 0.0 ? 1 : ceil(Int64, start_driving_time / 3600)
        driving_duration = 0 #seconds
        path_length = length(path)
        travel_speeds = load_travel_speed() #load the travel speed for each hour 
        for i in 2:path_length
            #get the frc and max speed
            edge_infos = props(manhaten_city_driving_graph, path[i-1], path[i])
            frc = edge_infos[:frc]
            max_speed = edge_infos[:maxSpeed]# km/s
            edge_length = edge_infos[:weight]# meters

            #as we have values that are not shown in the Insatnce file (xml) I rounded them up to solve the problem
            max_speed = ceil(max_speed / 10) * 10

            #get the speed for thi edge on this hour
            speed = travel_speeds[(frc, max_speed)][current_hour] # m/s
            driving_duration += edge_length / speed

            #count the new current hour
            current_hour = ceil(Int64, (start_driving_time + driving_duration) / 3600)

            if current_hour > 24
                current_hour %= 24
            end
        end
        return driving_duration / 60 # convert it ro minuts
    end
end

"""
    description : get the distance between two nodes ( meters)
"""
function get_trip_distance(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_driving_graph, [id_node1])))
    end
    distance = shortest_car_paths[id_node1].dists[id_node2]
    distance
end
"""
    get_trip_battery_consumption(id_node1::Int64, id_node2::Int64, Tcar::car_type)

get the battery consumtion of the trip between two nodes

    @inputs:
        id_node1::Int64 : the id of the origin node
        id_node2::Int64 : the id of the distination node
        Tcar::car_type : the type of the car

    @outputs:
        percentage_energy_consumption::Float64 : the percentage of the battery that will be consumed
"""
function get_trip_battery_consumption(id_node1::Int64, id_node2::Int64, Tcar::car_type)
    α = vehicle_specific_values[Tcar][:α]
    β = vehicle_specific_values[Tcar][:β]
    γ = vehicle_specific_values[Tcar][:γ]
    vij = fixed_driving_speed * 1000 / 3600
    battery_capacity = vehicle_specific_values[Tcar][:battery_capacity]

    distance = get_trip_distance(id_node1, id_node2)
    percentage_energy_consumption = distance * (α + β * vij^2 + γ / vij) / battery_capacity * 100

    percentage_energy_consumption
end

function get_potential_locations()
    global potential_locations
    if length(potential_locations) == 0
        potential_locations = collect(filter_vertices(manhaten_city_driving_graph, :type, 2))
    end
    potential_locations
end

function get_battery_level_needed(path)
    return 100
end

"""
    this function calculate the expected battery level at time @time (input)
"""
function get_expected_battery_levels(available_cars::DataFrame, time)
    expected_battery_level = copy(available_cars.last_battery_level)

    for i in 1:nrow(available_cars)
        current_car = available_cars[i, :]
        #count the expected battery level 
        if current_car.status in [CAR_PARKED, CAR_RESERVED]
            charging_time = (time - current_car.start_charging_time) * 60 * (work_with_time_slot ? time_slot_length : 1)

            battery_capacity = vehicle_specific_values[current_car.car_type][:battery_capacity]
            charging_rate = vehicle_specific_values[current_car.car_type][:fast_charging_rate]

            percentage_charging_amount = 100 * charging_time / (battery_capacity / charging_rate)

            expected_battery_level[i] = min(100, current_car.last_battery_level + percentage_charging_amount)
        end
    end
    expected_battery_level
end

"""
    this function recounts the battery level of the cars in the station 
"""
function refrech_battery_levels(station_id::Int64, current_time)
    for car in eachrow(stations[station_id].cars)

        if car.status in [CAR_PARKED, CAR_RESERVED] && car.last_battery_level < 100  # in both these cases the cars is supposed to be plugged-in the charger spot
            charging_time = (current_time - car.start_charging_time) * 60 * (work_with_time_slot ? time_slot_length : 1)

            battery_capacity = vehicle_specific_values[car.car_type][:battery_capacity]
            charging_rate = vehicle_specific_values[car.car_type][:fast_charging_rate]

            percentage_charging_amount = 100 * charging_time / (battery_capacity / charging_rate)

            car.last_battery_level = min(100, car.last_battery_level + percentage_charging_amount)
            car.start_charging_time = current_time
        end
    end
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
function book_trip(pickup_station_id, drop_off_station_id, car_id,
    parking_place_id, start_trip_time, expected_arriving_time)
    # get the car index inside the data frame
    car_indx = findfirst(x -> x == car_id, stations[pickup_station_id].cars.car_id)
    #= if pickup_station_id == 27 
        @show stations[pickup_station_id].cars
        @show car_indx
    end =#
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

    if stations[pickup_station_id].cars.status[car_indx] == CAR_PARKED
        stations[pickup_station_id].cars.status[car_indx] = CAR_RESERVED
    else # CAR_ON_WAY
        #just we increment the pending reservation 
        stations[pickup_station_id].cars.pending_reservation[car_indx] += 1
    end
    stations[pickup_station_id].cars.start_reservation_time[car_indx] = start_trip_time

    #decrease the battery level
    stations[pickup_station_id].cars.last_battery_level[car_indx] -= get_trip_battery_consumption(potential_locations[pickup_station_id], potential_locations[drop_off_station_id], stations[pickup_station_id].cars.car_type[car_indx])

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
    car.start_charging_time[1] = NaN
    append!(stations[drop_off_station_id].cars, car)
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

"""a function to set the walking time"""
function set_walking_time(new_walking_time::Int64)
    global maximum_walking_time = new_walking_time
end

"""function to set the cost factor"""
function set_cost_factor(new_cost_factor::Int64)
    global cost_factor = new_cost_factor
end

############################ solution functions ###################################

"""
    generate a random solution:
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
    sol

end

"""
    Check whether or not the solution is feasible according to constraints 2, 3, 7 and 8 in Hatice paper
    inputs:
        @sol: the solution
        @all_feasible_paths: the set all feasible paths given by the preprocessing procedure(get_ all_requests_feasible_paths)
    outputs:
        Feasible => Boolean which is true if the solution is feasible, false otherwise.
"""
function is_feasible_solution(sol::Solution)
    global failed
    #check if the dimention of the solution fields are correctly defined
    if length(get_potential_locations()) != length(sol.open_stations_state) != length(sol.initial_cars_number)
        print_simulation && printstyled(stdout, "the solution fields are not correctly defined\n", color=:light_red)
        failed = true
        return false
    end
    #check the initial number of cars (constraint 7 and 8)
    for i in eachindex(sol.open_stations_state)
        if sol.initial_cars_number[i] > (sol.open_stations_state[i] ? get_prop(manhaten_city_driving_graph, potential_locations[i], :max_number_of_charging_points) : 0)
            print_simulation && println("the initial number of cars in the station ", potential_locations[i], " is greater the the total allowed number (or a stations contains cars despite it is closed")
            failed = true
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
                failed = true
                return false
            end
            #check if each station in the selected trips is open (constraint 3)

            #get the oponed stations as dataframes (the number of node in the graph)
            open_stations_df = DataFrame(station_ids=get_potential_locations()[sol.open_stations_state]) # usuful for innerjoin

            #i used inner join to take select the opened stations used in the selected paths and then compared with the number of stations used
            if nrow(innerjoin(selected_paths, open_stations_df, on=:origin_station => :station_ids)) != nrow(selected_paths) ||
               nrow(innerjoin(selected_paths, open_stations_df, on=:destination_station => :station_ids)) != nrow(selected_paths)

                print_simulation && println("there is at least one closed station in the feasible paths")
                failed = true
                return false
            end
        end
        # constraint 4, 5 and 6  are to be checked in the simulation
    end

    return true
end

function get_stored_solution(sol_id=1)
    deserialize(project_path("Data/other/GIHH_sol.jls"))
end

function save_sol(sol::Solution, path)
    serialize("/Users/taki/Desktop/Preparation doctorat ERM/Projects/E_car_sharing/$(path)", sol)
end

function load_sol(path::String)
    sol = deserialize(path)
    sol
end
################################## heuristics ########################################
"""
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
    for sc_id in eachindex(scenario_list)
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
    This function retrun the paths where the station, defined by station_id, is origin or destination.
    in addition this function must ensure that the paths are for only requests that not served by another path
    inputs:
        @sol: the solution
        @scenario: the scenario that we want to have feasible paths for it
    outputs:
        Feasible_paths => list of feasible paths 
"""

function get_request_feasible_path(sol::Solution, stations_idx::Array{Int64,1}, scenario::Scenario)
    #get the stations nodes ids 
    station_nodes_idx = get_potential_locations()[stations_idx]

    #get the paths that are 
    potential_feasible_paths = filter(scenario.feasible_paths) do fp
        fp.origin_station in station_nodes_idx || fp.destination_station in station_nodes_idx
    end
end

"""
    This function clean up the solution to serve the requests using only 
    minimal number of cars
    it alters the initial number of cars in the solution
    inputs:
        @sol: the solution
    outputs:
        the new objective function
"""
function clean_up_cars_number!(sol::Solution)
    global used_cars

    #= df = DataFrame(stations = sol.open_stations_state, cars = sol.initial_cars_number)
    CSV.write("sol.csv", df) =#
    # the way that the cars are created are basically when we instantiate the stations
    E_carsharing_sim(sol)
    new_initial_number_of_cars = Int64[]
    car_id = 1

    for station_id in eachindex(sol.initial_cars_number)
        old_init_car_numbers = car_id:car_id+sol.initial_cars_number[station_id]-1
        car_id += sol.initial_cars_number[station_id]

        push!(new_initial_number_of_cars, count([x in used_cars for x in old_init_car_numbers]))
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
function clean_up_selected_paths!(sol::Solution)
    #println("************** cleaning up selected paths **************")
    global failed
    global current_scenario_id
    global trips_to_unselect
    #first the trivial case: unselect all the trips that contains a closed station
    for sc in eachindex(sol.selected_paths)
        for fp in eachindex(sol.selected_paths[sc])
            if sol.selected_paths[sc][fp]
                #get the trip information
                origin_station_id = findfirst(get_potential_locations() .== scenario_list[sc].feasible_paths.origin_station[fp])
                destination_station_id = findfirst(get_potential_locations() .== scenario_list[sc].feasible_paths.destination_station[fp])

                if !sol.open_stations_state[origin_station_id] || !sol.open_stations_state[destination_station_id]
                    sol.selected_paths[sc][fp] = false
                    #@info "unselecting trip $fp in scenario $sc"
                end
            end
        end
    end
    # second step is to run the simulation and see if everything is okay.
    # if there still trips that cause infeasible solution we unselect them
    failed = true
    fit_value = 10^16
    set_online_mode(false)
    while failed
        trips_to_unselect = Int64[]
        fit_value = E_carsharing_sim(sol)
        sol.selected_paths[current_scenario_id][trips_to_unselect] .= false
    end

    fit_value, sol.selected_paths
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
    global failed

    opend_stations_node_ids = get_potential_locations()[sol.open_stations_state]
    stations_to_increment_cars = Int64[]
    for sc_id in eachindex(requests_list)
        #loop over the requests
        for req_id in requests_list[sc_id]
            #get the feasible trips for the current request
            curr_req_feasible_trips_id = findall(x -> scenario_list[sc_id].feasible_paths.req[x] == req_id &&
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
                #= failed = false
                E_carsharing_sim(sol)

                if !failed
                    #@info "request $req_id is served by trip $trip_id"
                    #here we served the requests so we can break the loop
                    break
                else
                    #here the selected path results a non feasible solution

                    #check if we increment the number of cars in the origin station we can serve the request

                    origin_station_id = findfirst(get_potential_locations() .== scenario_list[sc_id].feasible_paths.origin_station[trip_id])
                    #@info "try to increment the cars number in the $origin_station_id of trip $trip_id"
                    sol.initial_cars_number[origin_station_id] += 1
                    #try to serve the request again
                    failed = false##
                    E_carsharing_sim(sol, sc_id)
                    if !failed
                        push!(stations_to_increment_cars, origin_station_id)
                        #@info "serve request after increment number of cars in station $origin_station_id"
                        break
                    else
                        sol.initial_cars_number[origin_station_id] -= 1
                        sol.selected_paths[sc_id][trip_id] = false
                    end
                    #sol.selected_paths[sc_id][trip_id] = false
                end =#
            end
        end
    end
    #return the new objective function
    E_carsharing_sim(sol), sol.selected_paths, stations_to_increment_cars
end

function get_trips_station(station_id, feasible_trips)
    station_node_id = get_potential_locations()[station_id]
    trips = filter(x -> x.origin_station == station_node_id || x.destination_station == station_node_id, feasible_trips)
    trips
end

"""
    can_use_trip(sc_id::Int64, trip::DataFrameRow, sol::Solution)

    check if we can serve the trip in the solution
"""
function can_use_trip(sc_id::Int64, trip::DataFrameRow, sol::Solution)
    #list of global variables
    global stations
    
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
