
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

########################### Simulation functions ###################################

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
function refrech_battery_levels!(cars::DataFrame, current_time)
    for car in eachrow(cars)

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

"""a function to set the walking time"""
function set_walking_time(new_walking_time::Int64)
    global maximum_walking_time = new_walking_time
end

"""function to set the cost factor"""
function set_cost_factor(new_cost_factor::Int64)
    global cost_factor = new_cost_factor
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
        Threads.lock(shortest_walking_paths_lock)
        try
            shortest_walking_paths[id_node1] = dijkstra_shortest_paths(manhaten_city_length_graph, [id_node1])
        finally
            Threads.unlock(shortest_walking_paths_lock)
        end
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
        Threads.lock(shortest_car_paths_lock)
        try
            shortest_car_paths[id_node1] = dijkstra_shortest_paths(manhaten_city_driving_graph, [id_node1])
        finally
            Threads.unlock(shortest_car_paths_lock)
        end
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


function get_trip_distance(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        Threads.lock(shortest_car_paths_lock)
        try
            shortest_car_paths[id_node1] = dijkstra_shortest_paths(manhaten_city_driving_graph, [id_node1])
        finally
            Threads.unlock(shortest_car_paths_lock)
        end
    end
    distance = shortest_car_paths[id_node1].dists[id_node2]
    distance
end

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

############################ solution functions ###################################
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

function get_station_selected_trips(sol::Solution, station_id::Int64)
    global scenario_list 
    station_node_id = get_potential_locations()[station_id]
    station_trips = [filter(x -> x.origin_station == station_node_id ||
                         x.destination_station == station_node_id ,
                         scenario.feasible_paths[sol.selected_paths[scenario.scenario_id], :])
                         for scenario in scenario_list]
                                  
    for i in eachindex(scenario_list)
        station_trips[i].scenario_id = i .* ones(Int, nrow(station_trips[i]))
    end

    return vcat(station_trips ...)
end

function get_station_feasible_trips(station_id::Int64)
    global scenario_list 
    station_node_id = get_potential_locations()[station_id]
    station_trips = [filter(x -> x.origin_station == station_node_id ||
                         x.destination_station == station_node_id ,
                         scenario.feasible_paths)
                         for scenario in scenario_list]
                         
    for i in eachindex(scenario_list)
        station_trips[i].scenario_id = i .* ones(Int, nrow(station_trips[i]))
    end

    return vcat(station_trips ...)
end

function set_trips_to_requets_var()
    global request_feasible_trips_ids = []

    for scenario in scenario_list
        fp = [Int[] for _ in eachindex(scenario.request_list.reqId)]
        for i in eachindex(scenario.feasible_paths.req)
            push!(fp[scenario.feasible_paths.req[i]], i)
        end
        push!(request_feasible_trips_ids, fp)
    end
end