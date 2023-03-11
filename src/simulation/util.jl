
######################### General functions ############################
#=
    description: create the graph of the city to solve the path problems using djisktra as an example
    inputs:
        @xml_path => the file where the data are stored (e.g., longitude, latitude, etc)
        @save_file (optional) => the path where the MetaDiGraph object will be stored
    outputs:
       MetaDiGraph object => the created network
=#
function create_graph_from_XML(xml_path::String; save_file::String="", GraphType=MetaGraph)
    @assert GraphType in [MetaGraph, MetaDiGraph] "the graph Type need to be Meta(Di)Graph "
    
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

        properties = Dict(:weight => length, :frc => frc, :maxSpeed => maxSpeed)
        # add the add_edge 
        if !add_edge!(graph, tail, head, properties)
            if has_edge(graph, tail, head)
                if properties[:weight] < get_prop(graph, tail, head, :weight) #= && properties[:weight] != 0 =#
                    set_prop!(graph, tail, head, :weight, properties[:weight])
                end

            end
        end
    end

    if save_file != ""
        savegraph(save_file, graph)
    end
    graph
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
######################## Scenario functions ###########################
#=
    inputs: 
        @scenario_path => the path of the file wich contain the scenarion details
    output:
        scenario_df => list of requests with ther details as dataFrame
=#
function requests_as_dataframe(scenario_path::String)
    # the scenario file contains only the id of the requests
    scenario_requests_list_df = CSV.read(scenario_path, DataFrame, header=false)
    rename!(scenario_requests_list_df, :Column1 => :reqId)

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

#=
    the preprocessing procedure
=#
function get_feasible_paths(requests_list::DataFrame, stations::Vector{T}, maximum_walking_time) where {T<:Int64}
    paths = DataFrame(req=Int64[], origin_station=Int64[], destination_station=Int64[], start_driving_time=[], arriving_time=[], Rev=[])
    β_w = (work_with_time_slot) ? maximum_walking_time / time_slot_length : maximum_walking_time
    for req in eachrow(requests_list)

        print_preprocessing && println("we are with req ", req.reqId)
        #check all possible starting station
        for origin_station_id in stations

            walking_time_to_origin_station = get_walking_time(req.ON, origin_station_id)
            print_preprocessing && println("--> walking duration from ", req.ON, " to the station $origin_station_id is $walking_time_to_origin_station")
            #check accessibilty to the origin station (walking time)
            if walking_time_to_origin_station <= β_w

                print_preprocessing && println("---->the origin station $origin_station_id is accesible")
                #check the ending station
                for destination_station_id in stations
                    if destination_station_id == origin_station_id
                        continue
                    end
                    walking_time_to_destination_node = get_walking_time(destination_station_id, req.DN)
                    print_preprocessing && println("---->walking duration from station $destination_station_id to distination node ", req.DN, " is $walking_time_to_destination_node")

                    if walking_time_to_destination_node <= β_w
                        print_preprocessing && println("------>the destination node is accesible from the station $destination_station_id")
                        #check the total length
                        trip_duration = get_trip_duration(origin_station_id, destination_station_id)
                        if trip_duration == Inf
                            #there is no path which link, the two stations
                            continue
                        end

                        total_time_for_trip = walking_time_to_destination_node + walking_time_to_origin_station + trip_duration
                        print_preprocessing && println("-------->the total trip duration considering start station $origin_station_id and destination station $destination_station_id is $total_time_for_trip")
                        # check the total length of the trip constraints

                        if total_time_for_trip <= 1.1 * get_trip_duration(req.ON, req.DN) &&
                            get_trip_battery_consumption(origin_station_id, destination_station_id, Smart_ED) <= get_battery_level_needed(())
                            
                            print_preprocessing && println("---------->the trip is accepted")

                            #=  
                                as the times in the scenario files represents the number of time slots so we have 
                                to do some special treatement 
                            =#

                            start_time = walking_time_to_origin_station + req.ST * (work_with_time_slot ? 1 : time_slot_length)
                            arriving_time = start_time + trip_duration

                            #= if work_with_time_slot
                                start_time = ceil(Int, start_time / time_slot_length)
                                arriving_time = start_time + ceil(Int, trip_duration / time_slot_length)
                            end =#

                            push!(paths, (req.reqId, origin_station_id, destination_station_id, start_time, arriving_time, req.Rev))
                        end
                    end
                end
            end

        end
    end
    paths
end

#=
    description: create the scenario object which gather the requests list and their feasible paths
    inputs:
        @scenario_path: the path to the file where the request od the scenario are stored
        @id: the id of the scenario
        @check_file: if true the function will check if the scenario is already serialized 
                      and if so it will return the serialized object
    outputs:
        scenario: the object of type scenario
 =#
function initialize_scenario(scenario_path::String, id::Int64=-1; check_file::Bool=true)
    #println("we are initializing scenarion number $id")
    file_path = serialized_scenarios_folder * "/scenario$id.sc"
    if check_file && isfile(file_path)
        deserialize(file_path)
    else
       
        # construct the requests lists 
        requests_list = requests_as_dataframe(scenario_path)

        # preprossesing procedure
        afp = get_feasible_paths(requests_list, get_potential_locations(), maximum_walking_time)

        # link the feasible paths to their corresponding requests
        feasible_paths_ranges = Array{Array{Int64,1}, 1}()
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
        id != -1 && serialize(file_path, sc)
        sc
    end
end

function initialize_scenarios(scenario_idx::Array{Int64,1})
    global scenarios_paths
    global scenario_list = [initialize_scenario(scenarios_paths[scenario_idx[i]], i) for i in eachindex(scenario_idx)]
end

########################### Simulation functions ###################################
#=
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
=#
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
            #one of the stations is close
            continue
        end

        battery_level_needed = get_battery_level_needed(path) # always 100%
        expected_start_riding_time = current_time + get_walking_time(req.ON, path.origin_station[1])

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
        if pickup_station_id != -1 && drop_off_station_id != -1 && selected_car_id != -1 && parking_place_id != -1
            # we need to memorize the selected requests in the online mode.
            global online_selected_paths[scenario.scenario_id][path_id] = true
            break
        end
    end
    return (pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)

end

function get_walking_time(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_walking_paths)
        merge!(shortest_walking_paths, Dict(id_node1 => dijkstra_shortest_paths(non_directed_manhaten_city_graph, [id_node1])))
    end

    walking_time = shortest_walking_paths[id_node1].dists[id_node2] / walking_speed / 60
    work_with_time_slot && walking_time != Inf && (walking_time = ceil(Int64, walking_time / time_slot_length))

    walking_time
end

#= 
    description : get the driving time between wo nodes
=#
function get_trip_duration(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_graph, [id_node1])))
    end

    driving_duration = shortest_car_paths[id_node1].dists[id_node2] / driving_speed / 1000 * 60 # convert it to minutes 
    work_with_time_slot && driving_duration != Inf && (driving_duration = ceil(Int64, driving_duration / time_slot_length))
    driving_duration
end

function get_trip_distance(id_node1, id_node2)
    if id_node1 ∉ keys(shortest_car_paths)
        merge!(shortest_car_paths, Dict(id_node1 => dijkstra_shortest_paths(manhaten_city_graph, [id_node1])))
    end
    distance = shortest_car_paths[id_node1].dists[id_node2]
    distance
end

function get_trip_battery_consumption(id_node1::Int64, id_node2::Int64, Tcar::car_type)
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
    global potential_locations
    if length(potential_locations) == 0
        potential_locations = collect(filter_vertices(manhaten_city_graph, :type, 2))
    end
    potential_locations
end

function get_battery_level_needed(path)
    return 100
end

#=
    this function calculate the expected battery level at time @time (input)
=#
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

#=
    this function recounts the battery level of the cars in the station 
=#
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

#=
    in this function we don't do any check about the paramaters. the necessary check is done befor

        inputs:
        @pickup_station_id: the station where the car to pick up is parked
        @drop_off_station_id: the station where the car will be parked after performing the trip
        @car_id: the id of the car to be taken
    outputs:
            Nothing
=#
function book_trip(pickup_station_id, drop_off_station_id, car_id,
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

#= a function to set the walking time =#
function set_walking_time(new_walking_time::Int64)
    global maximum_walking_time = new_walking_time
end

#= function to set the cost factor =#
function set_cost_factor(new_cost_factor::Int64)
    global cost_factor = new_cost_factor
end

############################ solution functions ###################################

#= 
    generate a random solution:
        1- decide how much station to open if it is not precised in @open_stations_number
        2- open random stations 
        3- set random initial number of cars for each station
    inputs:
        @open_stations_number (optional): the number of station to open
    outputs:
        sol:: Solution 
=#
function generate_random_solution(; open_stations_number=-1)
    sol = Solution()
    potential_locations = get_potential_locations()
    #decide the number of stations to open if it is not precised by the user
    if open_stations_number == -1
        open_stations_number = rand(1:length(potential_locations))
    end
    #randomly open stations
    sol.open_stations_state[sample(1:length(potential_locations), open_stations_number, replace=false)] .= true

    #set initial car number for each station
    for i in eachindex(sol.open_stations_state)
        max_number_of_charging_points = sol.open_stations_state[i] ? get_prop(manhaten_city_graph, get_potential_locations()[i], :max_number_of_charging_points) : 0
        sol.initial_cars_number[i] = rand(0:max_number_of_charging_points)
    end

    # initialize Randomly the selected paths
    for sc in scenario_list
        curr_sc_selected_paths = [rand(Bool) for _ in 1:nrow(sc.feasible_paths)]
        push!(sol.selected_paths, curr_sc_selected_paths)
    end
    sol
end

#=
    Check whether or not the solution is feasible according to constraints 2, 3, 7 and 8 in Hatice paper
    inputs:
        @sol: the solution
        @all_feasible_paths: the set all feasible paths given by the preprocessing procedure(get_ all_requests_feasible_paths)
    outputs:
        Feasible => Boolean which is true if the solution is feasible, false otherwise.
=#
function is_feasible_solution(sol::Solution)
    #check if the dimention of the solution fields are correctly defined
    if length(get_potential_locations()) != length(sol.open_stations_state) != length(sol.initial_cars_number)
        print_simulation && printstyled(stdout, "the solution fields are not correctly defined\n", color=:light_red)
    end
    #check the initial number of cars (constraint 7 and 8)
    for i in eachindex(sol.open_stations_state)
        if sol.initial_cars_number[i] > (sol.open_stations_state[i] ? get_prop(manhaten_city_graph, potential_locations[i], :max_number_of_charging_points) : 0)
            print_simulation && println("the initial number of cars in the station ", potential_locations[i], " is greater the the total allowed number (or a stations contains cars despite it is closed")
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
                return false
            end
            #check if each station in the selected trips is open (constraint 3)

            #get the oponed stations as dataframes (the number of node in the graph)
            open_stations_df = DataFrame(station_ids=get_potential_locations()[sol.open_stations_state]) # usuful for innerjoin

            #i used inner join to take select the opened stations used in the selected paths and then compared with the number of stations used
            if nrow(innerjoin(selected_paths, open_stations_df, on=:origin_station => :station_ids)) != nrow(selected_paths) ||
               nrow(innerjoin(selected_paths, open_stations_df, on=:destination_station => :station_ids)) != nrow(selected_paths)

                print_simulation && println("there is at least one closed station in the feasible paths")
                return false
            end
        end
        # constraint 4, 5 and 6  are to be checked in the simulation
    end

    return true
end

function get_stored_solution(sol_id=1)
    deserialize(project_path("Data/MIP/solutions/solution_$sol_id.jls"))
end

################################## heuristics ########################################
#=
    this function tries to serve new requests after opening new station for each scenario
    inputs:
        @sol: the solution
        @station_id: the id of the station that we recently oponed
    outputs:
        basicaly it return the new seletced paths ad assigne it as well to the onling_selected_paths so we 
        could get it from there as well
=#
function serve_requests_after_opening_station(sol::Solution, stations_idx::Array{Int64,1})
    @assert !online_request_serving && all(x -> x, sol.open_stations_state[stations_idx]) "we are in wrong mode or the station is closed"
    for sc_id in eachindex(scenario_list)
        sc = scenario_list[sc_id] # handle one scenario a time

        #get the ids of feasible paths that starts (or ends) from (or to) station_id
        potential_feasible_path = get_request_feasible_path(stations_idx, sc)
        #println("list of paths to be examined : $potential_feasible_path")
        if length(potential_feasible_path) > 0
            #index over the paths
            cur_fp_id = 1
            #the current request id allow as to skip paths if we succefly handle it with a previous path 
            cur_req_id::Int64 = sc.feasible_paths.req[potential_feasible_path[cur_fp_id]]
            #println("we are trying to handle request N° $cur_req_id")
            #iterate over the potential paths
            while cur_fp_id <= length(potential_feasible_path)

                #check if we are handling a new request
                if sc.feasible_paths.req[potential_feasible_path[cur_fp_id]] != cur_req_id
                    cur_req_id = sc.feasible_paths.req[potential_feasible_path[cur_fp_id]]
                    #println("**************************************************")
                    #println("we are trying to handle request N° $cur_req_id")
                end

                sol.selected_paths[sc_id][potential_feasible_path[cur_fp_id]] = true
                #println("we are trying serve request N° $cur_req_id with path N° $(potential_feasible_path[cur_fp_id])")

                global failed = false
                E_carsharing_sim(sol, sc_id)

                if !failed
                    #println("Success!! We could serve request N° $cur_req_id with the path N° $(potential_feasible_path[cur_fp_id])")
                    #here  the added path is succefuly excuted without any problem 
                    #se we need to skip all the paths corresponding to the current request
                    while cur_fp_id <= length(potential_feasible_path) && sc.feasible_paths.req[potential_feasible_path[cur_fp_id]] == cur_req_id
                        #println("we are skiping the paths for request $cur_req_id")
                        cur_fp_id += 1
                    end
                else
                    #println("failed!! We could not  serve request N° $cur_req_id with the path N° $(potential_feasible_path[cur_fp_id])")
                    #here the selected path result a non feasable solution
                    #se we need to reset the solution
                    sol.selected_paths[sc_id][potential_feasible_path[cur_fp_id]] = false
                    cur_fp_id += 1
                end
            end
        end
    end

    #save the solution 
    global online_selected_paths = sol.selected_paths
    E_carsharing_sim(sol)
end

#=
    This function retrun the paths where the station, defined by station_id, is origin or destination.
    in addition this function must ensure that the paths are for only requests that not served by another path
    inputs:
        @sol: the solution
        @scenario: the scenario that we want to have feasible paths for it
    outputs:
        Feasible_paths => list of feasible paths 
=#
function get_request_feasible_path(stations_idx::Array{Int64,1}, scenario::Scenario)
    station_nodes_idx = get_potential_locations()[stations_idx]

    fp_starting_from_req = findall(x -> x in station_nodes_idx, scenario.feasible_paths.origin_station)
    fp_ending_from_req = findall(x -> x in station_nodes_idx, scenario.feasible_paths.destination_station)

    union(fp_starting_from_req, fp_ending_from_req)
end

