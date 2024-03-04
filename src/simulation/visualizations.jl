export plot_solution, plot_stations, plot_stations_and_requests
function myLayout(g)

    y = [get_prop(g, v, :latitude) for v in vertices(g)]
    x = [get_prop(g, v, :longitude) for v in vertices(g)]

    return Point2f.(x, y)
end

function plot_manhatten()
    

    osm = graph_from_file("Data/other/manhatten_drive.json";
        graph_type=:light, # SimpleDiGraph
        weight_type=:distance
    )

    # use min and max latitude to calculate approximate aspect ratio for map projection
    #autolimitaspect = map_aspect(area.minlat, area.maxlat)

    # plot it
    fig, ax, plot = osmplot(osm#= ; axis=(; autolimitaspect) =#)
    return fig, ax, plot, osm
end

function plot_stations(; g=manhaten_city_driving_graph) #= , ax::Union{Nothing, Axis}=nothing =#
    #firt plot the map
    #f, ax, p, _ = plot_manhatten()
    ax = nothing
    # crate the stations graphs
    stations_graph = MetaGraph()
    locations = get_potential_locations()
    for v in vertices(g)
        if v in locations
            add_vertex!(stations_graph, props(g, v))
        end
    end

    #plot it
    if isnothing(ax)
        f = Figure()
        ax = Axis(f[1, 1])
    end
    
    p = graphplot!(ax, stations_graph,
        layout=myLayout(stations_graph),
        node_marker=[:rect for _ in 1:nv(stations_graph)],
        node_size=[20 for i in 1:nv(stations_graph)],
        node_color=[:black for i in 1:nv(stations_graph)],
        nlabels=["" for i in 1:nv(stations_graph)],
        nlabels_distance=100,
        nlabels_align=[(:center, :bottom) for i in 1:nv(stations_graph)],
        nlabels_fontsize=[32 for i in 1:nv(stations_graph)])


    #interactions
    deregister_interaction!(ax, :rectanglezoom)

    ## node hover functioon
    function node_hover_action(state, idx, event, axis)
        p.node_size[][idx] = state ? 300 : 20
        p.node_size[] = p.node_size[] # trigger observable


        p.nlabels[][idx] = state ? "$idx" : ""
        p.nlabels[] = p.nlabels[]
        #p.nlabels[][idx] = state ? "station" : ""
    end
    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)
    
    #display(GLMakie.Screen(), f)

    return f, p, ax, stations_graph
end


function plot_solution(sol::Solution, optimal_sol::Union{Nothing,Solution}=nothing)
    
    if isnothing(optimal_sol)
        optimal_sol = load_sol("Data/MIP/solutions/E_carsharing_mip_scenario_7_requests_1000_walking_time_5.jls")
    end

    f, p, ax, g = plot_stations()
    # change the colore of the stations
    # green the selected stations and black the register_interaction
    node_color = [sol.open_stations_state[i] ? :red : :black for i in 1:nv(g)]
    if !isnothing(optimal_sol)
        for i in 1:nv(g)
            if sol.open_stations_state[i] && optimal_sol.open_stations_state[i]
                node_color[i] = :green
            elseif sol.open_stations_state[i]
                node_color[i] = :red
            elseif optimal_sol.open_stations_state[i]
                node_color[i] = :blue
            else
                node_color[i] = :black
            end
        end
    end

    p.node_color[] = node_color

    # deregister the previous nhover interaction
    deregister_interaction!(ax, :nhover)
    function node_hover_action(state, idx, event, axis)
        # get the useful information
        label = "id = $idx"
        label = label * "\nCapacity = $(get_prop(g, idx, :max_number_of_charging_points))"
        if sol.open_stations_state[idx]
            current_number_of_cars = sol.initial_cars_number[idx]
            label = label * "\nsol.cars = $current_number_of_cars"
        end
        if !isnothing(optimal_sol) && optimal_sol.open_stations_state[idx]
            optimal_number_of_cars = optimal_sol.initial_cars_number[idx]
            label = label * "\noptimal.cars = $optimal_number_of_cars"
        end
        p.node_size[][idx] = state ? 300 : 20
        p.node_size[] = p.node_size[] # trigger observable

        p.nlabels[][idx] = state ? "$label" : ""
        p.nlabels[] = p.nlabels[]
        #p.nlabels[][idx] = state ? "station" : ""
    end

    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)

    # the legend traitement
    elem_1 = PolyElement(color=:green, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_2 = PolyElement(color=:blue, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_3 = PolyElement(color=:red, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_4 = PolyElement(color=:black, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_5 = PolyElement(color=:white, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_6 = PolyElement(color=:white, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])

    stations_sol = sum(sol.open_stations_state)
    stations_opt_sol = sum(optimal_sol.open_stations_state)
    cars_sol = sum(sol.initial_cars_number)
    cars_opt_sol = sum(optimal_sol.initial_cars_number)


    Legend(f[1, 2],
        [elem_1, elem_2, elem_3, elem_4, elem_5, elem_6],
        ["open in both sols", "open in optimal", "open in sol", "closed in both sols",
            "$stations_sol stations in sol with $cars_sol cars", "$stations_opt_sol stations in optimal with $cars_opt_sol cars"
        ],
        patchsize=(35, 35), rowgap=10)


    display(GLMakie.Screen(), f)
end

"""
    plot_stations_and_requests(requests_df::DataFrame; show_start_node::Bool = true, walking_time::Int = 5)

    plot the requests and the stations on the same graph. this plot can be used to use the ensity of requests near the stations
    and to see the number of requests that can be assigned to each station
    @requests_df: the requests dataframe (it can be a list requests from different scenarios)
    @show_start_node: if true the start node of the requests will be shown, else the destination node will be shown
    @walking_time: the walking time in minutes that the user can walk to reach the station (default = 5)
"""
function plot_stations_and_requests(scenario::Scenario; show_start_node::Bool=true, walking_time::Int=5, sol::Solution=Solution())
    #sol = load_sol("sol.jls")
    #scenario = scenario_list[1]
    requests_df = scenario.request_list[unique(scenario.feasible_paths.req), :]
    
    #f, ax, p, _ = plot_manhatten()
    
    f = Figure()
    ax = Axis(f[1, 1])
    g = manhaten_city_driving_graph

    # make the graph for the requests
    node_to_show = show_start_node ? :ON : :DN
    figure_title = show_start_node ? "requests (origine nodes) and stations" : "requests (destination nodes) and stations"
    #get the number of requests that can be assigned to each station
    locations = get_potential_locations()
    requests_per_station = Int64[]
    for st in locations
        push!(requests_per_station,
            nrow(filter(row -> get_walking_time(row[node_to_show], st) <= walking_time, requests_df)))
    end

    # crate the stations graphs
    requests_graph = MetaGraph()
    for node in requests_df[:, node_to_show]
        add_vertex!(requests_graph, props(g, node))
    end

    for node in locations
        add_vertex!(requests_graph, props(g, node))
    end

    #the edges
    feasible_requests = unique(scenario.feasible_paths.req)
    for trip in eachrow(scenario.feasible_paths)
        req_node = findfirst(feasible_requests .== trip.req)
        station_node = show_start_node ? locations_dict[trip.origin_station] + nrow(requests_df) : locations_dict[trip.destination_station] + nrow(requests_df)
        add_edge!(requests_graph, req_node, station_node)
    end
    
    node_marker_array = vcat(['•' for _ in 1:nrow(requests_df)],
        [:rect for _ in 1:length(locations)])
    node_color_array = vcat([:black for _ in 1:nrow(requests_df)],
        [sol.open_stations_state[i] ? :green : :black for i in 1:length(locations)])

    edge_color_array = [ :black for _ in edges(requests_graph)]
    req_origin = unique(scenario.feasible_paths[:, [:req, :origin_station]])
    selected_trips = scenario.feasible_paths[sol.selected_paths[1], :]
    for i in 1:nrow(selected_trips)
        edge_id = findfirst(req_origin.req .== selected_trips.req[i] .&& req_origin.origin_station .== selected_trips.origin_station[i])
        edge_color_array[edge_id] = :green
    end

    #plot it
    p = graphplot!(ax, requests_graph,
        layout=myLayout(requests_graph),
        node_marker=node_marker_array,
        node_size=[20 for i in 1:nv(requests_graph)],
        node_color=node_color_array,
        nlabels=["" for i in 1:nv(requests_graph)],
        nlabels_distance=100,
        nlabels_align=[(:center, :bottom) for i in 1:nv(requests_graph)],
        nlabels_fontsize=[32 for i in 1:nv(requests_graph)],
        edge_color=edge_color_array
    )

    #interactions
    deregister_interaction!(ax, :rectanglezoom)

    
    ## node hover functioon
    function node_hover_action(state, idx, event, axis)

        #p.node_size[][idx] = (state && idx > nrow(requests_df)) ? 300 : 20

        #p.node_size[][idx] = (state && idx > nrow(requests_df)) ? 300 : 20
        #p.node_size[] = p.node_size[] # trigger observable

        if state
           #nbr = idx > nrow(requests_df) ? requests_per_station[idx-nrow(requests_df)] : 0

            #p.nlabels[][idx] = (idx > nrow(requests_df)) ? "id = $(idx - nrow(requests_df)); ∑req = $nbr" : "$idx"
            p.nlabels[][idx] = (idx > nrow(requests_df)) ? "id = $(idx - nrow(requests_df))" : "$idx"
            #p.node_marker[][idx] = (idx > nrow(requests_df)) ? :circle : '+'
        else
            p.nlabels[][idx] = ""
            p.node_marker[][idx] = (idx <= nrow(requests_df)) ? '•' : :rect
        end
        p.nlabels[] = p.nlabels[]
        p.node_marker[] = p.node_marker[] # trigger observable

        p.node_size[] = p.node_size[] # trigger observable
        p.nlabels[] = p.nlabels[] # trigger observable

    end
    ax.title = figure_title
    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)


    # The legend traitement
    elem_1 = PolyElement(color=:green, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_2 = PolyElement(color=:black, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    
    Legend(f[1, 2],
        [elem_1, elem_2],
        ["station is open", "station is close"],
        patchsize=(35, 35), rowgap=10)

    display(GLMakie.Screen(), f)

end
