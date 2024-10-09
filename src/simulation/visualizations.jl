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
    fig, ax, plot = osmplot(osm) #= ; axis=(; autolimitaspect) =#
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

    #= if isnothing(optimal_sol)
        optimal_sol = load_sol("Data/MIP/solutions/E_carsharing_mip_scenario_7_requests_1000_walking_time_5.jls")
    end =#

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
    stations_opt_sol = isnothing(optimal_sol) ? 0 : sum(optimal_sol.open_stations_state)
    cars_sol = sum(sol.initial_cars_number)
    cars_opt_sol = isnothing(optimal_sol) ? 0 : sum(optimal_sol.initial_cars_number)


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
function plot_stations_and_requests(scenario::Scenario;
    show_start_node::Bool=true, walking_time::Int=5,
    sol::Solution=Solution(),
    opt_sol::Union{Nothing,Solution} = nothing)

   #=  sol = load_sol("sol.jls")
    scenario = scenario_list[1]
    show_start_node, walking_time = true, 10
    opt_sol = load_sol("Data/MIP/solutions/E_carsharing_mip_scenario_1_requests_1000_walking_time_$(walking_time).jls")
 =#
    requests_df = scenario.request_list[unique(scenario.feasible_paths.req), :]

    f = Figure()
    ax = Axis(f[1, 1])
    g = manhaten_city_driving_graph

    # make the graph for the requests
    node_to_show = show_start_node ? :ON : :DN
    figure_title = show_start_node ? "requests (origine nodes) and stations" : "requests (destination nodes) and stations"
    #get the number of requests that can be assigned to each station
    locations = get_potential_locations()

    # crate the stations graphs
    requests_graph = MetaGraph()
    for req in eachrow(requests_df)
        node = req[node_to_show]
        props_dict = props(g, node)
        props_dict[:node_type] = :request
        props_dict[:request_id] = req.reqId
        add_vertex!(requests_graph, props_dict)
    end

    for node in locations
        add_vertex!(requests_graph, props(g, node))
    end


    stations_color = [sol.open_stations_state[i] ? :red : :black for i in eachindex(sol.open_stations_state)]
    if !isnothing(opt_sol)
        for i in eachindex(sol.open_stations_state)
            if sol.open_stations_state[i] && opt_sol.open_stations_state[i]
                stations_color[i] = :green
            elseif sol.open_stations_state[i]
                stations_color[i] = :red
            elseif opt_sol.open_stations_state[i]
                stations_color[i] = :blue
            else
                stations_color[i] = :black
            end
        end
    end
    node_marker_array = vcat(['•' for _ in 1:nrow(requests_df)],
        [:rect for _ in 1:length(locations)])
    node_color_array = vcat([:black for _ in 1:nrow(requests_df)], stations_color)

    #plot it
    p = graphplot!(ax, requests_graph,
        layout=myLayout(requests_graph),
        node_marker=node_marker_array,
        node_size=[20 for i in 1:nv(requests_graph)],
        node_color=node_color_array,
        nlabels_distance=100,
        nlabels_align=[(:center, :bottom) for i in 1:nv(requests_graph)],
        nlabels_fontsize=[32 for i in 1:nv(requests_graph)]
    )

    #interactions
    deregister_interaction!(ax, :rectanglezoom)

    ## node hover functioon
    function node_click_action(idx, event, axis)

        if idx > nrow(requests_df)
            return
        end

        #get all the feasible paths of the request
        reqId = get_prop(requests_graph, idx, :request_id)
        @info "node $reqId"
        trips_ids = request_feasible_trips_ids[scenario.scenario_id][reqId]
        trips = scenario.feasible_paths[trips_ids, :]

        #delete all the edges
        [rem_edge!(requests_graph, e) for e in edges(requests_graph)]

        # add the edges to the graph
        for trip in eachrow(trips)
            #trip = trips[1, :]
            o_node = locations_dict[trip.origin_station] + nrow(requests_df)
            d_node = locations_dict[trip.destination_station] + nrow(requests_df)
            add_edge!(requests_graph, o_node, d_node)
        end

        # color the edges green for the optimal and the red for selected in sol
        selected_trip_id = findfirst(sol.selected_paths[scenario.scenario_id][trips_ids])
        opt_trip_id = findfirst(opt_sol.selected_paths[scenario.scenario_id][trips_ids])

        edges_color = [:black for e in edges(requests_graph)]
        if !isnothing(selected_trip_id)
            edges_color[selected_trip_id] = :red
        end
        if !isnothing(opt_trip_id)
            edges_color[opt_trip_id] = :blue
        end
        if !isnothing(selected_trip_id) && !isnothing(opt_trip_id)
            if selected_trip_id == opt_trip_id
                edges_color[selected_trip_id] = :green
            end
        end
        @info edges_color
        #reploting the graph
        empty!(axis)
        graphplot!(axis, requests_graph,
            layout=myLayout(requests_graph),
            node_marker=node_marker_array,
            node_size=[20 for i in 1:nv(requests_graph)],
            node_color=node_color_array,
            nlabels_distance=100,
            nlabels_align=[(:center, :bottom) for i in 1:nv(requests_graph)],
            nlabels_fontsize=[32 for i in 1:nv(requests_graph)],
            edge_color=edges_color
        )
        deregister_interaction!(axis, :nodeclick)
        nclick = NodeClickHandler(node_click_action)
        register_interaction!(axis, :nodeclick, nclick)

    end

    nclick = NodeClickHandler(node_click_action)
    register_interaction!(ax, :nodeclick, nclick)


    # The legend traitement
    elem_1 = PolyElement(color=:green, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_2 = PolyElement(color=:blue, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_3 = PolyElement(color=:red, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_4 = PolyElement(color=:black, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_5 = PolyElement(color=:white, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_6 = PolyElement(color=:white, points=Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])

    stations_sol = count(sol.open_stations_state)
    stations_opt_sol = count(opt_sol.open_stations_state)
    cars_sol = sum(sol.initial_cars_number)
    cars_opt_sol = sum(opt_sol.initial_cars_number)
    Legend(f[1, 2],
        [elem_1, elem_2, elem_3, elem_4, elem_5, elem_6],
        ["open in both sols", "open in optimal", "open in sol", "closed in both sols",
            "$stations_sol stations in sol with $cars_sol cars", "$stations_opt_sol stations in optimal with $cars_opt_sol cars"
        ],
        patchsize=(35, 35), rowgap=10)

    display(GLMakie.Screen(), f)

end

