using GeometryBasics

function myLayout(g)
    x = [get_prop(g, v, :latitude) for v in vertices(g)]
    y = [get_prop(g, v, :longitude) for v in vertices(g)]
    return Point2f.(x, y)
end

function plot_stations(;g=manhaten_city_driving_graph)
    
    # crate the stations graphs
    stations_graph = MetaGraph()
    locations = get_potential_locations() 
    for v in vertices(g)
        if v in locations
            add_vertex!(stations_graph, props(g, v))
        end
    end

    #plot it
    f, ax, p = graphplot(stations_graph,
        layout = myLayout(stations_graph), 
        node_marker = [:rect for _ in 1:nv(stations_graph)], 
        node_size = [40 for i in 1:nv(stations_graph)],
        node_color = [:black for i in 1:nv(stations_graph)],
        nlabels = ["" for i in 1:nv(stations_graph)],
        nlabels_distance = 100,
        nlabels_align= [(:center, :bottom) for i in 1:nv(stations_graph)],
        nlabels_fontsize =[32 for i in 1:nv(stations_graph)]);

    #interactions
    deregister_interaction!(ax, :rectanglezoom)

    ## node hover functioon
    function node_hover_action(state, idx, event, axis)
        p.node_size[][idx] = state ? 300 : 40
        p.node_size[] = p.node_size[] # trigger observable

       
        p.nlabels[][idx] = state ? "$idx" : ""
        p.nlabels[] = p.nlabels[]
        #p.nlabels[][idx] = state ? "station" : ""
    end
    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)


    display(GLMakie.Screen(), f)

    return f, ax, p, stations_graph
end

function plot_solution(sol::Solution; optimal_sol::Union{Solution, Nothing}=nothing)
    
    optimal_sol = load_sol("Data/other/scenario_1_opt_sol.jls")
    
    f, ax, p, g = plot_stations()
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
        p.node_size[][idx] = state ? 300 : 40
        p.node_size[] = p.node_size[] # trigger observable

        p.nlabels[][idx] = state ? "$label" : ""
        p.nlabels[] = p.nlabels[]
        #p.nlabels[][idx] = state ? "station" : ""
    end

    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)

    # the legend traitement
    elem_1 = PolyElement(color = :green, points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_2 = PolyElement(color = :blue, points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_3 = PolyElement(color = :red, points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_4 = PolyElement(color = :black, points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_5 = PolyElement(color = :white, points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    elem_6 = PolyElement(color = :white, points = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)])
    
    stations_sol = sum(sol.open_stations_state)
    stations_opt_sol = sum(optimal_sol.open_stations_state)
    cars_sol = sum(sol.initial_cars_number)
    cars_opt_sol = sum(optimal_sol.initial_cars_number)

   
    Legend(f[1, 2],
        [elem_1, elem_2, elem_3, elem_4, elem_5, elem_6],
        [ "open in both sols", "open in optimal", "open in sol", "closed in both sols",
         "$stations_sol stations in sol with $cars_sol cars", "$stations_opt_sol stations in optimal with $cars_opt_sol cars"
         ],
        patchsize = (35, 35), rowgap = 10)

   
end
