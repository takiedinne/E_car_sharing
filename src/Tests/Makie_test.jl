include("../E_car_sharing.jl")
using Main.E_car_sharing

using Graphs, MetaGraphs
using GLMakie, GraphMakie

const e = E_car_sharing

g = e.manhaten_city_driving_graph

# create the graph to be ploted 
function myLayout(g)
    x = [get_prop(g, v, :latitude) for v in vertices(g)]
    y = [get_prop(g, v, :longitude) for v in vertices(g)]
    return Point2f.(x, y)
end

# crate the stations graphs
stations_graph = MetaGraph()
locations = e.get_potential_locations() 
for v in vertices(g)
    if v in locations
        add_vertex!(stations_graph, props(g, v))
    end
end

#plot it
nodes_marker = [:rect for _ in 1:nv(stations_graph)]

f, ax, p = graphplot(stations_graph,
        layout = myLayout(stations_graph), 
        node_marker = nodes_marker, 
        node_size = [10 for i in 1:nv(stations_graph)],
        node_color = [:black for i in 1:nv(stations_graph)],
        nlabels = ["" for i in 1:nv(stations_graph)]);

#interactions
deregister_interaction!(ax, :rectanglezoom)

## node hover functioon
function node_hover_action(state, idx, event, axis)
    p.node_size[][idx] = state ? 100 : 10
    p.node_size[] = p.node_size[] # trigger observable

    p.node_color[][idx] = state ? :red : :black
    p.node_color[] = p.node_color[]

    p.nlabels[][idx] = state ? "$idx" : ""
    p.nlabels[] = p.nlabels[]
    #p.nlabels[][idx] = state ? "station" : ""
end
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)


display(GLMakie.Screen(), f)