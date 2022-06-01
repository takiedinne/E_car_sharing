include("../E_carsharing_sim.jl")

using Makie
using GLMakie
using GraphMakie
using Colors


################################################################################
########################### draw the graph #####################################

#cols = distinguishable_colors(nv(manhaten_city_graph), [RGB(1,1,1), RGB(0,0,0)], dropseed=true);
cols = []
node_size = Array{Real, 1}()
node_marker = [:circle for i in 0:nv(manhaten_city_graph)-1]

for i in 1:nv(manhaten_city_graph)
    if get_prop(manhaten_city_graph, i, :type) == 2
        push!(node_size, 40)
        push!(cols,  :blue )
        node_marker[i] = :rect
    else
        push!(node_size, 20)
        push!(cols, :black ) 
    end
end
req_colors= [:red, :yellow, :gray, :green]
for req in eachrow(scenario)
    if node_marker[req.ON ] == :circle
        node_marker[req.ON] = :star5
        node_size[req.ON] = 40
        cols[req.ON] = req_colors[req.reqId]
    end
    if node_marker[req.DN ] == :circle
        node_marker[req.DN ] = :diamond
        node_size[req.DN] = 40
        cols[req.DN] = req_colors[req.reqId]
    end
end

f, ax, p = graphplot(manhaten_city_graph, edge_width=[3 for i in 1:ne(manhaten_city_graph)],
                     node_size=node_size ,
                     node_color=cols,
                     nlabels = [string(i) for i in 0:nv(manhaten_city_graph)-1],
                     node_marker = node_marker)

deregister_interaction!(ax, :rectanglezoom)

register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p))



#####################################################################################
######################### test the simulation #######################################

# create the solution
sol = generate_random_solution( open_stations_number = 3)
# select the paths
sol.selected_paths = sol.selected_paths = [false, true, true, true, true]

initialize_sim(sol, scenario_path)

E_carsharing_sim(sol)

scenario