using GLMakie
using GraphMakie
using Colors
using Graphs, MetaGraphs
using DataFrames

small_graph_file = "Data/test_graph.mg"
small_graph = loadgraph(small_graph_file, MGFormat())
cols_small = distinguishable_colors(nv(small_graph), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true);

f, ax, p = graphplot(small_graph, edge_width=[3 for i in 1:ne(small_graph)],
    node_size=[20 for i in 1:nv(small_graph)],
    node_color=cols_small,
    nlabels=[string(i+1) for i in 0:nv(small_graph)-1]
    )

deregister_interaction!(ax, :rectanglezoom)

register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p))

props(small_graph, Edge(11, 12))

for a in edges(small_graph)
    println(a)
end

#= Manhatten_network_Metagraph_file = "Data/manhatten_graph.mg"
huge_graph = loadgraph(Manhatten_network_Metagraph_file, MGFormat())
cols_huge = distinguishable_colors(nv(huge_graph), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)

f, ax, p = graphplot(huge_graph, edge_width=[3 for i in 1:ne(g)],
    node_size=[20 for i in 1:nv(g)],
    node_color=cols,
    nlabels=[string(i) for i in 0:nv(g)-1])

deregister_interaction!(ax, :rectanglezoom)

register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p)) =#

