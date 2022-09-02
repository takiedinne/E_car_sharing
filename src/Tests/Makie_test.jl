using Makie
using GLMakie
using GraphMakie
using Colors


using Graphs, MetaGraphs


Manhatten_network_Metagraph_file = "manhatten_graph.mg"
small_graph_file = "test_graph.mg"

huge_graph = loadgraph(Manhatten_network_Metagraph_file, MGFormat())
small_graph = loadgraph(small_graph_file, MGFormat())

cols_huge = distinguishable_colors(nv(huge_graph), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
cols_small = distinguishable_colors(nv(small_graph), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

f, ax, p = graphplot(small_graph, edge_width=[3 for i in 1:ne(g)],
                     node_size=[20 for i in 1:nv(g)],
                     node_color=cols,
                     nlabels = [string(i) for i in 0:nv(g)-1])

deregister_interaction!(ax, :rectanglezoom)

register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p))

