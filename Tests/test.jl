include("../E_carsharing_sim.jl")

using BenchmarkTools
using Makie
using GLMakie
using GraphMakie
using Colors

################################################################################
########################### draw the graph #####################################
#= manhaten_city_graph = create_graph_from_XML("Tests/test_graph.xml", save_file ="Data/test_graph.mg")

global all_request_df = CSV.read(all_request_details_path, DataFrame)
scenario = scenario_as_dataframe(scenario_path)

draw_graph_and_scenario(manhaten_city_graph, scenario) =#

#####################################################################################
#####################################################################################
######################### test the simulation #######################################

include("../E_carsharing_sim.jl")
# create the solution
sol = generate_random_solution( open_stations_number = 30)

scenario = initilaize_scenario(scenario_path, sol)

using BenchmarkTools
@benchmark E_carsharing_sim(sol, scenario)


global time1 = 0
global time2 = 0
global time3 = 0
global time4 = 0
global time5 = 0




 
a1 = 100 * time1 / total_time
a2 = 100 * time2 / total_time
a3 = 100 * time3 / total_time
a4 = 100 * time4 / total_time
a5 = 100 * time5 / total_time

println("time $a1 %)")
println("time $a2 %)")
println("time $a3 %)")
println("time $a4 %)")
println("time $a5 %)")

get_each_requests_feasible_paths(scenario, sol)

a = get_all_requests_feasible_paths(scenario, sol)

@benchmark b = [filter( x-> x.req == i, a) for i in scenario.reqId]

filter( x-> x.req == 16738, a)

sol.selected_paths = c

scenario