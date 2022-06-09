include("../E_carsharing_sim.jl")

using Makie
using GLMakie
using GraphMakie
using Colors

################################################################################
########################### draw the graph #####################################
manhaten_city_graph = create_graph_from_XML("Tests/test_graph.xml", save_file ="Data/test_graph.mg")

global all_request_df = CSV.read(all_request_details_path, DataFrame)
scenario = scenario_as_dataframe(scenario_path)

draw_graph_and_scenario(manhaten_city_graph, scenario)

##################################################################################"
# test the number of paths in the preprocessing procedure

props(manhaten_city_graph, 19)

scenario = scenario_as_dataframe(scenario_path)

a = all_requests_feasible_paths(scenario, sol)


global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph

origine_node = 20
destination_node = 8
origin_station = 19
destination_station = 7

global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph


trips_total_duration = get_walking_time(origine_node,origin_station) + get_trip_duration(origin_station, destination_station) + get_walking_time(destination_station,destination_node)

threshold = 1.1 * get_trip_duration(origine_node, destination_node)

get_trip_duration(13, 18)
get_walking_time(1,2)
#####################################################################################
#####################################################################################
######################### test the simulation #######################################
include("../E_carsharing_sim.jl")
# create the solution
sol = generate_random_solution( open_stations_number = 30)

# create the solution manually
#= 
open_stations = [19, 41, 60]

initial_car_number = [2, 2, 2]
sol = Solution(open_stations, initial_car_number, [])

# select the paths
sol.selected_paths = [true, true, true,#= true, true,  false, true, true =#]

initialize_sim(sol, scenario_path)
 =#
using BenchmarkTools

@time E_carsharing_sim(sol)

