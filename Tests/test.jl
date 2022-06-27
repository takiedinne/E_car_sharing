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
using BenchmarkTools
# create the solution
sol = generate_random_solution( open_stations_number = 3)
sol.open_stations_ids = [7, 19, 24]
sol.initial_cars_number = [2, 0, 0]
sol.selected_paths = [true, true , true]

scenario = initilaize_scenario(scenario_path, sol)

E_carsharing_sim(sol, scenario)




prgrm = Model()
set_optimizer(prgrm, GLPK.Optimizer)

@variable(prgrm, 0<=x)
@variable(prgrm, 0<=y)

@constraint(prgrm, x <= 10)
@constraint(prgrm, y <= 7.5)

@objective(prgrm, Max, x+2y)

using MathOptInterface

write_to_file(prgrm, "test.mof.json")

print(prgrm)


model = read_from_file("test.mof.json")

println(model)