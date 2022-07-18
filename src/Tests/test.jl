include("../E_carsharing.jl")

using Main.E_carsharing


################################################################################
########################### draw the graph #####################################
#= manhaten_city_graph = create_graph_from_XML("Tests/test_graph.xml", save_file ="Data/test_graph.mg")

global all_request_df = CSV.read(all_request_details_path, DataFrame)
scenario = scenario_as_dataframe(scenario_path)

draw_graph_and_scenario(manhaten_city_graph, scenario) =#

#####################################################################################
#####################################################################################
######################### test the simulation #######################################

# create the solution
sol = generate_random_solution(open_stations_number=35)
sol.open_stations_ids = [7, 19, 24]
sol.initial_cars_number = [2, 0, 0]
sol.selected_paths = [true, true, true]

scenario = initilaize_scenario(E_carsharing.scenario_path, sol)

E_carsharing_sim(sol, scenario)

####################################################################################
####################################################################################
##############################validate the simulation###############################

#= 
station_node = 906
station_id = findfirst(sol.open_stations_ids .== station_node)
initial_car_number = sol.initial_cars_number[station_id]
station_capacity = get_prop(E_carsharing.manhaten_city_graph, station_node, :max_number_of_charging_points)

a = filter(x -> x.ST <= 65 && nrow(x.feasible_paths) > 0 && (x.feasible_paths.origin_station[1] == station_node || x.feasible_paths.destination_station[1] == station_node), scenario)
b = vcat(a.feasible_paths...)
draw_station_car_nbr(station_node, initial_car_number, b)

ids = findall(scenario.reqId .== 208)
scenario[ids, :].feasible_paths

findall(E_carsharing.all_feasible_paths.req .== 208)

# get the solution mixed integer program
revenue_sim = E_carsharing.revenues

stations_cost_sim = sum([station.charging_station_base_cost +
                     station.max_number_of_charging_points * station.charging_point_cost_fast
                     for station in E_carsharing.stations])

E_carsharing.stations_costs_mixed_integer < stations_cost_sim

E_carsharing.get_solutions() =#


include("../E_carsharing.jl")

using Main.E_carsharing
using Serialization

using DataFrames
using MetaGraphs
using Plots

function draw_station_car_nbr(station_node, initial_car_number, b)
    X = DataFrame(x=[0], y=[initial_car_number])
    y = initial_car_number

    df = DataFrame(time=[], signe=[])
    for i in 1:nrow(b)
        if b.origin_station[i] == station_node
            push!(df, (b.start_driving_time[i], -1))
        elseif b.destination_station[i] == station_node
            push!(df, (b.arriving_time[i], 1))
        end
    end

    sort!(df, [:time])
    for i in 1:nrow(df)
        y += df.signe[i]
        if df.time[i] in X.x
            X.y[findfirst(X.x .== df.time[i])] = y
        else
            push!(X, (df.time[i], y))
        end
    end
    plot(X.x, X.y, linetype=:steppre)

end

obj_values_mixed_integer = deserialize("Data/MIP/solutions/Objective_values.jls")
scenario_path_list = ["Data/Scenarios_1000_greaterthan2/Output1000_1.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_2.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_3.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_4.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_5.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_6.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_7.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_8.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_9.txt",
    "Data/Scenarios_1000_greaterthan2/Output1000_10.txt"]

# read solution   
sol_path_list = ["Data/MIP/solutions/solution_1.jls",
    "Data/MIP/solutions/solution_2.jls",
    "Data/MIP/solutions/solution_3.jls",
    "Data/MIP/solutions/solution_4.jls",
    "Data/MIP/solutions/solution_5.jls",
    "Data/MIP/solutions/solution_6.jls",
    "Data/MIP/solutions/solution_7.jls",
    "Data/MIP/solutions/solution_8.jls",
    "Data/MIP/solutions/solution_9.jls",
    "Data/MIP/solutions/solution_10.jls"]

obj_values_sim = []
for i in 1:length(sol_path_list)
    println("we are with solution number $i ...")
    sol = deserialize(sol_path_list[i])
    # read the scenario
    scenario = initilaize_scenario(scenario_path_list[i], sol)

    sim_obj_value = E_carsharing_sim(sol, scenario)
    push!(obj_values_sim, sim_obj_value)
end
round.(obj_values_mixed_integer, digits = 2) .== obj_values_sim

sol = deserialize(sol_path_list[6])
# read the scenario
scenario = initilaize_scenario(scenario_path_list[6], sol)

sim_obj_value = E_carsharing_sim(sol, scenario)
obj_values_sim[6]
obj_values_mixed_integer[6]


sol.open_stations_ids[39]
station_node = 2681
station_id = findfirst(sol.open_stations_ids .== station_node)
initial_car_number = sol.initial_cars_number[station_id]
station_capacity = get_prop(E_carsharing.manhaten_city_graph, station_node, :max_number_of_charging_points)

a = filter(x -> x.ST <= 230 && nrow(x.feasible_paths) > 0 && (x.feasible_paths.origin_station[1] == station_node || x.feasible_paths.destination_station[1] == station_node), scenario)
b = vcat(a.feasible_paths...)
draw_station_car_nbr(station_node, initial_car_number, b)

ids = findall(scenario.reqId .== 208)
scenario[ids, :].feasible_paths

findall(E_carsharing.all_feasible_paths.req .== 208)


#create solution
filter(x -> x.ST == 230 , scenario)


sim_obj_value = E_carsharing_sim(sol, scenario)