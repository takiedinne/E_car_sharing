include("../E_car_sharing.jl")
using Main.E_car_sharing
using Serialization
using DataFrames, Graphs, MetaGraphs
using CSV

const e = Main.E_car_sharing

e.preprocessing_experiment()
e.preprocessing_experiment2019()
e.generate_feasible_paths_for_generated_scenarios()
#= 
    #scenario scenarios_objects
    sc = e.scenario_list[1]
    #station id
    station_id = 60
    station_node_id = e.potential_locations[station_id]

    #get all the selected_requests and teir feasible paths
    selected_paths = sc.feasible_paths[sol.selected_paths[1], :]
    served_requests = sc.request_list[selected_paths.req, :]

    # get the feasible paths for the station
    selected_paths_starting_from_station = filter(x -> x.origin_station == station_node_id, selected_paths)
    selected_paths_ending_from_station = filter(x -> x.destination_station == station_node_id, selected_paths)
    interesting_selected_paths = sort(vcat(selected_paths_starting_from_station, selected_paths_ending_from_station), [:start_driving_time])
    interesting_selected_requests = filter(x-> x.reqId in interesting_selected_paths.req, served_requests)
    # initial cars number and capacity
    initial_cars_number = sol.initial_cars_number[station_id]
    capacity = get_prop(e.manhaten_city_graph, e.potential_locations[station_id], :max_number_of_charging_points)
 =#
# test the mixed integer programming


sc = initialize_scenario(e.scenarios_paths[1], 1, check_file=false)
length(sc.request_list.reqId)

push!(e.scenario_list, sc)


empty!(e.scenario_list)
e.scenario_list