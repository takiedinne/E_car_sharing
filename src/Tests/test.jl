using E_car_sharing
const e = E_car_sharing


initialize_scenarios([1])
opt_sol = e.load_sol("Data/other/scenario_1_opt_sol.jls")
sol = e.load_sol("Data/other/GIHH_sol.jls")
optimal_fitness = E_carsharing_sim(opt_sol)
e.plot_solution(sol)
current_fitness = E_carsharing_sim(sol)

# analysis:
file_path = "/Users/taki/Desktop/Preparation doctorat ERM/Projects/GIHH_V2.0/results/solve_single_scenario/1/GIHH__HCnt.csv"
plot_best_fitness_tracking(file_path, optimal_fitness = optimal_fitness)
get_nbr_improvement(file_path)
read_heuristics_performance(file_path)


station_id = 8
station_id1 = 50
station_node_id = e.get_potential_locations()[station_id]
station_node_id1 = e.get_potential_locations()[station_id1]

ofp = e.scenario_list[1].feasible_paths[opt_sol.selected_paths[1], :]
fp  = e.scenario_list[1].feasible_paths[sol.selected_paths[1], :]
afp = e.scenario_list[1].feasible_paths
trips = e.get_trips_station(station_id, ofp)
trips

#check if req 647 is served or not
findall(fp.req .== 647)

old_fit = E_carsharing_sim(sol)
sol.open_stations_state[47] = true
sol.initial_cars_number[47]

serve_requests_after_opening_station(sol, [station_id])

#get the new fitness after serving the requests
e.get_trips_station(47, fp)
e.get_trips_station(8, afp)

reqs =  e.get_trips_station(47, fp).req
for req in reqs
    i =findfirst(ofp.req .== req)
    println("req  $req:  ")
    println("$(ofp.origin_station[i]) => $(ofp.destination_station[i])")
    println("$(fp.origin_station[i]) => $(fp.destination_station[i])")
    println("##########################")
end
