using E_car_sharing
const e = E_car_sharing

opt_sol = e.load_sol("Data/other/scenario_1_opt_sol.jls")
sol = e.load_sol("Data/other/GIHH_sol.jls")

e.plot_solution(sol, optimal_sol=opt_sol)

initialize_scenarios([1])
sol_fp_details = e.scenario_list[1].feasible_paths[sol.selected_paths[1], :]

#get trips starts or ends from a station number 63
station = 16
sol = opt_sol
sol.open_stations_state[station] 
station_node_id = e.get_potential_locations()[station]

#filter the trips that starts with or end in station 63
trips = filter(row -> station_node_id in [row.origin_station, row.destination_station], sol_fp_details)

