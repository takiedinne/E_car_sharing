include("../E_car_sharing.jl")
using Main.E_car_sharing
using BenchmarkTools
using DataFrames
using Profile

const e = E_car_sharing

initialize_scenarios(collect(1:1))
n = length(e.get_potential_locations())
open_stations_state = falses(n)
initial_cars_number = zeros(n)
selected_paths = [falses(i) for i in [nrow(e.scenario_list[j].feasible_paths) for j in eachindex(e.scenario_list)]]

# treatement for the solution
open_stations_state[60] =  true
initial_cars_number[60] = 1

sol = Solution(open_stations_state, initial_cars_number, selected_paths)

sol.open_stations_state[9] = true
e.serve_requests_after_opening_station(sol, [9])

a = trues(10)
stations_idx = [1, 2]

all(x->x, [true, false])

scenario = e.scenario_list[1]

station_nodes_idx = e.get_potential_locations()[stations_idx]

fp_starting_from_req = findall(x->x in station_nodes_idx, scenario.feasible_paths.origin_station)

nbr = 0
for i in eachindex(fp_starting_from_req)
    path_id = fp_starting_from_req[i]
    if (scenario.feasible_paths.origin_station[path_id] in  station_nodes_idx)
        println("error with path number $path_id")
        nbr += 1
    end
end
nbr == length(fp_starting_from_req)

e.manhaten_city_graph

