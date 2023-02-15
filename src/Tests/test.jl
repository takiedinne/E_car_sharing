include("../E_car_sharing.jl")

using Main.E_car_sharing
using BenchmarkTools
using DataFrames

initialize_scenarios(collect(1:1))
#sol = Solution(Bool[true, true, true], Integer[2, 0, 0], [[true, false, true, false, false]])
sol = generate_random_solution()
set_online_mode(true)
@btime [E_carsharing_sim(sol) for i in 1:1806]

E_car_sharing.online_selected_paths

const e = E_car_sharing
e.get_trip_duration(2, 5)
100 - e.get_trip_battery_consumption(1, 5, e.Smart_ED)

e.stations[3].cars

e.refrech_battery_levels(3, 10)