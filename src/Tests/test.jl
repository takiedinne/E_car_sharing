
include("../E_car_sharing.jl")

using Main.E_car_sharing
using Serialization
using DataFrames

const E = Main.E_car_sharing
sol1= deserialize("Data/MIP/solutions/solution_1.jls")

initialize_scenarios([1,#= 2,3,5 =#])

sc1 = E.scenario_list[1]

E.all_feasible_paths = E.all_feasible_paths_scenario[1]
E.initialize_sim(sol1, sc1)

station = E.stations[4]

station.cars

total_cars_cost = 0
for i in 1:length(E.stations)
    if nrow(E.stations[i].cars) > 0
        total_cars_cost += sum([E.vehicle_specific_values[E.stations[i].cars.car_type[j]][:car_cost] for j in 1:nrow(E.stations[i].cars)])
    end
end
i=1
sum([sum([E.vehicle_specific_values[E.stations[i].cars.car_type[j]][:car_cost] for j in 1:nrow(E.stations[i].cars)]) for i in eachindex(E.stations)])

E.stations[4].cars

using BenchmarkTools
@benchmark sum([E.vehicle_specific_values[i][:car_cost] for i in vcat([E.stations[i].cars.car_type for i in eachindex(E.stations)]...)])


@benchmark begin
    total_cars_cost = 0
    for i in 1:length(E.stations)
        if nrow(E.stations[i].cars) > 0
            total_cars_cost += sum([E.vehicle_specific_values[E.stations[i].cars.car_type[j]][:car_cost] for j in 1:nrow(E.stations[i].cars)])
        end
    end
end
