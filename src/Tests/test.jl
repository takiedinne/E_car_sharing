include("../E_car_sharing.jl")

using Main.E_car_sharing
using BenchmarkTools
using DataFrames

initialize_scenarios(collect(1:1))
sol = Solution(Bool[true, true, true], Integer[2, 0, 0], [[true, true, true, false, true]])
#sol = generate_random_solution()
E_carsharing_sim(sol)

a= 10 ^16

typeof(a)

