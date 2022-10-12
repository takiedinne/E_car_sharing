include("../E_car_sharing.jl")
using Main.E_car_sharing
using Serialization
using BenchmarkTools

const e = Main.E_car_sharing
#read solution 
sol1 = deserialize("Data/MIP/solutions/solution_1.jls")
# read the scenario
e.initialize_scenarios([1])

# run the simulations 
f_obj = e.E_carsharing_sim(sol1, 1)

