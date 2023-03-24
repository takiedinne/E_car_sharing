include("../E_car_sharing.jl")
using Main.E_car_sharing
using DataFrames
const e = E_car_sharing

e.preprocessing_experiment2019()

using SparseArrayKit

a = SparseArray{Bool}(undef, 5120000, 85, 300)
