include("../E_car_sharing.jl")
using Main.E_car_sharing
const e = E_car_sharing

e.preprocessing_experiment()
e.preprocessing_experiment2019()

e.generate_feasible_paths_for_Ci()

