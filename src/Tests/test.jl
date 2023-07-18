include("../E_car_sharing.jl")
using Main.E_car_sharing

const e = E_car_sharing

sol = e.load_sol("Data/other/GIHH_sol.jls")
e.plot_solution(sol)