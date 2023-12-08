using E_car_sharing
const e = E_car_sharing

using BenchmarkTools

e.initialize_scenarios(collect(1:200));

sol = e.load_sol("sol.jls")
E_carsharing_sim(sol)
curr_fit, curr_sol = e.addCarsLS(sol, local_search_depth = 4)

sol = e.load_sol("sol_error.jls")

E_carsharing_sim(sol)
