using E_car_sharing
const e = E_car_sharing
using BenchmarkTools

initialize_scenarios([1])

opt_sol_path = "Data/other/scenario_1_opt_sol.jls"
file_path = "/Users/taki/Desktop/Preparation doctorat ERM/Projects/GIHH_V2.0/results/solve_multiple_scenario/1_scenarios/9_trial/GIHH__HCnt.csv"
opt_sol = e.load_sol(opt_sol_path)
E_carsharing_sim(opt_sol)

e.plot_best_fitness_tracking(file_path, optimal_fitness = -26817.4)

counter = 0
fit_sum = 0.
for i in 2:10
    file_path = "/Users/taki/Desktop/Preparation doctorat ERM/Projects/GIHH_V2.0/results/solve_multiple_scenario/1_scenarios/$(i-1)_trial/GIHH__HCnt.csv"
    counter += e.get_nbr_improvement(file_path)
    fit_sum += e.get_best_fitness(file_path)
end
counter/9
fit_sum/9
