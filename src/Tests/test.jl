using E_car_sharing
const e = E_car_sharing

using BenchmarkTools

e.initialize_scenarios(collect(1:200));

sol = generate_random_solution();
E_carsharing_sim(sol)