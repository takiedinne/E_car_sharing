
include("../E_car_sharing.jl")

using Main.E_car_sharing
using Serialization

const E = Main.E_car_sharing
sol1= deserialize("C:/Users/Folio/Desktop/Preparation doctorat ERM/Projects/E_car_sharing/Data/MIP/solutions/solution_1.jls")
sol2= deserialize("C:/Users/Folio/Desktop/Preparation doctorat ERM/Projects/E_car_sharing/Data/MIP/solutions/solution_2.jls")
initialize_scenarios([1,2#= ,3,5 =#])
E_carsharing_sim(sol2, 2)

sc1 = E.scenario_list[1]
sc2 =  E.scenario_list[2]

afp1 = Main.E_car_sharing.get_all_requests_feasible_paths(sc1, E.get_potential_locations(), E.maximum_walking_time)
afp2 = Main.E_car_sharing.get_all_requests_feasible_paths(sc2, E.get_potential_locations(), E.maximum_walking_time)

E.all_feasible_paths_scenario[2] == afp2
sol2.selected_paths

E.get_solutions()