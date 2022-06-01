#=
    this file regroups all the general parameters for the simulation
=#

# Different File paths
#global all_request_details_path = "Data/trips_ML_daily_greaterthan2_with_revenue&minutes_backup.txt"
global all_request_details_path = "Data/trips_ML_daily_greaterthan2_with_revenue&minutes.txt"

global Manhatten_network_details_file = "Data/Instances/manhattan-long-trips.xml" #the path of the file which contains the manhatten data
#global Manhatten_network_Metagraph_file = "Data/manhatten_graph.mg"
global Manhatten_network_Metagraph_file = "Data/test_graph.mg"

# the file which contains the scenario details (that's need to move from here)
#global scenario_path = "Data/Scenarios_1000_greaterthan2/Output1000_1.txt"
global scenario_path = "Data/Scenarios_1000_greaterthan2/1test_scenario.txt"

#general parameters for the simulation
global maximum_walking_time = 5 # in min
global walking_speed = 1/300 #( 1 link by 5 min) #1.34 # m/s
global driving_speed = 1/60 #(1 link by 1 min) # 50Km /have

global battery_capacity = 63360 # kWs
global charging_rate = 17.600  # kW

global car_cost = 20000 #euro
global fast_charging_unit = 22000 #euro
global cost_factor = 10^6 # or 10^5

global online_request_serving = false # true if we consider the requests by their arrival order
global penality = -10^6 # if the solution is infeasible so we return this penality

#different status for the car
global const CAR_PARKED = 1
global const CAR_RESERVED = 2
global const CAR_ON_WAY = 3

#different status for parking potential_selected_parked_cars_df
global const P_FREE = 1
global const P_RESERVED = 2
global const P_OCCUPIED = 3

# printing parameters
global print_preprocessing = false
global print_simulation = true