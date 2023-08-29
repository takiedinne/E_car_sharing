"""
    this file regroups all the general parameters for the simulation
"""

""" type of cars available in the simulation:"""
@enum car_type Smart_ED Mitsubishi_iMIEV Nissan_LEAF

"""the specification for the cars"""
global vehicle_specific_values = Dict{car_type, Any}(   Smart_ED => Dict(:battery_capacity => 63360000, :slow_charging_rate =>17600, 
                                                            :fast_charging_rate=>17600, :car_cost=>20000,
                                                            :α => 126.33, :β=>0.52, :γ => 1000), 
                                            Mitsubishi_iMIEV => Dict(:battery_capacity => 57600000, :slow_charging_rate =>22000, 
                                                            :fast_charging_rate=>30000, :car_cost=>24000,
                                                            :α => 126.33, :β=>0.52, :γ => 1000),# recalculate these values 
                                            Nissan_LEAF => Dict(:battery_capacity => 86400000, :slow_charging_rate =>22000, 
                                                            :fast_charging_rate=>38400, :car_cost=>29000,
                                                            :α => 126.33, :β=>0.52, :γ => 1000),# recalculate these values
                                                            )
"""
    project_path(parts...) = normpath(joinpath(@__DIR__, "../..", parts...))

    project_path give you the absulte path for each relative path in inputs
    inputs: parts... : the relative path (as parts)
    output: the absolute path
"""
project_path(parts...) = normpath(joinpath(@__DIR__, "../..", parts...))


# Different File paths
global all_request_details_path = project_path("Data/Instances/trips_ML_daily_greaterthan2_with_revenue&minutes.txt")
global Manhatten_network_details_file = project_path("Data/Instances/manhattan-long-trips.xml") #the path of the file which contains the manhatten data

global Manhatten_network_driving_graph_file = project_path("Data/other/manhatten_driving_graph.mg")
global Manhatten_network_walking_graph_file = project_path("Data/other/manhatten_walking_graph.mg")
global Manhatten_network_length_graph_file = project_path("Data/other/manhatten_length_graph.mg")

function custom_sort(filename)
    match_result = match(r"Output_(\d+)\.txt", filename)
    if match_result !== nothing
        return parse(Int, match_result.captures[1])
    else
        return Inf
    end
end
global scenarios_paths = project_path.("Data/Instances/C1_5000_500/scenario_txt_files/" .* filter!(x -> startswith(x, "Out"), readdir(project_path("Data/Instances/C1_5000_500/scenario_txt_files"))))
sort!(scenarios_paths, by=custom_sort)

#general parameters for the simulation:
global maximum_walking_time = 5 # in min
global walking_speed = 1.34 # m/s
global time_slot_length = 5 # min
global fixed_driving_speed =  50 #Km/h

# incase we want to use dynamic speeds:
global travel_speed_file_path = project_path("Data/other/travel_speeds.jls")

#different status for the car:
global const CAR_PARKED = 1
global const CAR_RESERVED = 2
global const CAR_ON_WAY = 3

#different status for parking potential_selected_parked_cars_df:
global const P_FREE = 1
global const P_RESERVED = 2
global const P_OCCUPIED = 3

# printing parameters:
global print_preprocessing = false
global print_simulation = false
global printMIP = false
# generale parameters:
global online_request_serving = false # true if we consider the requests by their arrival order
global penality = 10^16 # if the solution is infeasible so we return this penality
global work_with_time_slot = true
global cost_factor = 10^6 # or 10^5
global use_dynamic_speeds = false # to simulate the trafic congestion
global multiple_driving_speeds = true # take into account the max speed for each edge

"""
    set_online_mode(mode::Bool)

    assigne new value to the global variable online_request_serving from outside this module
"""
function set_online_mode(mode::Bool)
    global online_request_serving = mode
end

"""
    set_use_dynamic_speeds(use_dynamic_speed::Bool)

    assigne new value to the global variable use_dynamic_speeds from outside this module 
"""
function set_use_dynamic_speeds(use_dynamic_speed::Bool)
    global use_dynamic_speeds = use_dynamic_speed
end

