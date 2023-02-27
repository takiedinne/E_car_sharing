mutable struct Station
    cars::DataFrame 
    #= each row represents a car (car_id, car_type (ED) status,
            last_battery_level, start_charging_time,
            start_reservation_time, expected_arrival_time)
     status of the cars : 1=> parked,
                          2=> reserved (check the start reservation time), 
                          3=> in its way to come (check the excpected arrival time)
    =#
    parking_places::DataFrame 
    #= each row represents a parking place (id, 
        status (1=> free; 2=> reserved, 3=>occupied),
        cars (id) (if status == 3 ) 
    =#
    max_number_of_charging_points
    max_power
    charging_station_base_cost
    charging_point_cost_slow
    max_charging_rate_per_charging_point_slow
    charging_point_cost_fast
    max_charging_rate_per_charging_point_fast
end

Station(station_node_id::Int64, initial_cars_number::Int64, initial_id::Int64) = begin
    #get teh different information from the graph
    prop_dict = props(manhaten_city_graph, station_node_id)
    
    total_parking_places = prop_dict[:max_number_of_charging_points]
    max_power = prop_dict[:max_power]
    charging_station_base_cost = prop_dict[:charging_station_base_cost]
    charging_point_cost_slow = prop_dict[:charging_point_cost_slow]
    max_charging_rate_per_charging_point_slow = prop_dict[:max_charging_rate_per_charging_point_slow]
    charging_point_cost_fast = prop_dict[:charging_point_cost_fast]
    max_charging_rate_per_charging_point_fast = prop_dict[:max_charging_rate_per_charging_point_fast]
    
    id = initial_id
    
    cars = DataFrame(car_id = collect(id:id+initial_cars_number-1), car_type = [Smart_ED for i in 1:initial_cars_number],
                    status = CAR_PARKED .* ones(Int64, initial_cars_number),
                    last_battery_level = 100.0 .* ones(Float64, initial_cars_number), start_charging_time = zeros(initial_cars_number),
                    start_reservation_time = NaN * ones(initial_cars_number), pending_reservation = zeros(initial_cars_number), expected_arrival_time = NaN * ones(initial_cars_number))
    
    #create the parking places
    parking_places = DataFrame(p_id = collect(1:total_parking_places), 
                                status = vcat(P_OCCUPIED .* ones(Int64, initial_cars_number), P_FREE .* ones(Int64, total_parking_places - initial_cars_number)), 
                                cars = vcat(collect(id:id+initial_cars_number-1), -1 .* ones(Int64, total_parking_places - initial_cars_number)),
                                pending_reservation = zeros(Int64, total_parking_places) )
    
    
    Station(cars, parking_places, total_parking_places, max_power, charging_station_base_cost, charging_point_cost_slow,
             max_charging_rate_per_charging_point_slow, charging_point_cost_fast, max_charging_rate_per_charging_point_fast)
end

mutable struct Solution
    # decision variables
    open_stations_state::Array{Bool, 1}
    initial_cars_number::Array{Int64, 1}
    selected_paths::Array{Array{Bool,1}, 1}
end
Solution() = Solution(falses(length(get_potential_locations())), zeros(length(get_potential_locations())), []) 



mutable struct Scenario
    scenario_id::Int64
    request_list::DataFrame 
    feasible_paths::DataFrame
end

