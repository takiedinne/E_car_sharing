# data Structure:
mutable struct Ecar
    start_charging_time::Integer # time when the car is plugged-in the charging units
    last_battery_lavel::Float64
end

mutable struct Station
    cars::DataFrame 
    #= each row represents a car (car_id, status,
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
end

Station(total_parking_places::Integer, initial_cars_number::Integer, initial_id::Integer) = begin
    
    id = initial_id
    
    cars = DataFrame(car_id = collect(id:id+initial_cars_number-1), status = CAR_PARKED .* ones(Integer, initial_cars_number),
                    last_battery_level = 100.0 .* ones(Float64, initial_cars_number), start_charging_time = zeros(initial_cars_number),
                     start_reservation_time = zeros(initial_cars_number), expected_arrival_time = zeros(initial_cars_number))
    
    #create the parking places
    parking_places = DataFrame(p_id = collect(1:total_parking_places), 
                                status = vcat(P_OCCUPIED .* ones(Integer, initial_cars_number), P_FREE .* ones(Integer, total_parking_places - initial_cars_number)), 
                                cars = vcat(collect(id:id+initial_cars_number-1), -1 .* ones(Integer, total_parking_places - initial_cars_number)),
                                pending_reservation = zeros(Integer, total_parking_places) )
    
    
    Station(cars, parking_places)
end


mutable struct Solution
    # decision variables
    open_stations_ids::Array{Integer, 1}
    initial_cars_number::Array{Integer, 1}
    selected_paths::Array{Bool, 1}
    
    # other fields
    #stations::Array{Station, 1}
end
Solution() = Solution([], [], []) 

#= 
    generate a random solution:
        1- decide how much station to open if it is not precised in @open_stations_number
        2- open random stations 
        3- set random initila number of cars for each station
    inputs:
        @open_stations_number (optional): the number of station to open
    outputs:
        sol:: Solution 
=#
function generate_random_solution(; open_stations_number=-1)

    sol = Solution()
    potential_locations = get_potential_locations()
    #decide the number of stations to open if it is not precised by the user
    if open_stations_number == -1
        open_stations_number = rand(1:length(potential_locations))
    end
    #randomly open stations
    sol.open_stations_ids = sample(potential_locations, open_stations_number, replace=false)
    #sol.open_stations_ids = get_potential_locations()

    #set initial car number for each station
    for station in sol.open_stations_ids
        max_number_of_charging_points = get_prop(manhaten_city_graph, station, :max_number_of_charging_points)
        push!(sol.initial_cars_number, rand(1:max_number_of_charging_points))
    end

    sol
end

#=
    Check whether or not the solution is feasible according to constraints 2, 3, 7 and 8 in Hatiçe paper
    inputs:
        @sol: the solution
        @all_feasible_paths: the set all feasible paths given by the preprocessing procedure(all_requests_feasible_paths)
    outputs:
        Feasible => Boolean which is true if the solution is feasible, false otherwise.
=#
function is_feasible_solution(sol::Solution)
    global all_feasible_paths
    selected_paths = all_feasible_paths[sol.selected_paths, :]
    
    #check the initial number of cars (constraint 7 and 8)
    for i in 1:length(sol.open_stations_ids)
        if sol.initial_cars_number[i] > get_prop(manhaten_city_graph, sol.open_stations_ids[i], :max_number_of_charging_points)
            println("the initial number of cars in the station ", sol.open_stations_ids[i], " is greater the the total allowed number")
            return false
        end
    end

    # check if the customer is served by  at most one trip (constraint 2)
    if nrow(selected_paths) != length(unique(selected_paths.req))
        println("there is at least a customer served by more than one trip")
        return false
    end
    
    open_stations_df = DataFrame(station_ids = sol.open_stations_ids) # usuful for innerjoin
    
    #check if each station in the selected trips is open (constraint 3)
    if nrow(innerjoin(selected_paths, open_stations_df, on = :origin_station => :station_ids)) != nrow(selected_paths) ||
        nrow(innerjoin(selected_paths, open_stations_df, on = :destination_station => :station_ids)) != nrow(selected_paths)
        
        println("there is at least one closed station in the feasible paths")
        return false
    end

    # constraint 4, 5 and 6  are to be checked in the simulation

    return true
end


