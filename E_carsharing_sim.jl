using ResumableFunctions
using SimJulia

using DataFrames
using CSV

using Graphs, MetaGraphs
using EzXML

using Distances
using StatsBase

include("Solution.jl")
include("util.jl")
include("Vars.jl")

#create the graph
# manhaten_city_graph = create_graph_from_XML(Manhatten_network_details_file, save_file = "Data/manhatten_graph.mg")

#to load the graph use the following instruction ===> it is better to load the graph in terme of computational performance
global manhaten_city_graph = loadgraph(Manhatten_network_Metagraph_file, MGFormat())

global non_directed_manhaten_city_graph = MetaGraph(manhaten_city_graph) # for walking purposes we don't take into account the directed edges 

#read all requests details
global all_request_df = CSV.read(all_request_details_path, DataFrame)

#initialization
global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph

#result of the pre-processing step
global all_feasible_paths

global failed = false #for the offline mode it is needed to stop the simulation

global revenues = 0 #

global stations = Array{Station, 1}() 
global scenario # the requests list

@resumable function request_arrival_process(env::Environment, scenario::DataFrame, sol::Solution)

    #browse all the requests (the requests are sorted according to their arrival time)
    for req in eachrow(scenario)

        # waiting until a new request is arrived 
        sleeping_time = req.ST - now(env)
        @yield timeout(env, sleeping_time)
        print_simulation && println("Customer [", req.reqId, "]: request arrive at ", now(env))

        # get the trip information (pickup station, drop off station, the selected car id, parking place)
        # if the request can not be served, this function will return -1 in one of the information variables
        pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id = get_trip_info_for_request(req, sol, now(env)) # one row Dataframe

        # check the availabilty of paths to serve the request
        if pickup_station_id != -1 && drop_off_station_id != -1 && selected_car_id != -1 && parking_place_id != -1
            
            @process perform_the_trip_process(env, req, sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
        else #we couldn't serve the request there is no feasible path
            print_simulation && println("Customer [", req.reqId, "]: We can not serve the request there is no feasible path")

            #if we are working on the ofline mode so we need to return a penality 
            if !online_request_serving
                global failed = true
                break
            end
        end

    end

end


@resumable function perform_the_trip_process(env::Environment, req::DataFrameRow, sol::Solution, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
    #some information
    walking_duration = get_walking_time(req.ON, sol.open_stations_ids[pickup_station_id])
    start_walking_time = now(env)
    driving_duration = get_trip_duration(sol.open_stations_ids[pickup_station_id], sol.open_stations_ids[drop_off_station_id])
    
    #book the trip (the car + parking palce ,etc)
    book_trip(sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id, start_walking_time + walking_duration, start_walking_time + walking_duration + driving_duration)

    print_simulation && println("Customer [", req.reqId, "]:the request is accepted")
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station ", sol.open_stations_ids[pickup_station_id], " at ",now(env),"and he is taking the car number ", selected_car_id)
    
    print_simulation && println("Customer [", req.reqId, "]: start ridding at ", now(env))
    start_riding(pickup_station_id, selected_car_id)
    #simulate the riding time
    @yield timeout(env, driving_duration)
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", sol.open_stations_ids[drop_off_station_id], " at ", now(env))

    #drop the car
    drop_car(drop_off_station_id, selected_car_id, parking_place_id, now(env))

    #make the payment
    global revenues += req.Rev
end

function initialize_sim(sol::Solution, scenario_path)
    # construct the requests lists 
    global scenario = scenario_as_dataframe(scenario_path)
    global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
    global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph

    #result of the pre-processing step
    global all_feasible_paths
    
    global stations = Array{Station, 1}()
    #set the stations
    car_id = 1
    for i in 1:length(sol.open_stations_ids)
        total_parking_places = get_prop(manhaten_city_graph, sol.open_stations_ids[i], :max_number_of_charging_points)
        initial_cars_number = sol.initial_cars_number[i]

        push!(stations, Station(total_parking_places, initial_cars_number, car_id))
        car_id += initial_cars_number
    end
    
    # preprossesing 
    if !online_request_serving
        global all_feasible_paths = all_requests_feasible_paths(scenario, sol)
    end
end

function E_carsharing_sim(sol::Solution)
   
    initialize_sim(sol, scenario_path)
     # check the feasibilty of the solution
     is_feasible_solution(sol)
    
    sim = Simulation()
    @process request_arrival_process(sim, scenario, sol)
    run(sim)
end

