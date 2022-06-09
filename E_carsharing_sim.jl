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
# manhaten_city_graph = create_graph_from_XML(Manhatten_network_details_file, save_file = "Data/test_graph.mg")
# manhaten_city_graph = create_graph_from_XML("Tests/test_graph.xml", save_file = "Data/test_graph.mg")

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

global potential_refusal_requests = []

@resumable function request_arrival_process(env::Environment, scenario::DataFrame, sol::Solution)
    
    #browse all the requests (the requests are sorted according to their arrival time)
    for req in eachrow(scenario)
        
        # in the offline mode we check if we failed to serve a request so we stop the simulation
        if !online_request_serving && failed
            break
        end
        # waiting until a new request is arrived 
        sleeping_time = req.ST - now(env)
        @yield timeout(env, sleeping_time)
        print_simulation && println("Customer [", req.reqId, "]: request arrive at ", now(env))

        # get the trip information (pickup station, drop off station, the selected car id, parking place)
        # if the request can not be served, this function will return -1 in one of the information variables
        
        # in the offline mode we can check if the request is decided to be served
        if !online_request_serving && req.reqId ∉ all_feasible_paths[sol.selected_paths, :].req
            # the resuest is refused by the decision variables
            print_simulation && println("Customer [", req.reqId, "]: the requests is rejected according to the decision variables")
            continue
        end
        pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id = get_trip_info_for_request(req, sol, now(env)) # one row Dataframe

         
        # check the availabilty of paths to serve the request
        if pickup_station_id != -1 && drop_off_station_id != -1 && selected_car_id != -1 && parking_place_id != -1
            if pickup_station_id != drop_off_station_id
                @process perform_the_trip_process(env, req, sol, pickup_station_id, drop_off_station_id, selected_car_id, parking_place_id)
            end
        else
            #we couldn't serve the request there is no feasible path
            # the cause message
            if print_simulation
                if (pickup_station_id == -1 || drop_off_station_id == -1)
                    #there is no feasible path for the request
                    cause_message = "there is no feasible path"
                elseif print_simulation && selected_car_id == -1
                    cause_message = "there is no available car"
                elseif print_simulation && parking_place_id == -1
                    cause_message = "there is no parking place"
                end
            end
                 
            print_simulation && printstyled(stdout, "Customer [", req.reqId, "]: We can not serve the request there is no feasible path (",cause_message,")\n", color = :yellow)
           
            #if we are working on the ofline mode so we need to return a penality 
            if !online_request_serving
                # Here we add the request to the list so that we can check later if we are able to serve it.
                push!(potential_refusal_requests, req)
                #= global failed = true
                break =#
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
    
    #= if req.reqId == 11913
        @show stations[drop_off_station_id].cars
        @show stations[pickup_station_id].cars
    end =#
    print_simulation && println("Customer [", req.reqId, "]: the request is accepted")
    print_simulation && println("Customer [", req.reqId, "]: start walking from ", req.ON, " at ", now(env))

    # simulate the walking time
    @yield timeout(env, walking_duration)
    print_simulation && println("Customer [", req.reqId, "]: arrive at the station ", sol.open_stations_ids[pickup_station_id], " at ",now(env)," and he is taking the car number ", selected_car_id)
    
    print_simulation && println("Customer [", req.reqId, "]: start ridding at ", now(env))
    start_riding(pickup_station_id, selected_car_id)
    #simulate the riding time
    @yield timeout(env, driving_duration)
    print_simulation && println("Customer [", req.reqId, "]: drop the car off at the station ", sol.open_stations_ids[drop_off_station_id], " at ", now(env))

    #drop the car
    drop_car(drop_off_station_id, selected_car_id, parking_place_id, now(env))

    #make the payment
    global revenues += req.Rev

    # in the offline mode we check if we are able to serve any requests from the list of potential refusal requests
    if !online_request_serving
        serve_from_potential_refusal_requests()
    end
end

function initialize_sim(sol::Solution, scenario_path)
    global failed = false
    global revenues = 0
    # construct the requests lists 
    global scenario = scenario_as_dataframe(scenario_path)
    global shortest_car_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time
    global shortest_walking_paths = Dict{Integer,Any}() #the results of djisktra algorithms to avoid calling the algorithm each time based on non directed graph
    global potential_refusal_requests = []
    
    global stations = Array{Station, 1}()
    #set the stations
    car_id = 1
    for i in 1:length(sol.open_stations_ids)
        total_parking_places = get_prop(manhaten_city_graph, sol.open_stations_ids[i], :max_number_of_charging_points)
        initial_cars_number = sol.initial_cars_number[i]

        push!(stations, Station(sol.open_stations_ids[i], initial_cars_number, car_id))
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
    if failed || (!online_request_serving && !isempty(potential_refusal_requests))
        print_simulation && printstyled(stdout, "The simulation is stopped because a request couldn't be served\n", color = :light_red)
        return penality
    else
        # count the objective function
        print_simulation && println("counting the objective function")
        total_cars_cost = [[vehicle_specific_values[stations[i].cars.car_type[j]][:car_cost] for j in 1:nrow(stations[i].cars)] for i in 1:length(sol.open_stations_ids)] 

        total_station_cost = sum([station.charging_station_base_cost + 
                                    station.max_number_of_charging_points * station.charging_point_cost_fast
                                    for station in stations])
        return revenues - (total_cars_cost  + total_station_cost ) / cost_factor
    end
end

