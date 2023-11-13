#= using SimJulia, ResumableFunctions


    @resumable function start_proc(env::Environment)
        @process inner_fun(env)
    end

    @resumable function inner_fun(env::Environment)
        selected_car_id = 2
        x = selected_car_id 
        println("Car $selected_car_id is selected")
        @yield timeout(env, 1)
    end

    sim = Simulation()

    @process start_proc(sim)

    run(sim) 
=#

using E_car_sharing


set_number_of_requests_per_scenario(4000)
set_walking_time(15)

initialize_scenarios([1])
sol = generate_random_solution()
set_online_mode(true)
using BenchmarkTools

E_carsharing_sim(sol)

set_online_mode(false)
E_carsharing_sim(sol)