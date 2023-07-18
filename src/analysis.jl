using CSV
using DataFrames
using StatsPlots


"""
    plot_best_fitness_tracking(tracking_file_path::String; data_range::range=nothing)

    plot the progress of the best solution during the search process.
    @tracking_file_path: the path of the file containing the best fitness for each Iteration
    @data_range: the range of the data to be plotted (by default we plot all the data)
"""
function plot_best_fitness_tracking(tracking_file_path::String; data_range=nothing)
   data_df = CSV.read(tracking_file_path, DataFrame ,delim=';'; select = [:Iteration, :BestFitness], decimal = ',')
   
   #delete all the points wher fit equal to 10 ^ 16
   #filter!(row -> row.BestFitness != 1e16, data_df)

   #select the range of data
   data_df = isnothing(data_range) ? data_df : data_df[data_range, :]
   #plot the data
   @df data_df plot!(:Iteration, :BestFitness, xlabel="Iteration", ylabel="Fitness", title="Best Fitness")
   
end

"""
    read_heuristics_performance(performance_file_path::String; iteration_number::Int64)

    read the performance of the heuristics for a given iteration 
    it return a DataFrame of heuristicId the number of calls, the number of improvement,
    the number of equal and the number of worse solutions produced by heursitic.
    @performance_file_path: the path of the file containing the performance of the heuristics
    @iteration_number: the iteration number for which we want to read the performance 
                        (if isnothing it returns the last iteration)
    @save_file_path: the path of the file where the performance will be saved
"""
function read_heuristics_performance(performances_file_path::String; save_file_path::String = "", iter_num::Int64 = -1)
    
    df = CSV.read(performances_file_path, DataFrame, delim = ";")
    
    # get the number of heuristics
    heuristics_number = Int64((size(df, 2) - 4) / 6)

    # iteration number
    iter_num = (iter_num==-1) ? size(df, 1) : iter_num
    heuristics_performance_df = DataFrame(id = [], calls_number = [], best_found_number = [], improvement_number = [] , equal_number = [], worse_number = [], cpu_time = [])
    
    for i in 1:heuristics_number
        shift = 4+(i-1)*6
        push!(heuristics_performance_df,  (i, Tuple(df[iter_num, shift:(shift+5)])...) )
    end
    
    if save_file_path != ""
        CSV.write(save_file_path, heuristics_performance_df, decimal = ',')
    end
    heuristics_performance_df
end

"""
    get_nbr_improvement(tracking_file_path::String)

see the the number of times that the solver could improve the overall best solution
"""
function get_nbr_improvement(tracking_file_path::String)
    data_df = CSV.read(tracking_file_path, DataFrame ,delim=';'; select = [:Iteration, :BestFitness], decimal = ',')
   
    counter = 0
    for i in 2:nrow(data_df)
        if data_df[i, :BestFitness] < data_df[i-1, :BestFitness]
            counter += 1
        end
    end
    @info "[solver Improvement capability]: $(counter) improvement(s) out of $(data_df[end, :Iteration]) iteration(s)"
    
    return counter
end

tracking_file_path = "/Users/taki/Desktop/Preparation doctorat ERM/Projects/GIHH_V2.0/results/solve_single_scenario/1/GIHH__HCnt.csv"
plot_best_fitness_tracking(tracking_file_path#= , data_range = 1:10 =#)

get_nbr_improvement(tracking_file_path)

performances = read_heuristics_performance(tracking_file_path#= , save_file_path = "/Users/taki/Desktop/Preparation doctorat ERM/Projects/E_car_sharing/results/analysis/heuristics_performance.csv" =#);

sort!(performances, :improvement_number, rev = true)
