using JuMP
using Gurobi
#= 
    let try to  solve the following integer program:
        max f =3x1 +2x2
        s.t x1 + x2 <= 6
        x1, x2 >= 0, x1, x2 integer
 =#

#define the JuMP Model
m = Model()

# define the variables 
# @variable(model, nameOfVariable, Type)
@variable(m, x1, Int)
@variable(m, x2, Int)

# define the constraints
# @constarint(model, NameOfConstarint, expression)
@constraint(m, constraint1, x1 + x2 <= 6)
@constraint(m, constraint2, x1 >= 0)
@constraint(m, constraint3, x2 >= 0)

# define the objective function
@objective(m, Max, 3 * x1 +2 * x2)

#set the Optimizer
set_optimizer(m, Gurobi.Optimizer)

optimize!(m)

best_x1 , best_x2 = value(x1) , value(x2)
objValue = objective_value(m)

@variable(m, x, Int)
@constraint(m, con1, x <= 85)
@objective(m, Max, -x^2 + 85*x)

value(x)


x = collect(0:85)
f(x) = 85*x -x^2 

using Plots

plot(x, f)
42*43
