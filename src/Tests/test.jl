using E_car_sharing
using Serialization

sol =  deserialize("/Users/taki/Desktop/Preparation doctorat ERM/Projects/GIHH_V2.0/sol.jls")
initialize_scenarios([1])
E_carsharing_sim(sol)
sol = generate_random_solution()

a = 10 ^ 1000
b = 10 ^ 100 
a = typemax(Int128) 

b = a ^(1/(2+(10/ 1000)))
c = a ^(1/(2+(11/ 1000)))

b -c 
c -b 