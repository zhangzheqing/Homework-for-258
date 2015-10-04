
# Penalty funtion= l_2 norm



using MathProgBase
using PyPlot
using ECOS

using Convex
# Generate random problem data
m = 100;  n = 30
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

# The problem is to minimize ||Ax - b||^2 
problem = minimize(sum_squares(A * x - b), [x >= 0])

# Solve the problem by calling solve!
solve!(problem,ECOSSolver())

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval

# Get the residue values.
r=Variable(n)
r = A * x.value - b 

# hist plot using PyPlot.

fig = figure("pyplot_histogram",figsize=(10,4)) # Not strictly required
ax = axes() # Not strictly required
nbins = 50
h = plt[:hist](r, nbins) # Histogram, PyPlot.plt required to differentiate with conflicting hist command

grid("on")
xlabel("r")
ylabel("Counts")
axis([-5, 5, 0 ,40])
title("Histogram")

# plot penalty function
x = linspace(-5,5,1000); y = x.^2;
plot(x, y, color="red", linewidth=2.0, linestyle="--")



