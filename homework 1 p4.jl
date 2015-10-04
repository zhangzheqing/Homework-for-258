
using MathProgBase
using PyPlot
using ECOS
using SCS
using Convex
# Generate random problem data
m = 100;  n = 30
 A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

a = 1
# The problem is to minimize  sum (- a^2 * log(1 - (r/a)^2))

#### First I solve for the min(norm(Ax-b, inf)) to see if it is possible for us to have a starting point
#### If not I will reset A and b.
pre_problem = minimize(norm(A * x - b, Inf))

solve!(pre_problem, ECOSSolver())
A = A/(1.1*pre_problem.optval)
b = b/(1.1*pre_problem.optval)

function f(x)
    - a^2 * log(1- (x/a)^2)
end
    
problem = minimize(sum(-log(1-(A*x-b)) - log(1+(A*x-b)) ) )
#problem.constraints += [A*x-b <=0.999; A*x-b >= -0.999]

# Solve the problem by calling solve!
solve!(problem, ECOSSolver())

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
axis([-2, 2, 0 ,20])
title("Histogram")

# plot penalty function
x1 = linspace(-0.99999,0.99999,1000); y1 = -a^2 * log(1- (x1/a).^2);
plot(x1, y1, color="red", linewidth=2.0, linestyle="--")


