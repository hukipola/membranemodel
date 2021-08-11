using ModelingToolkit

@variables t x(t)  # independent and dependent variables
@parameters τ       # parameters
D = Differential(t) # define an operator for the differentiation w.r.t. time

# your first ODE, consisting of a single equation, indicated by ~
@named fol_model = ODESystem(D(x) ~ (1 - x)/τ)
      # Model fol_model with 1 equations
      # States (1):
      #   x(t)
      # Parameters (1):
      #   τ
      
using DifferentialEquations
using Plots

prob = ODEProblem(fol_model, [x => 0.0], (0.0,1.0), [τ => 1.0])
plot(solve(prob))
