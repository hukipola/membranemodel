#=
Start einer Umschreibung meiner Membranmodellierung mit Julia, wobei die Stoffwerte über cantera mit PyCall verwendet werden
=#

using PyCall
using Conda
using Plots
using RecursiveArrayTools, DiffEqBase, OrdinaryDiffEq
using LinearAlgebra

# Importieren der python Pakete
constants = pyimport("scipy.constants")
ct = pyimport("cantera")

# Hinzufügen der neuen Stoffdatenbank
ct.add_directory("/home/matthiskurth/Dokumente/Programmieren/python/Cantera/membran-reaktor/cantera/data")

# Aufrufen der Stoffdatenbank
gas = ct.Solution("surfaceNi.yaml")

# Zunächst die Referenzbedingunen einpflegen
include("ref.jl")
# Dann die Funktionen laden
include("funktionen.jl")


### Ab hier Beschreibung der Membran
# Die variablen Stoffwerte werden direkt in der Funktion aufgerufen

using ModelingToolkit, DifferentialEquations, Symbolics

# define the problem
# p(x) ist der Vektor der Partialdrücke
# T ist die Temperatur an jeder Stelle

function membranreaktor(du, u, params, x)
  T = u.x[1]
  
  wdot = kinetik()
  α_mem, α_heat = params
  du.x[1] .= 0
  du.x[2] .= ones(size(u.x[2]))
  #((wdot) * constants.R * gas.T * 1000)/ (u_0)  - (u.x[2]/ gas.T * du.x[1])
  nothing
end


u0 = ArrayPartition(T_0, gas.X)
params = [α_mem, α_heat]
xspan = (0, length)

problem = ODEProblem(membranreaktor, u0, xspan, params)
# für später
#D(p) ~ ((kinetik()- dgm_cantera() * Vo_Ar_fact * α_mem) * constants.R * gas.T * 1000)/ (u_0)  - (p(x)/ gas.T * D(x) )

soltn = solve(problem, Rosenbrock23())
plot(soltn)

