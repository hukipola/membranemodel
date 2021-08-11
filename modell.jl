#=
Start einer Umschreibung meiner Membranmodellierung mit Julia, wobei die Stoffwerte über cantera mit PyCall verwendet werden
=#

using PyCall
using Conda
using Plots
using DiffEqBase, OrdinaryDiffEq
using LinearAlgebra

# Importieren der python Pakete
constants = pyimport("scipy.constants")
ct = pyimport("cantera")
np = pyimport("numpy")

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

using DifferentialEquations

# define the problem

function membran!(du, u, params, x)
  T, p = u[1], u[2:end]
  println(T)
  gas.TPX = T, sum(p), p/sum(p)
  rho = gas.density_mole
  wdot = kinetik(T,p)
  du[1]= ( -(dot(gas.partial_molar_enthalpies, (wdot .* ν)) .+ (U .* Vo_Ar_fact .* ( T - T_w ) .* α_heat  ) )./ # J/(s*m³)
                  (rho .* gas.cp_mole) .* u_0) # J/(m³*K) --->  # Energiebilanz in K/m
  du[2:end] .= ((((wdot .* ν) .* constants.R .* T .* 1000)./u_0) - ((p./T) .* du[1]))[1:end]
end


u0 =vcat(T_0, gas.X*gas.P) 

params = [α_mem, α_heat]
xspan = (0, length)

prob = ODEProblem(membran!, u0, xspan, params)

sol = solve(prob, Tsit5())
plot(sol) 
