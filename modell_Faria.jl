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


function membran_faria!(du, u, params, x)
  leng = size(gas.X,1)
  T, dTdx, p_ges, c_u0 = u[1],u[2], u[3], u[4:end]
  wdot = kinetik(T, p_ges) 
  gas.TP = T, p_ges 
  println(T)
  gas.concentrations = c_u0/(u_0 * 1000)
  rho = gas.density_mole
  du[1] = dTdx  
  #du[2] = 0
  du[2] = (u_0 * rho * gas.cp_mass * length/λ_ax) * dTdx - (ρ_bed * length^2/λ_ax) *(-ΔH_R) - ((π*d_react * length^2/λ_ax) .* U .* ( T - T_w ) .* α_heat  )  
  du[3] = 0
  #du[4:end] .= 0
  du[4:end] .=  -1 * length * ν .* wdot 
end


u0 = vcat(T_0, 1, gas.P, gas.concentrations*1e3.*u_0) 
params = [α_mem, α_heat]
xspan = (0, length)

prob = ODEProblem(membran_faria!, u0, xspan, params)
sol = solve(prob)
plot(sol, vars = 1)
