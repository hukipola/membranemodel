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

# Ab hier Beschreibung der Membran
# Die variablen Stoffwerte werden direkt in der Funktion aufgerufen

using DifferentialEquations

# define the problem


function membran_faria!(du, u, params, x)
  leng = size(gas.X,1)
  T, dTdx, p_ges, c_u0 = u[1],u[2], u[3], u[4:end]
  
  println(gas.X)
  if gas.X[5] == 0
    wdot = 0
  else 
    wdot = kinetik(T, p_ges) * 1e3 / ρ_bed 
  end
  gas.TP = T, p_ges 
  #gas.concentrations = c_u0/(u_0 * 1e3)
  gas.concentrations = c_u0./(u_0 * 1e3)
  ρ_gas_mole = gas.density_mole 
  ρ_gas_mass = gas.density_mass

## Beginn der Beschreibung der Gleichungen
  du[1] = u[2] 
  #du[2] = 0
  du[2] = ((u_0 * ρ_gas_mass * gas.cp_mass * length/λ_ax) * du[1]) - # [K] 
  (π * d_react * length^2/(A_r * λ_ax) * U * ( gas.T - T_w ) * α_heat) + # [K]
  (ρ_bed * length^2/λ_ax) *(ΔH_R * wdot) #[K]
  du[3] = 0 
  #du[4:end] .= 0
  du[4:end] .=  ρ_bed * length * ν .* wdot 
end

T_in = T_0
dTdx_0 = - (u_0 * gas.density_mole*gas.cp_mass*length/λ_ax) * (T_in - T_0)
u0 = vcat(T_0, dTdx_0, gas.P, gas.concentrations*1e3.*u_0) 
params = [α_mem, α_heat]
xspan = (0, length)

prob = ODEProblem(membran_faria!, u0, xspan, params)
alg = AutoTsit5(Rosenbrock23())
sol = solve(prob,  reltol=1e-4, abstol=1e-4 )
plot(sol, vars = 1)
