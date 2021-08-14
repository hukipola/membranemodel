# Definition der Referenzbedingungen
T_ref  = 298.15
p_ref = ct.one_atm

# Startbedingungen
T_0 = 280 + 273.15

# Zusammensetzung
komponenten = ["CO2" "H2" "CH4"]
zusammensetzung = [1 4 2]
# Stöchiometrische Faktoren der Reaktion
ν = [1, 0, 0, -1, -4, 2]
comp_0 = Dict(zip(komponenten, zusammensetzung ))

P_ini = 10 * ct.one_atm

λ_ax = 1.403 # W/m2 *K
# Andere Eigenschaften des Reaktors
## Reaktor
length =  0.4 # *approximate* PFR length [m]
# u_0 = 5e-2  # inflow velocity [m/s]
GHVS = 1000 # Bei Standardbedingungen T_ref und p_ref
# todo: u_0 im funktionsaufruf selber berechnen, da es sich ändert

d_react = 15-3 # Reaktordurchmesser in m
A_r = π * (d_react/2)^2  # cross-sectional area [m**2]
volume = A_r * length
# Berechnung der Leerrohr-Geschwindigkeit bei Prozessbedingungen durch Verwendung des id. Gasgesetzes aus der GHVS
u_0 = GHVS * (T_0/T_ref * p_ref/P_ini) * volume /(A_r*3600) # 1/h * m³ * h/(m²*s) = m/s

mass_flow_rate = u_0 * gas.density * A_r # m/s * kg/m³ * m² = kg/s
volume_flow_rate = u_0 * A_r *60/(1000) # m/s * m² = m³/s = L/min

ε_bed = 0.39          # Porosität des Bettes des Reaktors/void fraction

ρ_bed =  1475 
ΔH_R = 165e3 # J/mol
#=
Membran
=#
mem_porosity = 0.38
x_out =  Dict( zip( gas.species_names, zeros(size(gas.species_names)) ) ) # Konzentration der Komponenten ausserhalb der Membran
d_0 = 20e-9 # pore size of membrane in m
L_mem = 350e-6 # Thickness of membrane in m
p_out = 1 # Druck ausserhalb der Membran in bar

#=
Reaktor
=#

# Wärmeübertragung
U = 20 # W/m²-K
T_w = 273.15 + 280 # Wandtemperaur in K
Vo_Ar_fact = 4/d_react

# Anpassungsfaktoren
α_mem = 1e-8
α_heat = 1e4


# Definiere die Startbedingungen des Gases
gas.TPX = T_0, P_ini, comp_0

