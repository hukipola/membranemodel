using Plots
using DataInterpolations
using PyCall
using Conda
using Plots

#
#
#
#Importieren der python Pakete
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

T_int = 273.15:5:700

komponenten = ["CO2" "H2" "H2O" "CH4" "CO"]
zusammensetzung = [rand(1) for i in 1:1:size(komponenten,2)]

comp_0 = Dict(zip(komponenten, zusammensetzung ))
gas.TPX = T[1], P_ini, comp_0 
kinetik()

Curvefit()
