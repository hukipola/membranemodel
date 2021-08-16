function kinetik(T,p_tot)

    # Methanisierung von CO2 Kinetik nach
    # Schlereth, David; Hinrichsen, Olaf (2014): A fixed-bed reactor modeling study on the methanation of CO2. In: Chemical Engineering Research and Design     92 (4), S. 702–712. DOI: 10.1016/j.cherd.2013.11.014.
    #     und
    # Methane steam reforming, methanation and water‐gas shift: I. Intrinsic kinetics - Xu - 1989 - AIChE Journal - Wiley Online Library. Online verfügbar unter https://onlinelibrary.wiley.com/doi/epdf/10.1002/aic.690350109, zuletzt geprüft am 07.02.2019.

        """
        Eine sehr einfache Darstellung der Methanisierung, die über einen einfachen Ansatz nach Koschany
        Eine Erweiterung wäre die Verwendung von Xu et al für die anderen komponenten.
        Gibt die Reaktionsrate zurück in (kmol/m³-s) """
        
        rho_kat = ρ_bed              # Dichte des Katalysators in kg/m3. Lässt sich aber auch anders berechnen
        # Wirkungsgrad des Katalysators
        nu_cat = 0.1

        # temp = gas.T # Temperatur in K
        temp = T # Temperatur in K
        # p_tot = gas.P # Druck in Pa
        # p_tot =  p # Druck in Pa
        mole_dict = Dict(zip( gas.species_names , gas.X )) # Molenbruch der einzelnen Komponenten als dict
        # Die Partialdrücke der einzelnen komponenten
        pH2 = get(mole_dict, "H2", 0) * p_tot/1e5 #[Pa -> bar]
        pCO = get(mole_dict, "CO", 0) * p_tot/1e5
        pCO2 =get(mole_dict, "CO2", 0)  * p_tot/1e5
        pCH4 = get(mole_dict, "CH4", 0)  * p_tot/1e5
        pH2O = get(mole_dict, "H20", 0)  * p_tot/1e5

        # Methanisierung -> Konstanten von kJ in J umgerechnet -> Faktor 1000 bzw e3
        Tref = 555 # K
        # Gleichgewichtskonstanten [Aparicio 1997]
        K3_eq = 137 * temp^(-3.998) * exp(158.7e3/(constants.R * temp))
        # Hier eine Alternative Beschreibung der Gleichgewichtskonstante nach Uhlmann, die ebenfalls verwendet werden kann.
        # K1 = (10^((-18.98177) + (-9520.569 /temp) + 10.91012 * log10(temp) + (-0.319556e-2)*temp + 0.3631894e-6 * temp^2) ) # [bar**2] CO + 3H2 <-> CH4 + H2O
        # K2 = (10^((-8.000980) + (2456.620 /temp) + 1.984098 * log10(temp) + (-0.3329441e-3)*temp + 0.563315e-7 * temp^2) ) # [1] CO + H2O <-> CO2 + H2
    # Umwandlung der Gleichgewichtslage in K_p zu K_x über die Umwandlung in Molenbrüche
        # K_uhl = 1/(K1 * K2)/((p_tot)^(-2))

          # # Reaktionsgeschwindigkeitskonstante
        k =  10^3 * 3.46e-4 * exp(77.5e3/constants.R * (1/Tref-1/gas.T))  #[mol/(bar*kg_cat*s)] -> inkl Umrechnung von g in kg_cat -> Vorfaktor 10**3

        # # Adsorptionskonstanten
        K_OH = 0.5 * exp(22.4e3/constants.R*(1/Tref-1/temp)) ##[bar**-05]
        K_H2 = 0.44 * exp(-6.2e3/constants.R*(1/Tref-1/temp)) ##[bar**-05]
        K_mix = 0.88 * exp(-10e3/constants.R*(1/Tref-1/temp)) ##[bar**-05]


        # gemeinsmer Nenner:
        DEN = (1 + K_OH*pH2O/(pH2^0.5) + K_H2*pH2^0.5 + K_mix*pCO2^0.5) #[-]

        r_I = ((k*pH2^0.5*pCO2^0.5*(1-(pCH4*pH2O^2)/(pCO2*pH2^4*K3_eq)))/DEN^2) * rho_kat / 1000 #[mol/(kg_cat*s)] # Änderung der GGW - Konstanten zu der Berechnung von          Uhlmanns # (mol/kg_cat-s) * rho_cat(kg/m3) * 1e-3 (kmol/mol) = (kmol/m³-s)

        # Erstellen eines Vektors aus Nullen
        r = zeros(gas.n_total_species)
        # Der Wert der einzelnen Komponenten entspricht wird mit der Reaktionsgeschwindigkeit überschrieben
        r[gas.species_index("CO2")+1] = -r_I
        r[gas.species_index("CH4")+1] = r_I
        r[gas.species_index("H2O")+1] = 2 * r_I
        r[gas.species_index("H2")+1] = - 4 * r_I

        # Verwendung von Katalysator void fraction mit u_s (z,t)= ε_bed∙u_(g,real) (z,t), da die Geschwindigkeit dadurch eine andere ist.

        # return r*nu_cat*(1-ε_bed)
        return r_I*nu_cat*(1-ε_bed)
end

function dgm_cantera()
        g = ct.DustyGas("surfaceNi.yaml")
        # Seite 1 ist Innenerhalb des Reaktors und Seite 2 ist Ausserhalb des Reaktors

        g.TPX = gas.T, gas.P , gas.X

        g.porosity = mem_porosity
        g.tortuosity = mem_porosity^(-0.5)

        g.mean_pore_radius = d_0
        #todo sehr fraglicher Parameter
        g.mean_particle_diameter = 10e-7  # lengths in meters

        # Berechne die Membraninnenseite
        T1, rho1, Y1 = g.TDY
        # Setze neue Daten ein und berechne ebenfalls die neuen Daten
        g.TP = gas.T, p_out * ct.one_atm
        T2, rho2 = g.TD

        return g.molar_fluxes(T1, T2, rho1, rho2, Y1, zeros(size(gas.species_names)), L_mem)
end

