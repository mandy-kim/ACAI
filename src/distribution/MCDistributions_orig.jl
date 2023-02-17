Base.@kwdef mutable struct MainDistributions <: MCDistributions 
    emis_unc_FB::ConstantDist
    emis_unc_NOx::ConstantDist
    avNOx2006::ConstantDist
    avfuel2005::ConstantDist
    lambda_O3::ConstantDist
    lambda_Contrail::Union{TriangDist,ConstantDist}
    lambda_Sulfates::ConstantDist
    lambda_Soot::ConstantDist
    lambda_H2O::ConstantDist
    lambda_CH4::ConstantDist
    lambda_Nitrate::ConstantDist
    RF_mod_unc::Union{UniformDist,ConstantDist}
    mass_to_conc_model_unc::Union{NormalDist,ConstantDist}
    SARF_adjust::Union{TriangDist,ConstantDist}
    NOx_Effect::Union{DiscreteDist,ConstantDist}
    RF_Contrail::Union{TriangDist,ConstantDist}
    RF_Sulfates::Union{TriangDist,ConstantDist}
    RF_Soot::Union{UniformDist,ConstantDist}
    RF_H2O::Union{UniformDist,ConstantDist}
    RF_Nitrate::Union{UniformDist,ConstantDist}
    RF2xCO2::TriangDist
    clim_sen::Union{RoeBakerDist,ConstantDist}
    spec_heat::Union{TriangDist,ConstantDist} #C1
    C2_in::Union{TriangDist,ConstantDist}
    advect::Union{TriangDist,ConstantDist} #F
    diffuse::Union{UniformDist,ConstantDist} #Kz
    z_in::Union{TriangDist,ConstantDist}
    a19::ConstantDist
    a29::Union{NormalDist,ConstantDist}
    a39::ConstantDist
    a19_US::ConstantDist
    a29_US::ConstantDist
    aUS_unc_fact::ConstantDist
    aUS_unc_add::ConstantDist
    aUS_tempAdjust::ConstantDist
    aUS_distrib::Union{NormalDist,ConstantDist}
end

Base.@kwdef mutable struct LCDistributions <: MCDistributions 
    emis_unc_lc_CO2::ConstantDist
    emis_unc_lc_CH4::ConstantDist
    emis_unc_lc_N2O::ConstantDist
    emis_unc_lc_ff_frac::ConstantDist
    co2init::NormalDist
    ch4init::NormalDist
    n2oinit::NormalDist
    RF_mod_unc_CO2::Union{NormalDist,ConstantDist}
    RF_mod_unc_CH4::Union{NormalDist,ConstantDist}
    RF_mod_unc_N2O::Union{NormalDist,ConstantDist}
    RF_mod_unc_CH4tropO3::Union{NormalDist,ConstantDist}
    RF_mod_unc_CH4stratH2O::Union{NormalDist,ConstantDist}
    RF_mod_unc_CH4ffCO2::Union{NormalDist,ConstantDist}
    tau_CH4::Union{NormalDist,ConstantDist}
    tau_N2O::Union{NormalDist,ConstantDist}
    tau_N2O_strat_delay::ConstantDist
    lambda_lc_CH4::ConstantDist
    lambda_lc_N2O::ConstantDist
    lambda_lc_CH4stratH2O::ConstantDist
    lambda_lc_CH4tropO3::ConstantDist
end

Base.@kwdef mutable struct CombustionDistributions <: MCDistributions
    #combustion multipliers
    H2O::ConstantDist #kg/kg
    SOx::ConstantDist #kg/kg
    AIC::ConstantDist #kg/kg
    CO2::ConstantDist #kg/kg
    NOx::Union{UniformDist,ConstantDist} #kg/kg
    BC::Union{UniformDist,ConstantDist} #kg/kg
    CO2_lcycle::ConstantDist #kg CO2/kg fuelburn
end


function MCDistributions(simulation::Simulations_orig; lc_analysis=false)
    config_file_path = simulation.config_file
    distrib_dict = YAML.load_file(config_file_path)
    nrun = simulation.nrun
    s = skip(SobolSeq(1),rand(1:1000000))
    S = [next!(s)[1] for i in 1:nrun]

    maindistribution =  MainDistributions(
        emis_unc_FB = generate_RVs(distrib_dict["emis_unc_FB"], shuffle(S)),
        emis_unc_NOx = generate_RVs(distrib_dict["emis_unc_NOx"], shuffle(S)),
        avNOx2006 = generate_RVs(distrib_dict["avNOx2006"], shuffle(S)),
        avfuel2005 = generate_RVs(distrib_dict["avfuel2005"], shuffle(S)),
        lambda_O3 = generate_RVs(distrib_dict["lambda_O3"], shuffle(S)),
        lambda_Contrail = generate_RVs(distrib_dict["lambda_Contrail"], shuffle(S)),
        lambda_Sulfates = generate_RVs(distrib_dict["lambda_Sulfates"], shuffle(S)),
        lambda_Soot = generate_RVs(distrib_dict["lambda_Soot"], shuffle(S)),
        lambda_H2O = generate_RVs(distrib_dict["lambda_H2O"], shuffle(S)),
        lambda_CH4 = generate_RVs(distrib_dict["lambda_CH4"], shuffle(S)),
        lambda_Nitrate = generate_RVs(distrib_dict["lambda_Nitrate"], shuffle(S)),
        RF_mod_unc = generate_RVs(distrib_dict["RF_mod_unc"], shuffle(S)),
        mass_to_conc_model_unc = generate_RVs(distrib_dict["mass_to_conc_model_unc"], shuffle(S)),
        SARF_adjust = generate_RVs(distrib_dict["SARF_adjust"], shuffle(S)),
        NOx_Effect = generate_RVs(distrib_dict["NOx_Effect"], shuffle(S)),
        RF_Contrail = generate_RVs(distrib_dict["RF_Contrail"], shuffle(S)),
        RF_Sulfates = generate_RVs(distrib_dict["RF_Sulfates"], shuffle(S)),
        RF_Soot = generate_RVs(distrib_dict["RF_Soot"], shuffle(S)),
        RF_H2O = generate_RVs(distrib_dict["RF_H2O"], shuffle(S)),
        RF_Nitrate = generate_RVs(distrib_dict["RF_Nitrate"], shuffle(S)),
        RF2xCO2 = generate_RVs(distrib_dict["RF2xCO2"], shuffle(S)),
        clim_sen = generate_RVs(distrib_dict["clim_sen"], shuffle(S)),
        spec_heat = generate_RVs(distrib_dict["spec_heat"], shuffle(S)),
        advect = generate_RVs(distrib_dict["advect"], shuffle(S)),
        C2_in = generate_RVs(distrib_dict["C2_in"], shuffle(S)),
        z_in = generate_RVs(distrib_dict["z_in"], shuffle(S)),
        diffuse = generate_RVs(distrib_dict["diffuse"], shuffle(S)),
        a19 = generate_RVs(distrib_dict["a19"], shuffle(S)),
        a29 = generate_RVs(distrib_dict["a29"], shuffle(S)),
        a39 = generate_RVs(distrib_dict["a39"], shuffle(S)),
        a19_US = generate_RVs(distrib_dict["a19_US"], shuffle(S)),
        a29_US = generate_RVs(distrib_dict["a29_US"], shuffle(S)),
        aUS_unc_fact = generate_RVs(distrib_dict["aUS_unc_fact"], shuffle(S)),
        aUS_unc_add = generate_RVs(distrib_dict["aUS_unc_add"], shuffle(S)),
        aUS_tempAdjust = generate_RVs(distrib_dict["aUS_tempAdjust"], shuffle(S)),
        aUS_distrib = generate_RVs(distrib_dict["aUS_distrib"], shuffle(S)),
    )

    #Combustion multipliers
    fueltype = simulation.fueltype
    Lcycle_GHG = [
        18.5657, 20.5199, -36.9515, -34.659, 43.4093, -46.6048, -53.6782,-57.2279, 
        109.5289, 45.5751, -40.3128, -31.8015, -19.5247, -8.9867, -31.57, 0.9412,
        ] .*1e-3 #kg per MJ
    energ_dens = 43.2 #MJ/kg conventional jet    
    
    if fueltype <= 4
        H2O_mult = ConstantDist(1)
        AIC_mult = ConstantDist(1)
        CO2_mult = ConstantDist(1)
        NOx_mult = ConstantDist(1)
        BC_mult = ConstantDist(1)
        if isodd(fueltype)
            #Conventional fuel
            SOx_mult = ConstantDist(1)
        else 
            #ULS fuel
            SOx_mult = ConstantDist(15/600); #ratio of ULS sulfur to conventional sulfur
        end

        CO2_lcycle_mult = ConstantDist(Lcycle_GHG[fueltype] * energ_dens)

    elseif fueltype >= 5
        #all SPK fuels
        H2O_mult = generate_RVs(distrib_dict["SPK_H2O"], shuffle(S))
        AIC_mult = generate_RVs(distrib_dict["AIC"], shuffle(S))
        CO2_mult = generate_RVs(distrib_dict["CO2"], shuffle(S))
        NOx_mult = generate_RVs(distrib_dict["NOx"], shuffle(S))
        BC_mult = generate_RVs(distrib_dict["Soot"], shuffle(S))
        SOx_mult = generate_RVs(distrib_dict["SOx"], shuffle(S))

        energ_dens = energ_dens*0.963 #MJ/kg SPK
        CO2_lcycle_mult = ConstantDist(Lcycle_GHG[fueltype] * energ_dens)
    end

    combustdistribution = CombustionDistributions(
        H2O = H2O_mult,
        SOx = SOx_mult,
        AIC = AIC_mult,
        CO2 = CO2_mult,
        NOx = NOx_mult,
        BC = BC_mult,
        CO2_lcycle = CO2_lcycle_mult,
    )
    
    if !lc_analysis
        return maindistribution, combustdistribution
    else
        #LC analysis variables
        lcdistribution = LCDistributions(
            emis_unc_lc_CO2 = generate_RVs(distrib_dict["emis_unc_lc_CO2"], shuffle(S)),
            emis_unc_lc_CH4 = generate_RVs(distrib_dict["emis_unc_lc_CH4"], shuffle(S)),
            emis_unc_lc_N2O = generate_RVs(distrib_dict["emis_unc_lc_N2O"], shuffle(S)),
            emis_unc_lc_ff_frac = generate_RVs(distrib_dict["emis_unc_lc_ff_frac"], shuffle(S)),
            co2init = generate_RVs(distrib_dict["co2init"], shuffle(S)),
            ch4init = generate_RVs(distrib_dict["ch4init"], shuffle(S)),
            n2oinit = generate_RVs(distrib_dict["n2oinit"], shuffle(S)),
            RF_mod_unc_CO2 = generate_RVs(distrib_dict["RF_mod_unc_CO2"], shuffle(S)),
            RF_mod_unc_CH4 = generate_RVs(distrib_dict["RF_mod_unc_CH4"], shuffle(S)),
            RF_mod_unc_N2O = generate_RVs(distrib_dict["RF_mod_unc_N2O"], shuffle(S)),
            RF_mod_unc_CH4tropO3 = generate_RVs(distrib_dict["RF_mod_unc_CH4tropO3"], shuffle(S)),
            RF_mod_unc_CH4stratH2O = generate_RVs(distrib_dict["RF_mod_unc_CH4stratH2O"], shuffle(S)),
            RF_mod_unc_CH4ffCO2 = generate_RVs(distrib_dict["RF_mod_unc_CH4ffCO2"], shuffle(S)),
            tau_CH4 = generate_RVs(distrib_dict["tau_CH4"], shuffle(S)),
            tau_N2O = generate_RVs(distrib_dict["tau_N2O"], shuffle(S)),
            tau_N2O_strat_delay = generate_RVs(distrib_dict["tau_N2O_strat_delay"], shuffle(S)),
            lambda_lc_CH4 = generate_RVs(distrib_dict["lambda_lc_CH4"], shuffle(S)),
            lambda_lc_N2O = generate_RVs(distrib_dict["lambda_lc_N2O"], shuffle(S)),
            lambda_lc_CH4stratH2O = generate_RVs(distrib_dict["lambda_lc_CH4stratH2O"], shuffle(S)),
            lambda_lc_CH4tropO3 = generate_RVs(distrib_dict["lambda_lc_CH4tropO3"], shuffle(S)),
        )
        return maindistribution, combustdistribution, lcdistribution
    end

end