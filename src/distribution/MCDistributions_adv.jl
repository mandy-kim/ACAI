Base.@kwdef mutable struct AQDistributions <: MCDistributions
    #VSL
    VSL_US_1990::Union{WeibullDist,ConstantDist}
    income_elasticity::Union{TriangDist,ConstantDist}
    #O3 LL CRF
    β_resp::Union{TriangDist,ConstantDist}
    β_AC::Union{TriangDist,ConstantDist}
    #PM2.5 LL CRF
    β_cardio::Union{TriangDist,ConstantDist}
    #PM2.5 GEMM CRFs
    θ_IHD::Union{NormalDist,ConstantDist}
    θ_stroke::Union{NormalDist,ConstantDist}
    θ_COPD::Union{NormalDist,ConstantDist}
    θ_LC::Union{NormalDist,ConstantDist}
    θ_LRI::Union{NormalDist,ConstantDist}
end

Base.@kwdef mutable struct ClimateDistributions <: MCDistributions 
    #CO2 RF
    mass_to_conc_model_unc::Union{NormalDist,ConstantDist}
    RF_mod_unc::UniformDist
    SARF_adjust::Union{TriangDist,ConstantDist}
    #Contrail RF
    lambda_Contrail::Union{TriangDist,ConstantDist} #ERF/RF ratio
    contrail_SATL_RF::Union{NormalDist,ConstantDist}
    contrail_CPAC_RF::Union{NormalDist,ConstantDist}
    contrail_CONUS_RF::Union{NormalDist,ConstantDist}
    contrail_MEAST_RF::Union{NormalDist,ConstantDist}
    contrail_NATL_RF::Union{NormalDist,ConstantDist}
    contrail_NPAC_RF::Union{NormalDist,ConstantDist}
    contrail_NAM_RF::Union{NormalDist,ConstantDist}
    contrail_RUS_RF::Union{NormalDist,ConstantDist}
    contrail_SEAS_RF::Union{NormalDist,ConstantDist}
    contrail_WEU_RF::Union{NormalDist,ConstantDist}
    contrail_global_RF::Union{NormalDist,ConstantDist}
    contrail_RF_unc::Union{TriangDist,ConstantDist}
    #Temperature response
    RF2xCO2::TriangDist
    clim_sen::Union{RoeBakerDist,ConstantDist} #ECS
    spec_heat::Union{TriangDist,ConstantDist} #C1
    C2_in::Union{TriangDist,ConstantDist}
    advect::Union{TriangDist,ConstantDist} #F
    diffuse::Union{UniformDist,ConstantDist} #Kz
    z_in::Union{TriangDist,ConstantDist}
    #Global Damage function
    a19::ConstantDist
    a29::Union{NormalDist,ConstantDist}
    a39::ConstantDist
end

function MCDistributions(simulation::Simulations_adv)
    config_file_path = simulation.config_file
    distrib_dict = YAML.load_file(config_file_path)
    nrun = simulation.nrun
    s = skip(SobolSeq(1),rand(1:1000000))
    S = [next!(s)[1] for i in 1:nrun]

    #Climate variables
    lambda_Contrail = generate_RVs(distrib_dict["lambda_Contrail"], shuffle(S))
    a19 = generate_RVs(distrib_dict["a19"], shuffle(S))
    a29 = generate_RVs(distrib_dict["a29"], shuffle(S))
    a39 = generate_RVs(distrib_dict["a39"], shuffle(S))
    RF2xCO2 = generate_RVs(distrib_dict["RF2xCO2"], shuffle(S))
    spec_heat = generate_RVs(distrib_dict["spec_heat"], shuffle(S))
    advect = generate_RVs(distrib_dict["advect"], shuffle(S))
    C2_in = generate_RVs(distrib_dict["C2_in"], shuffle(S))
    z_in = generate_RVs(distrib_dict["z_in"], shuffle(S))
    diffuse = generate_RVs(distrib_dict["diffuse"], shuffle(S))
    RF_mod_unc = generate_RVs(distrib_dict["RF_mod_unc"], shuffle(S))
    mass_to_conc_model_unc = generate_RVs(distrib_dict["mass_to_conc_model_unc"], shuffle(S))
    clim_sen = generate_RVs(distrib_dict["clim_sen"], shuffle(S))
    SARF_adjust = generate_RVs(distrib_dict["SARF_adjust"], shuffle(S))
    contrail_SATL_RF = generate_RVs(distrib_dict["contrail_SATL_RF"], shuffle(S))
    contrail_CPAC_RF = generate_RVs(distrib_dict["contrail_CPAC_RF"], shuffle(S))
    contrail_CONUS_RF = generate_RVs(distrib_dict["contrail_CONUS_RF"], shuffle(S))
    contrail_MEAST_RF = generate_RVs(distrib_dict["contrail_MEAST_RF"], shuffle(S))
    contrail_NATL_RF = generate_RVs(distrib_dict["contrail_NATL_RF"], shuffle(S))
    contrail_NPAC_RF = generate_RVs(distrib_dict["contrail_NPAC_RF"], shuffle(S))
    contrail_NAM_RF = generate_RVs(distrib_dict["contrail_NAM_RF"], shuffle(S))
    contrail_RUS_RF = generate_RVs(distrib_dict["contrail_RUS_RF"], shuffle(S))
    contrail_SEAS_RF = generate_RVs(distrib_dict["contrail_SEAS_RF"], shuffle(S))
    contrail_WEU_RF = generate_RVs(distrib_dict["contrail_WEU_RF"], shuffle(S))
    contrail_global_RF = generate_RVs(distrib_dict["contrail_global_RF"], shuffle(S))
    contrail_RF_unc = generate_RVs(distrib_dict["contrail_RF_unc"], shuffle(S))
    
    #Air quality variables
    β_cardio = generate_RVs(distrib_dict["β_cardio"], shuffle(S))
    β_resp = generate_RVs(distrib_dict["β_resp"], shuffle(S))
    β_AC = generate_RVs(distrib_dict["β_AC"], shuffle(S))
    θ_IHD = generate_RVs(distrib_dict["θ_IHD"], shuffle(S))
    θ_stroke = generate_RVs(distrib_dict["θ_stroke"], shuffle(S))
    θ_COPD = generate_RVs(distrib_dict["θ_COPD"], shuffle(S))
    θ_LC = generate_RVs(distrib_dict["θ_LC"], shuffle(S))
    θ_LRI = generate_RVs(distrib_dict["θ_LRI"], shuffle(S))
    VSL_US_1990 = generate_RVs(distrib_dict["VSL_US_1990"], shuffle(S))
    income_elasticity = generate_RVs(distrib_dict["income_elasticity"], shuffle(S))
    
    climatedistribution = ClimateDistributions(
        lambda_Contrail=lambda_Contrail,
        a19=a19,
        a29=a29,
        a39=a39,
        RF2xCO2=RF2xCO2,
        spec_heat=spec_heat,
        advect=advect,
        C2_in=C2_in,
        z_in=z_in,
        diffuse=diffuse,
        RF_mod_unc=RF_mod_unc,
        clim_sen=clim_sen,
        SARF_adjust=SARF_adjust,
        mass_to_conc_model_unc=mass_to_conc_model_unc,
        contrail_SATL_RF=contrail_SATL_RF,
        contrail_CPAC_RF=contrail_CPAC_RF,
        contrail_CONUS_RF=contrail_CONUS_RF,
        contrail_MEAST_RF=contrail_MEAST_RF,
        contrail_NATL_RF=contrail_NATL_RF,
        contrail_NPAC_RF=contrail_NPAC_RF,
        contrail_NAM_RF=contrail_NAM_RF,
        contrail_global_RF=contrail_global_RF,
        contrail_RUS_RF=contrail_RUS_RF,
        contrail_SEAS_RF=contrail_SEAS_RF,
        contrail_WEU_RF=contrail_WEU_RF,
        contrail_RF_unc=contrail_RF_unc,
    )

    aqdistribution = AQDistributions(
        VSL_US_1990=VSL_US_1990,
        income_elasticity=income_elasticity,
        β_resp=β_resp,
        β_cardio=β_cardio,
        β_AC=β_AC,
        θ_IHD=θ_IHD,
        θ_stroke=θ_stroke,
        θ_COPD=θ_COPD,
        θ_LC=θ_LC,
        θ_LRI=θ_LRI,
    )

    return climatedistribution, aqdistribution
end