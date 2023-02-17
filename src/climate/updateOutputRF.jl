function update_output_RFs(i::Int64, climateoutput::ClimateOutputs, simulation::Simulations_adv, emission::Emissions, contrail_region_mask::ContrailMaps, climatedistribution::ClimateDistributions, climateconsts::ClimateConstants)
    #Advanced approach
    #CO2
    CO2_inv = totalEmisSpecies("CO2", emission, simulation);
    RF_CO2av = calculate_RF_CO2(i, CO2_inv, climatedistribution, simulation, climateconsts);
    climateoutput.RF_CO2av[:,i] = RF_CO2av;

    #Contrail
    λ_Contrail = sample(climatedistribution.lambda_Contrail,i); 
    contrail_sens_map = update_contrail_map(i, contrail_region_mask, climatedistribution);
    RF_Contrail = calculateRF_Contrail(λ_Contrail, contrail_sens_map, simulation, emission);

    climateoutput.RF_Contrail[:,i] = RF_Contrail;
    return climateoutput
end

function update_output_RFs(i::Int64, mainoutput::MainOutputs, simulation::Simulations_orig, emission::Emissions_orig, combustdistribution::CombustionDistributions, maindistribution::MainDistributions, climateconsts::ClimateConstants)
    #Original approach
    #CO2
    CO2_inv = emission.CO2 .* sample(combustdistribution.CO2,i) .* sample(maindistribution.emis_unc_FB,i)
    RF_CO2av = calculate_RF_CO2(i, CO2_inv, maindistribution, simulation, climateconsts);
    mainoutput.RF_CO2av[:,i] = RF_CO2av;

    CO2_inv_lcycle = CO2_inv + emission.FB .* sample(combustdistribution.CO2_lcycle,i) .* sample(maindistribution.emis_unc_FB,i)
    RF_CO2av_lcycle = calculate_RF_CO2(i, CO2_inv_lcycle, maindistribution, simulation, climateconsts);
    mainoutput.RF_CO2av_lcycle[:,i] = RF_CO2av_lcycle;

    #NOx
    RF_Nitrate, RF_O3short, RF_O3long, RF_CH4 = calculate_RF_NOx(i, emission, maindistribution, simulation);
    mainoutput.RF_Nitrate[:,i] = RF_Nitrate;
    mainoutput.RF_O3short[:,i] = RF_O3short;
    mainoutput.RF_O3long[:,i] = RF_O3long;
    mainoutput.RF_CH4[:,i] = RF_CH4;

    #Short-term (non-CO2, non-NOx)
    RF_Contrail, RF_H2O, RF_BC, RF_S = calculate_RF_short(i, emission, maindistribution, simulation);
    mainoutput.RF_Contrail[:,i] = RF_Contrail;
    mainoutput.RF_H2O[:,i] = RF_H2O;
    mainoutput.RF_BC[:,i] = RF_BC;
    mainoutput.RF_S[:,i] = RF_S;

    return mainoutput
end

function update_output_RFs_lc(i::Int64, lcoutput::LCOutputs, simulation::Simulations_orig, emission::Emissions_orig, lcdistribution::LCDistributions, maindistribution::MainDistributions, climateconsts::ClimateConstants)
    RF_lc_CO2, RF_lc_N2O, RF_lc_CH4, RF_lc_CH4tropO3, RF_lc_CH4stratH2O, RF_lc_CH4_CO2 =calculate_RF_lc(i, emission, lcdistribution, maindistribution, simulation, climateconsts);
    lcoutput.RF_lc_CO2[:,i] = RF_lc_CO2
    lcoutput.RF_lc_N2O[:,i] = RF_lc_N2O
    lcoutput.RF_lc_CH4[:,i] = RF_lc_CH4
    lcoutput.RF_lc_CH4tropO3[:,i] = RF_lc_CH4tropO3
    lcoutput.RF_lc_CH4stratH2O[:,i] = RF_lc_CH4stratH2O
    lcoutput.RF_lc_CH4_CO2[:,i] = RF_lc_CH4_CO2

    return lcoutput
end
    
    