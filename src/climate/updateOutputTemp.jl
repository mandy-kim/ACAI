function update_output_ΔTs(i::Int64, NOx_impacts::TotalImpacts, climateoutput::ClimateOutputs, simulation::Simulations_adv, climatedistribution::ClimateDistributions, climateconsts::ClimateConstants)
    #Advanced approach
    ΔT_CO2av = deltaT(i, climateoutput.RF_CO2av[:,i], simulation, climatedistribution);
    climateoutput.ΔT_CO2av[:,i] = ΔT_CO2av;
    climateoutput.ΔT_Contrail[:,i] = deltaT(i, climateoutput.RF_Contrail[:,i], simulation, climatedistribution);
    climateoutput.ΔT_NOx[:,i] = deltaT(i, NOx_impacts.RF, simulation, climatedistribution);
    # output.ΔT_H2O[:,i] = 
    # output.ΔT_S[:,i] = 
    # output.ΔT_BC[:,i] = 

    #calculate background ΔT
    ΔT_background = loadBackgroundTemp(i,ΔT_CO2av,simulation,climatedistribution,climateconsts);
    climateoutput.ΔT_background1[:,i] = ΔT_background[:,1];
    climateoutput.ΔT_background2[:,i] = ΔT_background[:,2];

    return climateoutput
end

function update_output_ΔTs(i::Int64, mainoutput::MainOutputs, simulation::Simulations_orig, maindistribution::MainDistributions, climateconsts::ClimateConstants)
    #Original approach
    mainoutput.ΔT_CO2av[:,i] = deltaT(i, mainoutput.RF_CO2av[:,i], simulation, maindistribution);
    mainoutput.ΔT_CO2av_lcycle[:,i] = deltaT(i, mainoutput.RF_CO2av_lcycle[:,i], simulation, maindistribution);
    mainoutput.ΔT_Contrail[:,i] = deltaT(i, mainoutput.RF_Contrail[:,i], simulation, maindistribution);
    mainoutput.ΔT_Nitrate[:,i] = deltaT(i, mainoutput.RF_Nitrate[:,i], simulation, maindistribution);
    mainoutput.ΔT_O3short[:,i] = deltaT(i, mainoutput.RF_O3short[:,i], simulation, maindistribution);
    mainoutput.ΔT_O3long[:,i] = deltaT(i, mainoutput.RF_O3long[:,i], simulation, maindistribution);
    mainoutput.ΔT_CH4[:,i] = deltaT(i, mainoutput.RF_CH4[:,i], simulation, maindistribution);
    mainoutput.ΔT_H2O[:,i] = deltaT(i, mainoutput.RF_H2O[:,i], simulation, maindistribution);
    mainoutput.ΔT_BC[:,i] = deltaT(i, mainoutput.RF_BC[:,i], simulation, maindistribution);
    mainoutput.ΔT_S[:,i] = deltaT(i, mainoutput.RF_S[:,i], simulation, maindistribution);

    #calculate background ΔT
    ΔT_background = loadBackgroundTemp(i,mainoutput.ΔT_CO2av[:,i],simulation,maindistribution,climateconsts);
    ΔT_background_lcycle = loadBackgroundTemp(i,mainoutput.ΔT_CO2av_lcycle[:,i],simulation,maindistribution,climateconsts);
    mainoutput.ΔT_background1[:,i] = ΔT_background[:,1];
    mainoutput.ΔT_background2[:,i] = ΔT_background[:,2];
    #ΔT_background_lcycle[:,1] is the same as ΔT_background[:,1]
    mainoutput.ΔT_background2_lcycle[:,i] = ΔT_background_lcycle[:,2];

    return mainoutput
end

function update_output_ΔTs_lc(i::Int64, lcoutput::LCOutputs, simulation::Simulations_orig, maindistribution::MainDistributions, climateconsts::ClimateConstants)
    lcoutput.ΔT_lc_CO2[:,i] = deltaT(i, lcoutput.RF_lc_CO2[:,i], simulation, maindistribution);
    lcoutput.ΔT_lc_N2O[:,i] = deltaT(i, lcoutput.RF_lc_N2O[:,i], simulation, maindistribution);
    lcoutput.ΔT_lc_CH4[:,i] = deltaT(i, lcoutput.RF_lc_CH4[:,i], simulation, maindistribution);
    lcoutput.ΔT_lc_CH4tropO3[:,i] = deltaT(i, lcoutput.RF_lc_CH4tropO3[:,i], simulation, maindistribution);
    lcoutput.ΔT_lc_CH4stratH2O[:,i] = deltaT(i, lcoutput.RF_lc_CH4stratH2O[:,i], simulation, maindistribution);
    lcoutput.ΔT_lc_CH4_CO2[:,i] = deltaT(i, lcoutput.RF_lc_CH4_CO2[:,i], simulation, maindistribution);

    return lcoutput
end