function update_output_ΔT_costs(climateoutput::ClimateOutputs, simulation::Simulations_adv, climatedistribution::ClimateDistributions, climateconsts::ClimateConstants)
    #Advanced approach - calculate (only) global damages
    
    #load damage function
    damage_av_func = getDamageAvFunction(simulation, climatedistribution, climateoutput.ΔT_background1);
    GDP_matrix, _ = loadBackgroundGDPs(simulation, climatedistribution, climateconsts);

    #calculate costs
    climateoutput.cost_CO2av = GDP_matrix .* damage_av_func(climateoutput.ΔT_background1 .- climateoutput.ΔT_background2);
    climateoutput.cost_Contrail = GDP_matrix .* damage_av_func(climateoutput.ΔT_Contrail);
    climateoutput.cost_NOx = GDP_matrix .* damage_av_func(climateoutput.ΔT_NOx);
	# climateoutput.cost_H2O = ...
    # climateoutput.cost_S = 
    # climateoutput.cost_BC = 

    #discounted costs
    climateoutput.discountedcost_CO2av = discountedCostClimate(climateoutput.cost_CO2av, simulation);
    climateoutput.discountedcost_Contrail = discountedCostClimate(climateoutput.cost_Contrail, simulation);
    climateoutput.discountedcost_NOx = discountedCostClimate(climateoutput.cost_NOx, simulation);
	# climateoutput.discountedcost_H2O =  ...
    # climateoutput.discountedcost_S = 
    # climateoutput.discountedcost_BC = 

    #NPVs
    climateoutput.NPV_CO2av = NPV(climateoutput.discountedcost_CO2av);
    climateoutput.NPV_Contrail = NPV(climateoutput.discountedcost_Contrail);
    climateoutput.NPV_NOx = NPV(climateoutput.discountedcost_NOx);
    # climateoutput.NPV_H2O =  ...
    # climateoutput.NPV_BC = 
    # climateoutput.NPV_S = 

    return climateoutput
end

function update_output_ΔT_costs(mainoutput::MainOutputs, simulation::Simulations_orig, maindistribution::MainDistributions, climateconsts::ClimateConstants)
    #Original approach - calculate global and US damages

    #load damage functions
    damage_av_func = getDamageAvFunction(simulation, maindistribution, mainoutput.ΔT_background1);
    damage_av_func_US = getDamageAvFunction_US(simulation, maindistribution, mainoutput.ΔT_background1);
    GDP_matrix, GDP_matrix_US = loadBackgroundGDPs(simulation, maindistribution, climateconsts);

    #calculate costs
    mainoutput.cost_CO2av = GDP_matrix .* damage_av_func(mainoutput.ΔT_background1 .- mainoutput.ΔT_background2);
    mainoutput.cost_Contrail = GDP_matrix .* damage_av_func(mainoutput.ΔT_Contrail);
    mainoutput.cost_O3short = GDP_matrix .* damage_av_func(mainoutput.ΔT_O3short);
    mainoutput.cost_O3long = GDP_matrix .* damage_av_func(mainoutput.ΔT_O3long);
    mainoutput.cost_Nitrate = GDP_matrix .* damage_av_func(mainoutput.ΔT_Nitrate);
    mainoutput.cost_CH4 = GDP_matrix .* damage_av_func(mainoutput.ΔT_CH4);
    mainoutput.cost_H2O = GDP_matrix .* damage_av_func(mainoutput.ΔT_H2O);
    mainoutput.cost_BC = GDP_matrix .* damage_av_func(mainoutput.ΔT_BC);
    mainoutput.cost_S = GDP_matrix .* damage_av_func(mainoutput.ΔT_S);
    
    mainoutput.cost_US_CO2av = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_background1 .- mainoutput.ΔT_background2);
    mainoutput.cost_US_Contrail = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_Contrail);
    mainoutput.cost_US_O3short = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_O3short);
    mainoutput.cost_US_O3long = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_O3long);
    mainoutput.cost_US_Nitrate = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_Nitrate);
    mainoutput.cost_US_CH4 = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_CH4);
    mainoutput.cost_US_H2O = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_H2O);
    mainoutput.cost_US_BC = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_BC);
    mainoutput.cost_US_S = GDP_matrix_US .* damage_av_func_US(mainoutput.ΔT_S);
    

    #discounted costs
    mainoutput.discountedcost_CO2av = discountedCostClimate(mainoutput.cost_CO2av, simulation);
    mainoutput.discountedcost_Contrail = discountedCostClimate(mainoutput.cost_Contrail, simulation);
    mainoutput.discountedcost_O3short = discountedCostClimate(mainoutput.cost_O3short, simulation);
    mainoutput.discountedcost_O3long = discountedCostClimate(mainoutput.cost_O3long, simulation);
    mainoutput.discountedcost_Nitrate = discountedCostClimate(mainoutput.cost_Nitrate, simulation);
    mainoutput.discountedcost_CH4 = discountedCostClimate(mainoutput.cost_CH4, simulation);
    mainoutput.discountedcost_H2O = discountedCostClimate(mainoutput.cost_H2O, simulation);
    mainoutput.discountedcost_BC = discountedCostClimate(mainoutput.cost_BC, simulation);
    mainoutput.discountedcost_S = discountedCostClimate(mainoutput.cost_S, simulation);

    mainoutput.discountedcost_US_O3short = discountedCostClimate(mainoutput.cost_US_O3short, simulation);
    mainoutput.discountedcost_US_O3long = discountedCostClimate(mainoutput.cost_US_O3long, simulation);
    mainoutput.discountedcost_US_Nitrate = discountedCostClimate(mainoutput.cost_US_Nitrate, simulation);
    mainoutput.discountedcost_US_CH4 = discountedCostClimate(mainoutput.cost_US_CH4, simulation);
    mainoutput.discountedcost_US_H2O = discountedCostClimate(mainoutput.cost_US_H2O, simulation);
    mainoutput.discountedcost_US_BC = discountedCostClimate(mainoutput.cost_US_BC, simulation);
    mainoutput.discountedcost_US_S = discountedCostClimate(mainoutput.cost_US_S, simulation);


    #NPVs
    mainoutput.NPV_CO2av = NPV(mainoutput.discountedcost_CO2av);
    mainoutput.NPV_Contrail = NPV(mainoutput.discountedcost_Contrail);
    mainoutput.NPV_O3short = NPV(mainoutput.discountedcost_O3short);
    mainoutput.NPV_O3long = NPV(mainoutput.discountedcost_O3long);
    mainoutput.NPV_Nitrate = NPV(mainoutput.discountedcost_Nitrate);
    mainoutput.NPV_CH4 = NPV(mainoutput.discountedcost_CH4);
    mainoutput.NPV_H2O = NPV(mainoutput.discountedcost_H2O);
    mainoutput.NPV_BC = NPV(mainoutput.discountedcost_BC);
    mainoutput.NPV_S = NPV(mainoutput.discountedcost_S);

    mainoutput.NPV_US_CO2av = NPV(mainoutput.discountedcost_US_CO2av);
    mainoutput.NPV_US_Contrail = NPV(mainoutput.discountedcost_US_Contrail);
    mainoutput.NPV_US_O3short = NPV(mainoutput.discountedcost_US_O3short);
    mainoutput.NPV_US_O3long = NPV(mainoutput.discountedcost_US_O3long);
    mainoutput.NPV_US_Nitrate = NPV(mainoutput.discountedcost_US_Nitrate);
    mainoutput.NPV_US_CH4 = NPV(mainoutput.discountedcost_US_CH4);
    mainoutput.NPV_US_H2O = NPV(mainoutput.discountedcost_US_H2O);
    mainoutput.NPV_US_BC = NPV(mainoutput.discountedcost_US_BC);
    mainoutput.NPV_US_S = NPV(mainoutput.discountedcost_US_S);
    
    return mainoutput
end

function update_output_ΔT_costs_lc(lcoutput::LCOutputs, mainoutput::MainOutputs, simulation::Simulations_orig, maindistribution::MainDistributions, climateconsts::ClimateConstants)
    #Original approach - calculate global and US damages for lifecycle analysis

    #load damage functions
    damage_av_func = getDamageAvFunction(simulation, maindistribution, mainoutput.ΔT_background1);
    damage_av_func_US = getDamageAvFunction_US(simulation, maindistribution, mainoutput.ΔT_background1);
    GDP_matrix, GDP_matrix_US = loadBackgroundGDPs(simulation, maindistribution, climateconsts);

    #calculate costs
    lcoutput.cost_lc_CO2 = GDP_matrix .* damage_av_func(lcoutput.ΔT_lc_CO2)
    lcoutput.cost_lc_N2O = GDP_matrix .* damage_av_func(lcoutput.ΔT_lc_N2O)
    lcoutput.cost_lc_CH4 = GDP_matrix .* damage_av_func(lcoutput.ΔT_lc_CH4)
    lcoutput.cost_lc_CH4tropO3 = GDP_matrix .* damage_av_func(lcoutput.ΔT_lc_CH4tropO3)
    lcoutput.cost_lc_CH4stratH2O = GDP_matrix .* damage_av_func(lcoutput.ΔT_lc_CH4stratH2O)
    lcoutput.cost_lc_CH4_CO2 = GDP_matrix .* damage_av_func(lcoutput.ΔT_lc_CH4_CO2)

    lcoutput.cost_US_lc_CO2 = GDP_matrix_US .* damage_av_func_US(lcoutput.ΔT_lc_CO2)
    lcoutput.cost_US_lc_N2O = GDP_matrix_US .* damage_av_func_US(lcoutput.ΔT_lc_N2O)
    lcoutput.cost_US_lc_CH4 = GDP_matrix_US .* damage_av_func_US(lcoutput.ΔT_lc_CH4)
    lcoutput.cost_US_lc_CH4tropO3 = GDP_matrix_US .* damage_av_func_US(lcoutput.ΔT_lc_CH4tropO3)
    lcoutput.cost_US_lc_CH4stratH2O = GDP_matrix_US .* damage_av_func_US(lcoutput.ΔT_lc_CH4stratH2O)
    lcoutput.cost_US_lc_CH4_CO2 = GDP_matrix_US .* damage_av_func_US(lcoutput.ΔT_lc_CH4_CO2)

    #discounted costs
    lcoutput.discountedcost_lc_CO2 = discountedCostClimate(lcoutput.cost_lc_CO2, simulation);
    lcoutput.discountedcost_lc_N2O = discountedCostClimate(lcoutput.cost_lc_N2O , simulation)
    lcoutput.discountedcost_lc_CH4 = discountedCostClimate(lcoutput.cost_lc_CH4 , simulation)
    lcoutput.discountedcost_lc_CH4tropO3 = discountedCostClimate(lcoutput.cost_lc_CH4tropO3 , simulation)
    lcoutput.discountedcost_lc_CH4stratH2O = discountedCostClimate(lcoutput.cost_lc_CH4stratH2O , simulation)
    lcoutput.discountedcost_lc_CH4_CO2 = discountedCostClimate(lcoutput.cost_lc_CH4_CO2 , simulation)

    lcoutput.discountedcost_US_lc_CO2 = discountedCostClimate(lcoutput.cost_US_lc_CO2 , simulation)
    lcoutput.discountedcost_US_lc_N2O = discountedCostClimate(lcoutput.cost_US_lc_N2O , simulation)
    lcoutput.discountedcost_US_lc_CH4 = discountedCostClimate(lcoutput.cost_US_lc_CH4 , simulation)
    lcoutput.discountedcost_US_lc_CH4tropO3 = discountedCostClimate(lcoutput.cost_US_lc_CH4tropO3 , simulation)
    lcoutput.discountedcost_US_lc_CH4stratH2O = discountedCostClimate(lcoutput.cost_US_lc_CH4stratH2O , simulation)
    lcoutput.discountedcost_US_lc_CH4_CO2 = discountedCostClimate(lcoutput.cost_US_lc_CH4_CO2 , simulation)
        
    #NPVs
    lcoutput.NPV_lc_CO2 = NPV(lcoutput.discountedcost_lc_CO2)
    lcoutput.NPV_lc_N2O = NPV(lcoutput.discountedcost_lc_N2O)
    lcoutput.NPV_lc_CH4 = NPV(lcoutput.discountedcost_lc_CH4)
    lcoutput.NPV_lc_CH4tropO3 = NPV(lcoutput.discountedcost_lc_CH4tropO3)
    lcoutput.NPV_lc_CH4stratH2O = NPV(lcoutput.discountedcost_lc_CH4stratH2O)
    lcoutput.NPV_lc_CH4_CO2 = NPV(lcoutput.discountedcost_lc_CH4_CO2)
    
    lcoutput.NPV_US_lc_CO2 = NPV(lcoutput.discountedcost_US_lc_CO2)
    lcoutput.NPV_US_lc_N2O = NPV(lcoutput.discountedcost_US_lc_N2O)
    lcoutput.NPV_US_lc_CH4 = NPV(lcoutput.discountedcost_US_lc_CH4)
    lcoutput.NPV_US_lc_CH4tropO3 = NPV(lcoutput.discountedcost_US_lc_CH4tropO3)
    lcoutput.NPV_US_lc_CH4stratH2O = NPV(lcoutput.discountedcost_US_lc_CH4stratH2O)
    lcoutput.NPV_US_lc_CH4_CO2 = NPV(lcoutput.discountedcost_US_lc_CH4_CO2)

    return lcoutput
end