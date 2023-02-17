include("RF_GHG.jl")
include("RF_CO2.jl")

# @views
function calculate_RF_lc(i::Int64, emission::Emissions_orig , lcdistribution::LCDistributions, maindistribution::MainDistributions, simulation::Simulations_orig, climateconsts::ClimateConstants)
    totalimpact = simulation.totalimpact
    t0 = simulation.t0
    interp_years = simulation.interp_years

    #Define relevant variables
    c_over_co2_mass = 12/44
    ch4init = sample(lcdistribution.ch4init,i)

    mass_to_conc_model_unc = sample(maindistribution.mass_to_conc_model_unc,i)
    SARF_adjust = sample(maindistribution.SARF_adjust,i)

    tau_CH4 = sample(lcdistribution.tau_CH4,i)
    tau_N2O = sample(lcdistribution.tau_N2O,i)
    tau_N2O_strat_delay = sample(lcdistribution.tau_N2O_strat_delay,i)

    emis_unc_lc_CO2 = sample(lcdistribution.emis_unc_lc_CO2,i)
    emis_unc_lc_CH4 = sample(lcdistribution.emis_unc_lc_CH4,i)
    emis_unc_lc_N2O = sample(lcdistribution.emis_unc_lc_N2O,i)
    emis_unc_lc_ff_frac = sample(lcdistribution.emis_unc_lc_ff_frac,i)

    RF_mod_unc_CO2 = sample(lcdistribution.RF_mod_unc_CO2, i)
    RF_mod_unc_CH4 = sample(lcdistribution.RF_mod_unc_CH4, i)
    RF_mod_unc_N2O = sample(lcdistribution.RF_mod_unc_N2O, i)
    RF_mod_unc_CH4tropO3 = sample(lcdistribution.RF_mod_unc_CH4tropO3, i)
    RF_mod_unc_CH4stratH2O = sample(lcdistribution.RF_mod_unc_CH4stratH2O, i)
    RF_mod_unc_CH4ffCO2 = sample(lcdistribution.RF_mod_unc_CH4ffCO2, i)

    lambda_lc_CH4 = sample(lcdistribution.lambda_lc_CH4, i)
    lambda_lc_N2O = sample(lcdistribution.lambda_lc_N2O, i)
    lambda_lc_CH4stratH2O = sample(lcdistribution.lambda_lc_CH4stratH2O, i)
    lambda_lc_CH4tropO3 = sample(lcdistribution.lambda_lc_CH4tropO3, i)

    #Load background concentrations and RFs
    X_CO2_bg = climateconsts.co2_conc_SSP #CO2 concentration in ppm
    X_N2O_bg = climateconsts.n2o_conc_SSP #N2O concentration ppb
    X_CH4_bg = climateconsts.ch4_conc_SSP #CH4 concentration in ppb
    SARF_allCO2 = climateconsts.SARF_allCO2 #W/m^2 - stratospherically adjusted radiative forcing for all of anthropogenic CO2 emissions
    SARF_allN2O = climateconsts.SARF_allN2O #W/m^2 - stratospherically adjusted radiative forcing for all of anthropogenic N2O emissions
    SARF_allCH4 = climateconsts.SARF_allCH4 #W/m^2 - stratospherically adjusted radiative forcing for all of anthropogenic CH4 emissions

    #Emissions to concentrations
        # Calculate decay vectors
    CH4response = zeros(totalimpact, 1);
    N2Oresponse = zeros(totalimpact, 1);
    for yr in 1:totalimpact #years after emission
        CH4response[yr] = exp(-(yr-1)/tau_CH4) #fraction remaining
        N2Oresponse[yr] = exp(-(yr-1)/tau_N2O) #fraction remaining
    end

        # Add the N2O stratospheric decay delay used in MAGICC6
    N2Oresponse = cat(ones(Int(round(tau_N2O_strat_delay))), N2Oresponse[1:end-Int(tau_N2O_strat_delay)], dims=1)

        # Compute CO2 emissions due to fossil fuel CH4
    CH4_ff_response = 1 .- CH4response
    CH4_ff_response = cat([0], (CH4_ff_response[2:end] .- CH4_ff_response[1:end-1]), dims=1)
    CH4_ff_emis_mat = LowerTriangular(ones(totalimpact, totalimpact))
    for yr in 1:totalimpact
        CH4_ff_emis_mat[yr:end,yr] = CH4_ff_emis_mat[yr:end,yr] .* emission.lc_CH4_ff[yr]
    end
    kgCH4_kgCO2 = 44/16; 
    lc_CH4_ff_CO2_emis = vec(sum(CH4_ff_emis_mat,dims=2)) .* kgCH4_kgCO2 #kg CO2 emissions due to fossil fuel CH4 oxidized to CO2


        # Calculate matrix of linear IRF
    CH4response_mat = LowerTriangular(ones(totalimpact, totalimpact));
    CH4response_mat[:,1] = CH4response;
    N2Oresponse_mat = LowerTriangular(ones(totalimpact, totalimpact));
    N2Oresponse_mat[:,1] = N2Oresponse;
    for yr in 2:totalimpact
        CH4response_mat[yr:end,yr] = CH4response_mat[yr:end,yr] .* CH4response[1:end-yr+1]
        N2Oresponse_mat[yr:end,yr] = N2Oresponse_mat[yr:end,yr] .* N2Oresponse[1:end-yr+1]
    end

        # Conver emissions to concentrations
    ppb_over_kgCH4 = 1e-9 / 2.78 #ppb/kg CH4 (2.78 ppb to Tg CH4 from MAGICC6 - Meinshausen et al. 2011)
    ppb_over_kgN2O = 1e-9 * (28/4) / 4.81 #ppb/kg N2O (4.81 ppb to Tg N2O from MAGICC6 - Meinshausen et al. 2011)
    CH4response_mat = ppb_over_kgCH4 .* CH4response_mat #ppb/kg CH4 fraction remaining
    N2Oresponse_mat = ppb_over_kgN2O .* N2Oresponse_mat #ppb/kg N2O fraction remaining
    CO2response_mat = getco2response_mat(mass_to_conc_model_unc, climateconsts)

    lc_CO2emis = (emission.lc_CO2 .* emis_unc_lc_CO2) .* 1e-9 .* c_over_co2_mass #convert kgCO2 to TgC
    lc_CH4emis = (emission.lc_CH4 .* emis_unc_lc_CH4) #kg CH4
    lc_N2Oemis = (emission.lc_N2O .* emis_unc_lc_N2O) #kg N2O
    lc_CH4_ff_CO2_emis =  (lc_CH4_ff_CO2_emis.* emis_unc_lc_ff_frac) .* 1e-9 .* c_over_co2_mass #convert kgCO2 to TgC
        
        #Change in concentrations from lc emissions
    deltaX_lc_CO2 = vec(sum(lc_CO2emis[1:length(simulation.emissionyears)]' .* CO2response_mat, dims=2) ./ 2120) #convert TgC to ppm: (GtC)=2.120*Carbon (ppm), so we divide by 2120
    deltaX_lc_CH4 = vec(sum(lc_CH4emis[1:length(simulation.emissionyears)]' .* CH4response_mat, dims=2))
    deltaX_lc_N2O = vec(sum(lc_N2Oemis[1:length(simulation.emissionyears)]' .* N2Oresponse_mat, dims=2))
    deltaX_lc_CH4_CO2 = vec(sum(lc_CH4_ff_CO2_emis[1:length(simulation.emissionyears)]' .* CO2response_mat, dims=2) ./ 2120) 

        #Concentrations without lc emissions
    X_CO2_without_lc = X_CO2_bg .- deltaX_lc_CO2
    X_CH4_without_lc = X_CH4_bg .- deltaX_lc_CH4
    X_N2O_without_lc = X_N2O_bg .- deltaX_lc_N2O
    X_CH4_CO2_without_lc = X_CO2_bg .- deltaX_lc_CH4_CO2

    #Compute RFs
    SARF_CO2_without_lc = calculateCO2SARF.(X_CO2_without_lc, X_N2O_without_lc)
    SARF_CH4_without_lc = calculateCH4SARF.(X_N2O_without_lc, X_CH4_without_lc)
    SARF_N2O_without_lc = calculateN2OSARF.(X_CO2_without_lc, X_N2O_without_lc, X_CH4_without_lc)
    SARF_CH4_CO2_without_lc = calculateCO2SARF.(X_CH4_CO2_without_lc, X_N2O_without_lc)

    RF_lc_CO2 = (SARF_allCO2 .- SARF_CO2_without_lc) .* RF_mod_unc_CO2
    RF_lc_CH4 = (SARF_allCH4 .- SARF_CH4_without_lc) .* RF_mod_unc_CH4
    RF_lc_N2O = (SARF_allN2O .- SARF_N2O_without_lc) .* RF_mod_unc_N2O
    RF_lc_CH4_CO2 = (SARF_allCO2 .- SARF_CH4_CO2_without_lc) .* RF_mod_unc_CH4ffCO2

        #tropospheric ozone (Meinshausen et al. 2011, and comment from MAGICC6 code)
    RF_lc_CH4tropO3 = (RF_mod_unc_CH4tropO3 * 0.0335 * 5) .* ( log.(X_CH4_bg./ch4init) .- log.(X_CH4_without_lc./ch4init) )

        #stratospheric H2O (Meinshausen et al. 2011)
    β_CH4_H2O = 0.15;
    α_CH4 = 0.036;
    RF_lc_CH4stratH2O = (RF_mod_unc_CH4stratH2O * β_CH4_H2O * α_CH4) .* ( sqrt.(X_CH4_bg) .- sqrt.(X_CH4_without_lc) )
        
    return RF_lc_CO2, RF_lc_N2O, RF_lc_CH4, RF_lc_CH4tropO3, RF_lc_CH4stratH2O, RF_lc_CH4_CO2
end