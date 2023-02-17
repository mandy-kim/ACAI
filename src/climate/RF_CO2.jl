include("RF_GHG.jl")

function getco2response_mat(mass_to_conc_model_unc, climateconsts::ClimateConstants)
    co2response_mat = climateconsts.IRF_SSP
    ##Add uncertainty to co2response (IRF) matrix
    year_mat = climateconsts.year_mat
    co2response_mat = mass_to_conc_model_unc .* (1 .- exp.(-year_mat./10)) + co2response_mat
    ###Values should be in range[0,1]
    replace!(x -> x>1 ? 0 : x, co2response_mat)    
    replace!(x -> x<0 ? 0 : x, co2response_mat)

    return co2response_mat
end

"""
# RF_CO2 - compute RFs for total co2 and aviation co2 from co2 emissions
### Inputs
- i: iteration of for loop
- CO2_inv: CO2 emission inventory computed from totalCO2Emis function; length = simulation.totalimpact. CO2_inv[1] = CO2 emissions in year 1
- distribution: MCDistributions type variable to sample variables from relevant to CO2 RF calculations
- simulation: Simulation type variable containing fixed data from user input
- climateconsts: ClimateConstants type variable containing info for fixed data e.g., background CO2 concentration

### Outputs
- RF_avCO2: RF over years (length = simulation.totalimpact) from aviation's CO2 contribution
"""
function calculate_RF_CO2(i::Int64, CO2_inv::Array{Float64, 1}, distribution::MCDistributions, simulation::Simulations, climateconsts::ClimateConstants)
    totalimpact = simulation.totalimpact
    t0 = simulation.t0
    interp_years = simulation.interp_years

    #Define relevant variables
    c_over_co2_mass = 12/44
    RF_mod_unc = sample(distribution.RF_mod_unc,i)
    mass_to_conc_model_unc = sample(distribution.mass_to_conc_model_unc,i)
    SARF_adjust = sample(distribution.SARF_adjust,i)

    #Load relevant constants: CO2 IRF and background concentrations based on SSP scenario
    CO2response_mat = getco2response_mat(mass_to_conc_model_unc, climateconsts)
    X_CO2_bg = climateconsts.co2_conc_SSP #CO2 concentration in ppm
    X_N2O_bg = climateconsts.n2o_conc_SSP #N2O concentration ppb
    SARF_allCO2 = climateconsts.SARF_allCO2 #W/m^2 - stratospherically adjusted radiative forcing for all of anthropogenic CO2 emissions    

    #Aviation CO2 emissions -> aviation CO2 concentration
    co2emis_av = CO2_inv .* 1e-9 .* c_over_co2_mass #convert kgCO2 to TgC
    deltaX_CO2sum_av = vec(sum(co2emis_av[1:length(simulation.emissionyears)]' .* CO2response_mat, dims=2) ./ 2120) #convert TgC to ppm: (GtC)=2.120*Carbon (ppm), so we divide by 2120

    #CO2 concentration without aviation
    X_CO2_without_av = X_CO2_bg .- deltaX_CO2sum_av #ppm

    #Concentration -> RF [=] W/m^2
    SARF_CO2_without_av = calculateCO2SARF.(X_CO2_without_av, X_N2O_bg)

    #Impact(aviation CO2) = Impact(total anthropogenic CO2) - Impact(total anthropogenic CO2 - aviation CO2)
    SARF_avCO2 = (SARF_allCO2 - SARF_CO2_without_av) .* RF_mod_unc #add uniform 10% uncertainty (Etminan et al 2016)

    #Add adjustments for SARF -> ERF
    RF_avCO2 = SARF_avCO2 .* (1 + SARF_adjust) #tropospheric adjustment of the SARF from IPCC AR6 (7.3.2.1)

    return RF_avCO2
end