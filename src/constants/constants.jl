######### CLIMATE CONSTANTS

Base.@kwdef struct ClimateConstants
    year_mat::Array{Float64,2}
    iter_deltaT_SSP::Interpolations.Extrapolation
    IRF_SSP::Array{Float64,2}
    co2_conc_SSP::Array{Float64,1}
    n2o_conc_SSP::Array{Float64,1}
    ch4_conc_SSP::Array{Float64,1}
    SARF_allCO2::Array{Float64,1}
    SARF_allN2O::Array{Float64,1}
    SARF_allCH4::Array{Float64,1}
    GDP_matrix::Array{Float64,2} #units = billions USD
    GDP_matrix_US::Array{Float64,2} #units = billions USD
end


"""
# lowerDiagYearMat - create matrix of size totalimpact x nyears, where nyears = number of years with emissions, diagonals = 0. Used to add uncertainty to CO2 emissions in RF_CO2 function.
"""
function lowerDiagYearMat(totalimpact, nyears)
    A = LowerTriangular(ones(totalimpact,totalimpact))
    A[diagind(A)] .= 0
    year_mat = LowerTriangular(ones(totalimpact,totalimpact)) * A
    return year_mat[:,1:nyears]
end


@views function getClimateConstants(simulation::Simulations)
    year_mat = lowerDiagYearMat(simulation.totalimpact,length(simulation.emissionyears))

    ########## Load MAGICC6 deltaT background
    # assumes load_RCP_temp == 1. 
    # if load_RCP_temp == 0 (pre-v24 legacy), use CO2 RF to compute background temperature change. this would need to be put back in MC runs
    f =  jldopen("data/MAGICC6_SSP/deltaTsum_SSP_MAGICC6.jld2")
    deltaT_SSP_years = f["years"]   #rows
    deltaT_ECS_vals = f["ECS_vals"] #columns
    deltaT_SSP_array = f["deltaT_SSPs"][:,:,simulation.SSPscen]
    close(f)
    iter_deltaT_SSP = LinearInterpolation((deltaT_SSP_years, deltaT_ECS_vals), deltaT_SSP_array)

    ########## Load CO2, N2O, and CH4 background concentrations
    f = jldopen("data/MAGICC6_SSP/co2_conc_SSP.jld2")
    conc_years = f["years"]
    co2_conc_SSP_array = f["co2_conc_SSPs"] #values in ppm
    close(f)
    f = jldopen("data/MAGICC6_SSP/n2o_conc_SSP.jld2")
    n2o_conc_SSP_array = f["n2o_conc_SSPs"] #values in ppb
    close(f)
    f = jldopen("data/MAGICC6_SSP/ch4_conc_SSP.jld2")
    ch4_conc_SSP_array = f["ch4_conc_SSPs"] #values in ppb
    close(f)
    i = findall(x -> x == simulation.t0, conc_years)[1]
    j = findall(x-> x == simulation.interp_years[end], conc_years)[1]
    co2_conc_SSP = co2_conc_SSP_array[i:j, simulation.SSPscen] #background co2 concentration for years of simulation, ppm
    n2o_conc_SSP = n2o_conc_SSP_array[i:j, simulation.SSPscen] #background n2o concentration for years of simulation, ppb
    ch4_conc_SSP = ch4_conc_SSP_array[i:j, simulation.SSPscen] #background ch4 concentration for years of simulation, ppb

    # calculate SARF for all of anthropogenic GHGs
    SARF_allCO2 = calculateCO2SARF.(co2_conc_SSP, n2o_conc_SSP)
    SARF_allN2O = calculateN2OSARF.(co2_conc_SSP, n2o_conc_SSP, ch4_conc_SSP)
    SARF_allCH4 = calculateCH4SARF.(n2o_conc_SSP, ch4_conc_SSP)

    ########## Load CO2 impulse response functions
    f = jldopen("data/MAGICC6_SSP/IRFs_SSP_MAGICC6.jld2")
    IRF_years = f["years"]
    IRF_SSPs = f["IRF_SSPs"]
    close(f)
    i = simulation.t0 - 1999
    j = simulation.t0-2000+simulation.totalimpact
    yrs_of_emis = simulation.emissionyears[end] - simulation.emissionyears[1]
    IRF_SSP = IRF_SSPs[i:j, i:i+yrs_of_emis, simulation.SSPscen]


    ########## Load SSP global GDP forecast data ##########
    f = jldopen("data/GDP_data/GDP_forecast_SSP.jld2")
    SSP_years = f["years"]; #2000-2900
    GDP_global = f["GDP_global"]; #901x5 array
    GDP_US = f["GDP_US"] #901x5 matrix
    close(f)

    # SSP GDP values in 2005 USD, convert to tmoney USD
    f = jldopen("data/GDP_data/GDP_deflator.jld2")
    GDP_def_years = Int.(f["years"])
    GDP_def_vals = float.(f["values"])
    close(f)
    i_2005 = findall(x-> x==2005, GDP_def_years)[1]
    i_tmoney = findall(x-> x==simulation.tmoney, GDP_def_years)[1]
    conversion = GDP_def_vals[i_tmoney] / GDP_def_vals[i_2005]

    # find indices in SSP_years corresponding to interp_years years
    i = findall(x-> x == simulation.t0, SSP_years)[1]
    j = findall(x-> x == simulation.interp_years[end], SSP_years)[1]

    GDP_global = GDP_global[i:j,simulation.SSP] * conversion
    GDP_US = GDP_US[i:j,simulation.SSP] * conversion

    nrun=simulation.nrun
    GDP_matrix = repeat(GDP_global,1,nrun)
    GDP_matrix_US = repeat(GDP_US,1,nrun) 

    return ClimateConstants(
        year_mat=year_mat,
        iter_deltaT_SSP=iter_deltaT_SSP,
        IRF_SSP=IRF_SSP,
        co2_conc_SSP=co2_conc_SSP,
        n2o_conc_SSP=n2o_conc_SSP,
        ch4_conc_SSP=ch4_conc_SSP,
        SARF_allCO2=SARF_allCO2,
        SARF_allN2O=SARF_allN2O,
        SARF_allCH4=SARF_allCH4,
        GDP_matrix=GDP_matrix,
        GDP_matrix_US=GDP_matrix_US,
        )
end

function getCountryList()
    countries = [
        "AFG", "AGO", "ALB", "AND", "ARE", "ARG", "ARM", "ASM", "ATG", "AUS", "AUT", "AZE", 
        "BDI", "BEL", "BEN", "BFA", "BGD", "BGR", "BHR", "BHS", "BIH", "BLR", "BLZ", "BMU", "BOL", "BRA", "BRB", "BRN", "BTN", "BWA", 
        "CAF", "CAN", "CHE", "CHL", "CHN", "CIV", "CMR", "COD", "COG", "COK", "COL", "COM", "CPV", "CRI", "CUB", "CYP", "CZE", 
        "DEU", "DJI", "DMA", "DNK", "DOM", "DZA", 
        "ECU", "EGY", "ERI", "ESP", "EST", "ETH", 
        "FIN", "FJI", "FRA", "FSM", 
        "GAB", "GBR", "GEO", "GHA", "GIN", "GMB", "GNB", "GNQ", "GRC", "GRD", "GRL", "GTM", "GUM", "GUY", 
        "HND", "HRV", "HTI", "HUN", 
        "IDN", "IND", "IRL", "IRN", "IRQ", "ISL", "ISR", "ITA", 
        "JAM", "JOR", "JPN", 
        "KAZ", "KEN", "KGZ", "KHM", "KIR", "KNA", "KOR", "KWT", 
        "LAO", "LBN", "LBR", "LBY", "LCA", "LKA", "LSO", "LTU", "LUX", "LVA", 
        "MAR", "MCO", "MDA", "MDG", "MDV", "MEX", "MHL", "MKD", "MLI", "MLT", "MMR", "MNE", "MNG", "MNP", "MOZ", "MRT", "MUS", "MWI", "MYS", 
        "NAM", "NER", "NGA", "NIC", "NIU", "NLD", "NOR", "NPL", "NRU", "NZL", 
        "OMN", 
        "PAK", "PAN", "PER", "PHL", "PLW", "PNG", "POL", "PRI", "PRK", "PRT", "PRY", "PSE", 
        "QAT", 
        "ROU", "RUS", "RWA", 
        "SAU", "SDN", "SEN", "SGP", "SLB", "SLE", "SLV", "SMR", "SOM", "SRB", "SSD", "STP", "SUR", "SVK", "SVN", "SWE", "SWZ", "SYC", "SYR", 
        "TCD", "TGO", "THA", "TJK", "TKL", "TKM", "TLS", "TON", "TTO", "TUN", "TUR", "TUV", "TWN", "TZA", 
        "UGA", "UKR", "URY", "USA", "UZB", 
        "VCT", "VEN", "VIR", "VNM", "VUT", 
        "WSM", 
        "YEM", 
        "ZAF", "ZMB", "ZWE", 
        ]
    return countries
end

function getCountryFullNameList()
    #same order as countrylist (ISO codes), ordered by ISO codes
    countries = ["Afghanistan",
                "Angola",
                "Albania",
                "Andorra",
                "United Arab Emirates",
                "Argentina",
                "Armenia",
                "American Samoa",
                "Antigua and Barbuda",
                "Australia",
                "Austria",
                "Azerbaijan",
                "Burundi",
                "Belgium",
                "Benin",
                "Burkina Faso",
                "Bangladesh",
                "Bulgaria",
                "Bahrain",
                "Bahamas",
                "Bosnia and Herzegovina",
                "Belarus",
                "Belize",
                "Bermuda",
                "Bolivia",
                "Brazil",
                "Barbados",
                "Brunei Darussalam",
                "Bhutan",
                "Botswana",
                "Central African Republic",
                "Canada",
                "Switzerland",
                "Chile",
                "China",
                "Côte d'Ivoire",
                "Cameroon",
                "Dem Rep of Congo",
                "Congo",
                "Cook Islands",
                "Colombia",
                "Comoros",
                "Cabo Verde",
                "Costa Rica",
                "Cuba",
                "Cyprus",
                "Czechia",
                "Germany",
                "Djibouti",
                "Dominica",
                "Denmark",
                "Dominican Republic",
                "Algeria",
                "Ecuador",
                "Egypt",
                "Eritrea",
                "Spain",
                "Estonia",
                "Ethiopia",
                "Finland",
                "Fiji",
                "France",
                "Micronesia",
                "Gabon",
                "United Kingdom",
                "Georgia",
                "Ghana",
                "Guinea",
                "Gambia",
                "Guinea-Bissau",
                "Equatorial Guinea",
                "Greece",
                "Grenada",
                "Greenland",
                "Guatemala",
                "Guam",
                "Guyana",
                "Honduras",
                "Croatia",
                "Haiti",
                "Hungary",
                "Indonesia",
                "India",
                "Ireland",
                "Iran",
                "Iraq",
                "Iceland",
                "Israel",
                "Italy",
                "Jamaica",
                "Jordan",
                "Japan",
                "Kazakhstan",
                "Kenya",
                "Kyrgyzstan",
                "Cambodia",
                "Kiribati",
                "Saint Kitts and Nevis",
                "South Korea",
                "Kuwait",
                "Laos",
                "Lebanon",
                "Liberia",
                "Libya",
                "Saint Lucia",
                "Sri Lanka",
                "Lesotho",
                "Lithuania",
                "Luxembourg",
                "Latvia",
                "Morocco",
                "Monaco",
                "Republic of Moldova",
                "Madagascar",
                "Maldives",
                "Mexico",
                "Marshall Islands",
                "North Macedonia",
                "Mali",
                "Malta",
                "Myanmar",
                "Montenegro",
                "Mongolia",
                "Northern Mariana Islands",
                "Mozambique",
                "Mauritania",
                "Mauritius",
                "Malawi",
                "Malaysia",
                "Namibia",
                "Niger",
                "Nigeria",
                "Nicaragua",
                "Niue",
                "Netherlands",
                "Norway",
                "Nepal",
                "Nauru",
                "New Zealand",
                "Oman",
                "Pakistan",
                "Panama",
                "Peru",
                "Philippines",
                "Palau",
                "Papua New Guinea",
                "Poland",
                "Puerto Rico",
                "North Korea",
                "Portugal",
                "Paraguay",
                "Palestine",
                "Qatar",
                "Romania",
                "Russian Federation",
                "Rwanda",
                "Saudi Arabia",
                "Sudan",
                "Senegal",
                "Singapore",
                "Solomon Islands",
                "Sierra Leone",
                "El Salvador",
                "San Marino",
                "Somalia",
                "Serbia",
                "South Sudan",
                "Sao Tome and Principe",
                "Suriname",
                "Slovakia",
                "Slovenia",
                "Sweden",
                "Eswatini",
                "Seychelles",
                "Syrian Arab Republic",
                "Chad",
                "Togo",
                "Thailand",
                "Tajikistan",
                "Tokelau",
                "Turkmenistan",
                "Timor-Leste",
                "Tonga",
                "Trinidad and Tobago",
                "Tunisia",
                "Turkey",
                "Tuvalu",
                "Taiwan",
                "Tanzania",
                "Uganda",
                "Ukraine",
                "Uruguay",
                "USA",
                "Uzbekistan",
                "Saint Vincent and the Grenadines",
                "Venezuela",
                "U.S. Virgin Islands",
                "Vietnam",
                "Vanuatu",
                "Samoa",
                "Yemen",
                "South Africa",
                "Zambia",
                "Zimbabwe",
                ]
    return countries
end

function getPopulationGriddata(countries::Array{String,1})
    pop_dict = Dict{String,Array{Float64,2}}()
    pop_idx = Dict{String, Array{CartesianIndex{2},1}}()
    pop_totals = zeros(length(countries))
    pop_global = zeros(144,91)
    for i in 1:length(countries)
        country = countries[i]
        pop_file = NCDataset("data/population/grid/Population_2020_2x2.5_PC_DC_"*country*".nc");
        pop_grid = float.(pop_file["pop"][:,:,1]); #144x91 Array
        pop_global .+= pop_grid
        pop_dict[country] = pop_grid
        pop_totals[i] = sum(pop_grid)
        pop_idx[country] = findall(!iszero, pop_grid) #save non-zero indices only
    end
    pop_dict["global"] = pop_global

    return pop_dict, pop_idx, pop_totals
end


######### HEALTH CONSTANTS

Base.@kwdef struct HealthConstants
    n_exposure_years::Int64
    mortality_years::Array{Int64,1}
    VSLtype::String
    countries::Array{String,1}
    frac_25_plus::Array{Float64,2}
    IR_LC25::Array{Float64,2}
    IR_stroke25::Array{Float64,2}
    IR_IHD25::Array{Float64,2}
    IR_COPD25::Array{Float64,2}
    IR_LRI25::Array{Float64,2}
    frac_30_plus::Array{Float64,2}
    IR_cardio30::Array{Float64,2}
    IR_resp30::Array{Float64,2}
    IR_AC30::Array{Float64,2}
    X_PM25_base::Array{Float64,3}
    pop_dict::Dict{String,Array{Float64,2}}
    pop_totals::Array{Float64,1}
    pop_idx::Dict{String, Array{CartesianIndex{2},1}}
    pop_ratios_25::Array{Float64,2}
    pop_ratios_30::Array{Float64,2}
    GDP_deflator_ratio::Float64
    US_GDP_pc_ratio::Float64
    GDP_PPP_ratios::Array{Float64,2}
end


function getHealthConstants(simulation::Simulations)
    VSLtype = simulation.VSLtype

    n_exposure_years = length(simulation.emissionyears) + simulation.sensitivityAQyears - 1 #total number of years of exposures from sensitivities for all emission years
    exposure_years = collect(simulation.t0:1:simulation.emissionyears[end]+n_exposure_years-1)
    n_cessation_lag = 20 #number of years that cessation lag adds
    if simulation.cessationlag
        mortality_years = collect(simulation.t0:1:simulation.emissionyears[end]+simulation.sensitivityAQyears-1+n_cessation_lag-1) #years that mortalities occur in, including cessation lag
    else
        mortality_years = collect(simulation.t0:1:simulation.emissionyears[end]+simulation.sensitivityAQyears-1) #years that mortalities occur in, does NOT include cessation lag
    end

    #GBD countries 
    countries = getCountryList();

    # GBD info: fraction 25+ or 30+, incidence rates (IR)
    f = jldopen("data/GBD_data/GBD_mortalities.jld2")
    frac_25_plus = f["frac_25_plus"] #fraction of that country's population over 25 for the year 2019
    frac_30_plus = f["frac_30_plus"] #fraction of that country's population over 30 for the year 2019

    #PM2.5 - GEMM CRF from Burnett et al., applies to population over 25
    #Incidence rates for following diseases over age 25 from GBD2019, for year 2019 (latest year)
    IR_LC25 = f["IR_LC25"] ./ 100_000           #lung cancer
    IR_stroke25 = f["IR_stroke25"] ./ 100_000   #stroke
    IR_IHD25 = f["IR_IHD25"] ./ 100_000         #ischemic heart disease
    IR_COPD25 = f["IR_COPD25"] ./ 100_000       #chronic obstructive pulmonary disease
    IR_LRI25 = f["IR_LRI25"] ./ 100_000         #lower respiratory infections

    #PM2.5 - log-linear CRF using Hoek et al. (2013), applies to population over 30
    IR_cardio30 = f["IR_cardio30"] ./ 100_000         #Cardiovascular diseases

    #surface ozone - log-linear CRF using Turner et al., applies to population over 30: EITHER respiratory diseases or all-cause
    #Incidence rates for following diseases over age 30 from GBD2019, for year 2019 (latest year)
    IR_COPD30 = f["IR_COPD30"] ./ 100_000       #chronic obstructive pulmonary disease
    IR_LRI30 = f["IR_LRI30"] ./ 100_000         #lower respiratory infections
    IR_URI30 = f["IR_URI30"] ./ 100_000         #upper respiratory infections
    IR_asthma30 = f["IR_asthma30"] ./ 100_000   #asthma
    IR_pneum30 = f["IR_pneum30"] ./ 100_000     #pneumocosiosis
    IR_resp30 = IR_LRI30 .+ IR_URI30 .+ IR_COPD30 .+ IR_asthma30  #.+ IR_pneum30    #respiratory diseases
    IR_AC30 = f["IR_AC30"] ./ 100_000           #all-cause
    close(f)

    # Baseline PM2.5 concentration (needed for GEMM CRF)
    ds_X_PM25_base = NCDataset("data/airquality/PM_baseline_conc_2x25_2000_2011.nc")
    X_PM25_base = float.(ds_X_PM25_base["X_PM_base"][:]) #144 lon x 91 lat x 12 years (2000-2011)
    if n_exposure_years == 1
        X_PM25_base = mean(X_PM25_base,dims=3) #144 x 91 x 1 (use average of all 12 years)
    else
        X_PM25_base = X_PM25_base[:,:,1:n_exposure_years]
    end

    # population info
    ###SSP growth
    f = jldopen("data/population/population_SSP.jld2")
    SSP_years = 2010:5:2100
    idx_2020 = 3        #column index of year 2020 = 3
    interp_idx = LinearInterpolation(SSP_years, 1:19, extrapolation_bc=Line());

    population_25_SSP = f["population_25_SSPs"][:,:,simulation.SSP] #204 (ncountries) x 19 (years)
    population_30_SSP = f["population_30_SSPs"][:,:,simulation.SSP] 
    close(f)

    interp_pop_growth_25 = interpolate(population_25_SSP, BSpline(Linear()));
    interp_pop_growth_30 = interpolate(population_30_SSP, BSpline(Linear()));

    pop0_25 = interp_pop_growth_25[:, interp_idx[exposure_years]]
    pop_2020_25 = population_25_SSP[:,idx_2020] 
    pop_ratios_25 = pop0_25 ./ pop_2020_25 #country-specific ratios for age 25+ population compared to 2020 for every year in mortality years, ncountries x mortality years x SSP scenarios

    pop0_30 = interp_pop_growth_30[:, interp_idx[exposure_years]]
    pop_2020_30 = population_30_SSP[:,idx_2020,:] 
    pop_ratios_30 = pop0_30 ./ pop_2020_30 #country-specific ratios for age 30+ population compared to 2020 for every year in mortality years, ncountries x mortality years x SSP scenarios


    ###gridded pop data for year 2020, 2deg lat x 2.5 deg lon
    pop_dict, pop_idx, pop_totals = getPopulationGriddata(countries)

    # GDP info
    ##deflator ratio for tmoney
    f = jldopen("data/GDP_data/GDP_deflator.jld2")
    GDP_def_years = f["years"] #1990 - 2021
    GDP_def_ratios = f["ratios"]
    close(f)
    i = findall(x-> x==simulation.tmoney, GDP_def_years)[1]
    GDP_deflator_ratio = GDP_def_ratios[i]

    ##US GDP per capita for different years
    f = jldopen("data/GDP_data/GDP_US_pc.jld2")
    GDP_pc_years = f["years"] #1990 - 2020
    GDP_pc_ratios = f["ratios"]
    close(f)
    # get indices for ratio (GDP pc year 2020 / GDP pc year 1990) for 2020 - future projections will just use SSP GDP growth only
        # # if any year in mortality_years < 1990 or > 2050 (bounds), use closest value (1990 or 2050)
        # i = [findall(x-> y < years[1] ? x==years[1] : x==min(y,years[end]), years)[1] for y in mortality_years]
    US_GDP_pc_ratio = GDP_pc_ratios[end]

    ##GDP PPP for every country for different years
    f = jldopen("data/GDP_data/GDP_country_PPP.jld2")
    idx_US = findall(x->x=="USA",countries)[1] #row 193
        
    #GDP PPP pc projections, row = country, column = SSP_years (2010:5:2100)
    GDP_PPP = f["GDP_PPP_SSPs"][:,:,simulation.SSP] #204 (ncountries) x 19 (years)
    close(f)

    #PPP ratios = # GDP PPP_country_every year / GDP PPP_US_2020
    GDP_PPP_US_2020 = GDP_PPP[idx_US, idx_2020]
    GDP_PPP_ratios = GDP_PPP ./ GDP_PPP_US_2020
    if VSLtype == "US"
        GDP_PPP_ratios = repeat(GDP_PPP_ratios[idx_US:idx_US,:], 204,1)
    end

    #data only provides for every 5 years, need to interpolate for years in between
    interp_GDP_PPP_ratios = interpolate(GDP_PPP_ratios, BSpline(Linear()))

    #get indices for ratios (GDP PPP country i in year i/ GDP PPP USA in year i) - same method as above
    m_years = [y .< SSP_years[1] ? SSP_years[1] : min(y, SSP_years[end]) for y in mortality_years] #if mortality years are < 2010 or > 2100, use closest edge value (2010 or 2100)
    GDP_PPP_ratios = interp_GDP_PPP_ratios[:,interp_idx[m_years]] #204 x length(mortality years) array


    return HealthConstants(
                    n_exposure_years=n_exposure_years,
                    mortality_years=mortality_years,
                    VSLtype=VSLtype,
                    countries=countries,
                    frac_25_plus=frac_25_plus,
                    IR_LC25=IR_LC25,
                    IR_stroke25=IR_stroke25,
                    IR_IHD25=IR_IHD25,
                    IR_COPD25=IR_COPD25,
                    IR_LRI25=IR_LRI25,
                    frac_30_plus=frac_30_plus,
                    IR_cardio30=IR_cardio30,
                    IR_resp30=IR_resp30,
                    IR_AC30=IR_AC30,
                    X_PM25_base=X_PM25_base,
                    pop_dict=pop_dict,
                    pop_totals=pop_totals,
                    pop_idx=pop_idx,
                    pop_ratios_25=pop_ratios_25,
                    pop_ratios_30=pop_ratios_30,
                    GDP_deflator_ratio=GDP_deflator_ratio,
                    US_GDP_pc_ratio=US_GDP_pc_ratio,
                    GDP_PPP_ratios=GDP_PPP_ratios)  
end

