"""
# CRF_LL - Log-Linear CRF

### Inputs
- β: parameter estimate for pollutant (surface ozone from Turner et al. 2015)
- ΔX: gridded change in pollutant (surface ozone) exposure. size = 144 (lon) x 91 (lat) x nyears; ΔX[:,:,1] = change in exposure in the first year

### Output
- RR: relative risk at every gridcell, same array size as input ΔX
"""
function CRF_LL(β::Float64, ΔX::Union{Array{Float64},Float64})
    RR = exp.(β .* ΔX)
    return RR
end

"""
# CRF_GEMM - GEMM (Global Exposure Mortality Model) CRF from Burnett et al. 2018

### Inputs
- θ, α, μ, ν: parameter estimates for specific disease (IHD, stroke, COPD, LC, or LRI)
- X: PM₂.₅ exposure in μg/m³. size = 144 (lon) x 91 (lat) x nyears; X[:,:,1] = exposure in the first year

### Output
- RR: relative risk at every gridcell, same array size as input ΔX
"""
function CRF_GEMM(θ::Float64, α::Float64, μ::Float64, ν::Float64, X_PM::Union{Array{Float64},Float64})
    z = max.(X_PM.-2.4, 0.)

    f = log.(1 .+ (z ./ α))
    w = 1 ./ (1 .+ exp.( -(z .- μ) ./ ν) )
    T = f .* w 

    RR = exp.(θ .* T)
    return RR
end

"""
# cessationlag - add cessation lag to mortalities in each year 
### Input
- ΔM: total mortalities due to exposure in each year, size = ncountries x nyears, where nyears = years of different exposures 
- simulation: Simulation type variable containing fixed data from user input

### Output
- ΔM_yearly: ΔM with (or without) cessation lag included based on user input, size = ncountries x nyears or (nyears + 19)
"""
function cessationlag(ΔM::Array{Float64}, simulation::Simulations)
    if !simulation.cessationlag
        return ΔM
    else
        yearlyfraction = [0.3, 0.125, 0.125, 0.125, 0.125, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 
                        0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15, 0.2/15];
                        
        yearlyfraction = reshape(yearlyfraction, (1,20))
        
        ncountries, nyears = size(ΔM) #nyears = years of emissions
        ΔM_yearly = zeros(Float64,(ncountries,nyears+19))

        for i in 1:nyears
            ΔM_yeari = (ΔM[:,i]) * yearlyfraction
            ΔM_yearly[:,i:i+19] .+= ΔM_yeari
        end
        
        return ΔM_yearly
    end
end


"""
# healthimpact_O3 - calculate health impact (premature mortalities) from change in surface ozone exposure using log-linear CRF
### Inputs (MC variable):
- i: iteration of main Monte Carlo for loop
- ΔX: change in surface ozone exposure from aviation emissions, 144 lon x 91 lat x nyears. [:,:,1] = exposure in the first year
### Inputs (fixed):
- distribution: 
- healthconsts: contains data such as country list and GBD incidence rate data, etc.

### Output
- ΔM_yearly: array of size  ncountries x nyears or (nyears + 19) if cessation lag is included. ΔM[:,i] = mortalities from emissions in year i
"""
@views function healthimpact_O3(i::Int64, ΔX_O3::Array{Float64, 3}, distribution::MCDistributions, simulation::Simulations, healthconsts::HealthConstants)
    SSP = simulation.SSP #1-5
    if simulation.O3_risk == "respiratory"
        β = sample(distribution.β_resp, i) #respiratory diseases
        Iobs_O3 = healthconsts.IR_resp30 
    else
        β = sample(distribution.β_AC, i) #all-cause
        Iobs_O3 = healthconsts.IR_AC30
    end
    #Function will need to be updated if you want to use different beta (e.g. all-cause)
    #example: β = sample(distribution.β_AC, i) and Iobs_O3 below, etc.

    countries = healthconsts.countries
    frac_30_plus = healthconsts.frac_30_plus
    pop_dict = healthconsts.pop_dict
    pop_idx = healthconsts.pop_idx
    pop_ratios_30 = healthconsts.pop_ratios_30

    nyears = size(ΔX_O3)[end] # should match 2nd dimension of pop_ratios_30
    ncountries = length(countries)

    RR = CRF_LL(β, ΔX_O3)   #log linear CRF for ozone
    ΔIoverIobs = (RR .- 1) ./ RR

    ΔM_O3 = zeros(ncountries,nyears)
    for i in 1:ncountries 
        country = countries[i]
        idx_country = pop_idx[country]
        paff_grid = pop_dict[country] 
        
        ΔM_O3[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs[idx_country,:]), dims = 1);
    end

    ΔM_O3 = ΔM_O3 .* Iobs_O3 .* frac_30_plus .* pop_ratios_30
    ΔM_yearly_O3 = cessationlag(ΔM_O3, simulation)
    return ΔM_yearly_O3
end


"""
# healthimpact_PM - calculate health impact (premature mortalities) from change in PM2.5 exposure using GEMM CRF

### Inputs
- i: iteration of main Monte Carlo for loop
- ΔX_PM: change in PM2.5 exposure from aviation emissions, 144 lon x 91 lat x nyears. [:,:,1] = exposure in the first year

### Output
- ΔM_yearly: array of size  ncountries x nyears or (nyears + 19) if cessation lag is included. ΔM[:,i] = mortalities from emissions in year i
"""
function healthimpact_PM25_GEMM(i::Int64, ΔX_PM::Array{Float64,3}, distribution::MCDistributions, simulation::Simulations, healthconsts::HealthConstants)
    n_exposure_years = healthconsts.n_exposure_years
    X_PM25_base = healthconsts.X_PM25_base #baseline PM2.5 concentration for corresponding years
    X_PM = X_PM25_base .+ ΔX_PM

    #IHD
    Iobs_IHD = healthconsts.IR_IHD25
    θ_IHD = sample(distribution.θ_IHD,i)
    α_IHD = 1.9
    μ_IHD = 12.
    ν_IHD = 40.2
    RR_pol_IHD = CRF_GEMM(θ_IHD, α_IHD, μ_IHD, ν_IHD, X_PM)
    RR_base_IHD = CRF_GEMM(θ_IHD, α_IHD, μ_IHD, ν_IHD, X_PM25_base)
    ΔIoverIobs_IHD = RR_pol_IHD./RR_base_IHD .- 1.

    #stroke
    Iobs_stroke = healthconsts.IR_stroke25
    θ_stroke = sample(distribution.θ_stroke,i)
    α_stroke = 6.2
    μ_stroke = 16.7
    ν_stroke = 23.7
    RR_pol_stroke = CRF_GEMM(θ_stroke, α_stroke, μ_stroke, ν_stroke, X_PM)
    RR_base_stroke = CRF_GEMM(θ_stroke, α_stroke, μ_stroke, ν_stroke, X_PM25_base)
    ΔIoverIobs_stroke = RR_pol_stroke./RR_base_stroke .- 1.

    #COPD
    Iobs_COPD = healthconsts.IR_COPD25
    θ_COPD = sample(distribution.θ_COPD,i)
    α_COPD = 6.5
    μ_COPD = 2.5
    ν_COPD = 32.
    RR_pol_COPD = CRF_GEMM(θ_COPD, α_COPD, μ_COPD, ν_COPD, X_PM)
    RR_base_COPD = CRF_GEMM(θ_COPD, α_COPD, μ_COPD, ν_COPD, X_PM25_base)
    ΔIoverIobs_COPD = RR_pol_COPD./RR_base_COPD .- 1.

    #Lung cancer
    Iobs_LC = healthconsts.IR_LC25
    θ_LC = sample(distribution.θ_LC,i)
    α_LC = 6.2
    μ_LC = 9.3
    ν_LC = 29.8
    RR_pol_LC = CRF_GEMM(θ_LC, α_LC, μ_LC, ν_LC, X_PM)
    RR_base_LC = CRF_GEMM(θ_LC, α_LC, μ_LC, ν_LC, X_PM25_base)
    ΔIoverIobs_LC = RR_pol_LC./RR_base_LC .- 1.

    #Lower respiratory infections
    Iobs_LRI = healthconsts.IR_LRI25
    θ_LRI = sample(distribution.θ_LRI,i)
    α_LRI = 6.4
    μ_LRI = 5.7
    ν_LRI = 8.4
    RR_pol_LRI = CRF_GEMM(θ_LRI, α_LRI, μ_LRI, ν_LRI, X_PM)
    RR_base_LRI = CRF_GEMM(θ_LRI, α_LRI, μ_LRI, ν_LRI, X_PM25_base)
    ΔIoverIobs_LRI = RR_pol_LRI./RR_base_LRI .- 1.

    countries = healthconsts.countries
    ncountries = length(countries)
    frac_25_plus = healthconsts.frac_25_plus
    pop_dict = healthconsts.pop_dict
    pop_idx = healthconsts.pop_idx
    pop_ratios_25 = healthconsts.pop_ratios_25

    ΔM_IHD = zeros(ncountries,n_exposure_years)
    ΔM_stroke = zeros(ncountries,n_exposure_years)
    ΔM_COPD = zeros(ncountries,n_exposure_years)
    ΔM_LC = zeros(ncountries,n_exposure_years)
    ΔM_LRI = zeros(ncountries,n_exposure_years)
    
    for i in 1:ncountries 
        country = countries[i]
        idx_country = pop_idx[country]
        paff_grid = pop_dict[country] 
        
        ΔM_IHD[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs_IHD[idx_country,:]), dims = 1)
        ΔM_stroke[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs_stroke[idx_country,:]), dims = 1)
        ΔM_COPD[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs_COPD[idx_country,:]), dims = 1)
        ΔM_LC[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs_LC[idx_country,:]), dims = 1)
        ΔM_LRI[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs_LRI[idx_country,:]), dims = 1)
    end

    ΔM_IHD = ΔM_IHD .* Iobs_IHD
    ΔM_stroke = ΔM_stroke .* Iobs_stroke
    ΔM_COPD = ΔM_COPD .* Iobs_COPD
    ΔM_LC = ΔM_LC .* Iobs_LC
    ΔM_LRI = ΔM_LRI .* Iobs_LRI

    ΔM_PM = (ΔM_IHD .+ ΔM_stroke .+ ΔM_COPD .+ ΔM_LC .+ ΔM_LRI)   .* frac_25_plus .* pop_ratios_25
    ΔM_yearly_PM = cessationlag(ΔM_PM, simulation)
    return ΔM_yearly_PM
end

"""
Log-linear PM CRF using Hoek et al data
"""
@views function healthimpact_PM25_LL(i::Int64, ΔX_PM::Array{Float64,3}, distribution::MCDistributions, simulation::Simulations, healthconsts::HealthConstants)
    SSP = simulation.SSP #1-5
    β = sample(distribution.β_cardio, i) #cardiovascular 

    countries = healthconsts.countries
    frac_30_plus = healthconsts.frac_30_plus
    pop_dict = healthconsts.pop_dict
    pop_idx = healthconsts.pop_idx
    pop_ratios_30 = healthconsts.pop_ratios_30

    Iobs_PM = healthconsts.IR_cardio30 #ASSUMES RESPIRATORY INFECTIONS AS CAUSE

    nyears = size(ΔX_PM,3) # should match 2nd dimension of pop_ratios_30
    ncountries = length(countries)

    RR = CRF_LL(β, ΔX_PM)   #log linear CRF for ozone
    ΔIoverIobs = (RR .- 1) ./ RR

    ΔM_PM = zeros(ncountries,nyears)
    for i in 1:ncountries 
        country = countries[i]
        idx_country = pop_idx[country]
        paff_grid = pop_dict[country] 
        
        ΔM_PM[i,:] = sum((paff_grid[idx_country] .* ΔIoverIobs[idx_country,:]), dims = 1);
    end

    ΔM_PM = ΔM_PM .* Iobs_PM .* frac_30_plus .* pop_ratios_30
    ΔM_yearly_PM = cessationlag(ΔM_PM, simulation)
    return ΔM_yearly_PM
end
