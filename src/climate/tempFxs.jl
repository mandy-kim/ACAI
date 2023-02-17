"""
# deltaT - calculate temperature change over years of time horizon
### Inputs
- i: iteration of for loop
- RF: vector of length simulation.totalimpact
- simulation: Simulation type variable containing fixed data from user input
- distribution: MCDistributions type variable to sample variables from relevant to ΔT calculations from RF

### Output
- dT1: temperature change from RF, same size as RF input 
"""
 
function deltaT(i::Int, RF::Vector{Float64}, simulation::Simulations, distribution::MCDistributions)
	totalimpact = simulation.totalimpact
	clim_sen = sample(distribution.clim_sen,i); #climate sensitivity (ECS)
	RF2xCO2 = sample(distribution.RF2xCO2,i);

	spec_heat = sample(distribution.spec_heat,i); #J*K^-1*m^-2
	C1 = spec_heat * 0.71 # 71% percent ocean coverage
	C2 = sample(distribution.C2_in,i); #J*K^-1*kg^-1, heat capacity of deep ocean (3000m) 
	F = sample(distribution.advect,i) #kg*m^-2*s^-1,  advective mass flux of water from boundary layer
	Kz = sample(distribution.diffuse,i) #m^2*s^-1, diffusion coefficient for turbulent mixing
	z = sample(distribution.z_in,i) #m, mixing depth

	lambda = clim_sen / RF2xCO2 
	cw = 4.2*10^3; #J*K^-1*kg^-1, specific heat of liquid water
	rho = 1000; #kg/m^3, density of water
	
	tau = C1 * lambda
	alpha1 = cw/C1 * (F + (Kz * rho)/z)
	alpha2 = cw/C2 * (F + (Kz * rho)/z)
    sec_per_yr = 60*60*24*365.25 #seconds per year

	dT1 = zeros(totalimpact) # deltaT_species
    dT2 = zeros(totalimpact) # deltaT_species_deep

	for j = 1:totalimpact-1
		dT1[j+1] = sec_per_yr* (RF[j]/C1 - dT1[j]/tau - alpha1 * (dT1[j] - dT2[j])) + dT1[j]
		dT2[j+1] = sec_per_yr * alpha2 * (dT1[j] - dT2[j]) + dT2[j]
	end
	
	return dT1
end


"""
# backgroundTemp - linearly interpolate or extrapolate background temperature change from MAGICC6 based on RCP scenario and aviation temperature change
### Inputs
- i: iteration of for loop
- deltaTsum_av: deltaT from aviation CO2 (output of deltaT function)
- simulation: Simulation type variable containing fixed data from user input
- distribution: MCDistributions type variable to sample variables from relevant to ΔT calculations from RF
- climateconsts: ClimateConstants type variable for fixed data (MAGICC6 background temp lookup table)

### Output
- deltaT_background: background temperature change, used to calculate damages (see getDamageAvFunction function)
"""
function loadBackgroundTemp(i::Int,deltaTsum_av::Array{Float64},simulation::Simulations,distribution::MCDistributions,climateconsts::ClimateConstants)
	interp_years = simulation.interp_years
	clim_sen = sample(distribution.clim_sen,i)
	iter_deltaT_SSP = climateconsts.iter_deltaT_SSP

	#MAGICC6 background deltaT
	deltaTsum_b = iter_deltaT_SSP(interp_years, clim_sen)
	deltaT_background = hcat(deltaTsum_b, deltaTsum_b-deltaTsum_av)
	
	return deltaT_background
end
