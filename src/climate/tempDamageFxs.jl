"""
# loadBackgroundGDPs - 
### Inputs:
- 

"""
function loadBackgroundGDPs(simulation::Simulations, climatedistribution::MCDistributions, climateconsts::ClimateConstants)
	#climateconsts already selected relevant SSP scenario in getClimateConstants function and saved GDP_matrix
	GDP_matrix = climateconsts.GDP_matrix
	GDP_matrix_US = climateconsts.GDP_matrix_US
	return (GDP_matrix, GDP_matrix_US)
end

"""
# getDamageFunction - 
### Inputs:

"""
function getDamageFunction()
	damage_func = (a1,a2,a3,dT) -> 1 .- (1 ./ (1 .+ a1.*dT .+ a2.*dT.^a3));
	return damage_func
end

function getDamageAvFunction(simulation::Simulations, climatedistribution::MCDistributions, dT_background::Matrix{Float64})
	nrun = simulation.nrun
	totalimpact = simulation.totalimpact

	a19 = repeat([sample(climatedistribution.a19,i) for i in 1:nrun]', totalimpact,1);
	a29 = repeat([sample(climatedistribution.a29,i) for i in 1:nrun]', totalimpact,1);
	a39 = repeat([sample(climatedistribution.a39,i) for i in 1:nrun]', totalimpact,1);
	
	damage_func = getDamageFunction()
	damage_background = damage_func(a19, a29, a39, dT_background)
	damage_av_func = (dT_av) -> damage_background .- damage_func(a19, a29, a39, (dT_background .- dT_av))
	
	return damage_av_func
end

"""
# getDamageAvFunction_US - US Only Damages based on Hsiang et al. (2017) (Published in Science), directly from APMT-ICv24c
### Inputs

"""
function getDamageAvFunction_US(simulation::Simulations, climatedistribution::MCDistributions, dT_background::Matrix{Float64})
	# a19_US: Mid Value
	# a29_US: Mid Value
	# aUS_unc_fact: Uncertainty to be multiplied by Hsiang mid values
	# aUS_unc_add: Uncertainty to be added to by Hsiang mid values
	# aUS_distrib_rep:  Normal distribution with 0 mean and 1 std_dev
	# aUS_tempAdjust: Hsiang et al. (2017) computes damages relative to 1980-2010, while DICE computes relative to preindustrial.
	
	nrun = simulation.nrun
	totalimpact = simulation.totalimpact

	a19_US = repeat([sample(climatedistribution.a19_US,i) for i in 1:nrun]', totalimpact,1);
	a29_US = repeat([sample(climatedistribution.a29_US,i) for i in 1:nrun]', totalimpact,1);
	aUS_distrib = repeat([sample(climatedistribution.aUS_distrib,i) for i in 1:nrun]', totalimpact,1);
	aUS_tempAdjust = repeat([sample(climatedistribution.aUS_tempAdjust,i) for i in 1:nrun]', totalimpact,1)
	aUS_unc_add = repeat([sample(climatedistribution.aUS_unc_add,i) for i in 1:nrun]', totalimpact,1)
	aUS_unc_fact = repeat([sample(climatedistribution.aUS_unc_fact,i) for i in 1:nrun]', totalimpact,1)
	
	damage_func_US = (aUS_dist,dT) -> (1 .+ aUS_unc_fact .* aUS_dist) .* (a19_US .* (dT) .+ a29_US .* (dT).^2) .+ aUS_unc_add .* aUS_dist
	damage_background_US = damage_func_US(aUS_distrib, dT_background .- aUS_tempAdjust)
	damage_av_func_US = (dT_av) -> damage_background_US .- damage_func_US(aUS_distrib, dT_background .- dT_av .- aUS_tempAdjust)
	
	return damage_av_func_US
end

function discountedCostClimate(cost::Array{Float64,2}, simulation::Simulations)
	nrun = simulation.nrun
	DR = simulation.DR
	totalimpact = simulation.totalimpact

	discount_repeat = DR .* ones(totalimpact-1);
	discount_vector = vcat([1], exp.(-1 .* cumsum(discount_repeat)))
	discount_matrix = repeat(discount_vector,1,nrun)
	
	return cost .* discount_vector
end

function NPV(discounted_cost)
	NPV_value = sum(discounted_cost,dims=1)
	return NPV_value
end

