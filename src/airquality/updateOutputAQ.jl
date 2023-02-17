function update_output_ΔMs(i::Int64, NOx_impacts::TotalImpacts, aqoutput::AirQualityOutputs, distribution::MCDistributions, simulation::Simulations_adv, healthconsts::HealthConstants)
    aqoutput.ΔM_O3_NOx[:,:,i] = healthimpact_O3(i, NOx_impacts.Sfc_O3, distribution, simulation, healthconsts)
    if simulation.PM25_CRF == "GEMM"
        aqoutput.ΔM_PM25_NOx[:,:,i] = healthimpact_PM25_GEMM(i, NOx_impacts.PM25, distribution, simulation, healthconsts)
    elseif simulation.PM25_CRF == "LL"
        aqoutput.ΔM_PM25_NOx[:,:,i] = healthimpact_PM25_LL(i, NOx_impacts.PM25, distribution, simulation, healthconsts)
    end

    aqoutput.ΔM_total_global_NOx[:,:,i] = sum(aqoutput.ΔM_O3_NOx[:,:,i] .+  aqoutput.ΔM_PM25_NOx[:,:,i],dims=1)

    return aqoutput
end

function update_output_ΔM_costs(i::Int64, aqoutput::AirQualityOutputs, distribution::MCDistributions, simulation::Simulations_adv, healthconsts::HealthConstants)
    cost_ΔM_PM25_NOx = healthdamages(i, aqoutput.ΔM_PM25_NOx[:,:,i], distribution, simulation, healthconsts)
    cost_ΔM_O3_NOx = healthdamages(i, aqoutput.ΔM_O3_NOx[:,:,i], distribution, simulation, healthconsts)
    aqoutput.cost_ΔM_PM25_NOx[:,:,i] = cost_ΔM_PM25_NOx
    aqoutput.cost_ΔM_O3_NOx[:,:,i] = cost_ΔM_O3_NOx
    aqoutput.cost_ΔM_total_global_NOx[:,:,i] = sum(cost_ΔM_PM25_NOx .+ cost_ΔM_O3_NOx, dims=1)

    discountedcost_ΔM_PM25_NOx = discountedCost_ΔM(cost_ΔM_PM25_NOx, simulation)
    discountedcost_ΔM_O3_NOx = discountedCost_ΔM(cost_ΔM_O3_NOx, simulation)
    aqoutput.discountedcost_ΔM_PM25_NOx[:,:,i] = discountedcost_ΔM_PM25_NOx
    aqoutput.discountedcost_ΔM_O3_NOx[:,:,i] = discountedcost_ΔM_O3_NOx
    aqoutput.discountedcost_ΔM_total_global_NOx[:,:,i] = sum(discountedcost_ΔM_PM25_NOx .+ discountedcost_ΔM_O3_NOx, dims=1)

    NPV_ΔM_PM25_NOx = sum(discountedcost_ΔM_PM25_NOx, dims=2)
    NPV_ΔM_O3_NOx = sum(discountedcost_ΔM_O3_NOx, dims=2)
    aqoutput.NPV_ΔM_PM25_NOx[:,i] = NPV_ΔM_PM25_NOx
    aqoutput.NPV_ΔM_O3_NOx[:,i] = NPV_ΔM_O3_NOx
    aqoutput.NPV_ΔM_total_global_NOx[:,i] = sum(NPV_ΔM_PM25_NOx .+ NPV_ΔM_O3_NOx, dims=1)

    return aqoutput
end