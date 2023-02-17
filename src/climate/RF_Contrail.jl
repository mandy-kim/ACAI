Base.@kwdef struct ContrailMaps
    region_mask::Array{Float64,2}
    global_idx::Array{CartesianIndex{2},1}
    nam_idx::Array{CartesianIndex{2},1}
    conus_idx::Array{CartesianIndex{2},1}
    weu_idx::Array{CartesianIndex{2},1}
    rus_idx::Array{CartesianIndex{2},1}
    meast_idx::Array{CartesianIndex{2},1}
    seas_idx::Array{CartesianIndex{2},1}
    cpac_idx::Array{CartesianIndex{2},1}
    npac_idx::Array{CartesianIndex{2},1}
    natl_idx::Array{CartesianIndex{2},1}
    atl_idx::Array{CartesianIndex{2},1}
end

function load_region_map()
    f = jldopen("data/sensitivities/Contrail/contrail_mask.jld2")
    region_mask = f["region_map"]
    close(f)

    global_idx = findall(x->x == 0, region_mask)
    nam_idx = findall(x->x == 1, region_mask)
    conus_idx = findall(x->x == 2, region_mask)
    weu_idx = findall(x->x == 3, region_mask)
    rus_idx = findall(x->x == 4, region_mask)
    meast_idx = findall(x->x == 5, region_mask)
    seas_idx = findall(x->x == 6, region_mask)
    cpac_idx = findall(x->x == 7, region_mask)
    npac_idx = findall(x->x == 8, region_mask)
    natl_idx = findall(x->x == 9, region_mask)
    atl_idx = findall(x->x == 10, region_mask)

    contrail_map = ContrailMaps(
        region_mask=region_mask,
        global_idx=global_idx,
        nam_idx=nam_idx,
        conus_idx=conus_idx,
        weu_idx=weu_idx,
        rus_idx=rus_idx,
        meast_idx=meast_idx,
        seas_idx=seas_idx,
        cpac_idx=cpac_idx,
        npac_idx=npac_idx,
        natl_idx=natl_idx,
        atl_idx=atl_idx,
    )
    return contrail_map
end


function update_contrail_map(i::Int, contrail_map::ContrailMaps, distribution::MCDistributions)
    global_rf = sample(distribution.contrail_global_RF,i)
    nam_rf = sample(distribution.contrail_NAM_RF,i)
    conus_rf = sample(distribution.contrail_CONUS_RF,i)
    weu_rf = sample(distribution.contrail_WEU_RF,i)
    rus_rf = sample(distribution.contrail_RUS_RF,i)
    meast_rf = sample(distribution.contrail_MEAST_RF,i)
    seas_rf = sample(distribution.contrail_SEAS_RF,i)
    cpac_rf = sample(distribution.contrail_CPAC_RF,i)
    npac_rf = sample(distribution.contrail_NPAC_RF,i)
    natl_rf = sample(distribution.contrail_NATL_RF,i)
    atl_rf = sample(distribution.contrail_SATL_RF,i)
    contrail_rf_unc = sample(distribution.contrail_RF_unc,i)
    
    new_mask = zeros(size(contrail_map.region_mask))
    new_mask[contrail_map.global_idx] .= global_rf * (contrail_rf_unc)
    new_mask[contrail_map.nam_idx] .= nam_rf * (contrail_rf_unc)
    new_mask[contrail_map.conus_idx] .= conus_rf * (contrail_rf_unc)
    new_mask[contrail_map.weu_idx] .= weu_rf * (contrail_rf_unc) 
    new_mask[contrail_map.rus_idx] .= rus_rf * (contrail_rf_unc)
    new_mask[contrail_map.meast_idx] .= meast_rf * (contrail_rf_unc)
    new_mask[contrail_map.seas_idx] .= seas_rf * (contrail_rf_unc)
    new_mask[contrail_map.cpac_idx] .= cpac_rf * (contrail_rf_unc)
    new_mask[contrail_map.npac_idx] .= npac_rf * (contrail_rf_unc)
    new_mask[contrail_map.natl_idx] .= natl_rf * (contrail_rf_unc)
    new_mask[contrail_map.atl_idx] .= atl_rf * (contrail_rf_unc)
    
    return new_mask
end


function calculateRF_Contrail(lambda_Contrail::Float64, contrail_sens_map::Array{Float64,2}, simulation::Simulations, emission::Emissions)
    totalimpact = simulation.totalimpact
    emissionyears = simulation.emissionyears
    RF_CC_i = zeros(totalimpact)

    nyear = 0
    for yeari in emissionyears
        nyear+=1
        dist_i = emission.yearly_data[yeari].DIST; #lon x lat x lev, units = km
        if sum(dist_i) == 0
            println("Warning: no distances available, RF from contrails not calculated")
            RF_CC_i[nyear] = 0
            return RF_CC_i
        end

        #Display warning:contrail impacts not calculated for flights above 12km 
        #above 12km corresponding to lev 31 in GC vertical grids
        if size(dist_i)[3] > 31
            dist_i = @view dist_i[:,:,1:31]
            if sum(@view dist_i[:,:,32:end]) > 0
                println("Warning: RF from contrails not calculated for flights above 12km")
            end
        end
        RF_CC_i[nyear] = sum(contrail_sens_map .* dist_i) * lambda_Contrail
    end
    return RF_CC_i
end