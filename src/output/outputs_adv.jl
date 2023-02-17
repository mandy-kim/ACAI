Base.@kwdef mutable struct ClimateOutputs <: Outputs
    #Climate
    RF_CO2av::Array{Float64,2}      #size = (totalimpact, nrun)
    RF_Contrail::Array{Float64,2}   
    RF_NOx::Array{Float64,2}        #size = (totalimpact, 1) - doesn't change per MC run
    # RF_H2O::Array{Float64,2}        
    # RF_BC::Array{Float64,2}         
    # RF_S::Array{Float64,2}
    #
    ΔT_CO2av::Array{Float64,2}      #size = (totalimpact, nrun) 
    ΔT_Contrail::Array{Float64,2} 
    ΔT_NOx::Array{Float64,2} 
    # ΔT_H2O::Array{Float64,2} 
    # ΔT_BC::Array{Float64,2} 
    # ΔT_S::Array{Float64,2}
    ΔT_background1::Array{Float64,2} 
    ΔT_background2::Array{Float64,2}
    #
    cost_CO2av::Array{Float64,2} #size = (totalimpact, nrun) 
    cost_Contrail::Array{Float64,2} 
    cost_NOx::Array{Float64,2} 
    #
    discountedcost_CO2av::Array{Float64,2} #size = (totalimpact, nrun) 
    discountedcost_Contrail::Array{Float64,2} 
    discountedcost_NOx::Array{Float64,2} 
    #
    NPV_CO2av::Array{Float64,2}  #size = (1, nrun) 
    NPV_Contrail::Array{Float64,2} 
    NPV_NOx::Array{Float64,2} 
end

Base.@kwdef mutable struct AirQualityOutputs <: Outputs
    # Air quality
    #
    ΔX_PM25_NOx::Array{Float64,3}       #size = (lon, lat, exposure years) - doesn't change per MC run
    # ΔX_PM25_BC::Array{Float64}
    # ΔX_PM25_S::Array{Float64}
    # ΔX_PM25_CO::Array{Float64}
    # ΔX_PM25_HC::Array{Float64}
    # ΔX_PM25_OC::Array{Float64}
    #
    ΔX_O3_NOx::Array{Float64,3}       #size = (lon, lat, exposure years) - doesn't change per MC run
    #
    ΔX_PM25_avg_NOx::Array{Float64,2}  #size = (ncountries=ncountries, exposure years) - doesn't change per MC run
    ΔX_O3_avg_NOx::Array{Float64,2}
    #
    ΔM_PM25_NOx::Array{Float64,3}        #size = (ncountries=ncountries, mortality years, nrun)
    ΔM_O3_NOx::Array{Float64,3}
    ΔM_total_global_NOx::Array{Float64,3} #size = (1, mortality years, nrun) 
    #
    cost_ΔM_PM25_NOx::Array{Float64,3}   #size = (ncountries=ncountries, mortality years, nrun)
    cost_ΔM_O3_NOx::Array{Float64,3}
    cost_ΔM_total_global_NOx::Array{Float64,3} #size = (1, mortality years, nrun)
    #
    discountedcost_ΔM_PM25_NOx::Array{Float64,3}   #size = (ncountries=ncountries, mortality years, nrun)
    discountedcost_ΔM_O3_NOx::Array{Float64,3}
    discountedcost_ΔM_total_global_NOx::Array{Float64,3} #size = (1, mortality years, nrun)
    # 
    NPV_ΔM_PM25_NOx::Array{Float64,2}    #size = (ncountries=ncountries, nrun)
    NPV_ΔM_O3_NOx::Array{Float64,2}
    NPV_ΔM_total_global_NOx::Array{Float64,2} #size = (1, nrun)

    #O3 column
    Col_O3_NOx::Array{Float64,3}         #size = (lon, lat, O3 years) - doesn't change per MC run
    #Zonal O3
    Zonal_O3_NOx::Array{Float64,3}         #size = (alt, lat, O3 years) - doesn't change per MC run
end

function ClimateOutputs(simulation::Simulations_adv)
    nrun = simulation.nrun;
    totalimpact = simulation.totalimpact;

    climateoutput =  ClimateOutputs(
                #
                RF_CO2av=zeros(totalimpact,nrun),
                RF_Contrail=zeros(totalimpact,nrun),
                RF_NOx=zeros(totalimpact,1),
                #
                ΔT_CO2av=zeros(totalimpact,nrun),
                ΔT_NOx=zeros(totalimpact,nrun),
                ΔT_Contrail=zeros(totalimpact,nrun),
                ΔT_background1=zeros(totalimpact,nrun),
                ΔT_background2=zeros(totalimpact,nrun),
                #
                cost_CO2av=zeros(totalimpact,nrun),
                cost_NOx=zeros(totalimpact,nrun),
                cost_Contrail=zeros(totalimpact,nrun),
                #
                discountedcost_CO2av=zeros(totalimpact,nrun),
                discountedcost_NOx=zeros(totalimpact,nrun),
                discountedcost_Contrail=zeros(totalimpact,nrun),
                # 
                NPV_CO2av=zeros(1,nrun),
                NPV_NOx=zeros(1,nrun),
                NPV_Contrail=zeros(1,nrun),  
                )
    return climateoutput
end

function AirQualityOutputs(simulation::Simulations_adv, healthconsts::HealthConstants)
    nrun = simulation.nrun;
    totalimpact = simulation.totalimpact;
    ncountries = length(healthconsts.countries)

    aqoutput = AirQualityOutputs(
            #
            ΔX_PM25_NOx=zeros(144, 91, healthconsts.n_exposure_years),
            ΔX_O3_NOx=zeros(144, 91, healthconsts.n_exposure_years),
            ΔX_PM25_avg_NOx=zeros(ncountries,healthconsts.n_exposure_years),
            ΔX_O3_avg_NOx=zeros(ncountries,healthconsts.n_exposure_years),
            #
            ΔM_PM25_NOx=zeros(ncountries,length(healthconsts.mortality_years),nrun),
            ΔM_O3_NOx=zeros(ncountries,length(healthconsts.mortality_years),nrun),
            ΔM_total_global_NOx=zeros(1,length(healthconsts.mortality_years),nrun),
            #
            cost_ΔM_PM25_NOx=zeros(ncountries,length(healthconsts.mortality_years),nrun),
            cost_ΔM_O3_NOx=zeros(ncountries,length(healthconsts.mortality_years),nrun),
            cost_ΔM_total_global_NOx=zeros(1,length(healthconsts.mortality_years),nrun),
            #
            discountedcost_ΔM_PM25_NOx=zeros(ncountries,length(healthconsts.mortality_years),nrun),
            discountedcost_ΔM_O3_NOx=zeros(ncountries,length(healthconsts.mortality_years),nrun),
            discountedcost_ΔM_total_global_NOx=zeros(1,length(healthconsts.mortality_years),nrun),
            # 
            NPV_ΔM_PM25_NOx=zeros(ncountries,nrun),
            NPV_ΔM_O3_NOx=zeros(ncountries,nrun),
            NPV_ΔM_total_global_NOx=zeros(1,nrun),
            #
            Col_O3_NOx=zeros(144, 91, simulation.sensitivityO3years),
            Zonal_O3_NOx=zeros(72, 91, simulation.sensitivityO3years),
            );
    
    return aqoutput
end