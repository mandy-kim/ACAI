Base.@kwdef mutable struct MainOutputs <: Outputs
    RF_CO2av::Array{Float64,2}      #size = (totalimpact, nrun)
    RF_CO2av_lcycle::Array{Float64,2}
    RF_Contrail::Array{Float64,2}
    RF_O3short::Array{Float64,2}
    RF_O3long::Array{Float64,2}
    RF_Nitrate::Array{Float64,2}
    RF_CH4::Array{Float64,2}
    RF_H2O::Array{Float64,2}
    RF_BC::Array{Float64,2}
    RF_S::Array{Float64,2}
    #
    ΔT_CO2av::Array{Float64,2}      #size = (totalimpact, nrun) 
    ΔT_CO2av_lcycle::Array{Float64,2}
    ΔT_Contrail::Array{Float64,2} 
    ΔT_O3short::Array{Float64,2}        
    ΔT_O3long::Array{Float64,2}        
    ΔT_Nitrate::Array{Float64,2}        
    ΔT_CH4::Array{Float64,2}      
    ΔT_H2O::Array{Float64,2} 
    ΔT_BC::Array{Float64,2} 
    ΔT_S::Array{Float64,2} 
    ΔT_background1::Array{Float64,2} 
    ΔT_background2::Array{Float64,2}
    ΔT_background2_lcycle::Array{Float64,2}
    #
    cost_CO2av::Array{Float64,2} #size = (totalimpact, nrun) 
    cost_Contrail::Array{Float64,2} 
    cost_O3short::Array{Float64,2}        
    cost_O3long::Array{Float64,2}        
    cost_Nitrate::Array{Float64,2}        
    cost_CH4::Array{Float64,2}        
    cost_H2O::Array{Float64,2}        
    cost_BC::Array{Float64,2}         
    cost_S::Array{Float64,2} 

    cost_US_CO2av::Array{Float64,2} #size = (totalimpact, nrun) 
    cost_US_Contrail::Array{Float64,2} 
    cost_US_O3short::Array{Float64,2}        
    cost_US_O3long::Array{Float64,2}        
    cost_US_Nitrate::Array{Float64,2}        
    cost_US_CH4::Array{Float64,2}        
    cost_US_H2O::Array{Float64,2}        
    cost_US_BC::Array{Float64,2}         
    cost_US_S::Array{Float64,2} 
    #
    discountedcost_CO2av::Array{Float64,2} #size = (totalimpact, nrun) 
    discountedcost_Contrail::Array{Float64,2} 
    discountedcost_O3short::Array{Float64,2}        
    discountedcost_O3long::Array{Float64,2}        
    discountedcost_Nitrate::Array{Float64,2}        
    discountedcost_CH4::Array{Float64,2}        
    discountedcost_H2O::Array{Float64,2}        
    discountedcost_BC::Array{Float64,2}         
    discountedcost_S::Array{Float64,2}   
    
    discountedcost_US_CO2av::Array{Float64,2} #size = (totalimpact, nrun) 
    discountedcost_US_Contrail::Array{Float64,2} 
    discountedcost_US_O3short::Array{Float64,2}        
    discountedcost_US_O3long::Array{Float64,2}        
    discountedcost_US_Nitrate::Array{Float64,2}        
    discountedcost_US_CH4::Array{Float64,2}        
    discountedcost_US_H2O::Array{Float64,2}        
    discountedcost_US_BC::Array{Float64,2}         
    discountedcost_US_S::Array{Float64,2}    
    #
    NPV_CO2av::Array{Float64,2}  #size = (1, nrun) 
    NPV_Contrail::Array{Float64,2} 
    NPV_O3short::Array{Float64,2}        
    NPV_O3long::Array{Float64,2}        
    NPV_Nitrate::Array{Float64,2}        
    NPV_CH4::Array{Float64,2}        
    NPV_H2O::Array{Float64,2}        
    NPV_BC::Array{Float64,2}         
    NPV_S::Array{Float64,2}    

    NPV_US_CO2av::Array{Float64,2}  #size = (1, nrun) 
    NPV_US_Contrail::Array{Float64,2} 
    NPV_US_O3short::Array{Float64,2}        
    NPV_US_O3long::Array{Float64,2}        
    NPV_US_Nitrate::Array{Float64,2}        
    NPV_US_CH4::Array{Float64,2}        
    NPV_US_H2O::Array{Float64,2}        
    NPV_US_BC::Array{Float64,2}         
    NPV_US_S::Array{Float64,2} 
end

Base.@kwdef mutable struct LCOutputs <: Outputs
    RF_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    RF_lc_N2O::Array{Float64,2}
    RF_lc_CH4::Array{Float64,2}
    RF_lc_CH4tropO3::Array{Float64,2}
    RF_lc_CH4stratH2O::Array{Float64,2}
    RF_lc_CH4_CO2::Array{Float64,2}
    #
    ΔT_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    ΔT_lc_N2O::Array{Float64,2}
    ΔT_lc_CH4::Array{Float64,2}
    ΔT_lc_CH4tropO3::Array{Float64,2}
    ΔT_lc_CH4stratH2O::Array{Float64,2}
    ΔT_lc_CH4_CO2::Array{Float64,2}
    #
    cost_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    cost_lc_N2O::Array{Float64,2}
    cost_lc_CH4::Array{Float64,2}
    cost_lc_CH4tropO3::Array{Float64,2}
    cost_lc_CH4stratH2O::Array{Float64,2}
    cost_lc_CH4_CO2::Array{Float64,2}
    
    cost_US_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    cost_US_lc_N2O::Array{Float64,2}
    cost_US_lc_CH4::Array{Float64,2}
    cost_US_lc_CH4tropO3::Array{Float64,2}
    cost_US_lc_CH4stratH2O::Array{Float64,2}
    cost_US_lc_CH4_CO2::Array{Float64,2}
    #
    discountedcost_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    discountedcost_lc_N2O::Array{Float64,2}
    discountedcost_lc_CH4::Array{Float64,2}
    discountedcost_lc_CH4tropO3::Array{Float64,2}
    discountedcost_lc_CH4stratH2O::Array{Float64,2}
    discountedcost_lc_CH4_CO2::Array{Float64,2}
    
    discountedcost_US_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    discountedcost_US_lc_N2O::Array{Float64,2}
    discountedcost_US_lc_CH4::Array{Float64,2}
    discountedcost_US_lc_CH4tropO3::Array{Float64,2}
    discountedcost_US_lc_CH4stratH2O::Array{Float64,2}
    discountedcost_US_lc_CH4_CO2::Array{Float64,2}
    #
    NPV_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    NPV_lc_N2O::Array{Float64,2}
    NPV_lc_CH4::Array{Float64,2}
    NPV_lc_CH4tropO3::Array{Float64,2}
    NPV_lc_CH4stratH2O::Array{Float64,2}
    NPV_lc_CH4_CO2::Array{Float64,2}
    
    NPV_US_lc_CO2::Array{Float64,2}      #size = (totalimpact, nrun)
    NPV_US_lc_N2O::Array{Float64,2}
    NPV_US_lc_CH4::Array{Float64,2}
    NPV_US_lc_CH4tropO3::Array{Float64,2}
    NPV_US_lc_CH4stratH2O::Array{Float64,2}
    NPV_US_lc_CH4_CO2::Array{Float64,2}
    #
end


function MainOutputs(simulation::Simulations_orig)
    nrun = simulation.nrun;
    totalimpact = simulation.totalimpact;

    mainoutput =  MainOutputs(
                #
                RF_CO2av = zeros(totalimpact, nrun),
                RF_CO2av_lcycle = zeros(totalimpact, nrun),
                RF_Contrail = zeros(totalimpact, nrun), 
                RF_O3short = zeros(totalimpact, nrun),
                RF_O3long = zeros(totalimpact, nrun),
                RF_Nitrate = zeros(totalimpact, nrun),
                RF_CH4 = zeros(totalimpact, nrun),
                RF_H2O = zeros(totalimpact, nrun),
                RF_BC = zeros(totalimpact, nrun),
                RF_S = zeros(totalimpact, nrun),
                #
                ΔT_CO2av = zeros(totalimpact, nrun),
                ΔT_CO2av_lcycle = zeros(totalimpact, nrun),
                ΔT_Contrail = zeros(totalimpact, nrun),
                ΔT_O3short = zeros(totalimpact, nrun),
                ΔT_O3long = zeros(totalimpact, nrun),
                ΔT_Nitrate = zeros(totalimpact, nrun),
                ΔT_CH4 = zeros(totalimpact, nrun),
                ΔT_H2O = zeros(totalimpact, nrun),
                ΔT_BC = zeros(totalimpact, nrun),
                ΔT_S = zeros(totalimpact, nrun),
                ΔT_background1 = zeros(totalimpact, nrun),
                ΔT_background2 = zeros(totalimpact, nrun),
                ΔT_background2_lcycle = zeros(totalimpact, nrun),
                #
                cost_CO2av = zeros(totalimpact, nrun),
                cost_Contrail = zeros(totalimpact, nrun),
                cost_O3short = zeros(totalimpact, nrun),
                cost_O3long = zeros(totalimpact, nrun),
                cost_Nitrate = zeros(totalimpact, nrun),
                cost_CH4 = zeros(totalimpact, nrun),
                cost_H2O = zeros(totalimpact, nrun),
                cost_BC = zeros(totalimpact, nrun),
                cost_S = zeros(totalimpact, nrun),
                #
                cost_US_CO2av = zeros(totalimpact, nrun),
                cost_US_Contrail = zeros(totalimpact, nrun),
                cost_US_O3short = zeros(totalimpact, nrun),
                cost_US_O3long = zeros(totalimpact, nrun),
                cost_US_Nitrate = zeros(totalimpact, nrun),
                cost_US_CH4 = zeros(totalimpact, nrun),
                cost_US_H2O = zeros(totalimpact, nrun),
                cost_US_BC = zeros(totalimpact, nrun),
                cost_US_S = zeros(totalimpact, nrun),
                #
                discountedcost_CO2av = zeros(totalimpact, nrun),
                discountedcost_Contrail = zeros(totalimpact, nrun),
                discountedcost_O3short = zeros(totalimpact, nrun),
                discountedcost_O3long = zeros(totalimpact, nrun),
                discountedcost_Nitrate = zeros(totalimpact, nrun),
                discountedcost_CH4 = zeros(totalimpact, nrun),
                discountedcost_H2O = zeros(totalimpact, nrun),
                discountedcost_BC = zeros(totalimpact, nrun),
                discountedcost_S = zeros(totalimpact, nrun),
                #
                discountedcost_US_CO2av = zeros(totalimpact, nrun),
                discountedcost_US_Contrail = zeros(totalimpact, nrun),
                discountedcost_US_O3short = zeros(totalimpact, nrun),
                discountedcost_US_O3long = zeros(totalimpact, nrun),
                discountedcost_US_Nitrate = zeros(totalimpact, nrun),
                discountedcost_US_CH4 = zeros(totalimpact, nrun),
                discountedcost_US_H2O = zeros(totalimpact, nrun),
                discountedcost_US_BC = zeros(totalimpact, nrun),
                discountedcost_US_S = zeros(totalimpact, nrun),
                #
                NPV_CO2av = zeros(1, nrun),
                NPV_Contrail = zeros(1, nrun),
                NPV_O3short = zeros(1, nrun),
                NPV_O3long = zeros(1, nrun),
                NPV_Nitrate = zeros(1, nrun),
                NPV_CH4 = zeros(1, nrun),
                NPV_H2O = zeros(1, nrun),
                NPV_BC = zeros(1, nrun),
                NPV_S = zeros(1, nrun),
                #
                NPV_US_CO2av = zeros(1, nrun),
                NPV_US_Contrail = zeros(1, nrun),
                NPV_US_O3short = zeros(1, nrun),
                NPV_US_O3long = zeros(1, nrun),
                NPV_US_Nitrate = zeros(1, nrun),
                NPV_US_CH4 = zeros(1, nrun),
                NPV_US_H2O = zeros(1, nrun),
                NPV_US_BC = zeros(1, nrun),
                NPV_US_S = zeros(1, nrun),
                )
    return mainoutput
end

function LCOutputs(simulation::Simulations_orig)
    nrun = simulation.nrun;
    totalimpact = simulation.totalimpact;

    lcoutput =  LCOutputs(
        RF_lc_CO2 = zeros(totalimpact, nrun),
        RF_lc_N2O = zeros(totalimpact, nrun),
        RF_lc_CH4 = zeros(totalimpact, nrun),
        RF_lc_CH4tropO3 = zeros(totalimpact, nrun),
        RF_lc_CH4stratH2O = zeros(totalimpact, nrun),
        RF_lc_CH4_CO2 = zeros(totalimpact, nrun),
        #
        ΔT_lc_CO2 = zeros(totalimpact, nrun),
        ΔT_lc_N2O = zeros(totalimpact, nrun),
        ΔT_lc_CH4 = zeros(totalimpact, nrun),
        ΔT_lc_CH4tropO3 = zeros(totalimpact, nrun),
        ΔT_lc_CH4stratH2O = zeros(totalimpact, nrun),
        ΔT_lc_CH4_CO2 = zeros(totalimpact, nrun),
        #
        cost_lc_CO2 = zeros(totalimpact, nrun),
        cost_lc_N2O = zeros(totalimpact, nrun),
        cost_lc_CH4 = zeros(totalimpact, nrun),
        cost_lc_CH4tropO3 = zeros(totalimpact, nrun),
        cost_lc_CH4stratH2O = zeros(totalimpact, nrun),
        cost_lc_CH4_CO2 = zeros(totalimpact, nrun),

        cost_US_lc_CO2 = zeros(totalimpact, nrun),
        cost_US_lc_N2O = zeros(totalimpact, nrun),
        cost_US_lc_CH4 = zeros(totalimpact, nrun),
        cost_US_lc_CH4tropO3 = zeros(totalimpact, nrun),
        cost_US_lc_CH4stratH2O = zeros(totalimpact, nrun),
        cost_US_lc_CH4_CO2 = zeros(totalimpact, nrun),
        #
        discountedcost_lc_CO2 = zeros(totalimpact, nrun),
        discountedcost_lc_N2O = zeros(totalimpact, nrun),
        discountedcost_lc_CH4 = zeros(totalimpact, nrun),
        discountedcost_lc_CH4tropO3 = zeros(totalimpact, nrun),
        discountedcost_lc_CH4stratH2O = zeros(totalimpact, nrun),
        discountedcost_lc_CH4_CO2 = zeros(totalimpact, nrun),

        discountedcost_US_lc_CO2 = zeros(totalimpact, nrun),
        discountedcost_US_lc_N2O = zeros(totalimpact, nrun),
        discountedcost_US_lc_CH4 = zeros(totalimpact, nrun),
        discountedcost_US_lc_CH4tropO3 = zeros(totalimpact, nrun),
        discountedcost_US_lc_CH4stratH2O = zeros(totalimpact, nrun),
        discountedcost_US_lc_CH4_CO2 = zeros(totalimpact, nrun),
        #
        NPV_lc_CO2 = zeros(1, nrun),
        NPV_lc_N2O = zeros(1, nrun),
        NPV_lc_CH4 = zeros(1, nrun),
        NPV_lc_CH4tropO3 = zeros(1, nrun),
        NPV_lc_CH4stratH2O = zeros(1, nrun),
        NPV_lc_CH4_CO2 = zeros(1, nrun),
        
        NPV_US_lc_CO2 = zeros(1, nrun),
        NPV_US_lc_N2O = zeros(1, nrun),
        NPV_US_lc_CH4 = zeros(1, nrun),
        NPV_US_lc_CH4tropO3 = zeros(1, nrun),
        NPV_US_lc_CH4stratH2O = zeros(1, nrun),
        NPV_US_lc_CH4_CO2 = zeros(1, nrun),
    )
    return lcoutput
end