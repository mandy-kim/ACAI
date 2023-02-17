#RFs for NOx emission species
function calculate_RF_NOx(i::Int64, emission::Emissions_orig, maindistribution::MainDistributions, simulation::Simulations_orig)
    #Sample uncertain variables
    totalimpact = simulation.totalimpact; #years
    avNOx2006 = sample(maindistribution.avNOx2006,i); #Tg
    λ_Nitrate = sample(maindistribution.lambda_Nitrate,i);
    λ_O3 = sample(maindistribution.lambda_O3,i);
    λ_CH4 = sample(maindistribution.lambda_CH4,i);
    RF_Nitrate_i = sample(maindistribution.RF_Nitrate,i); #W/m2
    NOx_Effect = sample(maindistribution.NOx_Effect,i);

    #Short term Nitrate RF - scale by NOx inventory
    RF_Nitrate = λ_Nitrate .* RF_Nitrate_i .* (emission.NOx .* sample(maindistribution.emis_unc_NOx,i) ./ avNOx2006); #W/m2

    #NOx related effects: [1: Stevenson et al, 2: Wild, 3: Hoor et al] Wild et al values are shown in Stevenson et al.
    efold_matrix        = [11.53, 11.8, 10.7];        #Unit: year
    RFyr_CH4_matrix     = [-5.2, -5.7, -5.3]*1e-12;   #Unit: Wyr/m2/(kg NOx)
    RFyr_O3long_matrix  = [-0.95, -1.5, -1.8]*1e-12;  #Unit: Wyr/m2/(kg NOx)
    RFyr_O3short_matrix = [5.06,  7.9,  7.4]*1e-12;  #Unit: Wyr/m2/(kg NOx)

    #Short term O3 response
    RF_O3short = λ_O3 .* RFyr_O3short_matrix[NOx_Effect] .* emission.NOx;  #W/m2
    
    #Long term O3 and CH4 responses
    decay_constant = efold_matrix[NOx_Effect];
    decay_time = collect(1:totalimpact) .- 1;
    exp_decay = exp.(-decay_time ./decay_constant);

    RF_CH4_scale = λ_CH4 .* RFyr_CH4_matrix[NOx_Effect] ./ sum(exp_decay) .* exp_decay;
    RF_O3long_scale = λ_O3 .* RFyr_O3long_matrix[NOx_Effect] ./ sum(exp_decay) .* exp_decay;

    RF_CH4_temp = zeros(totalimpact, totalimpact);
    RF_O3long_temp = zeros(totalimpact, totalimpact);

    for yr in 1:totalimpact
        RF_CH4_temp[yr, yr:end] = RF_CH4_scale[yr] .* emission.NOx[1:totalimpact+1-yr]
        RF_O3long_temp[yr, yr:end] = RF_O3long_scale[yr] .* emission.NOx[1:totalimpact+1-yr]
    end

    RF_CH4 = sum(RF_CH4_temp,dims=2)
    RF_O3long = sum(RF_O3long_temp,dims=2)

    return RF_Nitrate, RF_O3short, RF_O3long, RF_CH4
end

#RFs for all non-CO2, non-NOx emission species
function calculate_RF_short(i::Int64, emission::Emissions_orig, maindistribution::MainDistributions, simulation::Simulations_orig)
    #Sample uncertain variables
    totalimpact = simulation.totalimpact; #years
    avfuel2005 = sample(maindistribution.avfuel2005,i); #Tg

    λ_Contrail = sample(maindistribution.lambda_Contrail,i);
    λ_H2O = sample(maindistribution.lambda_H2O,i);
    λ_BC = sample(maindistribution.lambda_Soot,i);
    λ_S = sample(maindistribution.lambda_Sulfates,i);

    RF_Contrail_i = sample(maindistribution.RF_Contrail,i);
    RF_H2O_i  = sample(maindistribution.RF_H2O,i);
    RF_BC_i = sample(maindistribution.RF_Soot,i);
    RF_S_i = sample(maindistribution.RF_Sulfates,i);

    #Calculate RFs based on sampled variables
    RF_Contrail = λ_Contrail * RF_Contrail_i .* emission.FB ./ (avfuel2005 * 1e9)
    RF_H2O = λ_H2O * RF_H2O_i .* emission.FB ./ (avfuel2005 * 1e9)
    RF_BC = λ_BC * RF_BC_i .* emission.FB ./ (avfuel2005 * 1e9)
    RF_S = λ_S * RF_S_i .* emission.FB ./ (avfuel2005 * 1e9)

    return RF_Contrail, RF_H2O, RF_BC, RF_S
end

