include("sensitivityStructs.jl")

#read NOx sensitivity files
function getNOxSensitivities(simulation::Simulations)
    n_case = simulation.sensitivitycases

    NOx_sensitivity_dict = Dict{Int64,SensitivityRegionImpacts}()
    NOx_folder = simulation.sensitivityfolders[1]

    for i in 1:n_case
        ds = NCDataset(NOx_folder*"/sensitivity_"*string(i)*".nc")
        RF = ds["RF"][:]
        PM25 = ds["PM25"][:]
        Sfc_O3 = ds["Sfc_O3"][:]
        Col_O3 = ds["Col_O3"][:]
        Zonal_O3 = ds["Zonal_O3"][:]
        mask = ds["mask"][:]
        lev_47 = ds["lev_47"][:]
        lev_72 = ds["lev_72"][:]
        
        NOx_sensitivity_dict[i] = SensitivityRegionImpacts(RF, PM25, Sfc_O3, Col_O3, Zonal_O3, mask, lev_47, lev_72)
    end

    NOx_sens = EmissionSensitivities(n_case, NOx_sensitivity_dict);
    return NOx_sens
end

#calculate NOx impacts based on saved sensitivity files and input emissions inventory
function calculateNOxImpacts(NOx_sens::EmissionSensitivities, climateoutput::ClimateOutputs, aqoutput::AirQualityOutputs, simulation::Simulations, emission::Emissions, healthconsts::HealthConstants)
    n_RFyears = simulation.sensitivityRFyears
    n_AQyears = simulation.sensitivityAQyears
    n_exposure_years = healthconsts.n_exposure_years
    n_O3years = simulation.sensitivityO3years
    n_cases = NOx_sens.n_case #46 for now: 1-40 = cruise, 41-46 = LTO (under 8km)

    multipliers = zeros(n_cases, length(simulation.emissionyears))
    RF_NOx = zeros(simulation.totalimpact);
    PM25_NOx = zeros(144,91,n_exposure_years);
    Sfc_O3_NOx = zeros(144,91,n_exposure_years);
    Col_O3_NOx = zeros(144,91,n_O3years + length(simulation.emissionyears) - 1);
    Zonal_O3_NOx = zeros(72,91,n_O3years + length(simulation.emissionyears) - 1);

    emissionyears = simulation.emissionyears;
    iyear = 0
    for yeari in emissionyears
        iyear+=1
        for i in 1:n_cases
            mask = NOx_sens.sensitivity_dict[i].mask
            if simulation.layers == 47
                lev = NOx_sens.sensitivity_dict[i].lev_47
            elseif simulation.layers == 72
                lev = NOx_sens.sensitivity_dict[i].lev_72
            end
            
            nox_temp = try emission.yearly_data[yeari].NOx[:,:,lev]
            catch
                0
            end
            
            multiplier = sum(nox_temp .* mask) * 1e-6
            multipliers[i, iyear] = multiplier
            RF_NOx[iyear:n_RFyears+iyear-1] .+= multiplier .* NOx_sens.sensitivity_dict[i].RF
            PM25_NOx[:,:,iyear:n_AQyears+iyear-1] .+= multiplier .* NOx_sens.sensitivity_dict[i].PM25
            Sfc_O3_NOx[:,:,iyear:n_AQyears+iyear-1] .+= multiplier .* NOx_sens.sensitivity_dict[i].Sfc_O3
            Col_O3_NOx[:,:,iyear:n_O3years+iyear-1] .+= multiplier .* NOx_sens.sensitivity_dict[i].Col_O3
            Zonal_O3_NOx[:,:,iyear:n_O3years+iyear-1] .+= multiplier .* NOx_sens.sensitivity_dict[i].Zonal_O3
        end
    end

    NOx_impacts = TotalImpacts(multipliers = multipliers, RF = RF_NOx, PM25 = PM25_NOx, Sfc_O3 = Sfc_O3_NOx, Col_O3 = Col_O3_NOx, Zonal_O3=Zonal_O3_NOx);
    #save in output variable, since impacts calculated based on sensitivity data don't change (for now - might want to incorporate uncertainties later)
    climateoutput.RF_NOx[:,1] = RF_NOx;
    aqoutput.ΔX_PM25_NOx = PM25_NOx;
    aqoutput.ΔX_O3_NOx = Sfc_O3_NOx;
    aqoutput.Col_O3_NOx = Col_O3_NOx;
    aqoutput.Zonal_O3_NOx = Zonal_O3_NOx;

    #calculate country-average ΔXs
    countries = healthconsts.countries
    ncountries = length(countries)
    pop_dict = healthconsts.pop_dict
    n_exposure_years = healthconsts.n_exposure_years
    PM25_NOx_countryavg = zeros(ncountries,n_exposure_years)
    O3_NOx_countryavg = zeros(ncountries,n_exposure_years)
    for i in 1:ncountries
        country = countries[i]
        paff_grid = pop_dict[country] 

        PM25_NOx_countryavg[i,:] = reshape(sum(PM25_NOx .* paff_grid, dims=(1,2)) / sum(paff_grid), (1, n_exposure_years))
        O3_NOx_countryavg[i,:] = reshape(sum(Sfc_O3_NOx .* paff_grid, dims=(1,2)) / sum(paff_grid), (1, n_exposure_years))
    end
    aqoutput.ΔX_PM25_avg_NOx = PM25_NOx_countryavg
    aqoutput.ΔX_O3_avg_NOx = O3_NOx_countryavg

    return NOx_impacts, climateoutput, aqoutput
end