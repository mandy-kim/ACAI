Base.@kwdef struct EmissionMultipliers 
    FB::Float64
    DIST::Float64
    NO::Float64
    HONO::Float64
    NO2::Float64
    CO::Float64
    OC::Float64
    BC::Float64
    HC::Float64
end

Base.@kwdef struct EmissionsSpecies
    FB::Array{Float64} #kg
    DIST::Array{Float64} #km
    CO2::Array{Float64} #kg 
    H2O::Array{Float64} #kg
    S::Array{Float64} #kg
    NOx::Array{Float64} #kg
    CO::Array{Float64} #kg
    OC::Array{Float64} #kg
    BC::Array{Float64} #kg
    HC::Array{Float64} #kg
end

## create EmissionsSpecies_Struct for original approach, ::Float64 NOT Array

Base.@kwdef struct Emissions
    lon::Array{Float64}
    lat::Array{Float64}
    alt::Array{Float64}
    area::Array{Float64}
    yearly_data::Dict{Int64,EmissionsSpecies} #Integer = year (e.g. 2020)
end

function createEmissionsSpecies(ds_transformed::NCDataset,approach)
    grid_areas = float.(ds_transformed["area"][:]);
    CO2 = float.(ds_transformed["CO2"][:] .* grid_areas);

    FB = float.(ds_transformed["FUELBURN"][:] .* grid_areas);
    DIST = float.(ds_transformed["DISTANCE"][:] .* grid_areas);
    H2O = float.(ds_transformed["H2O"][:] .* grid_areas);
    S = float.(ds_transformed["S"][:] .* grid_areas);
    NOx = float.(ds_transformed["NOx"][:] .* grid_areas);

    if approach == "original"
        #only need global sums
        FB = np.sum(FB)
        DIST = np.sum(DIST)
        H2O = np.sum(H2O)
        S = np.sum(S)
        NOx = np.sum(NOx)
        CO = 0
        OC = 0
        BC = 0
        HC = 0
    else
        # "advanced" approach, keep spatial information
        CO = ds_transformed["CO"][:] .* grid_areas;
        OC = ds_transformed["OC"][:] .* grid_areas;
        BC = ds_transformed["BC"][:] .* grid_areas;
        HC = ds_transformed["HC"][:] .* grid_areas;
    end
    
    #convert (kg/m2 or km/m2) to (kg or km)
    

    emissions_species = EmissionsSpecies(FB=FB, DIST=DIST, CO2=CO2, H2O=H2O, S=S, 
                                    NOx=NOx, CO=CO, OC=OC, BC=BC, HC=HC)
    return emissions_species
end

