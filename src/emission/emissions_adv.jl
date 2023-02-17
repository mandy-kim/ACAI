include("emissionStructs.jl")

"""
# totalEmisSpecies
### Inputs
- species: String e.g., "CO2", "DIST", or "NOx"
- emission: Emissions type variable
- simulation: Simulation type variable containing fixed data from user input
- region: optional input. 'global', 'NH' (Northern Hemisphere), or 'SH'.
- alt: optional input. 'total' (default), 'LTO' (under 8km), or 'cruise' (above 8km)

### Output
- total_emis_array: array of length simulation.totalimpact. total_emis_array[1] = emissions in year 1
"""
function totalEmisSpecies(species::String, emission::Emissions, simulation::Simulations_adv; region::String="global", alt::String="total")
    totalimpact = simulation.totalimpact
    total_emis_array = zeros(totalimpact)
    years = simulation.emissionyears
    ni = 0
    if region == "global"
        lats = 1:91
    elseif region == "NH"
        lats = 47:91
    elseif region == "SH"
        lats = 1:46
    end
    
    if alt == "total"
        alts = 1:length(emission.alt)
    elseif alt == "LTO"
        alts = 1:27
    elseif alt == "cruise"
        alts = 28:length(emission.alt)
    end
    
    for yeari in years
        ni+=1
        total_emis_array[ni] = sum(getfield(emission.yearly_data[yeari], Symbol(species))[:,lats,alts]) #kg
    end
    return total_emis_array
end

function getUnitMultipliers(ds::NCDataset, species::String, g_or_m::String)
    try
        if ds[species].attrib["units"] == "k"*"$g_or_m"*"/m2/s"
            return 1
        elseif ds[species].attrib["units"] == "$g_or_m"*"/m2/s"
            return 1e-3
        else
            error("units for "*"$species"*" must be in k"*"$g_or_m"*"/m2/s or "*"$g_or_m"*"/m2/s")
        end
    catch
        println("no units for "*"$species"*" specified, assumed k"*"$g_or_m"*"/m2/s")
        return 1
    end
end

function getEmissionMultipliers(folder::String)
    file_sample = readdir(folder)[1];
    ds = NCDataset(folder*file_sample);

    FB_mult = getUnitMultipliers(ds, "FUELBURN", "g")
    DIST_mult = getUnitMultipliers(ds, "DISTANCE", "m")
    NO_mult = getUnitMultipliers(ds, "NO", "g")
    HONO_mult = getUnitMultipliers(ds, "HONO", "g")
    NO2_mult = getUnitMultipliers(ds, "NO2", "g")
    CO_mult = getUnitMultipliers(ds, "CO", "g")
    BC_mult = getUnitMultipliers(ds, "BC", "g")
    HC_mult = getUnitMultipliers(ds, "HC", "g")
    OC_mult = getUnitMultipliers(ds, "OC", "g")

    emis_mult = EmissionMultipliers(
                    FB=FB_mult,
                    DIST=DIST_mult,
                    NO=NO_mult,
                    HONO=HONO_mult,
                    NO2=NO2_mult,
                    CO=CO_mult,
                    OC=OC_mult,
                    BC=BC_mult,
                    HC=HC_mult,
                    )

    return emis_mult
end

function getDaysArray(simulation::Simulations_adv)
    frequency = simulation.frequency;
    folders = simulation.emissionfolders;
    days = []
    if frequency == "monthly"
        for i in 1:length(folders)
            append!(days,[[31,28,31,30,31,30,31,31,30,31,30,31]])
        end
    elseif frequency == "daily"
        for i in 1:length(folders)
            append!(days,[ones(365)])
        end
    elseif frequency == "yearly"
        for i in 1:length(folders)
            append!(days,[[365]])
        end
    end
    return days
end

function getFileDims(folder::String)
    file_sample = readdir(folder)[1];
    ds_sample = NCDataset(folder*file_sample);
    dim_dict = ds_sample.dim; #lon x lat x lev x time
    long_dim = dim_dict["lon"];
    lat_dim = dim_dict["lat"];
    lev_dim = dim_dict["lev"];
    return (lon=long_dim, lat=lat_dim, lev=lev_dim)
end

function createOriginalDS(emi_folder_i::String,dims::NamedTuple,fb_array::AbstractArray,dist_array::AbstractArray,NOx_array::AbstractArray,CO_array::AbstractArray,OC_array::AbstractArray,BC_array::AbstractArray,HC_array::AbstractArray)
    EI_CO2 = 3.15
    EI_H2O = 1.26
    EI_S = 6e-4

    CO2_array = fb_array .* EI_CO2;
    H2O_array = fb_array .* EI_H2O;
    S_array = fb_array .* EI_S;

    #Create netcdf file of original resolution
    ds_original = NCDataset(emi_folder_i*"original.nc","c");
    defDim(ds_original,"lon",dims.lon)
    defDim(ds_original,"lat",dims.lat)
    defDim(ds_original,"lev",dims.lev)
    
    v = defVar(ds_original,"FUELBURN",Float64,("lon","lat","lev"));
    v[:,:,:] = fb_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"DISTANCE",Float64,("lon","lat","lev"));
    v[:,:,:] = dist_array;
    v.attrib["units"]="km/m2"

    v = defVar(ds_original,"NOx",Float64,("lon","lat","lev"));
    v[:,:,:] = NOx_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"CO",Float64,("lon","lat","lev"));
    v[:,:,:] = CO_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"OC",Float64,("lon","lat","lev"));
    v[:,:,:] = OC_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"BC",Float64,("lon","lat","lev"));
    v[:,:,:] = BC_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"HC",Float64,("lon","lat","lev"));
    v[:,:,:] = HC_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"CO2",Float64,("lon","lat","lev"));
    v[:,:,:] = CO2_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"H2O",Float64,("lon","lat","lev"));
    v[:,:,:] = H2O_array;
    v.attrib["units"]="kg/m2"

    v = defVar(ds_original,"S",Float64,("lon","lat","lev"));
    v[:,:,:] = S_array;
    v.attrib["units"]="kg/m2"

    close(ds_original)
end

py"""
import xarray as xr
import sparselt
import sparselt.xr
import sparselt.esmf
import sparselt.linear_transform

def transform_nc(original_path,regrid_weights,out_template_path,out_path):
    ds1 = xr.open_dataset(original_path)
    out_template = xr.open_dataset(out_template_path)

    transform = sparselt.esmf.load_weights(regrid_weights,
        input_dims=[('lat','lon'),(ds1.dims['lat'],ds1.dims['lon'])],
        output_dims=[('lat_b','lon_b'),(out_template.dims['lat'],out_template.dims['lon'])])
    
    ds2 = sparselt.xr.apply(transform,ds1,out_template)
    ds2.to_netcdf(out_path+"transformed.nc")

"""

function getAltArray(layers, len)
	if layers == 72
		alt = [78.14,74.594,71.812,69.440,67.243,65.115,63.004,60.902,58.816,56.752,54.717,52.716,
		50.754,48.835,46.962,45.134,43.355,41.627,39.951,38.328,36.759,35.244,33.782,32.372,
		31.015,29.701,28.423,27.180,25.971,24.794,23.648,22.531,21.438,20.364,19.309,18.269,
		17.243,16.222,15.198,14.170,13.134,12.086,11.021,9.936,8.846,7.943,7.237,6.585,5.980,
		5.413,4.879,4.375,3.896,3.439,3.074,2.792,2.517,2.249,1.988,1.759,1.584,1.436,1.290,
		1.146,1.004,0.864,0.726,0.589,0.454,0.320,0.189,0.058];
	elseif layers == 47
		alt = [72.180,63.053,54.834,47.135,40.166,34.024,28.654,25.307,23.020,20.836,18.727,
		17.243,16.222,15.198,14.170,13.134,12.086,11.021,9.936,8.846,7.943,7.237,6.585,5.980,
		5.413,4.879,4.375,3.896,3.439,3.074,2.792,2.517,2.249,1.988,1.759,1.584,1.436,1.290,
		1.146,1.004,0.864,0.726,0.589,0.454,0.320,0.189,0.058];
    else
        error("layers must be 47 or 72 GEOS-Chem vertical layers")
    end

    alt = alt[end:-1:1];
    alt = alt[1:len];

    return alt
end

function getYearlyEmis(simulation::Simulations_adv)
    foldername = simulation.emissionname;
    days = getDaysArray(simulation);
    sec_in_day = 24*3600
    
    grid_areas = missing;
    lon = missing;
    lat = missing;
    alt = missing;

    yearly_data = Dict();

    try mkdir("output/emissions/")
    catch
    end
    
    if isdir("output/emissions/"*"$foldername/")
        println("Warning: folders for emissions "*"(output/emissions/"*"$foldername/)"*" already exist, reading existing files")

        for year_i in simulation.emissionyears
            emi_folder_i = "output/emissions/"*"$foldername/"*"$year_i/"
            ds_transformed = NCDataset(emi_folder_i*"transformed.nc");

            if year_i == simulation.emissionyears[1]
                grid_areas = float.(ds_transformed["area"][:]);
                lon = float.(ds_transformed["lon"][:]);
                lat = float.(ds_transformed["lat"][:]);
                alt_len = Int.(ds_transformed.dim["lev"]);
                alt = getAltArray(simulation.layers, alt_len);
            end
            emissions_year_i = createEmissionsSpecies(ds_transformed, simulation.approach)
            yearly_data[year_i] = emissions_year_i 
        end
        emission = Emissions(lon, lat, alt, grid_areas, yearly_data)
        return emission

    end

    nfol = 0
    for folder_i in simulation.emissionfolders
        #each folder corresponds to year_i in years
        nfol += 1
        year_i = simulation.emissionyears[nfol]
        dims = getFileDims(folder_i)
        emi_mult = getEmissionMultipliers(folder_i)

        if nfol == 1
            mkdir("output/emissions/"*"$foldername/")
        end
        emi_folder_i = "output/emissions/"*"$foldername/"*"$year_i/"
        mkdir(emi_folder_i)

        All_d_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #km/m2/s
        All_fb_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #kg/m2/s
        All_NOx_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #kg/m2/s
        All_CO_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #kg/m2/s
        All_OC_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #kg/m2/s
        All_BC_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #kg/m2/s
        All_HC_per_m2 = zeros(dims.lon, dims.lat, dims.lev); #kg/m2/s

        # Aggregate data to total annual in km/m2 or kg/m2
        nfil = 0
        for file_i in readdir(folder_i)
            nfil += 1;
            ds = NCDataset(folder_i*file_i);
            try
                All_d_per_m2 = All_d_per_m2 .+ ds["DISTANCE"][:,:,:,1] .* emi_mult.DIST .* (days[nfol][nfil]*sec_in_day);
            catch
            end
            try
                All_fb_per_m2 = All_fb_per_m2 .+ ds["FUELBURN"][:,:,:,1] .* emi_mult.FB .* (days[nfol][nfil]*sec_in_day);
            catch
            end
            try
                All_NOx_per_m2 = All_NOx_per_m2 .+ (ds["NO"][:,:,:,1].*emi_mult.NO.*(46/30) .+ ds["NO2"][:,:,:,1].*emi_mult.NO2 .+ ds["HONO"][:,:,:,1].*emi_mult.HONO.*(46/47)) .* (days[nfol][nfil]*sec_in_day);
            catch
            end
            try
                All_CO_per_m2 = All_CO_per_m2 .+ ds["CO"][:,:,:,1] .*emi_mult.CO .* (days[nfol][nfil]*sec_in_day);
            catch
            end
            try
                All_OC_per_m2 = All_OC_per_m2 .+ ds["OC"][:,:,:,1] .*emi_mult.OC .* (days[nfol][nfil]*sec_in_day);
            catch
            end
            try
                All_BC_per_m2 = All_BC_per_m2 .+ ds["BC"][:,:,:,1] .*emi_mult.BC .* (days[nfol][nfil]*sec_in_day);
            catch
            end
            try
                All_HC_per_m2 = All_HC_per_m2 .+ ds["HC"][:,:,:,1] .*emi_mult.HC .* (days[nfol][nfil]*sec_in_day);
            catch
            end 
        end

        #Save aggregated annual data in original resolution
        createOriginalDS(emi_folder_i, dims,
                    All_fb_per_m2, All_d_per_m2,All_NOx_per_m2,
                    All_CO_per_m2,All_OC_per_m2,All_BC_per_m2,All_HC_per_m2);

        #Transform and save netcdf data to 2x2.5 (specified by "weights" input variable)
        out_template_path = "regrid_files/regular_lat_lon_91x144.nc";
        py"transform_nc"(emi_folder_i*"original.nc",simulation.weights,out_template_path,emi_folder_i);
        ds_transformed = NCDataset(emi_folder_i*"transformed.nc");
        if nfol == 1
            grid_areas = float.(ds_transformed["area"][:]);
            lon = float.(ds_transformed["lon"][:]);
            lat = float.(ds_transformed["lat"][:]);
            alt_len = Int.(ds_transformed.dim["lev"]);
            alt = getAltArray(simulation.layers, alt_len);
        end

        emissions_year_i = createEmissionsSpecies(ds_transformed, simulation.approach)
        yearly_data[year_i] = emissions_year_i        
    end

    emission = Emissions(lon, lat, alt, grid_areas, yearly_data)
    return emission
end

