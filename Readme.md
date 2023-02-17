# ACAI.jl
Aviation Climate and Air quality Impacts

Julia based model to evaluate the environmental impacts of aviation emissions. Includes two approaches:

1. original (legacy APMT-IC), and 
2. advanced. 

The original approach estimates climate impacts of subsonic aviation based on scaled consensus values. The advanced approach estimates climate, air quality, and ozone impacts of any emissions inventory based on GEOS-Chem emission sensitivities. 

The current model for the advanced approach only includes impacts from CO<sub>2</sub>, contrails, and NO<sub>x</sub>. Future versions to include sensitivities for all non-CO<sub>2</sub>, non-contrail species.

## Getting started

There are two example "run" files based on the approach to get started

1. `example_run_ACAI_adv.jl` : run advanced approach
2. `example_run_ACAI_orig.jl` : run original approach

These examples specify the input parameters saved in the input directory. 

The ACAI technical document, as well as more information about the input parameters is saved in the LAE research folder [here](https://mitprod.sharepoint.com/:f:/s/LAE/Egv3uzx0WQtDoNP6JEVzsl4BSV0wbLNBNnM8SheDej-nlw?e=BZ67wj).

The notebook in the postprocess folder shows examples of plotting the various outputs.

## Notes
The resolution of the input emissions inventory for the advanced approach is flexible since the model regrids the files into 2°x2.5° resolution. However, this is done within the model using the following python packages:

1. [xarray](https://docs.xarray.dev/en/stable/getting-started-guide/installing.html#)
2. [sparselt](https://github.com/LiamBindle/sparselt)

To run ACAI, make sure that these python modules (and their dependencies) are available in the python environment you are using (see the PyCall [README](https://github.com/JuliaPy/PyCall.jl)). 