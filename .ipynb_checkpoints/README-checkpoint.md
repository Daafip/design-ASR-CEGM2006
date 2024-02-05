# design of an Aquifer Storage Recharge (ASR) system TUDelft CEGM2006
Designing an ASR system for the cross over course Subsurface storage for climate, energy and water course (CEGM2006)
# Explanation of repo

- `bin` contains required moflow6 executable.
- `Model` contains all the code.
    - `1. basic_model_buoy_ml.ipynb` contains the initial model created with boyancy, the following are essentially copies with just a differen run.  
        - `2. optimise_nl_basic_model_buoy_ml.ipynb` runs code for different layers to check behavoir and run time. 
        - `3. vary_KH_buoy_ml.ipynb` is the start of the sensitivity analysis. Here different kh values are experimented with.
        - `4. worstcase_model_buoy_ml.ipynb` looks at a different scenario, not per sé _the_ worst case . 
    - `vary_params*.py` contains the sensitivity analysis for k and npor, MP = multiprocessing - runs loops at same time - be careful with memory and storage
    - `vary_alphaL*.py` contains the sensitivity analysis for alphaL, ""
    - `5. Vary_params_plotting_buoy_ml.ipynb` used to plot the results form the varying python script.
    - `6. adjust_injection_scheme_buoy_ml.ipynb` used to optimise injection scheme - run it to actually produce the desired volume. 
    - `7. Calculate_max_injection_extraction.ipynb` calculates max injection and extraction using emperical formulae. 
- `Sources` contains pdfs of papers used (please cite correctly).
- `Documents` contain initial workplan, final presentation and final report.

## running code
Advised to work in an anaconda environment with jupyter lab to run notebooks. 
Use `pip install flopy` to obtain the [FloPy](https://github.com/avaframe/FlowPy) package. 
The executable of modflow 6 is provided, strictly this should only be distribued by the USGS itself so get it [there](https://water.usgs.gov/water-resources/software/MODFLOW-6/). _Only in this repo to make sharing easier_. 



## Background
A drinking water company wants to build an ASR system to store drinking water for use during
the summer peak demand. The aquifer that they want to use consists of sand, is approximately
20 m thick, and the hydraulic conductivity is somewhere in the range 10 m/d till 40 m/d. Water in
the aquifer is salty with a concentration of 30 g/L. The drinking water company wants to store
enough drinking water such that they can extract a total of 40,000 m 3 of drinking water during the
summer months of July and August.
## Expectations
- Design an ASR system such that the drinking water company can extract a total of 40,000 m 3 of
drinking water during the summer months of July and August.
- Develop a schedule for injection, extraction and, if you wish, storage.
- Design the system such that the waste of injected drinking water is as small as possible.
- Make sure you meet the guidelines for injection pressure and maximum velocity during
extraction.
- Demonstrate that the ASR system can function for at least 10 years with possibly a start-up year
when the extracted volume is smaller.
- Discuss the uncertainties in your design.
- Report and present your findings.

