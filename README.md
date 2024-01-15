# design of an Aquifer Storage Recharge (ASR) system TUDelft CEGM2006
Designing an ASR system for the cross over course Subsurface storage for climate, energy and water course (CEGM2006)
# Explanation of repo

- `bin` contains required moflow6 executable.
- `Model` contains all the code.
    - `basic_model*.ipynb` contains the initial model created, the two versions contain slightly different versions - makes working together with git easier.
    - `changes_injection_scheme.ipynb` contains the code when injection scheme was altered to cater for the design requirements.
    - `vary_kh_test_injection.ipynb` is the start of the sensitivity analysis. As this was slower in notebook, this was also implemented in a `.py` files.
    - `vary_params.py` code to perform sensitivity analysis of hydraulic conductivity, porosity and coefficient of logitudinal dispersivity. 
    - `worstcase*.ipynb` contains a variation on the basic model but in the worst case as a result of the sensitivity analysis
    - `optimise_gw_model.py` used to optimise injection scheme the worst case, proved to be easier using seperate `.py` function. _Could_ optimise by adding this back in basic model. 
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

