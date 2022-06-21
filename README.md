# Teaching hydrological modelling: Illustrating model structure uncertainty with a ready-to-use teaching module
This repository contains a teaching module that can be used to showcase hydrologic model structure uncertainty (Knoben and Spieler, HESS, accepted).
The module contains data, example exercise sheets and example answer scripts.
Data are part of the CAMELS data set (Addor et al., 2017).
Models used are part of the MARRMoT framework (Knoben et al., 2019).
The module is split into two parts.
In the first, students obtain familiarity with the MARRMoT framework.
In the second, students are asked to calibrate two MARRMoT models for two CAMELS catchments and draw their own conclusions about model structure uncertainty.

A trial of this module at the Technische Universität Dresden (Germany) suggests that the module can completed in less than two afternoons, or a single full lecture day.

## Example of exercise sheets
Folders **exercise1_source** and **exercise2_source** contain `.pdf` and `.tex` files of handouts that can be given to students for the first and second part of the module. 
The first exercise sheet guides students through installing and using MARRMoT.
The second exercise sheet guides students through calibration and evaluation of two MARRMoT models for two CAMELS catchments, asking them to modify scripts they have been introduced to in the first part.
This document also contains questions to guide the four-way comparison of the performance of both models in both catchments, which leads to insights on model structure uncertainty.

## Catchment data
Files **Part 2 - catchment attributes.xlsx** and "Part 2 - catchment data.mat** contain data extracted from the CAMELS data set (Addor et al., 2017). 
The `.xlsx` file contains catchment descriptors, which can be used to gain a basic understanding of both catchments.
The `.mat` file contains time series of meteorological data and streamflow observations, which are needed for the calibration exercise.
How students shoud use both files is explained in the example handout for the second part of the exercise.

## Calibration outcomes
Files **calibrationScriptExample_MARRMoT_v1.m**,**calibrationScriptExample_MARRMoT_v2.m** and **calibrationResults.mat** contain the answers to the second part of the module.
The `.m` file shows how one of the scripts introduced in the first part of the exercise may be adapted to solve the second part of the exercise.
The `.mat` file contains calibrated parameter sets for the four combinations of models and catchments, and Kling-Gupta Efficiency (Gupta et al., 2009) scores obtained for calibration and evaluation time periods for all four cases.

### Note on MARRMoT versions
MARRMoT has been updated to version 2.x in 2022. One of the changes is inclusion of a dedicated calibration function contained within the main MARRMoT model class. 
To aid users of both the old (v1.x) and new (v2.x) versions, this repository contains example calibration code for both.
Calibration results are unchanged between both MARRMoT versions.
Overview:

- MARRMoT v1.x
	- Latest release: v1.4
	- Latest release DOI: https://dx.doi.org/10.5281/zenodo.6460624
	- Calibration code example: `calibrationScriptExample_MARRMoT_v1.m`
	- Paper: https://doi.org/10.5194/gmd-12-2463-2019
	- Paper supplements: https://gmd.copernicus.org/articles/12/2463/2019/gmd-12-2463-2019-supplement.pdf
- MARRMoT v2.x
	- Latest release: v2.1
	- Latest release DOI: https://doi.org/10.5281/zenodo.6484372
	- Calibration code example: `calibrationScriptExample_MARRMoT_v2.m`
	- Preprint: https://doi.org/10.5194/gmd-2022-135


## References
Addor, N., Newman, A. J., Mizukami, N., and Clark, M. P. (2017). _The CAMELS data set: catchment attributes and meteorology for large-sample studies_. Hydrology and Earth System Sciences, 21, 5293–5313, https://doi.org/10.5194/hess-2017-169/

Gupta,  H.  V.,  Kling,  H.,  Yilmaz,  K.  K.,  and  Martinez,  G.  F. (2009). _Decomposition  of  the  mean  squared  error  and  NSE  performance  criteria:Implications for improving hydrological modelling_. Journal of Hydrology, 377, 80–91, https://doi.org/10.1016/j.jhydrol.2009.08.003, http://dx.doi.org/10.1016/j.jhydrol.2009.08.003/

Knoben, W. J., Freer, J. E., Fowler, K. J., Peel, M. C., and Woods, R. A. (2019). _Modular Assessment of Rainfall-Runoff Models Toolbox (MARRMoT) v1.2: An open-source, extendable framework providing implementations of 46 conceptual hydrologic models as continuous state-space formulations_. Geoscientific Model Development, 12, 2463–2480, https://doi.org/10.5194/gmd-12-2463-2019/

Knoben, W. J. M. and Spieler, D. (accepted). _Teaching hydrological modelling: Illustrating model structure uncertainty with a ready-to-use teaching module_. HESS. Preprint DOI: https://doi.org/10.5194/hess-2021-30
