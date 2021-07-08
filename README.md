# Habitat Usage
Repo for code used in Selwyn, Usseglio & Hogan 20XX analysis of Coryphopterus Habitat Usage

Bash code was written to call SLURM scripts in sequence which subsequently call R scripts for analysis of DEM and orthomosaic data in relation to the location of <i>Coryphopterus personatus</i>/<i>hyalinus</i>.

Scripts in the `utils` folder are for post run analysis and figure creation for the manuscript.

Data is found here: **(WILL ADD ZENODO REPO post-publication)**

The order of scripts run is:
1. `SLURM_scripts/classify_habitat.slurm`
    - Fits multiple classification algorithms to identify sand/reef from RGB values in orthomosaics.
    - Use `utils/assess_classification.R` to determine best classification model.
2. `classification_selection.sh`
3. `topography.sh`
    - Calculates a variety of topographical metrics from DEM and classified Orthomosaics.
    - Use `utils/site_stats.R` to summarize topographic and habitat classification stats for each site.
4. `multi_inla.sh`
    - Fits multiple [inlabru](https://sites.google.com/inlabru.org/inlabru) models using randomized locations of individual fish around shoal centroid.
        - occassionally model will fail in ways I haven't been able to diagnose/replicate. As such may need to be run a second/third time to get the total number of model fits desired
    - After fitting all models post-processing of models will be run to output joint parameter estimates and diagnostics across all models
5. `shoal_composition.R`
    - Analyze shoal composition of genetically sequenced shoals based on shoal size and topography.
    - Use `utils/species_breakdown.R` to calculate individual shoal composition independent of topography.
6. Use `utils/parameter_plot.R` to plot parameter estimates of both the distribution and composition models.
    - Additionally calculates a variety of metrics useful for publication (e.g., posterior probability of positive slope)
