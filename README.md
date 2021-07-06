# Habitat Usage
Repo for code used in Selwyn, Usseglio & Hogan 20XX analysis of Coryphopterus Habitat Usage

Bash code was written to call SLURM scripts in sequence which subsequently call R scripts for analysis of DEM and orthomosaic data in relation to the location of <i>Coryphopterus personatus</i>/<i>hyalinus</i>.

Scripts in the `utils` folder are for post run analysis and figure creation for the manuscript.

The order of scripts run is:
1. `SLURM_scripts/classify_habitat.slurm`
    - Fits multiple classification algorithms to identify sand/reef from RGB values in orthomosaics
    - Use `utils/assess_classification.R` to determine best classification model
2. `classification_selection.sh`
3. `topography.sh`
    - Calculates a variety of topographical metrics from DEM and classified Orthomosaics
    - Use `utils/site_stats.R` to summarize topographic and habitat classification stats for each site
4. `multi_inla.sh`
    - Fits multiple [inlabru](https://sites.google.com/inlabru.org/inlabru) models using randomized locations of individual fish around shoal centroid.
