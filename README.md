# Habitat Usage
Repo for code used in Selwyn, Usseglio & Hogan 20XX analysis of Coryphopterus Habitat Usage

Bash code was written to call SLURM scripts in sequence which subsequently call R scripts for analysis of DEM and orthomosaic data in relation to the location of <i>Coryphopterus personatus</i>/<i>hyalinus</i>.

Scripts in the `utils` folder are for post run analysis and figure creation for the manuscript.

The order of scripts run is:
1. `SLURM_scripts/classify_habitat.slurm`
  -Use `utils/assess_classification.R` to determine best classification model
2. `classification_selection.sh`
3. `topography.sh`
4. `multi_inla.sh`
