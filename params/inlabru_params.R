#### Parameters file sourced into inlabru model ####

model_choice <- 'c5' #Choose habitat classification model to use (currently fixed at c5)
integration_strategy <- 'auto' #eb or auto (fast or slow)
number_samples <- 1000 #Number of times to sample posterior of posterior distributions 
pp_samples <- 100 #Number of posterior point processes to sample
number_workers <- 6 #number of parallel jobs at some steps - increase to go faster but risk running out of memory

#### Mesh Settings ####
mesh_setting <- list(
  max.edge = c(0.5, 5),
  min.angle = c(30, 21),
  max.n = c(48000, 16000), ## Safeguard against large meshes.
  max.n.strict = c(128000, 128000), ## Don't build a huge mesh!
  cutoff = 0.25, ## Filter away adjacent points.
  offset = c(2.5, 8), ## Offset for extra boundaries, if needed.
  crs = NA_character_
)

matern_priors <- list(
  prior.range = c(5, 0.5),
  prior.sigma = c(2, 0.05)
)

#### Various Formula to use for model fitting/prediction - make sure they all match ####
model_formula <- coordinates ~ 
  beta.boundaryDistance(main = continuous_covars(x, y, "boundaryDistance"), model = "linear", mean.linear = 0, prec.linear = 0.01) + 
  beta.complexity(main = continuous_covars(x, y, "complexity"), model = "linear", mean.linear = 0, prec.linear = 0.01) +
  beta.depth(main = continuous_covars(x, y, "depth"), model = "linear", mean.linear = 0, prec.linear = 0.01) + 
  beta.moran(main = continuous_covars(x, y, "moran"), model = "linear", mean.linear = 0, prec.linear = 0.01) + 
  beta.viewshed(main = continuous_covars(x, y, "viewshed"), model = "linear", mean.linear = 0, prec.linear = 0.01) + 
  spatialSmooth(main = coordinates, model = matern) + 
  
  habitat(main = habitat_type, model = "factor_contrast", hyper = list(prec = list(prior = 'loggamma', param = c(1, 5e-5)))) + 
  site(main = site_track, model = 'factor_contrast', hyper = list(prec = list(prior = 'gamma', param = c(2, 1)))) +
  
  Intercept(1, mean.linear = 0, prec.linear = 0.001)


pred_form <- ~ sum(weight * exp(Intercept + spatialSmooth + 
                                  beta.boundaryDistance +
                                  beta.complexity +
                                  beta.depth +
                                  beta.moran +
                                  beta.viewshed +
                                  habitat + 
                                  site))


spat_pred_form <- ~ exp(Intercept + spatialSmooth + 
                          beta.boundaryDistance + 
                          beta.complexity +
                          beta.depth +
                          beta.moran +
                          beta.viewshed +
                          habitat + 
                          site)

spat_pred_form_noSmooth <- ~ exp(Intercept + 
                                   beta.boundaryDistance + 
                                   beta.complexity +
                                   beta.depth +
                                   beta.moran +
                                   beta.viewshed +
                                   habitat + 
                                   site)
