# EEeWorm submodels (version 7.3) # STABLE

###########################################################                                                          
#                 PATCHES                                 #                                                          
#                                                         #
###########################################################

################## ENV FLUCTUATION #################################
# Function to set up environmental fluctuations
setup_env_fluctuations <- function(world,day,env_data=Rotham) {
  # hour: current hour of simulation (integer)
  # temperature, SWP, energy_content_food: existing worldMatrix objects (patch variables)
  # env_data: data.frame with columns temperature, swp, energy_content_food (22 rows)
  
  # Extract y-coords (pycor) for each patch
  pycor_vals <- temperature@pCoords[,2]
  
  if(hour == 1) {
    # Assign values from env_data based on pycor ranges (mimicking your NetLogo ifelse chains)
    
    # temperature update
    temperature_ref <<- NLset(temperature,
                         agents = patches(temperature),
                         val = ifelse(pycor_vals >= 89, env_data[day,1],
                                      ifelse(pycor_vals >= 79 & pycor_vals < 89, env_data[day,2],
                                             ifelse(pycor_vals >= 69 & pycor_vals < 79, env_data[day,3],
                                                    ifelse(pycor_vals >= 59 & pycor_vals < 69, env_data[day,4],
                                                           ifelse(pycor_vals >= 49 & pycor_vals < 59, env_data[day,5],
                                                                  ifelse(pycor_vals >= 39 & pycor_vals < 49, env_data[day,6],
                                                                         ifelse(pycor_vals < 39, env_data[day,7], NA))))))))
    
    # SWP update
    SWP_ref <<- NLset(SWP,
                 agents = patches(SWP),
                 val = ifelse(pycor_vals >= 89, env_data[day,8],
                              ifelse(pycor_vals >= 79 & pycor_vals < 89, env_data[day,9],
                                     ifelse(pycor_vals >= 69 & pycor_vals < 79, env_data[day,10],
                                            ifelse(pycor_vals >= 59 & pycor_vals < 69, env_data[day,11],
                                                   ifelse(pycor_vals >= 49 & pycor_vals < 59, env_data[day,12],
                                                          ifelse(pycor_vals >= 39 & pycor_vals < 49, env_data[day,13],
                                                                 ifelse(pycor_vals < 39, env_data[day,14], NA))))))))
    
    # energy_content_food update
    energy_content_food_ref <<- NLset(energy_content_food,
                                 agents = patches(energy_content_food),
                                 val = ifelse(pycor_vals >= 99, env_data[day,15],
                                              ifelse(pycor_vals >= 90 & pycor_vals < 99, env_data[day,16],
                                                     ifelse(pycor_vals >= 80 & pycor_vals < 90, env_data[day,17],
                                                            ifelse(pycor_vals >= 70 & pycor_vals < 80, env_data[day,18],
                                                                   ifelse(pycor_vals >= 60 & pycor_vals < 70, env_data[day,19],
                                                                          ifelse(pycor_vals >= 40 & pycor_vals < 60, env_data[day,20],
                                                                                 ifelse(pycor_vals >= 20 & pycor_vals < 40, env_data[day,21],
                                                                                        ifelse(pycor_vals < 20, env_data[day,22], NA)))))))))
  
    SWP <<- SWP_ref
    temperature <<- temperature_ref
    energy_content_food <<- energy_content_food_ref
  }else{
  SWP <<- SWP_ref + rnorm(length(SWP_ref), mean = 0, sd = 0.02)
  temperature <<- temperature_ref + rnorm(length(temperature_ref), mean = 0, sd = 0.5)
  energy_content_food <<- energy_content_food_ref + rnorm(length(energy_content_food_ref), mean = 0, sd = energy_content_food_ref/10)
  }
     # Return updated worldMatrices as a list
  return(list(temperature = temperature, SWP = SWP, energy_content_food = energy_content_food))
}

############# Setup patch quality ###################

# Function to set patch quality
setup_patch_quality <- function() {
  
  coords <- patches(temperature)
  pycors <- coords[, 2]
  
  # 1. patch_T
  temp_vals <- of(temperature, coords)
  patch_T_vals <- ifelse(temp_vals < 15,
                         0.065 * temp_vals,
                         1.75 - (0.05 * temp_vals))
  patch_T <<- NLset(world = patch_T, agents =  coords, val = patch_T_vals)
  
  # 2. patch_SWP
  swp_vals <- of(SWP, coords)
  patch_SWP_vals <- ifelse(swp_vals < 5,
                           0.80 * swp_vals,
                           5 - (0.20 * swp_vals))
  patch_SWP <<- NLset(world = patch_SWP,agents =  coords,val =  patch_SWP_vals)
  
  # 3. patch_Ex
  ex_vals <- of(energy_content_food, coords)
  patch_Ex_vals <- ifelse(ex_vals < 2.5,
                          0.80 * ex_vals,
                          1.8 + (0.12 * ex_vals))
  patch_Ex <<- NLset(world = patch_Ex,agents =  coords,val =  patch_Ex_vals)
  
  # 4. patch_quality = T + SWP + Ex
  quality_vals <- patch_T_vals + patch_SWP_vals + patch_Ex_vals
  
  # 5. Penalize if it's light and patch is in surface layer
  # 5. Penalize if it's light and patch is in surface layer
  if (!dark) {
    quality_vals[pycors > 99] <- quality_vals[pycors > 99] - 5
  }
  
  patch_quality <<- NLset(world = patch_quality, agents =  coords,val =  quality_vals)
  
  
  
  # 6. pcolor correction: soil, surface, burrow
  pcolor_vals <- pcolor@.Data
  
  for(idx in 1:nrow(coords)) {
    x <- pcolor@pCoords[idx,1]
    y <- pcolor@pCoords[idx,2]
    patch_type <- habitat_type@.Data[y, x]
    
    pcolor_vals[y, x] <- switch(
      patch_type,
      "soil" = 52 + quality_vals[idx],
      "soil-surface" = 33,     # couleur fixe pour surface
      "burrow" = 30,           # couleur fixe pour burrow
      52                        # fallback
    )
  }
  
  pcolor@.Data<<-pcolor_vals

  #plot(patch_quality)
  # Optional: Update pcolor logic (for visualization)
#  habitat_vals <- habitat_type
#  pcolor_vals <- ifelse(habitat_vals == "soil",
#                        52 + quality_vals,
#                        33)  # Only for soil
#  pcolor <<- NLset(world = pcolor, agents = coords,val =  pcolor_vals)
#  
  # Retpatch_SWP# Return updated world matrices
  return(list(patch_T = patch_T,
       patch_SWP = patch_SWP,
       patch_Ex = patch_Ex,
       patch_quality = patch_quality,
       pcolor = pcolor))
}

###########################################################                                                          
#                 TURTLES                                 #                                                          
#                                                         #
###########################################################

################## BACK GROUND MORTALITY ############################
calc_background_mortality <- function(worms) {
  
  # Extract worm properties
  ycor_vals  <- of(agents = worms, var = "ycor")
  breed_vals <- of(agents = worms, var = "breed")
  wormIDs    <- of(agents = worms, var = "who")
  
  dead_ids <- c()
  
  for (i in seq_along(wormIDs)) {
    who   <- wormIDs[i]
    ycor  <- ycor_vals[i]
    breed <- breed_vals[i]
    
    if (!dark && ycor > 95) {
      # High mortality at soil surface during daytime
      if (runif(1,min = 0,max = 100) < ((NLcount(worms) * 0.50)/100)) {
        dead_ids <- c(dead_ids, who)
      }
    } else {
      # Background mortality
      if (breed == "juveniles") {
        if (runif(n = 1,min = 0.0,max = 100) < ((NLcount(worms) * background_mortality) /100)) {
          dead_ids <- c(dead_ids, who)
        }
      } else if (breed == "adults") {
        if (runif(n = 1,min = 0.0,max=100) < ((NLcount(worms) * background_mortality * 5)/100 )) {
          dead_ids <- c(dead_ids, who)
        }
      }
    }
  }
  
  if (length(dead_ids) > 0) {
    worms <<- die(turtles = worms, who = dead_ids)
  }
  
  return(worms)
}

###########################################################                                                          
#                 LIFE STAGE TRANSITIONS                  #                                                          
#                                                         #
###########################################################

transform_cocoons <- function(focal) {
  # Select only cocoons
  cocoons <- NLwith(agents = focal, var = "breed", val = "cocoons")
  
  # Exit early if no cocoons
  if (NLcount(cocoons) == 0) return(cocoons)
  
  # Retrieve IDs and ages
  whococoons <- of(agents = cocoons, var = "who")
  ages <- of(agents = cocoons, var = "age")
  
  # Get temperature at cocoons' locations
  coords <- patchHere(world = temperature, turtles = cocoons)
  Temp <- of(world = temperature, agents = coords)
  Temp <- pmax(Temp+273.15, 0.1)  # Prevent division by zero
  
  # Calculate temperature-dependent incubation time
  incubation_time <- 90 * exp((-activation_energy / Boltz) * ((1 / reference_T) - (1 / Temp)))
  
  # Identify cocoons ready to transform
  to_transform <- whococoons[ages >= incubation_time]
  
  if (length(to_transform) > 0) {
    # Retrieve those cocoons as agents
    cocoonsToTransform <- turtle(turtles = focal, who = to_transform)
    
    # Set new properties: breed, color, size
    focal <- NLset(turtles = focal, agents = cocoonsToTransform, var = "breed", val = rep("juveniles", length(to_transform)))
    focal <- NLset(turtles = focal, agents = cocoonsToTransform, var = "color", val = rep("pink", length(to_transform)))
    focal <- NLset(turtles = focal, agents = cocoonsToTransform, var = "size", val = rep(2, length(to_transform)))  
  }
  
  return(focal)
}

transform_juveniles <- function(focal) {
  # Select juveniles
  juveniles <- NLwith(agents = focal, var = "breed", val = "juveniles")
  
  if (NLcount(juveniles) == 0) return(focal)
  
  # Retrieve who and mass
  whojuveniles <- of(agents = juveniles, var = "who")
  mass <- of(agents = juveniles, var = "mass")
  
  # Find juveniles to transform
  to_transform <- whojuveniles[mass >= mass_sexual_maturity]
  
  if (length(to_transform) > 0) {
    juvenilesToTransform <- turtle(turtles = focal, who = to_transform)
    
    focal <- NLset(
      turtles = focal,
      agents = juvenilesToTransform,
      var = c("breed", "color", "size"),
      val = cbind(
        rep("adults", length(to_transform)),
        rep("red", length(to_transform)),
        of(agents = juvenilesToTransform, var = "mass")  # set size = mass
      )
    )
    
    # === TESTING POST-TRANSFORMATION ===
    expect_true(all(of(agents = turtle(focal, who = to_transform),var =  "breed") == "adults"))
    expect_true(all(of(agents= turtle(focal, who = to_transform), var = "color") == "red"))
    expect_equal(of(agents = turtle(focal, who = to_transform), var = "size"),
                 of(agents = turtle(focal, who = to_transform), var = "mass"))
    
  }
  
  return(focal)
}

###########################################################                                                          
#                 BURROwING AND FORAGING BEHAVIOUR        #                                                          
#                                                         #
###########################################################

# check_burrow: sélectionne les groupes et appelle les sous-routines ciblées
check_burrow <- function(worms,focal) {
  if (is.null(focal)) {
    worms <- focal 
    }
  with_burrow <- NLwith(agents = focal, var = "burrowQ", val = TRUE)
  without_burrow <- NLwith(agents = focal, var = "burrowQ", val = FALSE)
  
  if (NLcount(with_burrow) > 0) {
    if (dark) {
      focal <- forage(worms = worms, focal = with_burrow)
    } else {
      focal <- move_in_burrow(worms = worms, focal = with_burrow)
    }
  }
  
  if (NLcount(without_burrow) > 0) {
    focal <- find_burrow_site(worms, focal = without_burrow)
  }
  
  return(focal)
}

# forage: applique la logique de forage sur le groupe 'focal' (ou sur tous si focal = NULL)
forage <- function(worms, focal = NULL) {
  if (is.null(focal)) focal <- worms
  who_vec <- of(agents = focal, var = "who")
  
  for (who in who_vec) {
    worm_i <- turtle(turtles = focal, who = who)
    
    # 1) move-to one-of patches with [pxcor = pxcor and habitat-type = "soil-surface"]
    current_patch <- patchHere(world = habitat_type, turtles = worm_i)
    if (!is.null(current_patch) && nrow(current_patch) > 0) {
      my_pxcor <- current_patch[, 1]  # extract pxcor
      # select all soil-surface patches at this pxcor
      soil_surface <- NLwith(world = habitat_type,
                             agents = patch(habitat_type, 
                                            x = rep(my_pxcor, nrow(habitat_type@pCoords)), 
                                            y = habitat_type@pCoords[,2]),
                             val = "soil-surface")
      if (NLcount(soil_surface) > 0) {
        target_patch <- oneOf(soil_surface)
        worm_i <- moveTo(turtles = worm_i, agents = target_patch)
        focal[focal@.Data[, "who"] == who, ] <- worm_i
        
      }
    }
    
    # 2) check nearby burrows in radius 50
    nearby <- inRadius(agents = worm_i,
                       agents2 = patches(habitat_type),
                       radius = 50,
                       world = habitat_type,
                       torus = FALSE)
    nearby_burrows <- NLwith(world = habitat_type, agents = nearby, val = "burrow")
    
    heading_val <- if (NLcount(nearby_burrows) > 0) 270 else 90
    worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "heading", val = heading_val)
    focal[focal@.Data[, "who"] == who, ] <- worm_i

    # 3) forward by random up to max_crawling_speed
    max_speed <- of(agents = worm_i, var = "max_crawling_speed")
    if (is.na(max_speed)) max_speed <- 0
    dist <- runif(1, 0, max_speed)
    worm_i <- fd(worm_i, dist = dist)
    focal[focal@.Data[, "who"] == who, ] <- worm_i
  }
    # 4) ingestion (updates worms globally)
    focal <- calc_ingestion_rate(worms = focal, focal = focal)

  
  return(focal)
}

# move_in_burrow: pour le groupe focal
move_in_burrow <- function(worms, focal = NULL) {
  if (is.null(focal)) focal <- worms
  who_vec <- of(agents = focal, var = "who")
  
  for (who in who_vec) {
    worm_i <- turtle(turtles = focal, who = who)
    
    # patch-here (gardes si absent)
    current_patch <- patchHere(world = habitat_type, turtles = worm_i)
    if (is.null(current_patch) || nrow(current_patch) == 0) next
    my_pxcor <- current_patch[, 1]
    my_pcolor <- of(world = pcolor, agents = current_patch)
    
    # candidate patches in same pxcor from patch_quality
    pcoords <- patch_quality@pCoords
    rows <- which(pcoords[, 1] == my_pxcor)
    if (length(rows) == 0) next
    ys <- pcoords[rows, 2]
    agents_candidates <- patch(patch_quality, x = rep(my_pxcor, length(ys)), y = ys, duplicate = TRUE)
    values <- of(world = patch_quality, agents = agents_candidates)
    
    # find favourite patch (max patch_quality)
    best_idx <- which.max(values)
    favourite_patch <- agents_candidates[best_idx, , drop = FALSE]
    
    # distance and burrowing speed
    dist <- NLdist(agents = current_patch, agents2 = favourite_patch, world = patch_quality)
    if (length(dist) > 1) dist <- dist[1]
    speed <- of(agents = worm_i, var = "burrowing_speed")
    if (is.na(speed)) speed <- 0
    
    if (dist > floor(speed)) {
      heading <- towards(agents = current_patch, agents2 = favourite_patch, world = patch_quality)
      worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "heading", val = heading)
      # move forward by burrowing_speed (int)
      worm_i <- fd(worm_i, dist = floor(speed))
    } else {
      worm_i <- moveTo(turtles = worm_i, agents = favourite_patch)
    }
    
    focal[focal@.Data[, "who"] == who, ] <- worm_i
    
    # If current patch is 'soil', convert to burrow and apply costs
    new_patch <- patchHere(world = habitat_type, turtles = worm_i)
    if (!is.null(new_patch) && nrow(new_patch) > 0) {
      patch_type <- of(world = habitat_type, agents = new_patch)
      if (!is.na(patch_type) && patch_type == "soil") {
        # update global patch layers
        habitat_type <<- NLset(world = habitat_type, agents = new_patch, val = "burrow")
        pcolor <<- NLset(world = pcolor, agents = new_patch, val = my_pcolor)
        # apply burrowing costs — fonction attendue acceptant focal (et optionnellement focal)
        worm_i <- calc_burrowing_costs(worms = focal, focal = worm_i ) 
        focal[focal@.Data[, "who"] == who, ] <- worm_i
      }
    }
  }
  
  return(focal)
}

# find_burrow_site: pour le groupe focal
find_burrow_site <- function(worms, focal = NULL) {
  if (is.null(focal)) focal <- worms
  who_vec <- of(agents = focal, var = "who")
  
  for (who in who_vec) {
    worm_i <- turtle(turtles = focal, who = who)
    current_patch <- patchHere(world = habitat_type, turtles = worm_i)
    if (is.null(current_patch) || nrow(current_patch) == 0) next
    
    # all burrow patches
    burrow_patches <- NLwith(world = habitat_type, agents = patches(habitat_type), val = "burrow")
    
    nearby_burrows <- inRadius(agents = current_patch, radius = 20,
                               agents2 = burrow_patches, world = habitat_type, torus = FALSE)
    
    if (nrow(nearby_burrows) > 0) {
      target_patch <- oneOf(nearby_burrows)
      worm_i <- moveTo(turtles = worm_i, agents = patch(habitat_type, x = target_patch[1], y = target_patch[2]))
      focal[focal@.Data[, "who"] == who, ] <- worm_i
      
      # set burrowQ, my_burrow
      focal <- NLset(turtles = focal, agents = worm_i, var = "burrowQ", val = TRUE)
      burrow_no_val <- of(world = burrow_no, agents = patch(habitat_type, x = target_patch[1], y = target_patch[2]))
      focal <- NLset(turtles = focal, agents = worm_i, var = "my_burrow", val = burrow_no_val)
      
      # then forage or move_in_burrow for *this* worm
      if (dark) {
        worm_i <- forage(worms = worms, focal = worm_i)
        focal[focal@.Data[, "who"] == who, ] <- worm_i
        } else {
        worm_i <- move_in_burrow(worms = worms, focal = worm_i)
        focal[focal@.Data[, "who"] == who, ] <- worm_i
        }
      
    } else {
      # if patch-here is soil -> maybe form new burrow
      patch_type <- of(world = habitat_type, agents = current_patch)
      if (!is.na(patch_type) && patch_type == "soil") {
        nearby2 <- inRadius(agents = current_patch, radius = 20,
                            agents2 = burrow_patches, world = habitat_type, torus = FALSE)
        if (nrow(nearby2) == 0) {
          worm_i <- form_new_burrow(worms = worms, focal = worm_i)
          focal[focal@.Data[, "who"] == who, ] <- worm_i
         
          }
      }
    }
  }
  
  return(focal)
}

# form_new_burrow: transforme la/les patches-here en nouveau burrow pour le groupe focal
form_new_burrow <- function(worms, focal = NULL) {
  if (is.null(focal)) focal <- worms
  who_vec <- of(agents = focal, var = "who")
  
  for (who in who_vec) {
    worm_i <- turtle(turtles = focal, who = who)
    here_patch <- patchHere(world = habitat_type, turtles = worm_i)
    if (is.null(here_patch) || nrow(here_patch) == 0) next
    
    burrow_patches <- NLwith(world = habitat_type, agents = patches(habitat_type), val = "burrow")
    if (NLcount(burrow_patches) == 0) {
      max_pcolor <- 0
      max_burrowno <- 0
    } else {
      max_pcolor <- max(of(world = pcolor, agents = burrow_patches), na.rm = TRUE)
      max_burrowno <- max(of(world = burrow_no, agents = burrow_patches), na.rm = TRUE)
    }
    
    new_pcolor <- max_pcolor + 1
    new_burrowno <- max_burrowno + 1
    
    # update global patch layers (globals)
    pcolor <<- NLset(world = pcolor, agents = here_patch, val = new_pcolor)
    burrow_no <<- NLset(world = burrow_no, agents = here_patch, val = new_burrowno)
    habitat_type <<- NLset(world = habitat_type, agents = here_patch, val = "burrow")
    
    # set my_burrow for the worm
    focal <- NLset(turtles = focal, agents = worm_i, var = "my_burrow", val = new_burrowno)
    focal <- NLset(turtles = focal, agents = worm_i, var = "burrowQ", val = TRUE)
    # apply burrowing costs (fonction attendue)
    worm_i <- calc_burrowing_costs(worms = focal, focal = worm_i)
    focal[focal@.Data[, "who"] == who, ] <- worm_i
  }
  
  return(focal)
}

# check_parental_burrow : pour chaque ver -> test mass > 1.5 puis logique parentale
check_parental_burrow <- function(focal) {
  who_vec <- of(agents = focal, var = "who")
  
  for (who in who_vec) {
    worm_i <- turtle(turtles = focal, who = who)
    mass <- of(agents = worm_i, var = "mass")
    
    if (!is.na(mass) && mass > 1.5) {
      my_burrow_val <- of(agents = worm_i, var = "my_burrow")
      focal <- NLset(turtles = focal, agents = worm_i, var = "last_burrow", val = my_burrow_val)
      
      same_burrow <- NLwith(agents = focal, var = "my_burrow", val = my_burrow_val)
      adults_in_same <- NLwith(agents = same_burrow, var = "breed", val = "adults")
      
      if (NLcount(adults_in_same) > 0) {
        # si adultes présents → exécuter check_burrow mais seulement pour le focal (on appelle les routines focales)
        if (dark) {
          worm_i <- forage(worms, focal = worm_i)
          focal[focal@.Data[, "who"] == who, ] <- worm_i
        } else {
          worm_i <- move_in_burrow(worms, focal = worm_i)
          focal[focal@.Data[, "who"] == who, ] <- worm_i
        }
      } else {
        worm_i <- find_burrow_site(worms, focal = worm_i)
        focal[focal@.Data[, "who"] == who, ] <- worm_i
      }
    }
  }
  
  return(focal)
}

###########################################################                                                          
#                 INDIVIDUAL ENERGY BUDGET                #                                                          
#                                                         #
###########################################################

calc_maintenance <- function(focal) {
  if (is.null(focal)) focal <- worms
  # Clamp coords inside world boundaries
  focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor), max_Pxcor)
  focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor), max_Pycor)
  
   who_vec <- of(agents = focal, var = "who")
  
  for (who in who_vec) {
    worm_i <- turtle(turtles = focal, who =  who)
    mass_i <- of(agents = worm_i, var = "mass")
    energy_assim_i <- of(agents = worm_i, var = "energy_assimilated")
    energy_reserve_i <- of(agents = worm_i, var = "energy_reserve")
    energy_reserve_max_i <- of(agents = worm_i, var = "energy_reserve_max")
    breed_i <- of(agents = worm_i, var = "breed")
    
    # --- Get local temperature ---
    patch_coords_T <- patchHere(world = temperature, turtles = worm_i)  
    T_here <- of(world = temperature, agents = patch_coords_T)
    T_i_kelvins    <- T_here + 273.15
    
    # --- Calculate BMR ---
    BMR_i <- B_0 * (mass_i^(3/4)) * exp(-activation_energy / (Boltz * T_i_kelvins))
    worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "BMR", val = BMR_i)
    
    # --- Energy consumption ---
    if (energy_assim_i >= BMR_i) {
      energy_assim_i <- energy_assim_i - BMR_i
    } else {
      # Not enough assimilated energy, draw from reserve
      energy_reserve_i <- energy_reserve_i + energy_assim_i - BMR_i
      energy_assim_i <- 0
    }
    
    # --- Update worm ---
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_assimilated", val = energy_assim_i)
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_reserve", val = energy_reserve_i)
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_reserve_max", val = energy_reserve_max_i)
    focal <- NLset(turtles = focal, agents = worm_i, var = "BMR", val = BMR_i)    
    # --- Starvation check ---
    if (energy_reserve_i < (energy_reserve_max_i * 0.5) && breed_i != "cocoons") {
      worm_i <- onset_starvation_strategy(focal = worm_i)
     if (is.na(worm_i["who"])){
      focal <- die(turtles = focal, who = who)
     } else {
      focal[focal@.Data[, "who"] == who, ] <- worm_i
     }
    }
    
    # Write back to global worms
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_assimilated", val = of(agents = worm_i, var = "energy_assimilated"))
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_reserve", val = of(agents = worm_i, var = "energy_reserve"))
    focal <- NLset(turtles = focal, agents = worm_i, var = "mass", val = of(agents = worm_i, var = "mass"))
    focal <- NLset(turtles = focal, agents = worm_i, var = "breed", val = of(agents = worm_i, var = "breed"))
   # points(focal, col="orange")
    }
  
 return(focal)
 
     }

onset_starvation_strategy <- function(focal) {
  if (is.null(focal)) focal <- worms
 
  # Clamp coords inside world boundaries
  focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor), max_Pxcor)
  focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor), max_Pycor)
  
   # Extract variables
  BMR <- of(agents = focal, var = "BMR")
  #energy_flesh <- of(agents = focal, var = "energy_flesh")
  #energy_synthesis <- of(agents = focal, var = "energy_synthesis")
  mass <- of(agents = focal, var = "mass")
  energy_reserve <- of(agents = focal, var = "energy_reserve")
  #mass_sexual_maturity <- of(agents = focal, var = "mass_sexual_maturity")
  #mass_birth <- of(agents = focal, var = "mass_birth")
  
  # Compute new values
  new_mass <- mass - (BMR / (energy_flesh + energy_synthesis))
  new_reserve <- energy_reserve + BMR
  
  # Update mass and reserve
  focal <- NLset(turtles = focal, agents = focal, var = "mass", val = new_mass)
  focal <- NLset(turtles = focal, agents = focal, var = "energy_reserve", val = new_reserve)
  
  # Identify individuals with mass < mass_sexual_maturity
  juveniles_inds <- which(new_mass < mass_sexual_maturity)
  if (length(juveniles_inds) > 0) {
    juveniles_who <- of(agents = focal, var = "who")[juveniles_inds]
    juveniles <- turtle(focal, who = juveniles_who)
    
    juveniles <- NLset(turtles = juveniles, agents = juveniles, var = "breed", val = rep("juveniles", length(juveniles_inds)))
    juveniles <- NLset(turtles = juveniles, agents = juveniles, var = "color", val = rep("pink", length(juveniles_inds)))
    juveniles <- NLset(turtles = juveniles, agents = juveniles, var = "size", val = rep(2, length(juveniles_inds)))
    focal[focal[,"breed"] == 'juveniles', ] <- juveniles
   
    }
  
  # Identify individuals with mass < mass_birth (vectorized)
  who_vals <- of(agents = focal, var = "who")
  mass_vals <- of(agents = focal, var = "mass")
  to_die <- who_vals[mass_vals < mass_birth]
  
  if (length(to_die) > 0) {
    focal <- die(turtles = focal, who = to_die)
    worms <<- die(turtles = worms, who = to_die)
    cat("Worms died of starvation: ", paste(to_die, collapse = ", "), "\n")
  }
  
  return(focal)
}

calc_burrowing_costs <- function(worms, focal = NULL) {
  if (is.null(focal)) focal <- worms
  for (i in seq_len(nrow(focal))) {
    # Clamp coords inside world boundaries
    focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor), max_Pxcor)
    focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor), max_Pycor)
    
      who <- focal@.Data[i, "who"]
    worm_i <- focal[focal@.Data[, "who"] == who, ]
    
    # Individual state
    burrowing_speed <- of(agents = worm_i, var = "burrowing_speed")
    mass            <- of(agents = worm_i, var = "mass")
    energy_assim    <- of(agents = worm_i, var = "energy_assimilated")
    energy_reserve  <- of(agents = worm_i, var = "energy_reserve")
    
    # --- Maintenance costs ---
    max_burrowing_costs <- burrowing_speed * (burrowing_costs * (mass ^ (2/5)))
    
    if (energy_assim > 0) {
      if (energy_assim > max_burrowing_costs) {
        energy_assim <- energy_assim - max_burrowing_costs
      } else {
        burrowing_speed <- (burrowing_speed / max_burrowing_costs) * energy_assim
        energy_assim <- 0
      }
    } else {
      burrowing_speed <- burrowing_speed / 4
      energy_reserve  <- energy_reserve - (max_burrowing_costs / 4)
    }
    
    # Update worm’s state
    focal <- NLset(turtles = focal, agents = worm_i, var = "burrowing_speed",     val = burrowing_speed)
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_assimilated", val = energy_assim)
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_reserve",     val = energy_reserve)
    
    # --- Only apply if patch-here == "burrow" ---
    patch_turt <- patchHere(world = habitat_type, turtles = worm_i)
    habitat_here <- of(world = habitat_type, agents = patch_turt)
    
    if (habitat_here == "burrow") {
      # Save current patch info
      last_pcolor   <- of(world = pcolor,    agents = patch_turt)
      last_burrowNo <- of(world = burrow_no, agents = patch_turt)
      
      # Determine heading
      neighs <- neighbors(world = patch_quality, agents = patch_turt, nNeighbors = 8)
      best   <- maxOneOf(world = patch_quality, agents = neighs)
     # points(worm_i, col="green")
      patch_ahead <- patchAhead(turtles = worm_i, dist = burrowing_speed, world = habitat_type)
      pycor_ahead <- patch_ahead[,2]
      if(is.na(pycor_ahead) ){pycor_ahead<-FALSE}
      if (pycor_ahead > 99) {
        new_heading <- 180
      } else {
        new_heading <- towards(agents = patch_turt, agents2 = best, world = patch_quality)
      }
      
      # Update heading and move
      dist <- NLdist(
        agents  = worm_i,
        agents2 = patch_ahead,
        world   = patch_quality
      )
     
      focal <- NLset(turtles = focal, agents = worm_i, var = "heading", val = new_heading)
      worm_i <- fd(turtles = turtle(turtles = focal, who = who), dist = dist)
      
      # Write updated turtle back into worms
      focal[focal@.Data[, "who"] == who, ] <- worm_i
      # Update patch after movement (only where this worm moved)
      new_patch <- patchHere(world = w, turtles = focal[focal@.Data[, "who"] == who, ])
      pcolor <<- NLset(world = pcolor, agents = new_patch,val = last_pcolor)
      burrow_no <<- NLset(world = burrow_no, agents = new_patch, val= last_burrowNo)
      habitat_type <<- NLset(world = habitat_type, agents = new_patch, val = "burrow")
    }
  }
 
  return(focal)
}

###########################################################                                                          
#                 INGESTION RATE                      #                                                          
#                                                         #
###########################################################
########################### Ingestion Rate 
# juveniles and adults calculate their ingestion rate (the amount of food ingested from the environment) which depends on the food density of the patch in which they are present and the mass dependent maximum ingestion rate.

# ------------ Function Definitions ------------

calc_ingestion_rate <- function(worms,focal = NULL) {
  if (is.null(focal)) focal <- worms
   who_vec <- of(agents = focal, var = "who")
  
  # Clamp coords inside world boundaries
  focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor), max_Pxcor)
  focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor), max_Pycor)
  
  for (who in who_vec) {
    worm_i <- turtle(focal, who)
    
    current_patch <- patchHere(world = food_density, turtles = worm_i)
    
    if (!is.null(current_patch) && nrow(current_patch) > 0) {
      px <- current_patch[, "pxcor"]
      py <- current_patch[, "pycor"]
      
      # Patch food density + soil water potential
      food_here <- of(world = food_density, agents = patch(world = food_density, x = px, y = py))
      swp_here  <- of(world = SWP, agents = patch(world = SWP, x = px, y = py))
      Arrhenius_here <- Arrhenius[current_patch[2], current_patch[1]]
      
      if (food_here > 0) {
        mass <- of(agents = worm_i, var = "mass")
        
        # Initial ingestion rate (before food limitation)
        ingestion_rate <- (max_ingestion_rate * Arrhenius_here) * (mass^(2/3)) * exp(-0.04 * swp_here)
        worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "ingestion_rate", val = ingestion_rate)
        focal<-NLset(turtles = focal, agents = worm_i, var = "ingestion_rate", val = ingestion_rate)
        # Update worm in global focal
                    
        # focal on this patch
        turtles_here <- focal[round(focal@.Data[, "xcor"],0) == px &
            round(focal@.Data[, "ycor"],0) == py,]
        total_IR <- sum(of(agents = turtles_here, var = "ingestion_rate"))
        
        if (total_IR > food_here) {
          n_here <- NLcount(turtles_here)
          share <- food_here / n_here
          
          for (who2 in of(agents = turtles_here, var = "who")) {
            focal <- NLset(turtles = focal, agents = turtle(focal, who2),
                           var = "ingestion_rate", val = share)
          }
          # deplete food
          food_density <- NLset(world = food_density, agents = patch(world = food_density, x = px, y = py), val = 0)
        } else {
          # subtract consumed food
          food_density <- NLset(world = food_density, agents = patch(world = food_density, x = px, y = py),
                                val = food_here - total_IR)
        }
      } else {
        worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "ingestion_rate", val = 0)
      }
      
      # Update worm in global focal
      focal <- NLset(turtles = focal, agents = worm_i, var = "ingestion_rate",
                     val = of(agents = worm_i, var = "ingestion_rate"))
    }
  }
  
  # Assimilation step
  focal <- calc_assimilation(focal = focal )
  
  return(focal)
}


calc_assimilation <- function(focal) {
  who_vec <- of(agents = focal, var = "who")
  
  # Clamp coords inside world boundaries
  focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor), max_Pxcor)
  focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor), max_Pycor)
  
  for (who in who_vec) {
    worm_i <- turtle(focal, who)
    
    # Patch where this worm is located
    current_patch <- patchHere(world = energy_content_food, turtles = worm_i)
  
        energy_content_here <- of(world = energy_content_food, agents = current_patch)
    
    ingestion_rate_i <- of(agents = worm_i, var = "ingestion_rate")
    
    # Assimilation
    energy_assimilated_i <- ingestion_rate_i * energy_content_here * assimilation_efficiency
 #   cat(" Ingestion for who =", who, "-> ingestion_rate:", ingestion_rate_i, "energy_content:", energy_content_here, "-> energy_assimilated:", energy_assimilated_i, "\n")
    # Update worm in global focal
    focal <- NLset(turtles = focal, agents = worm_i, var = "energy_assimilated", val = energy_assimilated_i)
    points(focal, col="blue")
  }
  
  return(focal)
}


###########################################################
#                     REPRODUCTION                        #
###########################################################
# ---------------- Mate Function ----------------
# reproduce: hatches one cocoon and returns updated worms + updated parent (as an agent subset)
reproduce <- function(reproducer_i, worms) {
  # get parent identifiers and values
  who_id     <- of(agents = reproducer_i, var = "who")
  hatchlings <- of(agents = reproducer_i, var = "hatchlings")
  R_old      <- of(agents = reproducer_i, var = "R")
  
  # save IDs before hatch
  ids_before <- of(agents = worms,var = "who")
  
  # hatch 1 (this returns updated worms with a new row)
  worms <- hatch(turtles = worms, who = who_id, n = 1, breed = "cocoons")
  
  # detect new IDs (robust even if parent fields were copied)
  ids_after <- of(agents = worms, var = "who")
  new_ids <- setdiff(ids_after, ids_before)
  
  if (length(new_ids) > 0) {
    new_cocoons <- NLwith(agents = worms, var = "who", val = new_ids)
    
    # set newborns exactly like NetLogo's hatch [ ... ] block
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "breed",              val = "cocoons")
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "color",              val = "white")
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "size",               val = 1)
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "energy_reserve",     val = mass_cocoon * energy_flesh)
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "energy_reserve_max", val = mass_cocoon * energy_flesh)
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "mass",               val = mass_birth)
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "energy_assimilated", val = 0)
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "age",                val = 0)
    worms <- NLset(turtles = worms, agents = new_cocoons, var = "hatchlings",         val = 0)
  }
  
  # update parent inside the worms object (so parent and worms remain consistent)
  parent_agent <- NLwith(agents = worms, var = "who", val = who_id)
  worms <- NLset(turtles = worms, agents = parent_agent, var = "hatchlings", val = hatchlings + 1)
  worms <- NLset(turtles = worms, agents = parent_agent, var = "R",          val = R_old  - (mass_cocoon * (energy_flesh + energy_synthesis)))
  
  # return both updated worms and the updated parent agent (subset)
  updated_parent <- NLwith(agents = worms, var = "who", val = who_id)
  cat("Hatching for parent who =", who_id, "-> new_ids:", paste(new_ids, collapse = ","), "\n")
  
  return(list(worms = worms, parent = updated_parent))
}

# calc_reproduction: handles energy allocation; when threshold reached calls reproduce()
calc_reproduction <- function(adult_i, worms) {
  mass_i             <- of(agents = adult_i, var = "mass")
  energy_assimilated <- of(agents = adult_i, var = "energy_assimilated")
  R                  <- of(agents = adult_i, var = "R")
  
  current_patch   <- patchHere(world = energy_content_food, turtles = adult_i)
  Arrhenius_here  <- Arrhenius[current_patch[2], current_patch[1]]
  max_R           <- (max_reproduction_rate * Arrhenius_here) * mass_i
  
  if (energy_assimilated >= max_R) {
    energy_assimilated <- energy_assimilated - max_R
    R <- R + max_R
  } else if (energy_assimilated > 0) {
    R <- R + energy_assimilated
    energy_assimilated <- 0
  }
  
  adult_i <- NLset(turtles = adult_i, agents = adult_i, var = "R", val = R)
  adult_i <- NLset(turtles = adult_i, agents = adult_i, var = "energy_assimilated", val = energy_assimilated)
  
 # cat(" R ", R >= mass_cocoon * (energy_flesh + energy_synthesis), " for adult who =", of(agents = adult_i, var = "who"), "\n")
  
  # reproduction step (call reproduce and accept both returns)
  if (R >= mass_cocoon * (energy_flesh + energy_synthesis)) {
    out <- reproduce(reproducer_i = adult_i, worms = worms)
    worms <- out$worms
    adult_i <- out$parent
    # get updated parent from the returned worms
   # who_id <- of(agents = adult_i, var = "who")
    #adult_i <- NLwith(agents = worms, var = "who", val = who_id)
    cat("Reproduce for adult who =", adult_i[,"who"], "\n")
  }
  
  # leftover energy -> growth (assumes calc_growth(adult_i) returns updated adult_i)
  if (energy_assimilated > 0) {
    adult_i <- calc_growth(focal = adult_i)
 # cat("Growth for adult who =", of(agents = adult_i, var = "who"), "\n")
  
    }
  
  return(list(adult_i = adult_i, worms = worms))
}

# mateQ: top-level loop over focal adults; receives and returns updated worms + focal
mateQ <- function(focal=worms, worms) {
  sperm_counter <- of(agents = focal, var = "sperm_counter")
  ready_to_mate <- which(sperm_counter > 2184 & dark == TRUE)
  
 # cat("Adults ready to mate (who):",ready_to_mate, "\n")
  
  for (i in seq_len(nrow(focal))) {
    who_i <- focal@.Data[i, "who"]
    adult_i <- focal[focal@.Data[, "who"] == who_i, , drop = FALSE]
    
   ## A VERIFIER  # mate test (same logic as NetLogo)
    who_ready <- of(agents = focal[ready_to_mate, , drop = FALSE], var = "who")
    if (who_i %in% who_ready) {
      nearby <- inRadius(
        agents  = adult_i,
        agents2 = focal[focal[, "breed"] == "adults", , drop = FALSE],
        radius  = 20,
        world   = habitat_type,
        torus   = FALSE
      )
      nearby_who <- setdiff(nearby[, "who"], who_i)
      if (length(nearby_who) > 0) {
        
        # reproduction bookkeeping: must accept both adult_i and worms returned
        out <- calc_reproduction(adult_i = adult_i, worms = worms)
        adult_i <- out$adult_i
        worms   <- out$worms
        # write updated adult back into focal (match by who)
        focal[focal@.Data[, "who"] == who_i, ] <- adult_i
        
        adult_i <- NLset(turtles = adult_i, agents = adult_i, var = "sperm_counter", val = 0)
        # also write this reset into worms (so parent in worms and focal stay consistent)
        worms <- NLset(turtles = worms, agents = NLwith(agents = worms, var = "who", val = who_i),
                       var = "sperm_counter", val = 0)
      }
    }else{
      # reproduction bookkeeping: must accept both adult_i and worms returned
      out <- calc_reproduction(adult_i = adult_i, worms = worms)
      adult_i <- out$adult_i
      worms   <- out$worms
      # write updated adult back into focal (match by who)
      focal[focal@.Data[, "who"] == who_i, ] <- adult_i
      
    }
  }
  
  return(list(focal = focal, worms = worms))
}


###########################################################
#                         GROWTH                          #
###########################################################

calc_growth <- function(focal) {
  # Worm-specific mass
  mass_i <- of(agents = focal, var = "mass")
  
  # Patch where the worm is located
  current_patch <- patchHere(world = SWP, turtles = focal)
  
  # Local Arrhenius from patch
  Arrhenius_here <- Arrhenius[current_patch[2], current_patch[1]]
  
  # Maximum growth rate based on von Bertalanffy
  max_growth_rate_focal <- growth_constant * Arrhenius_here * (mass_maximum^(1/3) * mass_i^(2/3) - mass_i)
  
  # Update focal's max_growth_rate
  focal <- NLset(turtles = focal, agents = focal, var = "max_growth_rate", val = max_growth_rate_focal)
  
  # If mass < mass_maximum, apply growth
  if (mass_i < mass_maximum) {
    focal <- grow(focal=focal)
  }
  
  return(focal)
}

grow <- function(focal) {
  # Worm-specific variables
  max_growth_rate_focal <- of(agents = focal, var = "max_growth_rate")
  energy_assimilated_focal <- of(agents = focal, var = "energy_assimilated")
  mass_focal <- of(agents = focal, var = "mass")
  
  # Energy required for max growth
  energy_growth_focal <- max_growth_rate_focal * (energy_flesh + energy_synthesis)
  
  # Full growth if enough energy
  if (energy_assimilated_focal >= energy_growth_focal) {
    growth_rate_focal <- max_growth_rate_focal
    mass_focal <- mass_focal + growth_rate_focal
    energy_assimilated_focal <- energy_assimilated_focal - energy_growth_focal
  } else {
    # Partial growth if insufficient energy
    growth_rate_focal <- (max_growth_rate_focal / energy_growth_focal) * energy_assimilated_focal
    mass_focal <- mass_focal + growth_rate_focal
    energy_assimilated_focal <- 0
  }
  
  # Update _focal in worms population
  focal <- NLset(turtles = focal, agents = focal, var = "growth_rate", val = growth_rate_focal)
  focal <- NLset(turtles = focal, agents = focal, var = "mass", val = mass_focal)
  focal <- NLset(turtles = focal, agents = focal, var = "energy_assimilated", val = energy_assimilated_focal)
  
  return(focal)
}

update_patches <- function(food_density) {
  food_density <- food_density
  coords_patches <- food_density@pCoords
  
  coords_turt <- cbind(worms@.Data[,"xcor"],worms@.Data[,"ycor"])
  ingestion_rate <- worms@.Data[,"ingestion_rate"]
  
  total_ingestion <- numeric(nrow(coords_patches))
  #to check if start from 0 : for(i in seq_len(nrow(coords_patches))-1)
  for(i in seq_len(nrow(coords_patches))) {
    idx <- which(coords_turt[,1] == coords_patches[i,1] & coords_turt[,2] == coords_patches[i,2])
    total_ingestion[i] <- sum(ingestion_rate[idx])
  }
  
  food_density_new <- ifelse(food_density <= 0, food_density, pmax(0, food_density - total_ingestion))
  
  food_density <- NLset(world  = food_density, agents = patches(food_density), val = food_density_new)
  return(food_density)
}

update_turtles <- function(worms) {
  # Only update if end of day
  if (hour == 24) {
    for (i in seq_len(NLcount(worms))) {
      worm_i <- turtle(turtles = worms, who = of(agents = worms, var = "who")[i])
      
      energy_assimilated <- of(agents = worm_i, var = "energy_assimilated")
      breed              <- of(agents = worm_i, var = "breed")
      mass               <- of(agents = worm_i, var = "mass")
      #mass_maximum       
      energy_reserve     <- of(agents = worm_i, var = "energy_reserve")
    #  energy_flesh       
     # energy_synthesis   
      
      # Energy update
      if (energy_assimilated > 0) {
        new_reserve <- energy_reserve +
          (energy_assimilated * (energy_flesh / (energy_flesh + energy_synthesis)))
        energy_assimilated <- 0
        
        # Adults and juveniles: update max reserve
        if (breed %in% c("adults", "juveniles")) {
          energy_reserve_max <- (mass / 2) * energy_flesh
          
          if (new_reserve > energy_reserve_max) {
            new_reserve <- energy_reserve_max
            mass <- mass + (energy_reserve_max / (energy_flesh + energy_synthesis))
            if (mass > mass_maximum) mass <- mass_maximum
          }
        }
        
        worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "energy_reserve", val = new_reserve)
        worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "energy_assimilated", val = energy_assimilated)
        worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "energy_reserve_max", val = energy_reserve_max)        
        if (breed %in% c("adults", "juveniles")) {
          worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "mass", val = mass)
        }
      }
      
      # Size update for visualization
      size_val <- ifelse(breed == "adults", 4,
                         ifelse(breed == "juveniles", 3,
                                ifelse(breed == "cocoons", 2, 1)))
      worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "size", val = size_val)
      
      # Adults: increment sperm counter
      if (breed == "adults") {
        sperm_counter <- of(agents = worm_i, var = "sperm_counter") + 1
        worm_i <- NLset(turtles = worm_i, agents = worm_i, var = "sperm_counter", val = sperm_counter)
      }
      
      # Write back updated worm_i
      worms[worms@.Data[, "who"] == of(agents = worm_i, var = "who"), ] <- worm_i
      
    }
  }
  
  # Ensure changes persist outside the function
  return(worms)
}


####### Update worms with  adults, juveniles and cocoons ######

upworms_with<-function (group, worms){
  worms <- NLset(turtles = worms, agents = group, var = "xcor", val = of(agents = group, var = "xcor"))
  worms <- NLset(turtles = worms, agents = group, var = "ycor", val = of(agents = group, var = "ycor"))
  worms <- NLset(turtles = worms, agents = group, var = "who", val = of(agents = group, var = "who"))
  worms <- NLset(turtles = worms, agents = group, var = "heading", val = of(agents = group, var = "heading"))
  worms <- NLset(turtles = worms, agents = group, var = "prevX", val = of(agents = group, var = "prevX"))
  worms <- NLset(turtles = worms, agents = group, var = "prevY", val = of(agents = group, var = "prevY"))
  worms <- NLset(turtles = worms, agents = group, var = "breed", val = of(agents = group, var = "breed"))
  worms <- NLset(turtles = worms, agents = group, var = "color", val = of(agents = group, var = "color"))
  worms <- NLset(turtles = worms, agents = group, var = "age", val = of(agents = group, var = "age"))
  worms <- NLset(turtles = worms, agents = group, var = "mass", val = of(agents = group, var = "mass"))
  worms <- NLset(turtles = worms, agents = group, var = "size", val = of(agents = group, var = "size"))
  worms <- NLset(turtles = worms, agents = group, var = "energy_reserve", val = of(agents = group, var = "energy_reserve"))
  worms <- NLset(turtles = worms, agents = group, var = "energy_assimilated", val = of(agents = group, var = "energy_assimilated"))
  worms <- NLset(turtles = worms, agents = group, var = "energy_reserve_max", val = of(agents = group, var = "energy_reserve_max"))
  worms <- NLset(turtles = worms, agents = group, var = "ingestion_rate", val = of(agents = group, var = "ingestion_rate"))
  worms <- NLset(turtles = worms, agents = group, var = "BMR", val = of(agents = group, var = "BMR"))
  worms <- NLset(turtles = worms, agents = group, var = "my_burrow", val = of(agents = group, var = "my_burrow"))
  worms <- NLset(turtles = worms, agents = group, var = "last_burrow", val = of(agents = group, var = "last_burrow"))
  worms <- NLset(turtles = worms, agents = group, var = "burrowing_speed", val = of(agents = group, var = "burrowing_speed"))
  worms <- NLset(turtles = worms, agents = group, var = "max_burrowing_costs", val = of(agents = group, var = "max_burrowing_costs"))
  worms <- NLset(turtles = worms, agents = group, var = "max_growth_rate", val = of(agents = group, var = "max_growth_rate"))
  worms <- NLset(turtles = worms, agents = group, var = "growth_rate", val = of(agents = group, var = "growth_rate"))
  worms <- NLset(turtles = worms, agents = group, var = "energy_growth", val = of(agents = group, var = "energy_growth"))
  worms <- NLset(turtles = worms, agents = group, var = "embryonic_development", val = of(agents = group, var = "embryonic_development"))
  worms <- NLset(turtles = worms, agents = group, var = "hatchlings", val = of(agents = group, var = "hatchlings"))
  worms <- NLset(turtles = worms, agents = group, var = "max_R", val = of(agents = group, var = "max_R"))
  worms <- NLset(turtles = worms, agents = group, var = "R", val = of(agents = group, var = "R"))
  worms <- NLset(turtles = worms, agents = group, var = "sperm_counter", val = of(agents = group, var = "sperm_counter"))
  worms <- NLset(turtles = worms, agents = group, var = "favourite_patch", val = of(agents = group, var = "favourite_patch"))
  worms <- NLset(turtles = worms, agents = group, var = "max_crawling_speed", val = of(agents = group, var = "max_crawling_speed"))
  worms <- NLset(turtles = worms, agents = group, var = "burrowQ", val = of(agents = group, var = "burrowQ"))
 return(worms) 
}

####### PLOT WORLD WITH BURROW AND WORMS  ######
library(ggplot2)
library(ggimage)

plot_world <- function(food_density, burrow_no, worms, hour, day, year, dark,
                       worm_icon = "C:/Users/qdevalloir/Pictures/worm.png") {
  
  # --- World dimensions ---
  nx <- food_density@maxPxcor - food_density@minPxcor + 1
  ny <- food_density@maxPycor - food_density@minPycor + 1
  pixel_size_x <- 1 / nx
  
  # --- Food density (background) ---
  mat_food <- food_density@.Data
  mat_food <- mat_food[nrow(mat_food):1, ]  # flip Y
  
  food_df <- expand.grid(
    x = seq(food_density@minPxcor, food_density@maxPxcor),
    y = seq(food_density@minPycor, food_density@maxPycor)
  )
  food_df$food <- as.vector(t(mat_food))
  
  # --- Burrows ---
  mat_burrow <- burrow_no@.Data
  mat_burrow <- mat_burrow[nrow(mat_burrow):1, ]
  
  burrow_df <- expand.grid(
    x = seq(burrow_no@minPxcor, burrow_no@maxPxcor),
    y = seq(burrow_no@minPycor, burrow_no@maxPycor)
  )
  burrow_df$burrow <- as.vector(t(mat_burrow))
  burrow_df <- burrow_df[burrow_df$burrow > 0, ]
  
  # --- Worms ---
  worms_df <- data.frame(
    x = worms@.Data[, "xcor"],
    y = worms@.Data[, "ycor"],
    size = of(agents = worms, var = "size"),
    color = of(agents = worms, var = "color"),
    icon = worm_icon
  )
  
  # --- Plot ---
  p <- ggplot() +
    # food resource background
    geom_tile(data = food_df, aes(x = x, y = y, fill = food)) +
    scale_fill_gradient(low = "lightyellow", high = "darkgreen",
                        name = "Food density") +
    
    # burrows
    geom_point(data = burrow_df, aes(x = x, y = y),
               shape = 15, color = "#654420",
               size = pixel_size_x * (nx/4), alpha = 0.8) +
    
    # worms as icons
    geom_image(
      data = worms_df,
      aes(x = x, y = y, image = icon),
      colour = worms_df$color,
      size = worms_df$size * (5 * pixel_size_x)   # scaling relative to cell size
    ) +
    
    coord_fixed() +
    theme_void() +
    ggtitle(paste0("hour = ", hour,
                   "; day = ", day,
                   "; year = ", year,
                   "; dark = ", dark))
    
  return(p)
}
