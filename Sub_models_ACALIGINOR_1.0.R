# =============================
# File: submodels_eeeworm_1.2.R
# -----------------------------

# Submodels for EEEWORM. These are functions called by main go().

# 1) setup_env_fluctuations: replicates the file-read based fluctuations used in NetLogo
  setup_env_fluctuations <- function(temperature,soil_moisture, env_data = NULL) {
  # world: stacked world or individual layers (we use temperature and SWP)
  # env_data: table of environmental inputs (rows ~ timesteps), column order expected to match the original weather file
  n_patches <- nrow(temperature@.Data)
  # pick the row from env_data (cycle if necessary)
  if (!is.null(env_data)) {
    ridx <- ((day - 1) %% nrow(env_data)) + 1
    rowv <- as.numeric(env_data[ridx, ])
    # Expect at least 14 columns (as in the original model) - map them by depth bands
    # If env_data has fewer columns, fall back to simple noise
    if (length(rowv) >= 14) {
      # assign temperature by depth bands (approx mapping used originally)
      # we set temperature layer directly as Celsius (NetLogo did T = file-read)
      py <- temperature@pCoords[,2]
      temperature_vals <- rep(NA_real_, length(py))
      temperature_vals[py >= 18] <- rnorm(sum(py >= 18), mean = rowv[1], sd = 0.2)
      temperature_vals[py >= 16 & py < 18] <- rnorm(sum(py >= 16 & py < 18), mean = rowv[2], sd = 0.2)
      temperature_vals[py >= 14 & py < 16] <- rnorm(sum(py >= 14 & py < 16), mean = rowv[3], sd = 0.2)
      temperature_vals[py >= 12 & py < 14] <- rnorm(sum(py >= 12 & py < 14), mean = rowv[4], sd = 0.2)
      temperature_vals[py >= 8 & py < 12] <- rnorm(sum(py >= 8 & py < 12), mean = rowv[5], sd = 0.2)
      temperature_vals[py >= 4 & py < 8] <- rnorm(sum(py >= 4 & py < 8), mean = rowv[6], sd = 0.2)
      temperature_vals[py < 4] <- rnorm(sum(py < 4), mean = rowv[7], sd = 0.2)
      temperature <<- NLset(world = temperature, agents = patches(temperature), val = temperature_vals)
      
      # SWP mapping (columns 8..14)
      soil_moisture_vals <- rep(NA_real_, length(py))
      soil_moisture_vals[py >= 18] <- rnorm(sum(py >= 18), mean = rowv[8], sd = rowv[8]/10)
      soil_moisture_vals[py >= 16 & py < 18] <- rnorm(sum(py >= 16 & py < 18), mean = rowv[9], sd = rowv[9]/10)
      soil_moisture_vals[py >= 14 & py < 16] <- rnorm(sum(py >= 14 & py < 16), mean = rowv[10], sd = rowv[10]/10)
      soil_moisture_vals[py >= 12 & py < 14] <- rnorm(sum(py >= 12 & py < 14), mean = rowv[11], sd = rowv[11]/10)
      soil_moisture_vals[py >= 8 & py < 12] <- rnorm(sum(py >= 8 & py < 12), mean = rowv[12], sd = rowv[12]/10)
      soil_moisture_vals[py >= 4 & py < 8] <- rnorm(sum(py >= 4 & py < 8), mean = rowv[13], sd = rowv[13]/10)
      soil_moisture_vals[py < 4] <- rnorm(sum(py < 4), mean = rowv[14], sd = rowv[14]/10)
      soil_moisture <<- NLset(world = SWP, agents = patches(SWP), val = soil_moisture_vals)
      
      return(invisible(TRUE))
    }
  }
}
  
  setup_movement_Ex <- function() {
    # Copy data
    ex_vals <- energy_content_food@.Data
    p_vals  <- pcolor@.Data
    
    # Apply rules
    p_vals[ex_vals >= 0.10 & ex_vals < 0.20] <- p_vals[ex_vals >= 0.10 & ex_vals < 0.20] + 1
    p_vals[ex_vals >= 0.20 & ex_vals < 0.30] <- p_vals[ex_vals >= 0.20 & ex_vals < 0.30] + 2
    p_vals[ex_vals >= 0.30 & ex_vals < 0.40] <- p_vals[ex_vals >= 0.30 & ex_vals < 0.40] + 3
    p_vals[ex_vals >= 0.40]                  <- p_vals[ex_vals >= 0.40] + 4
    
    # Update world
    pcolor@.Data <<- p_vals
    
    # Then call SWP version
    setup_movement_SWP()
  }
  
  setup_movement_SWP <- function() {
    # Copy data
    swp_vals <- SWP@.Data
    p_vals   <- pcolor@.Data
    
    # Apply rules
    p_vals[swp_vals >= 20]                <- p_vals[swp_vals >= 20] + 1
    p_vals[swp_vals >= 15 & swp_vals < 20] <- p_vals[swp_vals >= 15 & swp_vals < 20] + 2
    p_vals[swp_vals >= 7.5 & swp_vals < 15] <- p_vals[swp_vals >= 7.5 & swp_vals < 15] + 3
    p_vals[swp_vals >= 2.5 & swp_vals < 7.5] <- p_vals[swp_vals >= 2.5 & swp_vals < 7.5] + 4
    p_vals[swp_vals < 2.5]                 <- p_vals[swp_vals < 2.5] + 5
    
    # Update world
    pcolor@.Data <<- p_vals
  }
  
#############   movement 
  move <- function(focal) {
    # worms: AgentMatrix (turtles)
    # pcolor_world: WorldMatrix with pcolor values
    
    # 1. Get current patch color of each worm
    current_patches <- cbind(of(agents=focal, var = "who"),cbind(of(agents = focal, var = "xcor"),
                             of(agents = focal, var = "ycor")))
    
    current_pcolor <- cbind(of(agents=focal, var = "who"),patchHere(world = pcolor, turtles = focal) )
      
    current_pcolor_vals<-of(world = pcolor,agents = patchHere(world = pcolor, turtles = focal) )
    
    # 2. For each worm, pick one random neighbor and one max-color neighbor
    neighbors_list <- neighbors(world = pcolor, agents = focal, nNeighbors = 4)
    
    # --- Random neighbor (using oneOf) ---
    random_neighbors <- oneOf(agents = neighbors_list)
    
    best_neighbors <- matrix(NA, nrow = nrow(focal), ncol = 3)  # columns: pxcor, pycor
    
    # --- Best neighbor (using maxOf) ---
    for (i in 1:nrow(focal)){
    best_neighbors[i,] <- maxOneOf(world = pcolor, agents = neighbors_list[neighbors_list[,'id']==i,]) 
    }
    best_neighbors_pcolor_vals<-of(world = pcolor, agents = best_neighbors )
    
    # --- 5. Déterminer qui fait un mouvement aléatoire ---
    move_random <- current_pcolor_vals >= best_neighbors_pcolor_vals
    who_movers <- of(agents = focal,var = "who")[move_random]
    who_best <- of(agents = focal,var = "who")[!move_random]
    focal_movers <- NLwith(agents = focal,var = "who",val = who_movers)
    focal_best<- NLwith(agents = focal,var = "who",val = who_best)
   
     # --- 6. Mouvement aléatoire ---
    if (any(move_random,na.rm = TRUE)) {
     # worms_random <- turtle(turtles = focal, who = focal[move_random, "who"])
      n_random <- length(who_movers)
      #points(worms, col= 'green')
      worms_random <- right(turtles = focal_movers, angle = sample(0:15, n_random, replace = TRUE))
      worms_random <- left(turtles = focal_movers, angle = sample(0:15, n_random, replace = TRUE))
      worms_random <- fd(turtles = focal_movers, dist = 0.5, world = pcolor)
    
      focal<-NLset(turtles = focal, agents = worms_random, var = "xcor", val=of(agents = worms_random,var="xcor"))
      focal<-NLset(turtles = focal, agents = worms_random, var = "ycor", val=of(agents = worms_random,var="ycor"))      
  
      }
    
    # --- 7. Mouvement vers le meilleur patch ---
    if (any(!move_random)) {

    target_patches <- patch(world = pcolor,
                              x = best_neighbors[!move_random, 1],
                              y = best_neighbors[!move_random, 2],duplicate = TRUE)
      
      worms_best <- moveTo(turtles = focal_best, agents = target_patches)
      
      focal<-NLset(turtles = focal, agents = worms_best, var = "xcor", val=of(agents = worms_best,var="xcor"))
      focal<-NLset(turtles = focal, agents = worms_best, var = "ycor", val=of(agents = worms_best,var="ycor"))      
    }
    return(focal)
  }
    
############   Mortality Home made from EEEWORM
  calc_mortality <- function(worms){
    
    mort_prob <- 0.00096   # 0.096 %
    
    # juveniles
    juv <- which(worms[,"breed"] == "juveniles")
    if (length(juv) > 0 && runif(1) < mort_prob) {
      dead <- oneOf(worms[juv, ])
      worms <<- die(turtles = worms, who = dead@.Data[,"who"])
    }
    
    # adults
    adu <- which(worms[,"breed"] == "adults")
    if (length(adu) > 0 && runif(1) < mort_prob) {
      dead <- oneOf(worms[adu, ])
      worms <<- die(turtles = worms, who = dead@.Data[,"who"])
    }
    
    return(worms)
  }
  
####### Incubation perdiod ######
  #####!!!!!!!!!!!!! SET T in Kelvins !!!!!!!!!!
 calc_embryo_development <- function(focal) {
    # On ne l'applique qu’aux cocons
    cocoons <- NLwith(agents = focal, var = "breed", val =  "cocoons")
    current_embryo_dev<-of(agents = cocoons,var = "embryonic_development")
   # current_T <-cbind(of(agents = cocoons,var = "who"),patchHere(world = temperature,turtles = cocoons))
    current_T_vals<- of(world = temperature,agents = patchHere(world = temperature,turtles = cocoons))
    incubation_period_i<-c()
    for (i in 1:nrow(cocoons)){
    # Formule d’incubation (NetLogo : 62 * e ^ ((-Ea/B) * ((1/RefT) - (1/T))))
    incubation_period_i[i] <- 62 * exp((-activation_energy / Boltz) * ((1 / reference_T) - (1 / current_T_vals[i])))
    }
    focal <- NLset(turtles = focal, agents = cocoons,
                     var = "embryonic_development",
                     val= current_embryo_dev + (1 / incubation_period_i) * 100 )
   return(focal)
  }
  
######### Transform cocoons ####
  transform_cocoon0 <- function(focal) {
    cocoons <- NLwith(agents = focal, var = "breed", val = "cocoons")
    ready_cocoons <- cocoons[which(of(agents =cocoons,var = "age" )>=incubation_period),]
    to_transform_cocoons <- ready_cocoons[which(of(agents =ready_cocoons,var = "embryonic_development" )>=100),]
    if(length(to_transform_cocoons)>0){
    focal <- NLset(turtles = focal, agents = to_transform_cocoons, var = "breed", val="juveniles")
    focal <- NLset(turtles = focal, agents = to_transform_cocoons, var = "color", val="pink")
    focal <- NLset(turtles = focal, agents = to_transform_cocoons, var = "size", val="0.2")
    }
    return(focal)
  }
  
  transform_cocoon <- function(focal) {
    
    # Select cocoons
    cocoons <- NLwith(agents = focal, var = "breed", val = "cocoons")
    
    # If no cocoon, skip to avoid errors
    if (length(cocoons) == 0) return(focal)
    
    # Filter cocoons by age
    age_vec <- of(agents = cocoons, var = "age")
    ready_cocoons <- NLwith(agents = cocoons, var = "age", val = age_vec[age_vec >= incubation_period])
    
    if (length(ready_cocoons) == 0) return(focal)
    
    # Filter by embryonic development
    dev_vec <- of(agents = ready_cocoons, var = "embryonic_development")
    to_transform_cocoons <- NLwith(agents = ready_cocoons,
                                   var = "embryonic_development",
                                   val = dev_vec[dev_vec >= 100])
    
    if (length(to_transform_cocoons) == 0) return(focal)
    
    # Apply transformations
    focal <- NLset(turtles = focal, agents = to_transform_cocoons, var = "breed", val = "juveniles")
    focal <- NLset(turtles = focal, agents = to_transform_cocoons, var = "color", val = "pink")
    focal <- NLset(turtles = focal, agents = to_transform_cocoons, var = "size",  val = 0.2)
  
    return(focal)
  }
  
  
######### Transform juveniles #####
  transform_juvenile0 <- function(focal) {
    juveniles <- NLwith(agents = focal, var = "breed", val = "juveniles")
    to_transform_juveniles <- juveniles[which(of(agents =juveniles,var = "mass" )>=mass_sexual_maturity),]
    if(length(to_transform_juveniles)>0){
      focal <- NLset(turtles = focal, agents = to_transform_juveniles, var = "breed", val="adults")
      focal <- NLset(turtles = focal, agents = to_transform_juveniles, var = "color", val="red")
      focal <- NLset(turtles = focal, agents = to_transform_juveniles, var = "size", val="0.3")
    }
    return(focal)
  }
  
  transform_juvenile <- function(focal) {
    
    # Select juveniles
    juveniles <- NLwith(agents = focal, var = "breed", val = "juveniles")
    
    if (length(juveniles) == 0) return(focal)
    
    # Filter by mass maturity
    mass_vec <- of(agents = juveniles, var = "mass")
    to_transform_juveniles <- NLwith(agents = juveniles,
                                     var = "mass",
                                     val = mass_vec[mass_vec >= mass_sexual_maturity])
    
    if (length(to_transform_juveniles) == 0) return(focal)
    
    # Apply transformations
    focal <- NLset(turtles = focal, agents = to_transform_juveniles, var = "breed", val = "adults")
    focal <- NLset(turtles = focal, agents = to_transform_juveniles, var = "color", val = "red")
    focal <- NLset(turtles = focal, agents = to_transform_juveniles, var = "size",  val = 0.3)
    return(focal)
  }
  
###### Aestivating 
  aestivate <- function(focal) {
    
    ## -----------------------------------------------------------
    ## 1. Extract who values and SWP at current patches
    ## -----------------------------------------------------------
    who_vals <- of(agents = focal, var = "who")
    SWP_vals <- of(world = SWP, agents = patchHere(world = SWP, turtles = focal))
    
    ## -----------------------------------------------------------
    ## 2. Extract internal states
    ## -----------------------------------------------------------
    aestivating_vals      <- of(agents = focal, var = "aestivating")
    time_aestivating_vals <- of(agents = focal, var = "time_aestivating")
    
    ## -----------------------------------------------------------
    ## 3. Logical conditions (matching NetLogo)
    ## -----------------------------------------------------------
    cond_high_SWP         <- SWP_vals >= 25
    cond_keep_aestivate   <- aestivating_vals == TRUE  & time_aestivating_vals <= 60
    cond_start_aestivate  <- aestivating_vals == FALSE
    cond_awake            <- aestivating_vals == TRUE &
      SWP_vals <= 20 &
      time_aestivating_vals >= 14
    
    ## -----------------------------------------------------------
    ## 4. Convert conditions → who IDs
    ## -----------------------------------------------------------
    who_high_SWP        <- who_vals[cond_high_SWP]
    who_keep_aestivate  <- who_vals[cond_keep_aestivate]
    who_start_aestivate <- who_vals[cond_start_aestivate]
    who_awake           <- who_vals[cond_awake]
    
    ## -----------------------------------------------------------
    ## 5. Apply NetLogo logic with NLwith() and NLset()
    ## -----------------------------------------------------------
    
    ## -------------------- CASE 1: High SWP (≥25) ----------------
    if (any(cond_high_SWP, na.rm = TRUE)) {
      
      ## 1.1 Already aestivating (and ≤60 days)
      if (length(who_keep_aestivate) > 0) {
        focal_agents <- NLwith(agents = focal, var = 'who', val = who_keep_aestivate)
        
        focal <- NLset(
          turtles = focal,
          agents  = focal_agents,
          var     = "aestivating",
          val     = TRUE
        )
        
        focal_agents <- calc_metabolic_rate(focal = focal_agents)
        focal<-upworms_with(group = focal_agents,focal = focal)
       # worms<<-upworms_with(group = focal,focal = worms)
      }
      
      ## 1.2 Not yet aestivating → start aestivating
      if (length(who_start_aestivate) > 0) {
        focal_agents <- NLwith(agents = focal, var = 'who', val = who_start_aestivate)
        
        focal <- NLset(
          turtles = focal,
          agents  = focal_agents,
          var     = "aestivating",
          val     = TRUE
        )
        
        focal_agents <- calc_metabolic_rate(focal = focal_agents)
        focal<-upworms_with(group = focal_agents,focal = focal)
      }
      
      ## -------------------- CASE 2: Wake up ------------------------
    } else if (length(who_awake) > 0) {
      
      focal_agents <- NLwith(agents = focal, var = 'who', val = who_awake)
      
      focal <- NLset(
        turtles = focal,
        agents  = focal_agents,
        var     = "aestivating",
        val     = FALSE
      )
      
      focal <- NLset(
        turtles = focal,
        agents  = focal_agents,
        var     = "time_aestivating",
        val     = 1
      )
      
      focal <- NLset(
        turtles = focal,
        agents  = focal_agents,
        var     = "visible",
        val     = TRUE
      )
      
    }

    return(focal)
  }
  
###### calc_metabolic_rate
  calc_metabolic_rate0<-function (focal){
    Mass <- of( agents = focal, var = "mass")
    current_T <- of( world =  temperature, agents = patchHere(world =  temperature,turtles =  focal))
    time_aestivating_vals<- of (agents = focal, var = "time_aestivating")
    energy_reserve_vals <- of( agents = focal, var = "energy_reserve")
    BMR_vals <- of( agents = focal, var = "BMR")
    
    focal <- NLset(turtles = focal, agents = focal, var = "BMR",
                   val= ((B_0 * (Mass ^ (3 / 4)) * exp( - activation_energy / ( Boltz * current_T)) ) * exp (-0.0032 * time_aestivating_vals)))
    
    # 🧠 Cas 1 : énergie disponible > BMR → payer maintenance
    cond_pay <- energy_reserve_vals > BMR_vals
    if (any(cond_pay, na.rm = TRUE)) {
      energy_reserve_cp<-of(agents = focal[cond_pay], var="energy_reserve")
      BMR_cp<-of(agents = focal[cond_pay], var="BMR")
      focal <- NLset(turtles = focal, agents = focal[cond_pay],var =  "energy_reserve", val =  energy_reserve_cp - BMR_cp)
    }
    
    # 🧠 Cas 2 : énergie insuffisante → puiser dans les réserves (masse)
    if (any(!cond_pay,na.rm = TRUE)) {
      mass_dcp<-of(agents = focal[!cond_pay], var="mass")
      energy_reserve_dcp<-of(agents = focal[!cond_pay], var="energy_reserve")
      BMR_dcp<-of(agents = focal[!cond_pay], var="BMR")
      focal <- NLset(turtles = focal, agents = focal[!cond_pay],var =  "mass", val = mass_dcp - (BMR_dcp/(energy_flesh+energy_synthesis)))
      focal <- NLset(turtles = focal,agents = focal[!cond_pay], var = "energy_reserve", val = energy_reserve_dcp+BMR_dcp)
    }
    
        # 🧬 Régression de stade si perte de masse
    cond_rejuvenate <- Mass < mass_sexual_maturity
    if (any(cond_rejuvenate)) {
      focal <- NLset(turtles = focal,agents = focal[cond_rejuvenate],var =  "breed", val = "juveniles")
      focal <- NLset(turtles = focal,agents = focal[cond_rejuvenate],var =  "color",val =  "pink")
      focal <- NLset(turtles = focal,agents = focal[cond_rejuvenate],var =  "size",val =  0.2)
    }
    
    # 💀 Mort si masse < masse de naissance
    cond_die <- Mass < mass_birth
    if (any(cond_die)) {
     for(i in length(cond_die)){
       worms <<- die(turtles = worms, who = of(agents = worms[cond_die][i], var = "who"))
     }
           }
    
  return(focal)
  }  
  
  calc_metabolic_rate <- function(focal) {
    
    # --- PREP BASE VALUES ---
    Mass  <- of(agents = focal, var = "mass")
    current_T <- of(world = temperature,
                    agents = patchHere(world = temperature, turtles = focal))
    time_aest <- of(agents = focal, var = "time_aestivating")
    energy_reserve <- of(agents = focal, var = "energy_reserve")
    BMR_old <- of(agents = focal, var = "BMR")
    
    # --- 1. Compute NEW BMR ---
    new_BMR <- (B_0 * Mass^(3/4) *
                  exp(-activation_energy / (Boltz * current_T))) *
      exp(-0.0032 * time_aest)
    
    focal <- NLset(turtles = focal, agents = focal, var = "BMR", val = new_BMR)
    
    
    # --- 2. Pay maintenance or mobilize reserves ---
    cond_pay <- energy_reserve > BMR_old
    who_can_pay <- of(agents = focal, var = "who")[cond_pay]
    who_cannot_pay <- of(agents = focal, var = "who")[!cond_pay]
      focal_can_pay<- NLwith(agents = focal,var = "who",val = who_can_pay)
      focal_cannot_pay<- NLwith(agents = focal,var = "who",val = who_cannot_pay)
      cond_pay_idx <- which(cond_pay)
      cond_nopay_idx <- which(!cond_pay)
    
    # Case 1: Enough energy → pay BMR from reserves
    if (length(cond_pay_idx) > 0) {
      er_pay <- energy_reserve[cond_pay_idx]
      bmr_pay <- BMR_old[cond_pay_idx]
      
      focal <- NLset(
        turtles = focal,
        agents = focal_can_pay,
        var = "energy_reserve",
        val = er_pay - bmr_pay
      )
     
    }
    
    # Case 2: Not enough → use mass
    if (length(cond_nopay_idx) > 0) {
      mass_no <- Mass[cond_nopay_idx]
      er_no <- energy_reserve[cond_nopay_idx]
      bmr_no <- BMR_old[cond_nopay_idx]
      
      # Decrease mass
      focal <- NLset(
        turtles = focal,
        agents = focal_cannot_pay,
        var = "mass",
        val = mass_no - (bmr_no / (energy_flesh + energy_synthesis))
      )
      
      # Increase reserve (mobilizing)
      focal <- NLset(
        turtles = focal,
        agents = focal_cannot_pay,
        var = "energy_reserve",
        val = er_no + bmr_no
      )
      
      }
    
    
    # --- 3. Stage regression if below maturity threshold ---
    Mass_after <- of(agents = focal, var = "mass")
    cond_rejuvenate <- which(Mass_after < mass_sexual_maturity)
    who_rejuv <- of(agents = focal,var =  "who")[cond_rejuvenate]
    focal_rejuv<-NLwith(agents = focal,var = "who",val = who_rejuv)
    
    if (length(cond_rejuvenate) > 0) {
      focal <- NLset(turtles = focal,agents =  focal_rejuv,var =  "breed",val =  "juveniles")
      focal <- NLset(turtles = focal,agents =  focal_rejuv,var =  "color",val =  "pink")
      focal <- NLset(turtles = focal,agents =  focal_rejuv,var =  "size", val = 0.2)
      }
  
    # --- 4. Death if mass < newborn mass ---
    cond_die <- which(Mass_after < mass_birth)
    who_die<-of(agents = focal,var = "who")[cond_die]
    focal_die<-NLwith(agents = focal, var = "who",val = who_die)
    
    if (length(cond_die) > 0) {
      dying_who <- of(agents = focal, var = "who")[cond_die]
      
      for (id in dying_who) {
        focal <- die( turtles = focal, who= id)
        worms <<- die(turtles = worms, who = id)
        cat("-")
      }
    }
      return(focal)
  }
  
  # Ingestion Rate ;;;;;;;;;;;;;;;;;;;;
   #juveniles and adults calculate their ingestion rate (the amount of food ingested from the environment) which depends on the food density of the patch in which they are present and the mass dependent maximum ingestion rate.
  calc_ingestion_rate <- function(focal) {
    
    # Extraire les infos nécessaires
    Mass <- of( agents = focal, var = "mass")
    current_SWP_vals <- of( world =  SWP, agents = patchHere(world =  SWP,turtles =  focal))
    current_food_vals <- of( world =  food_density, agents = patchHere(world =  food_density,turtles =  focal))
    ingestion_rate_vals <- rep(0, NLcount(focal))  # initialisation
    Arrhenius_here<-of(agents = focal, var = "Arrhenius_here")
    # 🧠 Étape 1 : vers sur patch avec nourriture disponible
    
    cond_food <- current_food_vals > 0
    
    if (any(cond_food)) {
      # Cas SWP < 10 : pas d’effet
      cond_swp_low <- cond_food & (current_SWP_vals < 10)
      who_cond_swp_low <-of(agents = focal,var ="who")[cond_swp_low]
      focal_cond_swp_low <-NLwith(agents = focal, var = "who", val = who_cond_swp_low)
      ingestion_rate_vals_cond_swp_low <-(max_ingestion_rate * Arrhenius_here[cond_swp_low]) * (Mass[cond_swp_low]^(2/3))
      focal<-NLset(turtles = focal,agents = focal_cond_swp_low, var = 'ingestion_rate', val=ingestion_rate_vals_cond_swp_low )
      # Cas SWP >= 10 : effet exponentiel négatif
      who_cond_swp_high <-of(agents = focal,var ="who")[!cond_swp_low]
      focal_cond_swp_high <-NLwith(agents = focal, var = "who", val = who_cond_swp_high)
      ingestion_rate_vals_cond_swp_high <- (max_ingestion_rate * Arrhenius_here[!cond_swp_low]) * (Mass[!cond_swp_low]^(2/3)) * exp(-0.04 * current_SWP_vals[!cond_swp_low])
      focal<-NLset(turtles = focal,agents = focal_cond_swp_high, var = 'ingestion_rate', val=ingestion_rate_vals_cond_swp_high )
      
    }
    
    # 🧠 Étape 2 : compétition locale — si somme ingestion > nourriture dispo
    # pour chaque patch, on vérifie la somme des ingestion_rates
   
    # 🧭 1. Extract worms data (rounded patch coords)
    worm_coords <- round(coordinates(focal))
    ingestion   <- of(agents = focal,var =  "ingestion_rate")
    who <- of(agents = focal,var = "who")
    
    worms_df <- data.frame(
      x = worm_coords[, 1],
      y = worm_coords[, 2],
      ingestion_rate = ingestion,
      who=who
    )
    
    # 🧮 2. Aggregate total ingestion per patch
    worms_sum <- worms_df %>%
      group_by(x, y) %>%
      summarise(
        total_ingestion = sum(ingestion_rate, na.rm = TRUE),
        count_worms = n(),
        .groups = "drop"
      )
    
    # 🌍 3. Extract patch data
    patch_coords <- (patchHere(world = food_density,turtles = focal))
    food_vals    <- of(world = food_density, agents =  (patchHere(world = food_density,turtles = focal)))

    patches_df <- data.frame(
      x = patch_coords[, 1],
      y = patch_coords[, 2],
      food_density = food_vals
    )
    
    # 🔗 4. Merge worms and patches on (x, y)
    merged_df <- left_join(patches_df, worms_sum, by = c("x", "y"))
    
    # 🧪 5. Compare total_ingestion vs food_density
    merged_df <- merged_df %>%
      mutate(
        ingestion_exceeds_food = total_ingestion > food_density
      )
    
    if(any(duplicated(patch_coords))==TRUE){
    focal<-NLset(turtles = focal, agents = focal, var= "ingestion_rate", val= merged_df$food_density/merged_df$count_worms )
    food_density<- NLset(world = food_density,agents = (patch_coords), val= 0)
    }
      
    if(any(food_vals<=0)){
      no_food <- NLwith(world = food_density, agents = patch_coords, val=0)
      worm_no_food<-turtlesOn(world = food_density, turtles = focal, agents = no_food,simplify = TRUE)
      focal<-NLset(turtles = focal, agents = worm_no_food, var = "ingestion_rate", val=0)    
    } 
    
    # 🧬 Assimilation step (depends on ingestion)
    focal <- calc_assimilation(focal = focal,energy_content_food = energy_content_food )
    focal<-NLset(turtles = focal,agents = focal,var = "energy_assimilated",
                 val= of(agents = focal,var = "energy_assimilated"))
    
    return(focal)
  }
  
  calc_ingestion_rate1 <- function(focal) {
    
    # --- Extract basic values ---
    Mass    <- of(agents = focal, var = "mass")
    Arrh    <- of(agents = focal, var = "Arrhenius_here")
    
    current_SWP <- of(world = SWP,
                      agents = patchHere(world = SWP, turtles = focal))
    current_food <- of(world = food_density,
                       agents = patchHere(world = food_density, turtles = focal))
    
    n <- NLcount(focal)
    ingestion_rate_vals <- rep(0, n)
    
    who_all <- of(agents = focal, var = "who")
    
    # ============================================================
    # 🧠 STEP 1 — Local ingestion rate before competition
    # ============================================================
    cond_food <- which(current_food > 0)
    if (length(cond_food) > 0) {
      # CAS 1: SWP < 10
      cond_low_swp <- cond_food[current_SWP[cond_food] < 10]
      if (length(cond_low_swp) > 0) {
        idx <- cond_low_swp
        ingestion_rate_vals[idx] <-
          (max_ingestion_rate * Arrh[idx]) * (Mass[idx]^(2/3))
        
        focal <- NLset(
          turtles = focal,
          agents = NLwith(agents = focal,var =  "who",val =  who_all[idx]),
          var = "ingestion_rate",
          val = ingestion_rate_vals[idx]
        )
      }
      
      # CAS 2: SWP >= 10 (negative exponential)
      cond_high_swp <- cond_food[current_SWP[cond_food] >= 10]
      if (length(cond_high_swp) > 0) {
        idx <- cond_high_swp
        ingestion_rate_vals[idx] <-
          (max_ingestion_rate * Arrh[idx]) * (Mass[idx]^(2/3)) *
          exp(-0.04 * current_SWP[idx])
        
        focal <- NLset(
          turtles = focal,
          agents = NLwith(agents = focal, var = "who", val = who_all[idx]),
          var = "ingestion_rate",
          val = ingestion_rate_vals[idx]
        )
      }
    }
    
    
    # Refresh ingestion after update
    ingestion <- of(agents = focal, var = "ingestion_rate")
    # ============================================================
    # 🧭 STEP 2 — Patch-level competition 
    # ============================================================
    
    # Worm coords
    worm_xy <- round(coordinates(focal))
    xw <- worm_xy[, 1]
    yw <- worm_xy[, 2]
    
    # Patch food
    patch_xy <- patchHere(world = food_density, turtles = focal)
    food_vals <- of(world = food_density, agents = patch_xy)
    
    # Build worm table
    worms_df <- data.frame(
      who = who_all,
      x = xw,
      y = yw,
      ingestion = ingestion
    )
    
    patch_df <- data.frame(
      x = patch_xy[,1],
      y = patch_xy[,2],
      food = food_vals
    )
    
    # Aggregate ingestion by patch
    agg <- worms_df %>%
      group_by(x, y) %>%
      summarise(
        total_ing = sum(ingestion),
        n_worms   = n(),
        .groups = "drop"
      )
    
    merged <- left_join(patch_df, agg, by = c("x","y"))
    
    # ============================================================
    # 🧪 STEP 3 — If ingestion > food, ration equally
    # ============================================================
    
    exceed_idx <- which(merged$total_ing > merged$food)
    
    if (length(exceed_idx) > 0) {
      for (j in exceed_idx) {
        px <- merged$x[j]
        py <- merged$y[j]
        
        # worms on that patch
        worms_on_patch <- NLwith(agents = focal, var="who", 
                                 val = worms_df$who[worms_df$x==px & worms_df$y==py])
        
        if (NLcount(worms_on_patch) > 0) {
          ration <- merged$food[j] / merged$n_worms[j]
          
          focal <- NLset(
            turtles = focal,
            agents = worms_on_patch,
            var = "ingestion_rate",
            val = ration
          )
          
          # Patch is fully consumed
          food_density <<- NLset(
            world = food_density,
            agents = NLwith(world = food_density,
                            var = "pxcor", val = px) %>%
            NLwith(world = food_density,
                   var = "pycor", val = py),
            val = 0
          )
        }
      }
    }
    
    
    # ============================================================
    # 🥀 STEP 4 — If food ≤ 0 → ingestion = 0
    # ============================================================
    
    no_food_idx <- which(food_vals <= 0)
    if (length(no_food_idx) > 0) {
      who_no_food <- who_all[no_food_idx]
      focal <- NLset(
        turtles = focal,
        agents = NLwith(agents = focal,var =  "who", val = who_no_food),
        var = "ingestion_rate",
        val = 0
      )
    }
    
    
    # ============================================================
    # 🧬 STEP 5 — Assimilation
    # ============================================================
    
    focal <- calc_assimilation(focal, energy_content_food)
    return(focal)
  }
  
  # calc_assimilation ;;;;;;;;;;;;;;;;;;;;
 calc_assimilation0<-function(focal,energy_content_food){
   ingestion_rate_vals<-of(agents = focal,var = "ingestion_rate")
    pycor<-energy_content_food@pCoords[,2]
    
          if (day > 335 || day <= 32 ){
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor>=18,], val =  rnorm(n = length(energy_content_food[pycor>=18]) ,mean = 0.38,sd =  0.038))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 16 & pycor < 18), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 16 & pycor < 18) ]) ,mean = 0.36,sd =  0.036))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 14 & pycor < 16), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 14 & pycor < 16) ]) ,mean = 0.35,sd =  0.035))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 12 & pycor < 14), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 12 & pycor < 14) ]) ,mean = 0.33,sd =  0.033))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 8 & pycor < 12) ,], val =  rnorm(n = length(energy_content_food[ (pycor >= 8 & pycor < 12) ]) ,mean = 0.32,sd =  0.032))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 4 & pycor < 8) ,], val =  rnorm(n = length(energy_content_food[ (pycor >= 4 & pycor < 8) ]) ,mean = 0.29,sd =  0.029))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor<4,], val =  rnorm(n = length(energy_content_food[pycor<4]) ,mean = 0.26,sd =  0.026))
  }
      if (day > 32 && day <= 121 ){
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor>=18,], val =  rnorm(n = length(energy_content_food[pycor>=18]) ,mean = 0.44,sd =  0.044))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 16 & pycor < 18), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 16 & pycor < 18) ]) ,mean = 0.42,sd =  0.042))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 14 & pycor < 16), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 14 & pycor < 16) ]) ,mean = 0.4,sd =  0.04))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 12 & pycor < 14), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 12 & pycor < 14) ]) ,mean = 0.38,sd =  0.038))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 8 & pycor < 12) ,], val =  rnorm(n = length(energy_content_food[ (pycor >= 8 & pycor < 12) ]) ,mean = 0.36,sd =  0.036))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 4 & pycor < 8) ,], val =  rnorm(n = length(energy_content_food[ (pycor >= 4 & pycor < 8) ]) ,mean = 0.33,sd =  0.033))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor<4,], val =  rnorm(n = length(energy_content_food[pycor<4]) ,mean = 0.3,sd =  0.03))
  }  
      if (day > 121 && day <= 213){
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor>=18,], val =  rnorm(n = length(energy_content_food[pycor>=18]) ,mean = 0.49,sd =  0.049))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 16 & pycor < 18), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 16 & pycor < 18) ]) ,mean = 0.46,sd =  0.046))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 14 & pycor < 16), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 14 & pycor < 16) ]) ,mean = 0.44,sd =  0.044))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 12 & pycor < 14), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 12 & pycor < 14) ]) ,mean = 0.41,sd =  0.041))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 8 & pycor < 12) ,], val =  rnorm(n = length(energy_content_food[ (pycor >= 8 & pycor < 12) ]) ,mean = 0.4,sd =  0.04))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 4 & pycor < 8), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 4 & pycor < 8) ]) ,mean = 0.36,sd =  0.036))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor<4,], val =  rnorm(n = length(energy_content_food[pycor<4]) ,mean = 0.33,sd =  0.033))
  } 
      if (day > 213 && day <= 335){
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[pycor>=18,], val =  rnorm(n = length(energy_content_food[pycor>=18]) ,mean = 0.41,sd =  0.041))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 16 & pycor < 18), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 16 & pycor < 18) ]) ,mean = 0.37,sd =  0.037))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 14 & pycor < 16),], val =  rnorm(n = length(energy_content_food[ (pycor >= 14 & pycor < 16) ]) ,mean = 0.37,sd =  0.037))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 12 & pycor < 14), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 12 & pycor < 14) ]) ,mean = 0.35,sd =  0.035))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 8 & pycor < 12), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 8 & pycor < 12) ]) ,mean = 0.33,sd =  0.033))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[ (pycor >= 4 & pycor < 8), ], val =  rnorm(n = length(energy_content_food[ (pycor >= 4 & pycor < 8) ]) ,mean = 0.3,sd =  0.03))
    energy_content_food<<-NLset(world = energy_content_food, agents = patches(energy_content_food)[(pycor<4),], val =  rnorm(n = length(energy_content_food[pycor<4]) ,mean = 0.28,sd =  0.028))
      } 
    
    ingestion_rate_vals<-of(agents = focal, var = "ingestion_rate")
    energy_content_food_here<-of(world = energy_content_food,agents = patchHere(world = energy_content_food,turtles = focal))
    focal<-NLset(turtles = focal,agents = focal,var = "energy_assimilated",
          val= ingestion_rate_vals* (energy_content_food_here*assimilation_efficiency))
   
   return(focal) 
  }
 calc_assimilation <- function(focal, energy_content_food) {
   
   # --- 0. quick exit if no agents
   if (NLcount(focal) == 0) return(focal)
   
   # --- 1. Prepare patch-level energy-content vector (seasonal by day)
   pycor <- energy_content_food@pCoords[, 2]
   npatches <- length(pycor)
   energy_vals <- numeric(npatches)
   
   # seasonal means & sds for depth bins
   # Define helper to fill energy_vals for an index mask
   fill_vals <- function(mask, mean_v, sd_v) {
     if (any(mask)) {
       energy_vals[mask] <<- stats::rnorm(sum(mask), mean = mean_v, sd = sd_v)
     }
   }
   
   if (day > 335 || day <= 32) {
     fill_vals(pycor >= 18,     mean_v = 0.38, sd_v = 0.038)
     fill_vals(pycor >= 16 & pycor < 18, mean_v = 0.36, sd_v = 0.036)
     fill_vals(pycor >= 14 & pycor < 16, mean_v = 0.35, sd_v = 0.035)
     fill_vals(pycor >= 12 & pycor < 14, mean_v = 0.33, sd_v = 0.033)
     fill_vals(pycor >= 8  & pycor < 12, mean_v = 0.32, sd_v = 0.032)
     fill_vals(pycor >= 4  & pycor < 8,  mean_v = 0.29, sd_v = 0.029)
     fill_vals(pycor < 4,      mean_v = 0.26, sd_v = 0.026)
   } else if (day > 32 && day <= 121) {
     fill_vals(pycor >= 18,     mean_v = 0.44, sd_v = 0.044)
     fill_vals(pycor >= 16 & pycor < 18, mean_v = 0.42, sd_v = 0.042)
     fill_vals(pycor >= 14 & pycor < 16, mean_v = 0.40, sd_v = 0.040)
     fill_vals(pycor >= 12 & pycor < 14, mean_v = 0.38, sd_v = 0.038)
     fill_vals(pycor >= 8  & pycor < 12, mean_v = 0.36, sd_v = 0.036)
     fill_vals(pycor >= 4  & pycor < 8,  mean_v = 0.33, sd_v = 0.033)
     fill_vals(pycor < 4,      mean_v = 0.30, sd_v = 0.030)
   } else if (day > 121 && day <= 213) {
     fill_vals(pycor >= 18,     mean_v = 0.49, sd_v = 0.049)
     fill_vals(pycor >= 16 & pycor < 18, mean_v = 0.46, sd_v = 0.046)
     fill_vals(pycor >= 14 & pycor < 16, mean_v = 0.44, sd_v = 0.044)
     fill_vals(pycor >= 12 & pycor < 14, mean_v = 0.41, sd_v = 0.041)
     fill_vals(pycor >= 8  & pycor < 12, mean_v = 0.40, sd_v = 0.040)
     fill_vals(pycor >= 4  & pycor < 8,  mean_v = 0.36, sd_v = 0.036)
     fill_vals(pycor < 4,      mean_v = 0.33, sd_v = 0.033)
   } else if (day > 213 && day <= 335) {
     fill_vals(pycor >= 18,     mean_v = 0.41, sd_v = 0.041)
     fill_vals(pycor >= 16 & pycor < 18, mean_v = 0.37, sd_v = 0.037)
     fill_vals(pycor >= 14 & pycor < 16, mean_v = 0.37, sd_v = 0.037)
     fill_vals(pycor >= 12 & pycor < 14, mean_v = 0.35, sd_v = 0.035)
     fill_vals(pycor >= 8  & pycor < 12, mean_v = 0.33, sd_v = 0.033)
     fill_vals(pycor >= 4  & pycor < 8,  mean_v = 0.30, sd_v = 0.030)
     fill_vals(pycor < 4,      mean_v = 0.28, sd_v = 0.028)
   }
   
   # write the patch-level energy content in a single call
   energy_content_food <<- NLset(
     world = energy_content_food,
     agents = patches(energy_content_food),
     val = energy_vals
   )
   
   # --- 2. Compute assimilation per turtle ---
   ingestion_rate_vals <- of(agents = focal, var = "ingestion_rate")
   energy_content_here <- of(world = energy_content_food,
                             agents = patchHere(world = energy_content_food, turtles = focal))
   
   # energy_assimilated = ingestion_rate * (energy_content_here * assimilation_efficiency)
   assimilated_vals <- ingestion_rate_vals * (energy_content_here * assimilation_efficiency)
   
   focal <- NLset(
     turtles = focal,
     agents = focal,
     var = "energy_assimilated",
     val = assimilated_vals
   )
   return(focal)
 }
 
 #   CHANGE IN FOOD DENSITY
 ###==============================
 calc_change_food_density <- function(worms, change_food_density) {
   # For each patch, compute the change in food density per timestep
   if(any(food_density <= 0)){
   change_food_density <- NLset(world = change_food_density,agents = patches(change_food_density),val= 0)
   }else{
     # 🧭 1. Extract worms data (rounded patch coords)
     worm_coords <- round(coordinates(worms))
     ingestion   <- of(agents = worms,var =  "ingestion_rate")
     who <- of(agents = worms,var = "who")
     
     worms_df <- data.frame(
       x = worm_coords[, 1],
       y = worm_coords[, 2],
       ingestion_rate = ingestion,
       who=who
     )
     
     # 🧮 2. Aggregate total ingestion per patch
     worms_sum <- worms_df %>%
       group_by(x, y) %>%
       summarise(
         total_ingestion = sum(ingestion_rate, na.rm = TRUE),
         count_worms = n(),
         .groups = "drop"
       )
     
     # 🌍 3. Extract patch data
     patch_coords <- (patchHere(world = change_food_density,turtles = worms))
     food_vals    <- of(world = change_food_density, agents =  (patchHere(world = change_food_density,turtles = worms)))
     
     patches_df <- data.frame(
       x = patch_coords[, 1],
       y = patch_coords[, 2],
       food_density = food_vals
     )
     
     # 🔗 4. Merge worms and patches on (x, y)
     merged_df <- left_join(patches_df, worms_sum, by = c("x", "y"))
     
     # 🧪 5. Compare total_ingestion vs food_density
     merged_df <- merged_df %>%
       mutate(
         ingestion_exceeds_food = total_ingestion > food_density
       )
     ingestion_rate_sum_vals <-merged_df$ingestion_exceeds_food
     
     change_food_density <- NLset(world = change_food_density,agents = worm_coords,val= ingestion_rate_sum_vals)
   }
  return(change_food_density)   
 }
 
 ###==============================
 ###   SOMATIC MAINTENANCE
 ###==============================
 calc_maintenance0 <- function(focal) {
   # Clamp coords inside world boundaries
   focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor_move), max_Pxcor_move)
   focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor_move), max_Pycor_move)
   
    temperatureh_here_vals<-of(world = temperature,agents = patchHere(world = temperature,turtles = focal))
    mass_vals<-of(agents = focal, var = "mass")
    energy_assimilated_vals<-of(agents = focal, var = "energy_assimilated")
    energy_reserve_vals<-of(agents = focal, var = "energy_reserve")
    energy_reserve_max_vals<-of(agents = focal, var = "energy_reserve_max")
   
    # Basal Metabolic Rate (BMR)
    BMR_vals <- B_0 * (mass_vals^(3/4)) * exp(-activation_energy / (Boltz * temperatureh_here_vals))
    focal<-NLset(turtles = focal, agents = focal,var = "BMR",val = BMR_vals)
    
    ea_above<-(energy_assimilated_vals > BMR_vals)
      # Energy balance logic
     focal<-NLset(turtles = focal,agents = focal[ea_above] ,var = "energy_assimilated",val = energy_assimilated_vals[ea_above] - BMR_vals[ea_above])
     focal<-NLset(turtles = focal,agents = focal[!ea_above],var = "energy_reserve",val = energy_reserve_vals[!ea_above] + energy_assimilated_vals[!ea_above])
     focal<-NLset(turtles = focal,agents = focal[!ea_above],var = "energy_reserve",val = energy_reserve_vals[!ea_above] - BMR_vals[!ea_above])
     focal<-NLset(turtles = focal,agents = focal[!ea_above],var = "energy_assimilated",val = 0)

   # Starvation conditions
   if (any(energy_reserve_vals < (energy_reserve_max_vals * 0.5))) {
     w_adults<-NLwith(agents = focal[(energy_reserve_vals < (energy_reserve_max_vals * 0.5)) ],var = 'breed',val = 'adults')
     w_adults <- onset_starvation_strategy(focal= w_adults)
     focal<-upworms_with(group = w_adults,focal  = focal)
   }
    if (any(energy_reserve_vals < (energy_reserve_max_vals * 0.5))) {
     w_juveniles<-NLwith(agents = focal[(energy_reserve_vals < (energy_reserve_max_vals * 0.5))],var = 'breed',val = 'juveniles')
     w_juveniles <- onset_starvation_strategy(focal=w_juveniles)
     focal<-upworms_with(group = w_juveniles, focal  = focal)
   }
   
   return(focal)
 }
 calc_maintenance <- function(focal) {
   
   # =============================
   # 1. SAFELY CLAMP COORDINATES
   # =============================
   # We must use NLset — direct slot modification is illegal in RNetLogo.
   x_vals <- of(agents = focal, var = "xcor")
   y_vals <- of(agents = focal, var = "ycor")
   
   x_vals <- pmin(pmax(x_vals, min_Pxcor_move), max_Pxcor_move)
   y_vals <- pmin(pmax(y_vals, min_Pycor_move), max_Pycor_move)
   
   focal <- NLset(turtles = focal, agents = focal, var = "xcor", val = x_vals)
   focal <- NLset(turtles = focal, agents = focal, var = "ycor", val = y_vals)
   
   # =============================
   # 2. RETRIEVE VALUES
   # =============================
   temp_here <- of(world = temperature, agents = patchHere(world = temperature, turtles = focal))
   mass_vals <- of(agents = focal, var = "mass")
   EA_vals   <- of(agents = focal, var = "energy_assimilated")
   ER_vals   <- of(agents = focal, var = "energy_reserve")
   ERmax_vals <- of(agents = focal, var = "energy_reserve_max")
   breeds     <- of(agents = focal, var = "breed")
   
   # =============================
   # 3. Compute BMR and store
   # =============================
   BMR_vals <- B_0 * (mass_vals^(3/4)) * exp(-activation_energy / (Boltz * temp_here))
   focal <- NLset(turtles = focal, agents = focal, var = "BMR", val = BMR_vals)
   
   # =============================
   # 4. ENERGY BALANCE
   # =============================
   ea_above <- EA_vals > BMR_vals
   
   # Case 1: enough assimilation to pay BMR
   if (any(ea_above)) {
     focal <- NLset(
       turtles = focal,
       agents = focal[ea_above],
       var = "energy_assimilated",
       val = EA_vals[ea_above] - BMR_vals[ea_above]
     )
   }
   
   # Case 2: not enough → draw from reserves
   if (any(!ea_above)) {
     # Update reserves first: ER = ER + EA - BMR
     new_reserve <- ER_vals[!ea_above] + EA_vals[!ea_above] - BMR_vals[!ea_above]
     
     focal <- NLset(
       turtles = focal,
       agents = focal[!ea_above],
       var = "energy_reserve",
       val = new_reserve
     )
     
     # assimilation reset to 0
     focal <- NLset(
       turtles = focal,
       agents = focal[!ea_above],
       var = "energy_assimilated",
       val = 0
     )
   }
   
   # REFRESH after update
   ER_vals <- of(agents = focal, var = "energy_reserve")
   
   # =============================
   # 5. STARVATION LOGIC
   # =============================
   starving <- ER_vals < (0.5 * ERmax_vals)
   if (any(starving)) {
     
     # --- adults
     adult_ids <- which(starving & breeds == "adults")
     if (length(adult_ids) > 0) {
       adults <- focal[adult_ids]
       adults <- onset_starvation_strategy(focal = adults)
       focal  <- upworms_with(group = adults, focal = focal)
     }
     
     # --- juveniles
     juv_ids <- which(starving & breeds == "juveniles")
     if (length(juv_ids) > 0) {
       juveniles <- focal[juv_ids]
       juveniles <- onset_starvation_strategy(focal = juveniles)
       focal     <- upworms_with(group = juveniles, focal = focal)
     }
   }
   
   return(focal)
 }
 
 ###==============================
 ###   STARVATION STRATEGY
 ###==============================
 onset_starvation_strategy0 <- function(focal) {
   mass_vals<-of(agents = focal, var = "mass")
   energy_reserve_vals<-of(agents = focal, var = "energy_reserve")
   temperatureh_here_vals<-of(world = temperature,agents = patchHere(world = temperature,turtles = focal))
   
   BMR_vals <- B_0 * (mass_vals^(3/4)) * exp(-activation_energy / (Boltz * temperatureh_here_vals))
   
    focal<-NLset(turtles = focal, agents = focal, var = "mass", val= mass_vals - (BMR_vals/(energy_flesh+energy_synthesis)))
    focal<-NLset(turtles = focal, agents = focal, var = "energy_reserve", val= energy_reserve_vals+BMR_vals)
 
   mass_vals <- mass_vals - (BMR_vals / (energy_flesh + energy_synthesis))
   energy_reserve_vals <- energy_reserve_vals + BMR_vals
   
   # Regression to juvenile if mass < puberty mass
   if (any(mass_vals < mass_sexual_maturity)) {
     focal<-NLset(turtles = focal, agents = focal[(mass_vals < mass_sexual_maturity)], var="breed", val = "juveniles")
     focal<-NLset(turtles = focal, agents = focal[(mass_vals < mass_sexual_maturity)], var="color", val = "pink")
     focal<-NLset(turtles = focal, agents = focal[(mass_vals < mass_sexual_maturity)], var= "size", val = 0.2)
   }
   
   # Mortality if mass < birth mass
   if (any(mass_vals < mass_birth)) {
     focal_die<-NLwith(agents = focal, var = "mass",val = mass_vals[which(mass_vals < mass_birth)])
    who_die<- of(agents = focal_die,var = "who") 
   worms<<-die(turtles = worms,who = who_die)
   }
   
   return(focal)
 }
 onset_starvation_strategy <- function(focal) {
   
   # Extract current values
   mass_vals  <- of(agents = focal, var = "mass")
   er_vals    <- of(agents = focal, var = "energy_reserve")
   temp_vals  <- of(world = temperature,
                    agents = patchHere(world = temperature, turtles = focal))
   
   # Compute BMR
   BMR_vals <- B_0 * (mass_vals^(3/4)) *
     exp(-activation_energy / (Boltz * temp_vals))
   
   ### -------------------------------------------
   ### UPDATE MASS AND ENERGY FIRST
   ### -------------------------------------------
   new_mass <- mass_vals - (BMR_vals / (energy_flesh + energy_synthesis))
   new_er   <- er_vals + BMR_vals
   
   focal <- NLset(turtles = focal, agents = focal,
                  var = "mass", val = new_mass)
   
   focal <- NLset(turtles = focal, agents = focal,
                  var = "energy_reserve", val = new_er)
   
   ### -------------------------------------------
   ### RECHECK UPDATED VALUES
   ### -------------------------------------------
   mass_vals <- new_mass      # refresh
   er_vals   <- new_er
   
   ### -------------------------------------------
   ### REGRESSION TO JUVENILE
   ### -------------------------------------------
   idx_regress <- which(mass_vals < mass_sexual_maturity)
   
   if (length(idx_regress) > 0) {
     focal_regress <- focal[idx_regress]
     
     focal <- NLset(turtles = focal, agents = focal_regress,
                    var="breed", val="juveniles")
     focal <- NLset(turtles = focal, agents = focal_regress,
                    var="color", val="pink")
     focal <- NLset(turtles = focal, agents = focal_regress,
                    var="size",  val=0.2)
   }
   
   ### -------------------------------------------
   ### MORTALITY (mass < birth mass)
   ### -------------------------------------------
   idx_die <- which(mass_vals < mass_birth)
   
   if (length(idx_die) > 0) {
     focal_die <- focal[idx_die]
     who_die   <- of(agents = focal_die, var = "who")
     
     # kill worms in global population
     focal <- die(turtles = focal, who = who_die)
     worms <<- die(turtles = worms, who = who_die)
   }
   
   return(focal)
 }
 
 ### ==============================
 ###   REPRODUCTION
 ### ==============================
 mate <- function(focal) {
   if (any(of(agents = focal,var = "mass") >= mass_sexual_maturity) &&
       any(duplicated(patchHere(world = food_density,turtles = focal)))) {
     
     # 🧭 1. Coordonnées et ID individuels
     worm_coords <- round(coordinates(focal))
     who <- of(agents = focal, var = "who")
     
     worms_df <- data.frame(
       x = worm_coords[, 1],
       y = worm_coords[, 2],
       who = who
     )
     
     # 🧮 2. Identifier les groupes de vers par patch (x, y)
     worms_groups <- worms_df %>%
       group_by(x, y) %>%
       mutate(
         count_worms = n()
       ) %>%
       ungroup()
     
     # 🎯 3. Ne garder que les groupes avec ≥ 2 vers
     paired_groups <- worms_groups %>%
       filter(count_worms >= 2)%>% group_by(x, y) %>% mutate(  group_id = cur_group_id())
     
     # 🚀 4. Pour chaque groupe, choisir un ver et le déplacer vers un autre ver du même groupe
     if (nrow(paired_groups) > 0) {
       for (gid in unique(paired_groups$group_id)) {
         group_members <- paired_groups[paired_groups$group_id == gid, ]
         
         if (nrow(group_members) >= 2) {
           # Sélectionner aléatoirement un ver "mâle" et une cible du même groupe
           mover <- sample(group_members$who, 1)
           if( (nrow(group_members) == 2)){
             target <- setdiff(group_members$who, mover)  
           }else{
             target <- sample(setdiff(group_members$who, mover),1)
           }
           
           
           # Mouvement : le ver se déplace vers la cible
           focal_mover <- moveTo(turtles = focal[which(focal[,"who"]==mover)], agents = focal[which(focal[,"who"]==target)])
           focal <- NLset(turtles = focal ,agents =focal_mover ,var = "xcor" ,val = of(agents = focal_mover,var = "xcor"))
           focal <- NLset(turtles = focal ,agents =focal_mover ,var = "ycor" ,val = of(agents = focal_mover,var = "ycor"))
         }
       }
     }
     
     # 5. Reproduction pour les individus matures
     focal <- calc_reproduction(focal)
   }
   
   return(focal)
 }
 
 calc_reproduction <- function(focal) {
   mass_vals<-of(agents = focal,var = "mass")
   Arrhenius_here<-of(agents = focal, var = "Arrhenius_here")
   #Arrhenius_here<-patchHere(world = Arrhenius,turtles = focal)
   focal <- NLset(turtles = focal, agents = focal,var = "max_R",val =  (max_reproduction_rate * Arrhenius_here) * mass_vals)
   energy_assimilated_vals <- of(agents = focal,var = "energy_assimilated")
   max_R_vals<-of(agents = focal,var = "max_R")
   R_vals<-of(agents = focal,var = "R")
   # Energy allocation to reproduction
   if (any(energy_assimilated_vals >= max_R_vals)) {
     ea_above_maxR<-(energy_assimilated_vals >= max_R_vals)
     who_above_max_R <- of(agents = focal,var = "who")[which(energy_assimilated_vals >= max_R_vals)]
     focal_above_max_R <- NLwith(agents = focal,var = "who",val = who_above_max_R)
     focal<- NLset(turtles= focal, agents= focal_above_max_R, var="energy_assimilated", val= energy_assimilated_vals[ea_above_maxR] - max_R_vals[ea_above_maxR])
     focal<- NLset(turtles= focal, agents= focal_above_max_R, var="R",val= R_vals[ea_above_maxR] + max_R_vals[ea_above_maxR])
   } else if (any(energy_assimilated_vals > 0)) {
     ea_above_zero<-(energy_assimilated_vals > 0)
     who_above_zero<-of(agents = focal,var = "who")[which(energy_assimilated_vals > 0)]
     focal_above_zero<- NLwith(agents = focal,var = "who",val = who_above_zero)
     focal<- NLset(turtles= focal,agents= focal_above_zero,var="R",val= R_vals[ea_above_zero] + max_R_vals[ea_above_zero])
     # focal<- NLset(turtles= focal,agents= focal_above_zero,var="energy_assimlated",val= rep(0,nrow(focal_above_zero)))     
   }
   
   if (any(energy_assimilated_vals > 0)) { 
     focal<-calc_growth(focal) 
   }
   
   mass_cocoon<- -0.01768 + 0.04035 * mass_vals
   focal<-NLset(turtles = focal, agents = focal, var = "mass_cocoons", val= mass_cocoon)
   if (any(R_vals >= mass_cocoon * (energy_flesh + energy_synthesis))) {
     who_above_R<-of(agents = focal,var = "who")[which(R_vals >= mass_cocoon * (energy_flesh + energy_synthesis))]
     focal_above_zero<- NLwith(agents = focal,var = "who",val = who_above_R)
     focal<-reproduce(focal, focal_above_zero)
     #focal<-upworms_with(group = focal_above_zero,focal = focal)
   }
   
   return(focal)
 }
 
 reproduce <- function(focal, repro_agents) {
   who_reproduce <- of(agents = repro_agents,var = "who")
   ids_before <- of(agents = focal, var = "who")
     focal <- hatch(turtles = focal, who = who_reproduce, n = 1, breed = "cocoons")
     cat("new cocoons", nrow(focal)-length(ids_before))
   
   # detect new IDs (robust even if parent fields were copied)
   ids_after <- of(agents = focal, var = "who")
   new_ids <- setdiff(y = ids_before,x = ids_after)
   
   if (length(new_ids) > 0) {
     new_cocoons <- NLwith(agents = focal, var = "who", val = new_ids)
     
     # set newborns exactly like NetLogo's hatch [ ... ] block
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "color",              val = "white")
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "size",               val = 0.1)
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "energy_reserve",     val = of(agents = new_cocoons,var = "mass_cocoons") * energy_flesh)
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "mass",               val = -0.00152 + 0.794 * of(agents = new_cocoons,var = "mass_cocoons"))
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "energy_assimilated", val = 0)
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "age",                val = 0)
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "hatchlings",         val = 0)
   }
   
   # update parent inside the worms object (so parent and worms remain consistent)
   parent_agent <- NLwith(agents = focal, var = "who", val = who_reproduce)
   focal <- NLset(turtles = focal, agents = parent_agent, var = "hatchlings", val = of(agents = parent_agent,var = "hatchlings") + 1)
   focal <- NLset(turtles = focal, agents = parent_agent, var = "R",          val = of(agents = parent_agent,var = "R")  - (of(agents = parent_agent,var = "mass_cocoons") * (energy_flesh + energy_synthesis)))
   
   return(focal)
 }
 
  ###==============================
 ###   REPRODUCTION0
 ###==============================
 mate0 <- function(focal) {
   if (any(of(agents = focal,var = "mass") >= mass_sexual_maturity) &&
       any(duplicated(patchHere(world = food_density,turtles = focal)))) {
         
         # 🧭 1. Coordonnées et ID individuels
         worm_coords <- round(coordinates(focal))
         who <- of(agents = focal, var = "who")
         
         worms_df <- data.frame(
           x = worm_coords[, 1],
           y = worm_coords[, 2],
           who = who
         )
         
         # 🧮 2. Identifier les groupes de vers par patch (x, y)
         worms_groups <- worms_df %>%
           group_by(x, y) %>%
           mutate(
             count_worms = n()
           ) %>%
           ungroup()
         
         # 🎯 3. Ne garder que les groupes avec ≥ 2 vers
         paired_groups <- worms_groups %>%
           filter(count_worms >= 2)%>% group_by(x, y) %>% mutate(  group_id = cur_group_id())
         
         # 🚀 4. Pour chaque groupe, choisir un ver et le déplacer vers un autre ver du même groupe
         if (nrow(paired_groups) > 0) {
           for (gid in unique(paired_groups$group_id)) {
             group_members <- paired_groups[paired_groups$group_id == gid, ]
             
             if (nrow(group_members) >= 2) {
               # Sélectionner aléatoirement un ver "mâle" et une cible du même groupe
               mover <- sample(group_members$who, 1)
               if( (nrow(group_members) == 2)){
                 target <- setdiff(group_members$who, mover)  
               }else{
                 target <- sample(setdiff(group_members$who, mover),1)
               }
               
               
               # Mouvement : le ver se déplace vers la cible
               focal_mover <- moveTo(turtles = focal[which(focal[,"who"]==mover)], agents = focal[which(focal[,"who"]==target)])
                focal <- NLset(turtles = focal ,agents =focal_mover ,var = "xcor" ,val = of(agents = focal_mover,var = "xcor"))
                focal <- NLset(turtles = focal ,agents =focal_mover ,var = "ycor" ,val = of(agents = focal_mover,var = "ycor"))
               }
           }
         }
         
         # 5. Reproduction pour les individus matures
         focal <- calc_reproduction(focal)
       }
       
       return(focal)
     }
 
 calc_reproduction0 <- function(focal) {
   mass_vals<-of(agents = focal,var = "mass")
   Arrhenius_here<-of(agents = focal, var = "Arrhenius_here")
   #Arrhenius_here<-patchHere(world = Arrhenius,turtles = focal)
   focal <- NLset(turtles = focal, agents = focal,var = "max_R",val =  (max_reproduction_rate * Arrhenius_here) * mass_vals)
   energy_assimilated_vals <- of(agents = focal,var = "energy_assimilated")
   max_R_vals<-of(agents = focal,var = "max_R")
   R_vals<-of(agents = focal,var = "R")
   # Energy allocation to reproduction
   if (any(energy_assimilated_vals >= max_R_vals)) {
     ea_above_maxR<-(energy_assimilated_vals >= max_R_vals)
     who_above_max_R <- of(agents = focal,var = "who")[which(energy_assimilated_vals >= max_R_vals)]
     focal_above_max_R <- NLwith(agents = focal,var = "who",val = who_above_max_R)
  focal<- NLset(turtles= focal, agents= focal_above_max_R, var="energy_assimilated", val= energy_assimilated_vals[ea_above_maxR] - max_R_vals[ea_above_maxR])
  focal<- NLset(turtles= focal, agents= focal_above_max_R, var="R",val= R_vals[ea_above_maxR] + max_R_vals[ea_above_maxR])
   } else if (any(energy_assimilated_vals > 0)) {
     ea_above_zero<-(energy_assimilated_vals > 0)
     who_above_zero<-of(agents = focal,var = "who")[which(energy_assimilated_vals > 0)]
     focal_above_zero<- NLwith(agents = focal,var = "who",val = who_above_zero)
     focal<- NLset(turtles= focal,agents= focal_above_zero,var="R",val= R_vals[ea_above_zero] + max_R_vals[ea_above_zero])
    # focal<- NLset(turtles= focal,agents= focal_above_zero,var="energy_assimlated",val= rep(0,nrow(focal_above_zero)))     
   }
   
   if (any(energy_assimilated_vals > 0)) { 
     focal<-calc_growth(focal) 
     }
   
   mass_cocoon<- -0.01768 + 0.04035 * mass_vals
   focal<-NLset(turtles = focal, agents = focal, var = "mass_cocoons", val= mass_cocoon)
   if (any(R_vals >= mass_cocoon * (energy_flesh + energy_synthesis))) {
     who_above_R<-of(agents = focal,var = "who")[which(R_vals >= mass_cocoon * (energy_flesh + energy_synthesis))]
     focal_above_zero<- NLwith(agents = focal,var = "who",val = who_above_R)
     focal_above_zero<-reproduce(focal = focal_above_zero)
     focal<-upworms_with(group = focal_above_zero,focal = focal)
   }
   
   return(focal)
 }
 
  reproduce0 <- function(focal) {
   who_reproduce <- of(agents = focal,var = "who")
   ids_before <- of(agents = focal, var = "who")
   for(i in length(who_reproduce)){
   worms <<- hatch(turtles = worms, who = who_reproduce[i], n = 1, breed = "cocoons")
   cat("worms +1")
   }
   # detect new IDs (robust even if parent fields were copied)
   ids_after <- of(agents = focal, var = "who")
   new_ids <- setdiff(ids_after, who_reproduce)
   
   if (length(new_ids) > 0) {
     new_cocoons <- NLwith(agents = focal, var = "who", val = new_ids)
     
     # set newborns exactly like NetLogo's hatch [ ... ] block
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "color",              val = "white")
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "size",               val = 0.1)
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "energy_reserve",     val = of(agents = new_cocoons,var = "mass_cocoons") * energy_flesh)
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "mass",               val = -0.00152 + 0.794 * of(agents = new_cocoons,var = "mass_cocoons"))
     focal <- NLset(turtles = focal, agents = new_cocoons, var = "energy_assimilated", val = 0)
     worms <- NLset(turtles = worms, agents = new_cocoons, var = "age",                val = 0)
     worms <- NLset(turtles = worms, agents = new_cocoons, var = "hatchlings",         val = 0)
   }
   
   # update parent inside the worms object (so parent and worms remain consistent)
   parent_agent <- NLwith(agents = focal, var = "who", val = who_reproduce)
   focal <- NLset(turtles = focal, agents = parent_agent, var = "hatchlings", val = of(agents = parent_agent,var = "hatchlings") + 1)
   focal <- NLset(turtles = focal, agents = parent_agent, var = "R",          val = of(agents = parent_agent,var = "R")  - (of(agents = parent_agent,var = "mass_cocoons") * (energy_flesh + energy_synthesis)))
   
   return(focal)
 }

  mate1 <- function(focal) {
    # quick return if no agents
    if (NLcount(focal) == 0) return(focal)
    
    # need at least one mature worm and at least one patch with >1 worm
    if (!any(of(agents = focal, var = "mass") >= mass_sexual_maturity)) return(focal)
    if (!any(duplicated(patchHere(world = food_density, turtles = focal)))) return(focal)
    
    # 1. coordinates and who
    coords <- round(coordinates(focal))
    who_all <- of(agents = focal, var = "who")
    
    worms_df <- data.frame(
      x = coords[,1],
      y = coords[,2],
      who = who_all,
      stringsAsFactors = FALSE
    )
    
    # 2. group worms by patch and keep patches with >= 2 worms
    groups <- worms_df %>%
      group_by(x, y) %>%
      mutate(count_worms = n()) %>%
      ungroup()
    
    paired_groups <- groups %>% filter(count_worms >= 2)
    
    if (nrow(paired_groups) == 0) return(focal)
    
    paired_groups <- paired_groups %>% group_by(x,y) %>% mutate(group_id = cur_group_id()) %>% ungroup()
    
    # 3. For each group, choose a mover and a target and move the mover to the target
    for (gid in unique(paired_groups$group_id)) {
      members <- paired_groups %>% filter(group_id == gid)
      if (nrow(members) < 2) next
      
      mover_id <- sample(members$who, 1)
      if (nrow(members) == 2) {
        target_id <- setdiff(members$who, mover_id)
      } else {
        target_id <- sample(setdiff(members$who, mover_id), 1)
      }
      
      mover_agent  <- NLwith(agents = focal, var = "who", val = mover_id)
      target_agent <- NLwith(agents = focal, var = "who", val = target_id)
      
      # call your moveTo function (should accept valid agentsets)
      moved_agent <- tryCatch(
        moveTo(turtles = mover_agent, agents = target_agent),
        error = function(e) mover_agent
      )
      
      # write changes back into focal (and ensure upworms_with merges any changes)
      focal <- upworms_with(group = moved_agent, focal = focal)
    }
    
    # 4. After all movers moved, perform reproduction calculations
    focal <- calc_reproduction(focal)
    
    return(focal)
  }
  calc_reproduction1 <- function(focal) {
    if (NLcount(focal) == 0) return(focal)
    
    mass_vals <- of(agents = focal, var = "mass")
    Arrh_vals <- of(agents = focal, var = "Arrhenius_here")
    who_all   <- of(agents = focal, var = "who")
    
    # set max_R
    max_R_vals <- (max_reproduction_rate * Arrh_vals) * mass_vals
    focal <- NLset(turtles = focal, agents = focal, var = "max_R", val = max_R_vals)
    
    energy_assimilated_vals <- of(agents = focal, var = "energy_assimilated")
    R_vals <- of(agents = focal, var = "R")
    
    # Case: energy_assimilated >= max_R -> allocate max_R to R and subtract from assimilated
    idx_full <- which(energy_assimilated_vals >= max_R_vals)
    if (length(idx_full) > 0) {
      who_full <- who_all[idx_full]
      agents_full <- NLwith(agents = focal, var = "who", val = who_full)
      
      focal <- NLset(turtles = focal, agents = agents_full,
                     var = "energy_assimilated",
                     val = energy_assimilated_vals[idx_full] - max_R_vals[idx_full])
      
      focal <- NLset(turtles = focal, agents = agents_full,
                     var = "R",
                     val = R_vals[idx_full] + max_R_vals[idx_full])
    }
    
    # Case: some assimilation but < max_R -> allocate whatever assimilation exists to R
    idx_partial <- which(energy_assimilated_vals > 0 & energy_assimilated_vals < max_R_vals)
    if (length(idx_partial) > 0) {
      who_part <- who_all[idx_partial]
      agents_part <- NLwith(agents = focal, var = "who", val = who_part)
      
      focal <- NLset(turtles = focal, agents = agents_part,
                     var = "R",
                     val = R_vals[idx_partial] + energy_assimilated_vals[idx_partial])
      
      focal <- NLset(turtles = focal, agents = agents_part,
                     var = "energy_assimilated",
                     val = 0)
    }
    
    # if any energy was allocated, call growth
    if (any(energy_assimilated_vals > 0)) {
      focal <- calc_growth(focal)
    }
    
    # compute mass required to produce a cocoon (mass_cocoon formula you used)
    mass_cocoon <- (-0.01768 + 0.04035) * mass_vals
    focal <- NLset(turtles = focal, agents = focal, var = "mass_cocoons", val = mass_cocoon)
    
    # find who have enough R to reproduce
    R_vals <- of(agents = focal, var = "R")
    req_energy <- mass_cocoon * (energy_flesh + energy_synthesis)
    idx_repro <- which(R_vals >= req_energy)
    
    if (length(idx_repro) > 0) {
      who_repro <- who_all[idx_repro]
      repro_agents <- NLwith(agents = focal, var = "who", val = who_repro)
      
      # call reproduce for the selected agents
      focal <- reproduce(focal, repro_agents)
      
    }
    
    return(focal)
  }
  reproduce1 <- function(focal, repro_agents) {
    # 'focal' here should be the parents agentset (subset of worms) that will reproduce
    if (NLcount(focal) == 0) return(focal)
    
    parents_who <- of(agents = repro_agents, var = "who")
    
    # snapshot current global worm ids
    global_before <- of(agents = focal, var = "who")
    
    # hatch one cocoon per parent
    for (parent_id in parents_who) {
      # call hatch on global worms: this returns an updated global worms object
      focal <- hatch(turtles = focal, who = parent_id, n = 1, breed = "cocoons")
      cat(".")
    }
    #if(length(parents_who)>0) {cat(" +1 cocoons" ) }
    # find newly created ids by comparing global who lists
    global_after <- of(agents = focal, var = "who")
    new_ids <- setdiff(global_after, global_before)
    
    new_cocoons <- NULL
    if (length(new_ids) > 0) {
      # agentset of new cocoons in the GLOBAL worms object
      new_cocoons <- NLwith(agents = focal, var = "who", val = new_ids)
      
      # initialize newborns on the global worms object
      focal <<- NLset(turtles = focal, agents = new_cocoons, var = "color", val = "white")
      focal <<- NLset(turtles = focal, agents = new_cocoons, var = "size",  val = 0.1)
      
      # set mass_cocoons value on new cocoons from parent's mass_cocoons
      # For simplicity, assume each parent produced exactly one cocoon,
      # and order corresponds to the order of parents_who. We'll map by parent->new cocoon.
      # Build parent->new mapping robustly: find cocoons with age==0 and breed=="cocoons"
      candidates <- NLwith(agents = focal, var = "age", val = of(agents = focal, var = "age")[of(agents = focal, var = "age") == 0])
      candidates <- NLwith(agents = candidates, var = "breed", val = of(agents = candidates, var = "breed")[of(agents = candidates, var = "breed") == "cocoons"])
      # intersect with new_ids to be safe
      new_cocoons <- NLwith(agents = worms, var = "who", val = new_ids)
      
      # set energy_reserve and mass for new cocoons from their parent's mass_cocoons
      # We'll find parents' mass_cocoons and assign same value to their respective cocoon(s)
      # Get parents' mass_cocoons values (from focal)
      parents_mass_cocoon <- of(agents = focal, var = "mass_cocoons")
      # If number of new cocoons matches parents, assign in order
      if (length(parents_mass_cocoon) == length(new_ids)) {
        focal <- NLset(turtles = focal, agents = new_cocoons,
                       var = "energy_reserve",
                       val = parents_mass_cocoon * energy_flesh)
        focal <- NLset(turtles = focal, agents = new_cocoons,
                       var = "mass",
                       val = -0.00152 + 0.794 * (parents_mass_cocoon))
      } else {
        # fallback: use parent's mean mass_cocoon for all new cocoons
        mean_mc <- mean(parents_mass_cocoon, na.rm = TRUE)
        focal <- NLset(turtles = focal, agents = new_cocoons,
                       var = "energy_reserve",
                       val = mean_mc * energy_flesh)
        focal <- NLset(turtles = focal, agents = new_cocoons,
                       var = "mass",
                       val = -0.00152 + 0.794 * mean_mc)
      }
      
      # initialize other newborn fields on global worms
      focal <- NLset(turtles = focal, agents = new_cocoons, var = "energy_assimilated", val = 0)
      focal <- NLset(turtles = focal, agents = new_cocoons, var = "age", val = 0)
      focal <- NLset(turtles = focal, agents = new_cocoons, var = "hatchlings", val = 0)
      focal <- NLset(turtles = focal, agents = new_cocoons, var = "breed", val = "cocoons")
    }
    
    # Update parents (both in global worms and in returned focal)
    if (length(parents_who) > 0) {
      parent_agents_global <- NLwith(agents = focal, var = "who", val = parents_who)
      # increment hatchlings and deduct R from parent in global worms
      parents_mass_cocoon <- of(agents = focal, var = "mass_cocoons")
      deduction <- parents_mass_cocoon * (energy_flesh + energy_synthesis)
      
      # update global worms parents
      focal <- NLset(turtles = focal, agents = parent_agents_global,
                     var = "hatchlings",
                     val = of(agents = parent_agents_global, var = "hatchlings") + 1)
      
      focal <- NLset(turtles = focal, agents = parent_agents_global,
                     var = "R",
                     val = of(agents = parent_agents_global, var = "R") - deduction)
      
      # also update the focal agentset representing parents to keep it consistent
      parent_agents_focal <- NLwith(agents = focal, var = "who", val = parents_who)
      focal <- NLset(turtles = focal, agents = parent_agents_focal,
                     var = "hatchlings",
                     val = of(agents = parent_agents_global, var = "hatchlings"))
      focal <- NLset(turtles = focal, agents = parent_agents_focal,
                     var = "R",
                     val = of(agents = parent_agents_global, var = "R"))
    }
    
    # Return an agentset of the newly created cocoons so caller can merge if needed.
    return(focal)
  }
   
###==============================
###   GROWTH
###==============================

calc_growth0 <- function(focal) {
  # Extract values
  mass_vals <- of(agents = focal, var = "mass")
  Arrhenius_here <- of(agents = focal, var = "Arrhenius_here")
 
  # Calculate the maximum growth rate (von Bertalanffy type)
  max_growth_rate <- (growth_constant * Arrhenius_here) * 
    ((mass_maximum^(1/3) * mass_vals^(2/3)) - mass_vals)
  
  # Only individuals under their maximum mass can grow
  can_grow <- which(mass_vals < mass_maximum)
  if (length(can_grow) > 0) {
    focal <- grow(focal = focal[can_grow], max_growth_rate = max_growth_rate)
  }
  
  return(focal)
}
  
calc_growth <- function(focal) {
    # Extract values
    mass_vals <- of(agents = focal, var = "mass")
    Arrhenius_here <- of(agents = focal, var = "Arrhenius_here")
    
    # Calculate the maximum growth rate (von Bertalanffy type)
    max_growth_rate <- (growth_constant * Arrhenius_here) *
      ((mass_maximum^(1/3) * mass_vals^(2/3)) - mass_vals)
    
    # Filter individuals under their maximum mass
    can_grow <- which(mass_vals < mass_maximum)
    who_can_grow <- of(agents = focal[can_grow], var = "who")
      focal_can_grow <- NLwith(agents = focal,var = "who",val = who_can_grow)
      #max_growth_rate_who <- 
    if (length(can_grow) > 0) {
      focal_can_grow <- grow(
        focal = focal_can_grow,
        max_growth_rate = max_growth_rate[can_grow]
      )
    focal <- upworms_with(group = focal_can_grow,focal = focal)
        }
    
    return(focal)
  }
  
###==============================
###   GROW FUNCTION
###==============================
grow0 <- function(focal, max_growth_rate) {
  energy_assimilated_vals <- of(agents = focal, var = "energy_assimilated")
  mass_vals <- of(agents = focal, var = "mass")
  growth_rate_vals <- of(agents = focal, var = "growth_rate")
  # Energy cost of growth
  energy_growth <- max_growth_rate * (energy_flesh + energy_synthesis)
  
  for (i in seq_along(mass_vals)) {
   
      if (energy_assimilated_vals[i] >= energy_growth[i]) {
        mass_vals[i] <- mass_vals[i] + max_growth_rate[i]
        energy_assimilated_vals[i] <- energy_assimilated_vals[i] - energy_growth[i]
      } else {
        growth_rate_vals[i] <- (max_growth_rate[i] / energy_growth[i]) * energy_assimilated_vals[i]
        mass_vals[i] <- mass_vals[i] + growth_rate_vals[i]
        energy_assimilated_vals[i] <- 0
      }
    }
  
  focal <- NLset(turtles = focal, agents = focal, var = "mass", val = mass_vals)
  focal <- NLset(turtles = focal, agents = focal, var = "energy_assimilated", val = energy_assimilated_vals)
  focal <- NLset(turtles = focal, agents = focal, var = "mass", val = growth_rate_vals)
  
  return(focal)
}

grow <- function(focal, max_growth_rate) {
    
    energy_assimilated_vals <- of(agents = focal, var = "energy_assimilated")
    mass_vals <- of(agents = focal, var = "mass")
    growth_rate_vals <- of(agents = focal, var = "growth_rate")
    
    # Energy cost of growth
    energy_growth <- max_growth_rate * (energy_flesh + energy_synthesis)
    
    for (i in seq_along(mass_vals)) {
      
      if (energy_assimilated_vals[i] >= energy_growth[i]) {
        # Full growth possible
        mass_vals[i] <- mass_vals[i] + max_growth_rate[i]
        energy_assimilated_vals[i] <- energy_assimilated_vals[i] - energy_growth[i]
        growth_rate_vals[i] <- max_growth_rate[i]
        
      } else {
        # Partial growth
        growth_fraction <- energy_assimilated_vals[i] / energy_growth[i]
        
        achieved_growth <- max_growth_rate[i] * growth_fraction
        
        mass_vals[i] <- mass_vals[i] + achieved_growth
        growth_rate_vals[i] <- achieved_growth
        
        energy_assimilated_vals[i] <- 0
      }
    }
    
    # Write changes back to turtles
    focal <- NLset(turtles = focal, agents = focal, var = "mass", val = mass_vals)
    focal <- NLset(turtles = focal, agents = focal, var = "energy_assimilated", val = energy_assimilated_vals)
    focal <- NLset(turtles = focal, agents = focal, var = "growth_rate", val = growth_rate_vals)
   
    return(focal)
  }
  
###==============================
###   PATCH UPDATES
###==============================
update_patches <- function(food_density) { ### TAKE CARE IF < 0 no changes because I suppose it will be 0 to check 
  # Update food density (after consumption)
  calc_change_food_density(worms, food_density)
  food_vals <- of(world = food_density, agents = patches(food_density))
  change_food_vals <- of(world = change_food_density, agents = patches(change_food_density))
  
  # 3️⃣ Update food: food_density = max(0, food_density - change_food_density)
  new_food_vals <- pmax(0, food_vals - change_food_vals)
  
  # 4️⃣ Write back updated food values
  food_density <- NLset(world = food_density, agents = patches(food_density), val = new_food_vals)
  
  return(food_density)
}

  ###==============================
###   TURTLE UPDATES
###==============================
update_turtles0 <- function(worms) {
  # Clamp coords inside world boundaries
  worms@.Data[,"xcor"] <- pmin(pmax(worms@.Data[,"xcor"], min_Pxcor_move), max_Pxcor_move)
  worms@.Data[,"ycor"] <- pmin(pmax(worms@.Data[,"ycor"], min_Pycor_move), max_Pycor_move)
  
   # Extract key variables
  breed_vals              <- of(agents = worms, var = "breed")
  energy_assimilated_vals <- of(agents = worms, var = "energy_assimilated")
  energy_reserve_vals     <- of(agents = worms, var = "energy_reserve")
  mass_vals               <- of(agents = worms, var = "mass")
  aestivating_vals        <- of(agents = worms, var = "aestivating")
  time_aestivating_vals   <- of(agents = worms, var = "time_aestivating")
  
  # 1️⃣ Energy reserve capacity for adults and juveniles
  adults_or_juv <- which(breed_vals %in% c("adults", "juveniles"))
  if (length(adults_or_juv) > 0) {
    energy_reserve_max_vals <- (mass_vals[adults_or_juv] / 2) * energy_flesh
    worms <- NLset(turtles = worms,
                   agents = worms[adults_or_juv],
                   var = "energy_reserve_max",
                   val = energy_reserve_max_vals)
  }
  
  # 2️⃣ Handle aestivation (pause activity)
  aestivating_idx <- which(aestivating_vals == TRUE)
  if (length(aestivating_idx) > 0) {
    # Increment aestivation time
    time_aestivating_vals[aestivating_idx] <- time_aestivating_vals[aestivating_idx] + 1
    # Hide aestivating worms (optional, depending on visualization)
    worms <- NLset(turtles = worms,
                   agents = worms[aestivating_idx],
                   var = "visible",
                   val = rep(FALSE,length(aestivating_idx)))
  }
  
  # 3️⃣ Assimilated energy → energy reserves
  energy_reserve_vals <- energy_reserve_vals + (
    energy_assimilated_vals * (energy_flesh / (energy_flesh + energy_synthesis))
  )
  
  # 4️⃣ Cap reserves at max, and convert excess energy to body mass
  energy_reserve_max_vals <- of(agents = worms, var = "energy_reserve_max")
  over_max <- (energy_reserve_vals > energy_reserve_max_vals)
  if (length(over_max) > 0) {
    energy_reserve_vals[over_max] <- energy_reserve_max_vals[over_max]
    mass_vals[over_max] <- mass_vals[over_max] + (energy_reserve_vals[over_max] / 10.6)
  }
  
  # 6️⃣ Write updates back to worms
  worms <- NLset(turtles = worms, agents = worms, var = "energy_reserve", val = energy_reserve_vals)
  worms <- NLset(turtles = worms, agents = worms, var = "mass", val = mass_vals)
  worms <- NLset(turtles = worms, agents = worms, var = "time_aestivating", val = time_aestivating_vals)
  
  return(worms)
}
    update_turtles1 <- function(worms) {
    # Clamp coords inside world boundaries
    worms@.Data[,"xcor"] <- pmin(pmax(worms@.Data[,"xcor"], min_Pxcor_move), max_Pxcor_move)
    worms@.Data[,"ycor"] <- pmin(pmax(worms@.Data[,"ycor"], min_Pycor_move), max_Pycor_move)
    
    # Extract key variables
    breed_vals              <- of(agents = worms, var = "breed")
    energy_assimilated_vals <- of(agents = worms, var = "energy_assimilated")
    energy_reserve_vals     <- of(agents = worms, var = "energy_reserve")
    mass_vals               <- of(agents = worms, var = "mass")
    aestivating_vals        <- of(agents = worms, var = "aestivating")
    time_aestivating_vals   <- of(agents = worms, var = "time_aestivating")
    
    # 1️⃣ Energy reserve capacity for adults and juveniles
    adults_or_juv <- which(breed_vals %in% c("adults", "juveniles"))
    if (length(adults_or_juv) > 0) {
      energy_reserve_max_vals <- (mass_vals[adults_or_juv] / 2) * energy_flesh
      worms <- NLset(
        turtles = worms,
        agents = worms[adults_or_juv],
        var = "energy_reserve_max",
        val = energy_reserve_max_vals
      )
    }
    
    # 2️⃣ Handle aestivation (pause activity)
    aestivating_idx <- which(aestivating_vals)
    if (length(aestivating_idx) > 0) {
      # Increment aestivation time
      time_aestivating_vals[aestivating_idx] <- time_aestivating_vals[aestivating_idx] + 1
      # Optional: hide aestivating worms
      worms <- NLset(
        turtles = worms,
        agents = worms[aestivating_idx],
        var = "visible",
        val = rep(FALSE, length(aestivating_idx))
      )
    }
    
    # 3️⃣ Assimilated energy → energy reserves
    energy_reserve_vals <- energy_reserve_vals + (
      energy_assimilated_vals * (energy_flesh / (energy_flesh + energy_synthesis))
    )
    
    # 4️⃣ Cap reserves at max, convert excess energy to body mass
    energy_reserve_max_vals <- of(agents = worms, var = "energy_reserve_max")
    over_max <- energy_reserve_vals > energy_reserve_max_vals
    if (any(over_max)) {
      excess_energy <- energy_reserve_vals[over_max] - energy_reserve_max_vals[over_max]
      mass_vals[over_max] <- mass_vals[over_max] + (excess_energy / 10.6)
      energy_reserve_vals[over_max] <- energy_reserve_max_vals[over_max]
    }
    
    # 5️⃣ Write updates back to worms
    worms <- NLset(turtles = worms, agents = worms, var = "energy_reserve", val = energy_reserve_vals)
    worms <- NLset(turtles = worms, agents = worms, var = "mass", val = mass_vals)
    worms <- NLset(turtles = worms, agents = worms, var = "time_aestivating", val = time_aestivating_vals)
    
    return(worms)
  }
  update_turtles2 <- function(worms) {
    # -------------------------------
    # 1️⃣ Clamp coordinates inside world boundaries
    # -------------------------------
    worms@.Data[,"xcor"] <- pmin(pmax(worms@.Data[,"xcor"], min_Pxcor_move), max_Pxcor_move)
    worms@.Data[,"ycor"] <- pmin(pmax(worms@.Data[,"ycor"], min_Pycor_move), max_Pycor_move)
    
    # -------------------------------
    # 2️⃣ Extract key variables
    # -------------------------------
    breed_vals              <- of(agents = worms, var = "breed")
    energy_assimilated_vals <- of(agents = worms, var = "energy_assimilated")
    energy_reserve_vals     <- of(agents = worms, var = "energy_reserve")
    mass_vals               <- of(agents = worms, var = "mass")
    aestivating_vals        <- of(agents = worms, var = "aestivating")
    time_aestivating_vals   <- of(agents = worms, var = "time_aestivating")
    energy_reserve_max_vals <- of(agents = worms, var = "energy_reserve_max")
    
    # -------------------------------
    # 3️⃣ Energy reserve for adults & juveniles
    # -------------------------------
    adults_or_juv <- which(breed_vals %in% c("adults", "juveniles"))
    if (length(adults_or_juv) > 0) {
      energy_reserve_max_vals[adults_or_juv] <- (mass_vals[adults_or_juv] / 2) * energy_flesh
    }
    
    # -------------------------------
    # 4️⃣ Handle aestivation
    # -------------------------------
    aestivating_idx <- which(aestivating_vals == TRUE)
    if (length(aestivating_idx) > 0) {
      # Increment aestivation time
      time_aestivating_vals[aestivating_idx] <- time_aestivating_vals[aestivating_idx] + 1
      # Optionally hide aestivating worms
      worms <- NLset(turtles = worms,
                     agents = worms[aestivating_idx],
                     var = "visible",
                     val = rep(FALSE, length(aestivating_idx)))
    }
    
    # -------------------------------
    # 5️⃣ Assimilated energy → energy reserves
    # -------------------------------
    energy_reserve_vals <- energy_reserve_vals + (
      energy_assimilated_vals * (energy_flesh / (energy_flesh + energy_synthesis))
    )
    
    # -------------------------------
    # 6️⃣ Cap reserves at max and convert excess energy to mass
    # -------------------------------
    over_max <- which(energy_reserve_vals > energy_reserve_max_vals)
    if (length(over_max) > 0) {
      excess_energy <- energy_reserve_vals[over_max] - energy_reserve_max_vals[over_max]
      energy_reserve_vals[over_max] <- energy_reserve_max_vals[over_max]
      mass_vals[over_max] <- mass_vals[over_max] + excess_energy / 10.6
    }
    
    # -------------------------------
    # 7️⃣ Handle worm mortality (mass <= 0)
    # -------------------------------
    dead_idx <- which(mass_vals <= 0)
    if (length(dead_idx) > 0) {
      cat("⚠️ Removing", length(dead_idx), "worms due to zero/negative mass.\n")
      who_die <- of(agents = worms[dead_idx], var = "who")
      worms <<- die(turtles = worms, who = who_die)  # Remove from global worms object
      # Remove dead worms from local vectors
      energy_reserve_vals <- energy_reserve_vals[-dead_idx]
      mass_vals <- mass_vals[-dead_idx]
      time_aestivating_vals <- time_aestivating_vals[-dead_idx]
      energy_reserve_max_vals <- energy_reserve_max_vals[-dead_idx]
    }
    
    # -------------------------------
    # 8️⃣ Write updates back to worms
    # -------------------------------
    worms <- NLset(turtles = worms, agents = worms, var = "energy_reserve", val = energy_reserve_vals)
    worms <- NLset(turtles = worms, agents = worms, var = "mass", val = mass_vals)
    worms <- NLset(turtles = worms, agents = worms, var = "time_aestivating", val = time_aestivating_vals)
    worms <- NLset(turtles = worms, agents = worms, var = "energy_reserve_max", val = energy_reserve_max_vals)
    
    return(worms)
  }
  
  update_turtles <- function(worms) {
    
    # -------------------------------
    # 1️⃣ Clamp coordinates inside boundaries
    # -------------------------------
    worms@.Data[,"xcor"] <- pmin(pmax(worms@.Data[,"xcor"], min_Pxcor_move), max_Pxcor_move)
    worms@.Data[,"ycor"] <- pmin(pmax(worms@.Data[,"ycor"], min_Pycor_move), max_Pycor_move)
    
    # -------------------------------
    # 2️⃣ Extract variables
    # -------------------------------
    breed_vals              <- of(agents = worms,var =  "breed")
    energy_assimilated_vals <- of(agents = worms,var =  "energy_assimilated")
    energy_reserve_vals     <- of(agents = worms,var =  "energy_reserve")
    energy_reserve_max_vals <- of(agents = worms,var =  "energy_reserve_max")
    mass_vals               <- of(agents = worms,var =  "mass")
    aestivating_vals        <- of(agents = worms,var =  "aestivating")
    time_aestivating_vals   <- of(agents = worms,var =  "time_aestivating")
    #experienced_conc_vals   <- of(agents = worms,var =  "experienced_conc")
    #patch_conc_vals         <- of(agents = worms,var =  "patch_conc")
    
    # -------------------------------
    # 3️⃣ Adults & juveniles: reserve max = 50% mass * energy_flesh
    # -------------------------------
    adults_or_juv <- which(breed_vals %in% c("adults", "juveniles"))
    if (length(adults_or_juv) > 0) {
      energy_reserve_max_vals[adults_or_juv] <- (mass_vals[adults_or_juv] / 2) * energy_flesh
    }
    
    # -------------------------------
    # 4️⃣ Handle aestivation
    # -------------------------------
    aestivating_idx <- which(aestivating_vals == TRUE)
    if (length(aestivating_idx) > 0) {
      time_aestivating_vals[aestivating_idx] <- time_aestivating_vals[aestivating_idx] + 1
      worms <- NLset(
        turtles = worms,
        agents  = worms[aestivating_idx],
        var     = "visible",
        val     = rep(FALSE, length(aestivating_idx))
      )
    }
    
    # -------------------------------
    # 5️⃣ energy_assimilated → energy_reserve  
    #    (with synthesis cost included)
    # -------------------------------
    add_to_reserve <- energy_assimilated_vals * (energy_flesh / (energy_flesh + energy_synthesis))
    energy_reserve_vals <- energy_reserve_vals + add_to_reserve
    
    # -------------------------------
    # 6️⃣ Cap reserves at max & convert excess to mass
    #    (adults + juveniles only)
    # -------------------------------
    over_max <- which(energy_reserve_vals > energy_reserve_max_vals &
                        breed_vals %in% c("adults", "juveniles"))
    
    if (length(over_max) > 0) {
      excess <- energy_reserve_vals[over_max] - energy_reserve_max_vals[over_max]
      energy_reserve_vals[over_max] <- energy_reserve_max_vals[over_max]
      mass_vals[over_max] <- mass_vals[over_max] + excess / 10.6
    }
    
    # -------------------------------
    # 7️⃣ Update experienced concentration  
    #    (experienced_conc + patch_conc) / 2
    # -------------------------------
  #  experienced_conc_vals <- (experienced_conc_vals + patch_conc_vals) / 2
    
    # -------------------------------
    # 8️⃣ Write updates
    # -------------------------------
    worms <- NLset(turtles = worms,agents =  worms, var = "energy_reserve",      val =  energy_reserve_vals)
    worms <- NLset(turtles = worms,agents =  worms, var = "energy_reserve_max",  val =  energy_reserve_max_vals)
    worms <- NLset(turtles = worms,agents =  worms, var = "mass",                val =  mass_vals)
    worms <- NLset(turtles = worms,agents =  worms, var = "time_aestivating",    val =  time_aestivating_vals)
   # worms <- NLset(turtles = worms,agents =  worms, var = "experienced_conc",    val =  experienced_conc_vals)
    
    return(worms)
  }
  
  
  
upworms_with<-function (group, focal){
  who_ids <- of(agents = group, var = "who")
  group_who <- NLwith(agents = focal,var="who",val=who_ids)
  focal <- NLset(turtles = focal, agents = group_who, var = "xcor", val = of(agents = group, var = "xcor"))
  focal <- NLset(turtles = focal, agents = group_who, var = "ycor", val = of(agents = group, var = "ycor"))
  focal <- NLset(turtles = focal, agents = group_who, var = "who", val = of(agents = group, var = "who"))
  focal <- NLset(turtles = focal, agents = group_who, var = "heading", val = of(agents = group, var = "heading"))
  focal <- NLset(turtles = focal, agents = group_who, var = "prevX", val = of(agents = group, var = "prevX"))
  focal <- NLset(turtles = focal, agents = group_who, var = "prevY", val = of(agents = group, var = "prevY"))
  focal <- NLset(turtles = focal, agents = group_who, var = "breed", val = of(agents = group, var = "breed"))
  focal <- NLset(turtles = focal, agents = group_who, var = "color", val = of(agents = group, var = "color"))
  focal <- NLset(turtles = focal, agents = group_who, var = "age", val = of(agents = group, var = "age"))
  focal <- NLset(turtles = focal, agents = group_who, var = "mass", val = of(agents = group, var = "mass"))
  focal <- NLset(turtles = focal, agents = group_who, var = "energy_reserve", val = of(agents = group, var = "energy_reserve"))
  focal <- NLset(turtles = focal, agents = group_who, var = "size", val = of(agents = group, var = "size"))
  focal <- NLset(turtles = focal, agents = group_who, var = "energy_assimilated", val = of(agents = group, var = "energy_assimilated"))
  focal <- NLset(turtles = focal, agents = group_who, var = "energy_reserve_max", val = of(agents = group, var = "energy_reserve_max"))
  focal <- NLset(turtles = focal, agents = group_who, var = "ingestion_rate", val = of(agents = group, var = "ingestion_rate"))
  focal <- NLset(turtles = focal, agents = group_who, var = "BMR", val = of(agents = group, var = "BMR"))
  focal <- NLset(turtles = focal, agents = group_who, var = "mortality", val = of(agents = group, var = "mortality"))
  focal <- NLset(turtles = focal, agents = group_who, var = "max_growth_rate", val = of(agents = group, var = "max_growth_rate"))
  focal <- NLset(turtles = focal, agents = group_who, var = "growth_rate", val = of(agents = group, var = "growth_rate"))
  focal <- NLset(turtles = focal, agents = group_who, var = "energy_growth", val = of(agents = group, var = "energy_growth"))
  focal <- NLset(turtles = focal, agents = group_who, var = "embryonic_development", val = of(agents = group, var = "embryonic_development"))
  focal <- NLset(turtles = focal, agents = group_who, var = "hatchlings", val = of(agents = group, var = "hatchlings"))
  focal <- NLset(turtles = focal, agents = group_who, var = "max_R", val = of(agents = group, var = "max_R"))
  focal <- NLset(turtles = focal, agents = group_who, var = "R", val = of(agents = group, var = "R"))
  focal <- NLset(turtles = focal, agents = group_who, var = "aestivating", val = of(agents = group, var = "aestivating"))
  focal <- NLset(turtles = focal, agents = group_who, var = "time_aestivating", val = of(agents = group, var = "time_aestivating"))
  focal <- NLset(turtles = focal, agents = group_who, var = "mass_cocoons", val = of(agents = group, var = "mass_cocoons"))
  focal <- NLset(turtles = focal, agents = group_who, var = "Arrhenius_here", val = of(agents = group, var = "Arrhenius_here"))
  focal <- NLset(turtles = focal, agents = group_who, var = "visible", val = of(agents = group, var = "visible"))
 
     return(focal) 
}

# plot_world: simple ggplot visualization wrapper (returns ggplot)
plot_world <- function(food_density, burrow_no, worms, hour, day, year, dark) {
  pd <- data.frame(x = habitat_type@pCoords[,1], y = habitat_type@pCoords[,2], habitat = habitat_type@.Data)
  worms_df <- as.data.frame(worms@.Data)
  p <- ggplot(pd, aes(x = x, y = y)) + geom_raster(aes(fill = habitat)) +
    geom_point(data = worms_df, aes(x = xcor, y = ycor), colour = "red", size = 2) +
    ggtitle(sprintf("hour=%02d day=%03d year=%d dark=%s", hour, day, year, dark)) + coord_equal()
  p
}


###########################################################
#              PESTICIDE FLUCTUATIONS                     #
###########################################################
# Function to setup pesticide fluctuations  
setup_pesticide_fluctuations <- function(pesticide, base_swp) {
  pycor_vals <- base_swp@pCoords[,2]
  
  if(day == 1) {
    # Set baseline pesticide values using SWP as template, plus randomness
    pesticide_vals <- numeric(nrow(base_swp@.Data))
    
    # 7 vertical bands (like temperature/SWP)
    bands <- list(
      which(pycor_vals >= 89),
      which(pycor_vals >= 79 & pycor_vals < 89),
      which(pycor_vals >= 69 & pycor_vals < 79),
      which(pycor_vals >= 59 & pycor_vals < 69),
      which(pycor_vals >= 49 & pycor_vals < 59),
      which(pycor_vals >= 39 & pycor_vals < 49),
      which(pycor_vals < 39)
    )
    
    for(i in seq_along(bands)) {
      idx <- bands[[i]]
      if(length(idx) > 0) {
        # get patch coordinates as a 2-column matrix
        patch_coords <- patches(base_swp)[idx, , drop = FALSE]
        mean_conc <- mean(of(world = base_swp, agents = patch_coords))
        pesticide_vals[idx] <- rnorm(length(idx), mean = mean_conc, sd = 0.2 * mean_conc)
      }
    }
    
    pesticide <- NLset(world = pesticide, agents = patches(pesticide), val = pesticide_vals)
    pesticide_ref <<- pesticide   # global reference for later fluctuations
  } else {
    # Small random fluctuation around reference
    patch_coords <- patches(pesticide_ref)
    pesticide <- NLset(world = pesticide, agents = patch_coords,
                       val = of(world = pesticide_ref, agents = patch_coords) +
                         rnorm(length(patch_coords[,1]), mean = 0, 
                               sd = 0.05 * of(world = pesticide_ref, agents = patch_coords)))
  }
  
  return(pesticide)
}

setup_pesticide_fluctuations_hotspot <- function(pesticide,
                                                 base_swp,
                                                 base_spacing,        # distance entre les ceps (hotspots)
                                                 hotspot_strength,       # facteur d’intensité
                                                 hotspot_radius) {      # rayon de diffusion (en patches)
  # hotspot_strength : facteur multiplicatif au centre des hotspots
  # hotspot_radius : rayon d'influence d’un hotspot (en unités de patchs)
  # base_spacing : espacement entre les ceps (en m ou en unités de patchs)
  
  pycor_vals <- base_swp@pCoords[, 2]
  pxcor_vals <- base_swp@pCoords[, 1]
  
  if (day == 1) {
    
    # --- 1️⃣ Baseline pesticide values (comme avant) ---
    pesticide_vals <- numeric(nrow(base_swp@.Data))
    
    bands <- list(
      which(pycor_vals >= 89),
      which(pycor_vals >= 79 & pycor_vals < 89),
      which(pycor_vals >= 69 & pycor_vals < 79),
      which(pycor_vals >= 59 & pycor_vals < 69),
      which(pycor_vals >= 49 & pycor_vals < 59),
      which(pycor_vals >= 39 & pycor_vals < 49),
      which(pycor_vals < 39)
    )
    
    for (i in seq_along(bands)) {
      idx <- bands[[i]]
      if (length(idx) > 0) {
        patch_coords <- patches(base_swp)[idx, , drop = FALSE]
        mean_conc <- mean(of(world = base_swp, agents = patch_coords))
        pesticide_vals[idx] <- rnorm(length(idx),
                                     mean = mean_conc,
                                     sd = 0.2 * mean_conc)
      }
    }
    
    # --- 2️⃣ Créer plusieurs hotspots alignés tous les 2.5 m en haut de la carte ---
    max_x <- max(pxcor_vals)
    min_x <- min(pxcor_vals)
    max_y <- max(pycor_vals)
    
    hotspot_x <- seq(from = min_x, to = max_x, by = base_spacing)
    hotspot_y <- rep(max_y, length(hotspot_x))
    
    # --- 3️⃣ Calculer la contribution de chaque hotspot ---
    # Initialiser le facteur total d'influence à 1 (pas d'effet)
    hotspot_factor_total <- rep(1, length(pxcor_vals))
    
    for (i in seq_along(hotspot_x)) {
      dist_to_hotspot <- sqrt((pxcor_vals - hotspot_x[i])^2 + (pycor_vals - hotspot_y[i])^2)
      hotspot_factor_i <- exp(- (dist_to_hotspot^2) / (2 * hotspot_radius^2))
      hotspot_factor_total <- hotspot_factor_total + hotspot_strength * hotspot_factor_i
    }
    
    # --- 4️⃣ Appliquer les effets cumulés ---
    pesticide_vals <- pesticide_vals * hotspot_factor_total
    
    # --- 5️⃣ Mettre à jour l’objet NetLogoR ---
    pesticide <- NLset(world = pesticide,
                       agents = patches(pesticide),
                       val = pesticide_vals)
    pesticide_ref <<- pesticide
    
  } else {
    # --- Fluctuations légères après l’heure 1 ---
    patch_coords <- patches(pesticide_ref)
    pesticide <- NLset(world = pesticide,
                       agents = patch_coords,
                       val = of(world = pesticide_ref, agents = patch_coords) +
                         rnorm(length(patch_coords[, 1]), mean = 0,
                               sd = 0.05 * of(world = pesticide_ref, agents = patch_coords)))
  }
  
  return(pesticide)
}



########################## PESTICIDE #################################
calc_PPP_exposure <- function(focal, pesticide_layer) {
  who_vec <- of(agents = focal, var = "who")
  
  # Clamp coordinates inside world boundaries
  focal@.Data[,"xcor"] <- pmin(pmax(focal@.Data[,"xcor"], min_Pxcor), max_Pxcor)
  focal@.Data[,"ycor"] <- pmin(pmax(focal@.Data[,"ycor"], min_Pycor), max_Pycor)
  
  for (who in who_vec) {
    worm_i <- turtle(focal, who)
    
    # Patch where this worm is located
    current_patch <- patchHere(world = pesticide_layer, turtles = worm_i)
    
    # Pesticide concentration at this patch
    pesticide_here <- of(world = pesticide_layer, agents = current_patch)
    
    # Increment worm's turtles-own variable instead of overwriting
    current_value <- of(agents = worm_i, var = "PPP_external_concentration")
    new_value <- current_value + pesticide_here
    
    focal <- NLset(turtles = focal, agents = worm_i, var = "PPP_external_concentration", val = new_value)
  }
  
  return(focal)
}

setup_env_fluctuations_EEE <- function(world,day,env_data=Rotham) {
  # hour: current hour of simulation (integer)
  # temperature, SWP, energy_content_food: existing worldMatrix objects (patch variables)
  # env_data: data.frame with columns temperature, swp, energy_content_food (22 rows)
  
  # Extract y-coords (pycor) for each patch
  pycor_vals <- temperature@pCoords[,2]
  
  if(day == 1) {
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
    
    SWP_EEE <- SWP_ref
    temperature_EEE <- temperature_ref
    energy_content_food_EEE <- energy_content_food_ref
  }else{
    SWP_EEE <- SWP_ref + rnorm(length(SWP_ref), mean = 0, sd = 0.02)
    temperature_EEE <- temperature_ref + rnorm(length(temperature_ref), mean = 0, sd = 0.5)
    energy_content_food_EEE <- energy_content_food_ref + rnorm(length(energy_content_food_ref), mean = 0, sd = energy_content_food_ref/10)
  }
  # Return updated worldMatrices as a list
  return(SWP_EEE)
}

######################################################################
check_worms <- function(worms, step = "unknown") {
  
  if (length(worms) == 0) {
    warning(paste("No worms left after:", step))
    return(FALSE)
  }
  
  masses <- of(agents = worms, var = "mass")
  energy <- of(agents = worms, var = "energy_reserve")
  breed  <- of(agents = worms, var = "breed")
  
  # ---- BASIC CHECKS ----
  if (any(is.na(masses))) warning(paste(step,"-> NA mass detected"))
  if (any(is.infinite(masses))) warning(paste(step,"-> Inf mass detected"))
  
  if (any(masses < 0)){
    warning(paste(step,"-> NEGATIVE mass detected!"))
  }
  
  if (any(energy < 0)){
    warning(paste(step,"-> NEGATIVE energy detected!"))
  }
  
  # ---- STAGE CONSISTENCY ----
  if (any(breed == "adults" & masses < mass_sexual_maturity)) {
    warning(paste(step,"-> adults below maturity mass"))
  }
  
  if (any(breed == "juveniles" & masses >= mass_sexual_maturity)) {
    warning(paste(step,"-> juveniles above maturity mass"))
  }
  
  # ---- POSITION CHECK ----
  if (any(worms@.Data[,"xcor"] < min_Pxcor_move |
          worms@.Data[,"xcor"] > max_Pxcor_move)) {
    warning(paste(step,"-> X out of bounds"))
  }
  
  if (any(worms@.Data[,"ycor"] < min_Pycor_move |
          worms@.Data[,"ycor"] > max_Pycor_move)) {
    warning(paste(step,"-> Y out of bounds"))
  }
  
  cat("✔", step, "- n =",
      NLcount(worms),
      " | min mass = ", round(min(masses),3),
      " | max mass = ", round(max(masses),3), "\n")
  
  return(TRUE)
}

#######################################

aestivate00 <- function(focal=worms_juv) {
  
  current_SWP <-cbind(of(agents = focal,var = "who"),patchHere(world = SWP,turtles = focal))
  current_SWP_vals<- of(world = SWP,agents = patchHere(world = SWP,turtles = focal))
  current_SWP<-data.frame(cbind(current_SWP,current_SWP_vals))
  
  # Extraire états internes
  aestivating_vals     <- of(agents = focal, var = "aestivating")
  time_aestivating_vals <- of(agents = focal, var = "time_aestivating")
  
  # Condition 1: SWP >= 25
  cond_high_SWP <- current_SWP_vals >= 25
  # Condition 1.1 and 1.2: aestivate T or F
  cond_aestivate <- aestivating_vals==TRUE & (time_aestivating_vals <= 60)
  cond_No_aestivate <- aestivating_vals==FALSE 
  # COndititon 2: aestivate = T & SWP <= 20 & time_aestivate >=14
  cond_awake <- aestivating_vals==TRUE & (current_SWP_vals <= 20)  &(time_aestivating_vals >= 14)
  
  if (any(cond_high_SWP,na.rm = TRUE)) {
    if (any(cond_aestivate)){
      focal <- NLset(turtles = focal,  agents = focal[cond_aestivate],var = "aestivating",val = TRUE)
      focal <- calc_metabolic_rate(focal = focal[cond_aestivate])  # à adapter selon ton modèle
    }
    if (any(cond_No_aestivate)){
      focal <- NLset(turtles = focal,  agents = focal[cond_No_aestivate],var = "aestivating",val = TRUE)
      focal <- calc_metabolic_rate(focal[cond_No_aestivate])  # à adapter selon ton modèle
    }
  }else if (any(cond_awake)){
    focal <- NLset(turtles = focal, agents = focal[cond_awake],var = "aestivating",val = FALSE)
    focal <- NLset(turtles = focal, agents = focal[cond_awake],var = "time_aestivating",val = 1)
    focal <- NLset(turtles = focal,
                   agents = focal[cond_awake],
                   var = "visible",
                   val = TRUE)        
  }
  return(focal)
}