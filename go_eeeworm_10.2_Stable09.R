
# EEEWORM model - Main function 'go' - Version 10 #STABLE
# Author: Quentin DEVALLOIR
# Contact: @QUENTIN2VA
# Version: 10.2
# Note: Please contact me before using it, this is still a trial version.

source(file ='submodels_eeeworm_7.4_STABLE09.R')
source(file ='init_eeeworm_2_STABLE.R')

# test file Rotham
Rotham <- read.table(file = "YOUR FILE", header = FALSE)
Rotham
# define globally, once
hour <<- 0
day  <<- 1
year <<- 1

dark <<- ifelse(hour <= 11, FALSE, TRUE)


go <- function() {
  if (nrow(worms) == 0) return(invisible())
  
  # --- Time management ---
  hour <<- hour + 1
  dark <<- ifelse(hour <= 11, FALSE, TRUE)
  
  # --- Environment fluctuations ---
  setup_env_fluctuations(world = w,day = day, env_data = Rotham)
  
  pycor_vals <- w@pCoords[,2]
  food_density_vals <- ifelse(
    pycor_vals >= 80, rnorm(length(pycor_vals), mean = 3, sd = 0.30),
    ifelse(pycor_vals >= 40 & pycor_vals < 80, rnorm(length(pycor_vals), mean = 2.5, sd = 0.25),
           rnorm(length(pycor_vals), mean = 2, sd = 0.2))
  )
  food_density <<- NLset(food_density, agents = patches(food_density), val = food_density_vals)
  
  # Reset burrows where habitat-type != "burrow"
  new_burrow_no <- ifelse(habitat_type@.Data != "burrow", 0, burrow_no@.Data)
  burrow_no@.Data <- new_burrow_no@.Data
  
  # Temperature & patch quality
  T_vals <<- 273.15 + temperature@.Data
  setup_patch_quality()
  Arrhenius <<- exp((-activation_energy / Boltz) * ((1 / T_vals) - (1 / reference_T)))
  
  # --- Turtle Initialization (hour 1) ---
  if (hour == 1 && day == 1) {
    burrow_patches <- NLwith(agents = patches(habitat_type), world = habitat_type, val = "burrow")
    rand_ids <- sample(1:nrow(burrow_patches), size = nrow(worms), replace = TRUE)
    random_burrows <- burrow_patches[rand_ids, ]
    worms <<- moveTo(turtles = worms, agents = random_burrows)
    my_burrow_vals <- burrow_no@.Data[cbind(random_burrows[,2], random_burrows[,1])]
    worms <<- NLset(turtles = worms, agents = worms, var = "my_burrow", val = my_burrow_vals)
  }
  
  # --- Neighbors check ---
  neighbors_set <- inRadius(agents = worms, radius = 5,
                            agents2 = patches(habitat_type), world = habitat_type)
  neigh_types <- of(world = habitat_type, agents = neighbors_set[, c("pxcor","pycor")])
  burrow_flag <- tapply(neigh_types == "burrow", INDEX = neighbors_set[, "id"], FUN = any)
  worms <<- NLset(turtles = worms, agents = worms, var = "burrowQ",
                  val = as.logical(burrow_flag[as.character(worms@.Data[, "who"]+1)]))
  
  # If not in burrow → my_burrow = 0
  no_burrow_idx <- which(worms[, "burrowQ"] == FALSE)
  if (length(no_burrow_idx) > 0) {
    worms <<- NLset(turtles = worms, agents = worms[no_burrow_idx, ], var = "my_burrow", val = 0)
  }
  
  # --- Daytime movement ---
  if (!dark) {
    current_patch <- patchHere(world = habitat_type, turtles = worms)
    current_type  <- of(world = habitat_type, agents = current_patch)
    
    move_idx <- which(current_type != "burrow")
    
    if (length(move_idx) > 0) {
      burrow_patches <- NLwith(
        agents = patches(habitat_type),
        world  = habitat_type,
        val    = "burrow"
      )
      
      if (nrow(burrow_patches) > 0) {
        for (i in move_idx) {
          # pick 1 burrow patch at random
          rand_id <- sample(1:nrow(burrow_patches), size = 1)
          target  <- burrow_patches[rand_id, , drop = FALSE]
          
          # move this worm to that patch
          worms[i, ] <- moveTo(
            turtles = worms[i, ], 
            agents  = patch(world = habitat_type,
                            x = target[1],
                            y = target[2])
          )
          
          # update my_burrow
          my_burrow_val <- burrow_no@.Data[target[2], target[1]]
          worms <- NLset(
            turtles = worms,
            agents  = worms[i, ],
            var     = "my_burrow",
            val     = my_burrow_val
          )
        }
      }
    }
  }
  
  
  # --- Update speeds ---
  mass_vals <- of(agents = worms, var = "mass")
  worms <<- NLset(turtles = worms, agents = worms, var = "max_crawling_speed", val = 18.81 * (mass_vals ^ 0.0006))
  worms <<- NLset(turtles = worms, agents = worms, var = "burrowing_speed", val = 0.97 * (mass_vals ^ 0.06))
  
  # --- Mortality ---
  worms <<- calc_background_mortality(worms = worms)
  
  # ===================================================
  # === Adults ========================================
  # ===================================================
  adults_idx <- which(as.character(of(agents = worms, var = "breed")) == "adults")
  if (length(adults_idx) > 0) {
    adults <- worms[adults_idx, , drop = FALSE]
    
    # processes
    adults <- calc_maintenance(focal = adults)
    adults <- check_burrow(worms = worms, focal = adults)
    
    # reproduction pipeline
    out <- mateQ(focal = adults, worms = worms)
    adults <- out$focal
    worms  <- out$worms
   
    # update adults back
    worms<- upworms_with(group=adults, worms=worms)
  }
  
  # ===================================================
  # === Juveniles =====================================
  # ===================================================
  juveniles_idx <- which(as.character(of(agents = worms, var = "breed")) == "juveniles")
  if (length(juveniles_idx) > 0) {
    juveniles <- worms[juveniles_idx, , drop = FALSE]
    
    juveniles <- transform_juveniles(focal = juveniles)
    juveniles <- calc_maintenance(focal = juveniles)
    juveniles <- check_burrow(worms = worms, focal = juveniles)
    
    # growth
    growers_idx <- which(of(agents = juveniles, var = "energy_assimilated") > 0)
    
    if (length(growers_idx) > 0) {
      for (i in growers_idx) {
        grower_i <- juveniles[i, , drop = FALSE]    # extract this juvenile
        grower_i <- calc_growth(focal = grower_i)   # apply growth to one turtle
        juveniles[i, ] <- grower_i                  # put it back
      }
    }
    
    
    # update juveniles back
    worms<- upworms_with(group=juveniles, worms=worms)
    
  }
  # ===================================================
  # === Cocoons =======================================
  # ===================================================
  cocoons_idx <- which(as.character(of(agents = worms, var = "breed")) == "cocoons")
  if (length(cocoons_idx) > 0) {
    cocoons <- worms[cocoons_idx, , drop = FALSE]
    
    cocoons <- transform_cocoons(focal = cocoons)
    cocoons <- calc_maintenance(focal = cocoons)
    
    # update cocoons back
    worms<- upworms_with(group=cocoons, worms=worms)
    
}
  
  # --- Update patches & turtles ---
  food_density <- update_patches(food_density = food_density)
  worms<-update_turtles(worms = worms)
 
  # --- Stack world ---
  w <- stackWorlds(burrow_no, patch_quality, patch_T, patch_SWP, patch_Ex,
                   habitat_type, food_density, energy_content_food,
                   temperature, SWP, pcolor)
  
# --- Update age ---
  if(hour > 1){
    worms <<- NLset(turtles = worms, agents = worms, var = "age",
                    val = of(agents = worms,var = "age")+1)
  }
  
  # --- Daily cycle ---
  if (hour >= 24) {
    hour <<- 0
    day <<- day + 1
    worms <<- NLset(turtles = worms, agents = worms, var = "age",
                    val = of(agents = worms,var = "age")+1)
  }
  
  # --- Yearly cycle ---
  if (day == 365) {
    day <<- 1
    year <<- year + 1
  }
  
  # --- Diagnostics ---
  if (day %% 5 == 0 & hour == 14) {
    cat("Day", day, ":", "\n")
    cat("  mean mass =", mean(of(agents = worms, var = "mass")), "\n")
    cat("  mean energy_reserve =", mean(of(agents = worms, var = "energy_reserve")), "\n")
    cat("  mean energy_assimilated =", mean(of(agents = worms, var = "energy_assimilated")), "\n")
    cat("  mean ingestion_rate =", mean(of(agents = worms, var = "ingestion_rate")), "\n")
    cat("  mean maintenance_cost =", mean(of(agents = worms, var = "energy_growth")), "\n")
    cat("  adults =", mean(of(agents = worms, var = "age")), "\n")
    cat("  cocoons =", sum(of(agents = worms,var =  "breed") == "cocoons"), "\n")
    cat("  sperm_counter =", max(of(agents = worms,var =  "sperm_counter")), "\n")
    cat("  patches with food >", 0, "=", sum(of(world = food_density, agents = patches(world = food_density)) > 0), "\n\n")
  }
}
# Initial plot
plot(pcolor, main = paste0("hour = ",0, "; day = ",1,"; year = ",1, "; dark = ", dark) )
points(worms, pch = 16,
       col= (worms@.Data[,"color"]*2), cex = worms@.Data[,11]/10)
go()
plot(pcolor, main = paste0("hour = ",hour, "; day = ", day, "; year = ", year,"; dark = ", dark))
points(worms, pch = 16, col =(worms@.Data[,"color"]*2), cex = worms@.Data[,11]/5)


# --- Daily worm density monitor ---
worm_density <<- data.frame(
  day = numeric(0),
  density_surface = numeric(0),
  density_20cm = numeric(0),
  adults_n = numeric(0),
  juveniles_n = numeric(0)
)


pycor_20 <- habitat_type@pCoords[,2]

for(i in 1:1200) {
  go()
  # --- Worm density monitor (record at hour 8) ---
  if (hour == 8) {
    j<-day
    worms_a <- NLwith(agents = worms,var = "breed",val = c('adults','juveniles'))
    All_ind <- patchHere(world = habitat_type,turtles = worms_a)
    Ai20 <- NLwith(agents = worms_a, var = "ycor", val = All_ind[which(All_ind[,2]>80),2])
    Ais <- NLwith(agents = worms_a, var = "ycor", val = All_ind[which(All_ind[,2]>98),2])
    in_20cm<-NLcount(agents = Ai20)
   in_surface<-NLcount(agents = Ais)
   in_20_adults <- NLcount(NLwith(agents = worms_a, var = "breed", val= 'adults'))
   in_20_juveniles <- NLcount(NLwith(agents = worms_a, var = "breed", val= 'juveniles'))

   worm_density[j,1]<-day
   worm_density[j,2]<-in_surface
   worm_density[j,3]<-in_20cm
   worm_density[j,4]<-in_20_adults
   worm_density[j,5]<-in_20_juveniles
  }
  worm_density
  }

#### test speed #################
all.equal(worms_serial, worms_parallel)

worms_serial<-system.time(for(i in 1:10) go())

worms_parallel<-system.time(foreach(i=1:10) %do% go())

library(profvis)
profvis({
  for (i in 1:24) go()
})


######## 

library(doParallel)
library(foreach)

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

worm_density_list <- foreach(run = 1:10, .packages = c("NetLogoR")) %dopar% {
  source('submodels_eeeworm_7.4_STABLE09.R')
  source('init_eeeworm_2_STABLE.R')
  hour <<- 0; day <<- 1; year <<- 1
  worm_density <- data.frame(day = numeric(0), density_surface = numeric(0),
                             density_20cm = numeric(0), adults_n = numeric(0),
                             juveniles_n = numeric(0))
  for (i in 1:24) {
    go()
    # record as before
  }
  return(worm_density)
}

stopCluster(cl)
