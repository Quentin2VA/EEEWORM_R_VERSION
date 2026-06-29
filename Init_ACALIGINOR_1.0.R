# =============================
# File: init_eeeworm_1.2.R
# EEEWORM initialization script (version 1.2)
# Author: converted for R by ChatGPT (based on Quentin Devalloir)
# Date: 2025-10
# -----------------------------

# Libraries
if (!requireNamespace("NetLogoR", quietly = TRUE)) stop("Install the NetLogoR package first")
if (!requireNamespace("SpaDES.tools", quietly = TRUE)) warning("SpaDES.tools recommended for advanced utilities")
library(NetLogoR)
library(SpaDES.tools)
library(ggplot2)
library(dplyr)
library(sf)

# -----------------------------
# 1. Global parameters (defaults)
# -----------------------------
assimilation_efficiency <- 0.19       # no units
B_0 <- 967                            # baseline B0
activation_energy <- 0.25             # eV
energy_flesh <- 7                     # kJ/g
energy_synthesis <- 3.6               # kJ/g
max_ingestion_rate <- 0.805            # g/hour/g^(2/3) 
half_saturation_coef <- NA       # half saturation coefficient, calculated by fitting a Holling Type II response curve to food ration vs. ingestion rate data.
mass_birth <- 0.025                   # g                  # g
mass_sexual_maturity <- 0.5           # g
mass_maximum <- 2                   # g
mass_cocoon <- 0
growth_constant <- 0.049             # /hour
max_reproduction_rate <- 0.054      # kJ/g/hour
incubation_period <- 62               # days
reference_T <- 288.15                 # Kelvins
Arrhenius <- NA
Boltz <- (8.62*(10^-5))                      # eV K^-1

T_kelvins <- 283.15                   # default patch fallback if NA

# Time globals (created in main script but can be set here if sourced alone)

day  <<- 0
year <<- 1


# -----------------------------
# 2. World / patch layers
# -----------------------------
min_Pycor_move <- 1
max_Pycor_move <- 20
min_Pxcor_move <- 1
max_Pxcor_move <- 60

min_Pycor <- 1
max_Pycor <- 20
min_Pxcor <- 1
max_Pxcor <- 60

# create blank worlds (NetLogoR)
w <- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                 minPycor = min_Pycor, maxPycor = max_Pycor, data = NA_integer_)

change_food_density <- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                             minPycor = min_Pycor, maxPycor = max_Pycor, data = 0)
food_density<- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                           minPycor = min_Pycor, maxPycor = max_Pycor, data = 3.5)
energy_content_food<- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                           minPycor = min_Pycor, maxPycor = max_Pycor, data = 0)
temperature<- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                          minPycor = min_Pycor, maxPycor = max_Pycor, data = 0)
soil_moisture <-createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                            minPycor = min_Pycor, maxPycor = max_Pycor, data = 0)
SWP<- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                  minPycor = min_Pycor, maxPycor = max_Pycor, data = 0)
pcolor<- createWorld(minPxcor = min_Pxcor, maxPxcor = max_Pxcor,
                  minPycor = min_Pycor, maxPycor = max_Pycor, data = 50)
pesticide_concentration <- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0) # for future use
pesticide_ref <- pesticide_concentration
# SETUP patches
#food_density <- NLset(world = food_density, agents = patches(food_density), val = 3.5)

# Stack world layers (master stacked world for quick plotting/inspection)

# -----------------------------
# 3. Turtles (worms) initialization
# -----------------------------
n_worms_adults <- 100
n_worms_juveniles <- 100
n_worms_cocoons <- 100
n_worms<- n_worms_adults + n_worms_juveniles + n_worms_cocoons

coords <- cbind(runif(n_worms, min_Pxcor, max_Pxcor), runif(n_worms, min_Pycor, max_Pycor))

worms <- createTurtles(n = n_worms, coords = coords,
                       breed = c(rep("adults", n_worms_adults),rep("juveniles", n_worms_juveniles),rep("cocoons", n_worms_cocoons)),
                       heading = runif(n_worms, 0, 360),
                       color = c(rep("red", n_worms_adults),rep("pink", n_worms_juveniles),rep("white", n_worms_cocoons)))

# Add turtles-own variables (NetLogoR's turtlesOwn helper)
worms <- turtlesOwn(turtles = worms, tVar = "age", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "mass", tVal = rep(mass_sexual_maturity, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "energy_reserve", tVal = (of(agents = worms, var = "mass") / 2) * energy_flesh)
worms <- turtlesOwn(turtles = worms, tVar = "size", tVal = rep(mass_sexual_maturity, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "energy_assimilated", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "energy_reserve_max", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "ingestion_rate", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "BMR", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "mortality", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "max_growth_rate", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "growth_rate", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "energy_growth", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "embryonic_development", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "hatchlings", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "max_R", tVal = rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "R", tVal = rep(0,n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "aestivating", tVal = rep(FALSE, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "time_aestivating",tVal =  rep(1, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "mass_cocoons",tVal =  rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "Arrhenius_here",tVal =  rep(0, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "visible",tVal =  rep(TRUE, n_worms))
worms <- turtlesOwn(turtles = worms, tVar = "die",tVal =  rep(FALSE, n_worms))
adults<-NLwith(agents = worms,var = 'breed',val='adults')
juveniles<-NLwith(agents = worms,var = 'breed',val='juveniles')
cocoons<-NLwith(agents = worms,var = 'breed',val='cocoons')


# Add turtles-own variable for pesticide exposure
worms <- turtlesOwn(worms, "PPP_external_concentration", rep(0, n_worms))
worms <- turtlesOwn(worms, "PPP_internal_concentration", rep(0, n_worms))


## SETUP turtles
#adults
worms<-NLset(turtles = worms, agents = adults,var = "size", val= 0.3)
worms<-NLset(turtles = worms, agents = adults,var = "aestivating", val= FALSE)
worms<-NLset(turtles = worms, agents = adults,var = "mass", val= mass_sexual_maturity)
worms<-NLset(turtles = worms, agents = adults,var = "time_aestivating", val= 1)
# juveniles
worms<-NLset(turtles = worms, agents = juveniles,var = "size", val= 0.2)
worms<-NLset(turtles = worms, agents = juveniles,var = "mass", val= mass_birth)
juveniles<-NLset(turtles = juveniles, agents = juveniles,var = "mass", val= mass_birth)
worms <- NLset(turtles = worms, agents = juveniles, var = "energy_reserve", val = (of(agents = juveniles, var = "mass") / 2) * energy_flesh)
worms<-NLset(turtles = worms, agents = juveniles,var = "aestivating", val= FALSE)
worms<-NLset(turtles = worms, agents = juveniles,var = "time_aestivating", val= 1)
# cocoons
worms<-NLset(turtles = worms, agents = cocoons,var = "size", val= 0.1)
worms<-NLset(turtles = worms, agents = cocoons,var = "mass", val= mass_birth)
worms<-NLset(turtles = worms, agents = cocoons,var = "age", val= 0)
worms<-NLset(turtles = worms, agents = cocoons,var = "energy_reserve", val= mass_cocoon*energy_flesh)

adults<-NLwith(agents = worms,var = 'breed',val='adults')
juveniles<-NLwith(agents = worms,var = 'breed',val='juveniles')
cocoons<-NLwith(agents = worms,var = 'breed',val='cocoons')



message("Initialization complete: init_eeeworm_1.2.R sourced. Variables created: w, worms, food_density, habitat_type, burrow_no, pcolor, etc.")

