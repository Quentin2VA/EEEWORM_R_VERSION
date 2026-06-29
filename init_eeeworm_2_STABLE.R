
# EEEWORM model -  'init' - Version 2 #STABLE
# Author: Quentin DEVALLOIR
# Contact: @QUENTIN2VA
# Version: 2.0
# Note: Please contact me before using it, this is still a trial version.

# Earthworm individual-based model in R using NetLogoR package

library(NetLogoR)
library(SpaDES.tools)

#1. Initialize Global Parameters in R

  assimilation_efficiency = 0.55       # no units
  B_0 = 360                            # no units
  activation_energy = 0.32             # eV
  energy_flesh = 7                     # kJ/g
  energy_synthesis = 3.6               # kJ/g
  max_ingestion_rate = 0.37            # g/hour/g^(2/3)
  burrowing_costs = 0.0103             # kJ/cm/hour
  mass_birth = 0.053                   # g
  mass_cocoon = 0.061                  # g
  mass_sexual_maturity = 4.2           # g
  mass_maximum = 8.5                   # g
  growth_constant = 0.0023             # /hour
  max_reproduction_rate = 0.00021      # kJ/g/hour
  incubation_period = 90               # days
  reference_T = 288.15                 # Kelvins
  Boltz = 8.62e-5                     # eV K^-1
  background_mortality = 0.004         # %/hour
  T_kelvins = 283.15                           # Temperature in Kelvins


#2. Define Patch Variables

# Create the patch world
  #World dimensions
  min_Pxcor <- 1
  max_Pxcor <- 301
  min_Pycor <- 1
  max_Pycor <- 106
  
  
w <- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= NA_integer_)
burrow_no <- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
patch_quality <- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
patch_T<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
patch_SWP<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
patch_Ex<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
habitat_type<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= "soil")
food_density<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
energy_content_food<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
temperature<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
SWP<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 0)
pcolor<- createWorld(min_Pxcor, max_Pxcor, min_Pycor, max_Pycor, data= 52) # default color for soil = 52

# Initialize patches-own variables as layers in a stacked world
 #burrow_no:  unique burrow IDs
 #patch_quality:  initialize with 0
 #patch_T:  temperature component
 #patch_SWP:  soil water potential
 #patch_Ex:  energy content of food
 #habitat_type:  habitat category, text
 #food_density:  grams
 #energy_content_food:  kJ/g
 #temperature:  degrees Celsius
 #SWP:  soil water potential (-kPa)


# Assign soil-surface vs soil
surface_coords <- which(habitat_type@pCoords[,2] > 99)
habitat_type <- NLset(world = habitat_type, agents = patches(habitat_type)[surface_coords,], val = "soil-surface")
pcolor <- NLset(pcolor, agents = patches(pcolor)[surface_coords,], val = 37)

### 3. Setup Patches (NetLogo translation)
# Burrow positions as per NetLogo setup

## === BURROW COLUMNS (override baseline) ===
px <- c(  9,  34,  59,  84, 109, 135, 159, 185, 209, 229, 254, 295)
off<- c( 20,  30,  10,  40,  20,  30,  40,  10,  20,  20,  30,  10)
bno<- c(  1,   2,   3,   4,   5,   6,   7,   6,   7,   8,   9,  10)
pc <- c( 38,  37,  36,  34,  33,  37,  36,  37,  36,  38,  37,  36)

for (i in seq_along(px)) {
  y_start <- off[i] + sample(0:20, 1) + 1     # strictly greater than (offset + random 20)
  if (y_start <= min(99, max_Pycor)) {
    coords <- cbind(pxcor = px[i], pycor = y_start:min(99, max_Pycor))
    habitat_type <- NLset(world = habitat_type, agents = coords, val = "burrow")
    pcolor       <- NLset(world = pcolor,       agents = coords, val = pc[i])
    burrow_no    <- NLset(world = burrow_no,    agents = coords, val = bno[i])
  }
}

### 4. Stack world from setup
# Create stacked world
w <- stackWorlds(w, burrow_no,
                 patch_quality,
                 patch_T,
                 patch_SWP,
                 patch_Ex,
                 habitat_type,
                 food_density,
                 energy_content_food,
                 temperature,
                 SWP,
                 pcolor)

# 5. Define Turtles Variables (turtles-own)
###  Setup worms
n_worms <- 10
coords <- cbind(runif(n_worms, min_Pxcor, max_Pxcor), rep(98, n_worms)) # like setxy random-xcor 98

worms <- createTurtles(
  n = n_worms,
  coords = coords,
  breed = rep("adults", n_worms),
  heading = runif(n_worms, 0, 360),
  color = rep("red", n_worms)
)


# Add and set turtles-own variables with default values
worms <- turtlesOwn(worms, "age", rep(0, n_worms))
worms <- turtlesOwn(worms, "mass", rep(mass_sexual_maturity, n_worms))          # start at birth mass
worms <- turtlesOwn(worms, "size",rep(mass_sexual_maturity, n_worms))   # size = mass
worms <- turtlesOwn(worms, "energy_reserve", (of(agents = worms,var = "mass") / 2) * energy_flesh)

worms <- turtlesOwn(worms, "energy_assimilated", rep(0, n_worms))
worms <- turtlesOwn(worms, "energy_reserve_max", rep(0, n_worms))           # calculate based on mass later
worms <- turtlesOwn(worms, "ingestion_rate", rep(0, n_worms))
worms <- turtlesOwn(worms, "BMR", rep(0, n_worms))

worms <- turtlesOwn(worms, "my_burrow", rep(NA_integer_, n_worms))
worms <- turtlesOwn(worms, "last_burrow", rep(NA_integer_, n_worms))
worms <- turtlesOwn(worms, "burrowing_speed", rep(0, n_worms))
worms <- turtlesOwn(worms, "max_burrowing_costs", rep(0, n_worms))

worms <- turtlesOwn(worms, "max_growth_rate", rep(0, n_worms))
worms <- turtlesOwn(worms, "growth_rate", rep(0, n_worms))
worms <- turtlesOwn(worms, "energy_growth", rep(0, n_worms))

worms <- turtlesOwn(worms, "embryonic_development", rep(NA, n_worms))
worms <- turtlesOwn(worms, "hatchlings", rep(0, n_worms))
worms <- turtlesOwn(worms, "max_R", rep(NA, n_worms))
worms <- turtlesOwn(worms, "R", runif(n_worms, 0, 0.4))
worms <- turtlesOwn(worms, "sperm_counter", sample(0:2184, n_worms, replace = TRUE))
worms <- turtlesOwn(worms, "favourite_patch", rep(NA_integer_, n_worms))
worms <- turtlesOwn(worms, "max_crawling_speed", rep(0, n_worms))
worms <- turtlesOwn(worms, "burrowQ", rep(FALSE, n_worms))

worms
###################

# Visualize soil:

#library(rgdal);library(png);library(sp) 
# load icons in PNG format
#icon1 <- readPNG('C:/Users/qdevalloir/Pictures/icon1.png')
# need to offset the x/y location for the icon (depends on desired icon size)
#offset <- 15
#rasterImage(icon1, worms@.Data[,1]-offset, worms@.Data[,2]-offset, worms@.Data[,1]+offset, worms@.Data[,2]+offset)

plot(habitat_type)
points(worms, pch= 16, col ="red",cex=worms@.Data[,11]/10)
