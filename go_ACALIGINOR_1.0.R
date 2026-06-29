

###PROBLEM WITH JUVENIL MAYBE LINKED TO WORMS AT START BECAUSE NL IS COMPILED ? OR MAINTENANCE ? 

# =============================
# File: main_eeeworm_1.2.R (go)
# -----------------------------
# main go() file - sources init and submodels and runs hourly update

# When you source this main file it expects init_eeeworm_1.2.R and submodels_eeeworm_1.2.R to be in the same directory.

# Example usage (from interactive R):
source('Init_ACALIGINOR_1.0.R')
source('Sub_models_ACALIGINOR_2.1_PPP_cleaned.R')
 
Rotham_Acaliginosa<-read.table(file = "Rotham_Acaliginosa.txt", header = FALSE)
Rotham<-read.table(file = "C:/Users/qdevalloir/Downloads/Rothamsted.txt", header = FALSE)
# If the user prefers, the main loop below can be executed manually after sourcing the three files.

go <- function() {
  if (exists("worms") && nrow(worms@.Data) == 0) return(invisible())
  
  # --- Environment fluctuations ---
  #if (exists("Rotham")) setup_env_fluctuations(world = w, day = day, env_data = Rotham) else setup_env_fluctuations(world = w, day = day, env_data = NULL)
  
  if (length(worms) == 0) return(invisible(NULL))
    
    # Increment day
    day <<- day + 1
    
    # 1️⃣ Age and mortality update
    worms <- NLset(turtles = worms,
                    agents = worms,
                    var = "age",
                    val = of(agents = worms, var = "age") + 1)
    
    alive_at_start <- of(agents = worms, var = "who")
    Worms_at_start<- NLwith(agents = worms,var = "who",val = alive_at_start)
    #worms <- Worms_at_start
   # dead_worms<-NLwith(agents = worms,var = "die",val = TRUE)
   # dead_who<-of(agents = dead_worms,var = "who")
   # worms<-die(turtles = worms,who = dead_who)
      
    # ---- CHECKER A: mass after mortality ----
    if (length(worms) > 0) {
      m <- of(agents = worms, var = "mass")
      zero_or_neg <- which(m <= 0 | is.na(m))
      
      if (length(zero_or_neg) > 0) {
        cat("\n⚠️ MASS ISSUE detected immediately after mortality at day", day, "\n")
        cat("Worm IDs:", zero_or_neg, "\n")
        cat("Mass values:", m[zero_or_neg], "\n\n")
      }
    }
    # Get patch coordinates
    py <- food_density@pCoords[,2]  # get patch y-coordinates
    
   
    # --- 1. Food density by soil depth ---
    food_vals <- numeric(length(py))
    food_vals[py >= 15] <- rnorm(sum(py >= 15), mean = 4, sd = 0.4)
    food_vals[py >= 7 & py < 15] <- rnorm(sum(py >= 7 & py < 15), mean = 3.5, sd = 0.35)
    food_vals[py < 7] <- rnorm(sum(py < 7), mean = 3, sd = 0.3)
    food_density <<- NLset(world = food_density,
                           agents = patches(food_density),
                           val = food_vals)
    
    # --- 2. Temperature per patch ---
    temp_vals <-  of(world = temperature, agents = patches(temperature))
    temperature <<- NLset(world = temperature,
                          agents = patches(temperature),
                          val = temp_vals)
    T_val <- 273.15+ temperature
    
      # --- 3. Soil Water Potential (SWP) by depth ---
      soil_m_vals <- of(world = soil_moisture, agents = patches(soil_moisture))
      swp_vals <- numeric(length(py))
      
      # silty loam
      idx <- py >= 15
      swp_vals[idx] <- exp(13.31 - 97.38 * soil_m_vals[idx] + 298.80 * (soil_m_vals[idx]^2) - 349.70 * (soil_m_vals[idx]^3))
      
      # silty clay loam
      idx <- py >= 7 & py <15
      swp_vals[idx] <- exp(14.92 - 87.10 * soil_m_vals[idx] + 209.30 * (soil_m_vals[idx]^2) - 188.30 * (soil_m_vals[idx]^3))
      
      # clay
      idx <- py < 7
      swp_vals[idx] <- exp(50.75 - 339.1 * soil_m_vals[idx] + 861.8 * (soil_m_vals[idx]^2) - 766.2 * (soil_m_vals[idx]^3))
      
      SWP <<- NLset(world = SWP,
                    agents = patches(SWP),
                    val = swp_vals)
      
      # --- 4. Patch color (for visualization) ---
      pcolor_vals <- rep(50, length(py))
      pcolor <<- NLset(world = pcolor,
                       agents = patches(pcolor),
                       val = pcolor_vals)   
      
      # --- 5. Movement setup and food update ---
      setup_movement_Ex()   # assuming this function operates globally
      update_patches(food_density)   # your previously defined update_patches function
   # plot(food_density) 
      # 3️⃣ Environmental fluctuations (optional)
        setup_env_fluctuations(temperature = temperature,soil_moisture = soil_moisture,env_data = Rotham_Acaliginosa)
      
    # 4️⃣ Arrhenius factor
    T_vals <- of(world = temperature, agents =  patchHere(world = temperature,turtles = Worms_at_start))+273.15
    Arrhenius_vals <- exp((-activation_energy / Boltz) * ((1 / T_vals) - (1 / reference_T)))
        Worms_at_start<-NLset(turtles = Worms_at_start,agents = Worms_at_start,var = "Arrhenius_here", val= Arrhenius_vals)
    # ---- CHECKER B: Arrhenius sanity ----
    if (any(!is.finite(Arrhenius_vals))) {
      cat("\n⚠️ NON-FINITE ARRHENIUS at day", day, "\n")
      cat("T values:", T_vals[!is.finite(Arrhenius_vals)], "\n")
      cat("T values:", T_vals[(Arrhenius_vals)==0], "\n")
      cat("Arrhenius:", Arrhenius_vals[!is.finite(Arrhenius_vals)], "\n\n")
    }
    
    # 5️⃣ Turtle processes
       # A<- summary(worms[which(worms[,'breed']=='adults'),]@.Data)%>%data.frame
       # A[c(165,158,137,109,102, 88, 74),]
        
        adults_idx <- which(of(agents = Worms_at_start, var = "breed") == "adults")
    if (length(adults_idx) > 0) {
      #### PARAMETRIZING #####
      #worms_tick<- NLwith(agents = Worms_at_start,var = "who",val = alive_at_start)
      worms_adults <- NLwith(agents = worms,var = "breed",val='adults')
      
      #worms_who_before<-of(agents = worms, var = "who")
     # D<-summary(worms_adults@.Data)%>%data.frame
     # D[c(165,158,137,109,102, 88, 74),]
      
             #### INGESTION ########
      worms_adults <- calc_ingestion_rate(focal=worms_adults, all= Worms_at_start)
      
      #### MAINTENANCE ########
      worms_adults <- calc_maintenance(focal=worms_adults)
      dead_adu_worms_mtn<-NLwith(agents = worms_adults,var = "die",val = 1)
      dead_adu_who_mtn<-of(agents = dead_adu_worms_mtn,var = "who")
      cat("Adults dead: ", length(dead_adu_who_mtn),"\n")
      worms_adults<-die(turtles = worms_adults,who = dead_adu_who_mtn)
      worms<-die(turtles = worms,who = dead_adu_who_mtn)
      worms <- upworms_with(group = worms_adults, focal = worms)
      # Reset rejuvenators
      worms_adults <- NLwith(agents = worms,var = "breed",val='adults')
      
      #### MATE ###########
      # 1. Record which adults we start with
      worms_adults_who_before <- of(agents = worms_adults,var =  "who")
      worms_who_before <- of(agents = worms,var =  "who")
      before_mating<-of(agents = worms_adults, var = "hatchlings")
      # 2.
      worms_mated <- mate(focal = worms_adults)    # mate() no longer writes to global worms
      #cat("worms after", NLcount(agents = worms_adults_pups))
      
      # 3. Recover adults using original WHO (never indices)
      worms_adults_maters <- NLwith(agents = worms_mated, var = "breed", val = 'adults')
      
      # 4. Extract new cocoons (those with breed="cocoons")
      new_cocoons <- NLwith(agents = worms_mated, var    = "breed", val = "cocoons")
      
      # 5. Hatch new cocoons into global worms
      # (each adult producing 1 cocoon when hatchlings increase)
      hatchlings_post_mate <- of(agents = worms_adults_maters,var =  "hatchlings")
      worms_adults<-NLset(turtles = worms_adults,agents = worms_adults_maters,var = "hatchlings", val=hatchlings_post_mate)
      after_mating <- of(agents = worms_adults,var =  "hatchlings")
      breeders_id  <- which(after_mating > before_mating)
      breeders_who <- of(agents = worms_adults,var =  "who")[breeders_id]
      
      worms <- hatch(turtles = worms, who = breeders_who, n = 1, breed = "cocoons")
      
      # 6. Find the newly created WHO values
      worms_who_after <- of(agents = worms,var =  "who")
      new_who <- setdiff(worms_who_after, worms_who_before)
      
      # 7. Assign new WHO to the cocoons in worms_mated
      worms_mated <- NLset(turtles = worms_mated, agents = new_cocoons,
        var = "who", val  = new_who)
      
      # 8. Merge updated adults + cocoons back into worms
      worms_adults <- upworms_with( group = worms_adults_maters,
        focal = worms_adults)
     
      worms <- upworms_with( group = worms_mated, focal = worms)
    #  tail(worms)
      # 9. Safety check
      if (any(of(agents = worms_adults,var =  "breed") != "adults")) {
        stop("ERROR: Cocoons leaked into adults")
      }
     
      ####### MOVE ###########
      worms_adults <- move(focal=worms_adults)
     # worms_adults[,'breed']
      ###### AESTIVATE ########
      #of(agents = worms_adults,var = "mass")%>%mean(na.rm = TRUE)
      worms_adults <- aestivate(worms_adults)
      dead_adu_worms<-NLwith(agents = worms_adults,var = "die",val = 1)
      dead_adu_who<-of(agents = dead_adu_worms,var = "who")
      cat("Adults dead: ", length(dead_adu_who),"\n")
      worms_adults<-die(turtles = worms_adults,who = dead_adu_who)
      worms<-die(turtles = worms,who = dead_adu_who)
      ####### UPDATE ##########  
      worms <- upworms_with(group = worms_adults, focal = worms)
  #   points(worms_adults)
    }
        #B<- summary(worms[which(worms[,'breed']=='juveniles'),]@.Data)%>%data.frame
        #B[c(165,158,137,109,102, 88, 74),]
        
        juveniles_idx <- which(of(agents = Worms_at_start, var = "breed") == "juveniles")
        if (length(juveniles_idx) > 0) {
          #worms_tick<- NLwith(agents = Worms_at_start,var = "who",val = alive_at_start)
          worms_juv <- NLwith(agents = worms,var = "breed",val='juveniles')
      
      D<-summary(worms_juv@.Data)%>%data.frame
      D[c(165,158,137,109,102, 88, 74),]
      
      worms_juv <- calc_ingestion_rate(focal = worms_juv,all= Worms_at_start)
      worms_juv <- calc_maintenance(focal =  worms_juv)
      dead_juv_worms_mtn<-NLwith(agents = worms_juv,var = "die",val = 1)
      dead_juv_who_mtn<-of(agents = dead_juv_worms_mtn,var = "who")
      cat("Juveniles dead: ", length(dead_juv_who_mtn),"\n")
      worms_juv<-die(turtles = worms_juv,who = dead_juv_who_mtn)
      worms<-die(turtles = worms,who = dead_juv_who_mtn)
      worms <- upworms_with(group = worms_juv, focal = worms)
      # Reset juveniles
      worms_juv <- NLwith(agents = worms,var = "breed",val='juveniles')
      
      worms_juv <- calc_growth(focal = worms_juv)
      worms_juv <- transform_juvenile(focal  = worms_juv)
      worms <- upworms_with(group = worms_juv, focal = worms)
      worms_juv <- NLwith(agents = worms,var = "breed",val='juveniles')
      worms_juv <- move(focal  = worms_juv)
       worms_juv <- aestivate(focal =  worms_juv)
       dead_juv_worms<-NLwith(agents = worms_juv,var = "die",val = 1)
       dead_juv_who<-of(agents = dead_juv_worms,var = "who")
       cat("Juv dead: ", length(dead_juv_who),"\n")
       worms_juv<-die(turtles = worms_juv,who = dead_juv_who)
       worms<-die(turtles = worms,who = dead_juv_who)
      worms<- upworms_with(group=worms_juv, focal=worms)
    }
    # ---- CHECKER C: mass after growth / maintenance / movement ----
    
    #cat("\nDay", day, 
    #    "| min mass:", round(min(of(agents = worms, var = "mass")), 4),
    #    "| mean mass:", round(mean(of(agents = worms, var = "mass")), 4),
    #    "| dead/zero:", sum(of(agents = worms, var = "mass") <= 0),
    #    "\n")
      
            
    cocoons_idx <- which(of(agents = Worms_at_start, var = "breed") == "cocoons")
    if (length(cocoons_idx) > 0) {
     # worms_tick<- NLwith(agents = Worms_at_start,var = "who",val = alive_at_start)
      worms_cocoons <- NLwith(agents = worms,var = "breed",val='cocoons')
      
      ## EMBRYPO DEV
      worms_cocoons <- calc_embryo_development(focal = worms_cocoons)
      
      ### MAINTENANCE
      worms_cocoons <- calc_maintenance(focal = worms_cocoons)
     # dead_cocoons_worms_mtn<-NLwith(agents = worms_cocoons,var = "die",val = 1)
     # dead_cocoons_who_mtn<-of(agents = dead_cocoons_worms_mtn,var = "who")
     # cat("cocoons dead: ", length(dead_cocoons_who_mtn),"\n")
     # worms_cocoons<-die(turtles = worms_cocoons,who = dead_cocoons_who_mtn)
     # worms<-die(turtles = worms,who = dead_cocoons_who_mtn)
     # worms <- upworms_with(group = worms_cocoons, focal = worms)
     # worms_cocoons <- NLwith(agents = worms,var = "breed",val='cocoons')
      ## TRANSFORM
      worms_cocoons <- transform_cocoon(focal = worms_cocoons)
      worms <-upworms_with(group=worms_cocoons, focal=worms)
    }
         
    worms <<- update_turtles(worms)
    #worms <<- calc_PPP_exposure(focal = worms,pesticide_layer =  pesticide_concentration)
  
    
    # 7️⃣ Year/day rollover
    if (day > 365) {
      day <- 1
      year <- year + 1
    }
     
   #  if (any(of(agents = worms, var = "mass") <= 0, na.rm = TRUE)) {
   #    stop("Simulation stopped: At least one worm reached zero or negative mass")
#     }
}

# NOT MOVING
for(i in 1:30){
go()
day
colsoil<-colorRampPalette(c("springgreen","green2","green4"))
colworm<-colorRampPalette(c("red","pink","white"))
plot(energy_content_food, main = paste("Nadults: ", NLcount(worms[which(worms[,'breed']=='adults'),]),"\n",
                                     "Njuv: ", NLcount(worms[which(worms[,'breed']=='juveniles'),]),"\n",
                                     "Ncocoons: ", NLcount(worms[which(worms[,'breed']=='cocoons'),]),"\n",
                                     "day: ", day))#sum(worms@.Data[which(worms[,'breed']=='cocoons'),'aestivating']),"\n"))#, col=colsoil(50))
points(worms, col=colworm(worms[,'color']), pch=18,cex=0.7,)

}

NLwith(agents = worms,var = 'breed',val = 'juveniles')
(worms[which(worms[,'breed']=='juveniles'),])
mass_birth

test <- patchHere(world = food_density, turtles = worms)

cat("NAs in patchHere:", sum(is.na(test)), "\n")


plot(x=0,y=0,xlim=c(0,200), ylim=c(0,700))
z<-c() 
y<-c()
x<-c()

a <- c()
b <- c()
c <- c()
d <- c()
e <- c()
f <- c()
g <- c()

PPPext <- c()
PPPint <- c()

365+28
for(i in 1:200){
go()  
# points(x=worms@.Data[1,1],y=worms@.Data[1,2]) 

  x<-rbind(x,NLcount(worms[which(worms[,"breed"]=='adults'),]))#,var = "ingestion_rate"),na.rm = TRUE))
  y<-rbind(y,NLcount(worms[which(worms[,"breed"]=='juveniles'),]))#,var = "mass"),na.rm = TRUE))
  z<-rbind(z,NLcount(worms[which(worms[,"breed"]=='cocoons'),]))#,var = "energy_assimilated"),na.rm = TRUE))
  
  a<-rbind(a,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "mass"),na.rm = TRUE))
  b<-rbind(b,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "R"),na.rm = TRUE))
  c<-rbind(c,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "max_R"),na.rm = TRUE))
  d<-rbind(d,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "energy_reserve"),na.rm = TRUE))
  e<-rbind(e,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "energy_assimilated"),na.rm = TRUE))
  f<-rbind(f,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "BMR"),na.rm = TRUE))
  g<-rbind(g,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "ingestion_rate"),na.rm = TRUE)) 
  h<-rbind(c,mean(of(agents= worms[which(worms[,"breed"]=='juveniles'),], var = "energy_reserve_max"),na.rm = TRUE))
  #PPPext <- cbind(PPPext,of(agents= worms[which(worms[,"breed"]=='adults'),], var = "PPP_external_concentration"))
  #PPPint <- cbind(PPPint,of(agents= worms[which(worms[,"breed"]=='adults'),], var = "PPP_internal_concentration"))
  
  points(y=x[i], day,col="red",cex=0.2)
  points(y=y[i], day,col="blue",cex=0.2)
  points(y=z[i], day,col="green",cex=0.2)
  
i+1

cat("Day", day, "Total cocoons:", sum(of(agents = worms, var = "breed") == "cocoons"), "\n",
    "Total adults:", sum(of(agents = worms, var = "breed") == "adults"), "\n")

}

#write.csv(data.frame(PPPext), file = "C:/Users/qdevalloir/Documents/ACALIGINOR/TEST_FILE/strat_uni_28D_ADULTS01_.CSV")
#write.csv(data.frame(PPPint), file = "C:/Users/qdevalloir/Documents/ACALIGINOR/TEST_FILE/strat_uni_28D_ADULTS01_INT.CSV")

plot(x=0,y=0,xlim=c(0,365), ylim=c(0,700),xlab= "Days", ylab="Population density")
lines(y=x, x= 1:nrow(x),col="red")
lines(y=y, x= 1:nrow(y),col="pink")
lines(y=z, x= 1:nrow(z),col="blue")
lines(y=x+y+z, x= 1:nrow(z),col="green")


plot(x=0,y=0,xlim=c(0,400), ylim=c(0,2),xlab= "Days", ylab="Population density")
lines(y=a, x= 1:nrow(a),col="red")
lines(y=b, x= 1:nrow(b),col="blue")
lines(y=c, x= 1:nrow(c),col="yellow")
lines(y=d, x= 1:nrow(d),col="pink")
lines(y=e, x= 1:nrow(e),col="grey")
lines(y=f, x= 1:nrow(f),col="orange")
lines(y=g, x= 1:nrow(g),col="green")
lines(y=h, x= 1:nrow(h),col="black")


# 
go()   # run first step (HOUR = 0 + 1 = 1)
# 2nd plot (HOUR = 1 , DAY = 1, Year = 1)       
plot(pcolor, main = paste0("day = ", day, "; year = ", year))
points(worms, pch = 16, col =(worms@.Data[,"color"]*2), cex = worms@.Data[,11]/50)
worms@.Data[,'ycor']
