R18_VSMOW <<- 0.0020052 # ratio
R2_VSMOW <<- 0.00015576 # ratio
diffratio_18 <<- 1/0.9723 # oxygen D/D* for diffusion of water vapor through air.  Merlivat, 1978
diffratio_2 <<- 1/0.9755 # hydrogen D/D* for diffusion of water vapor through air. Merlivat, 1978

alpha_eq = function(temp){
  alpha18_l_v <<- exp(-2.0667 * 10^-3 - 0.4156/(temp+273.15) + (1.137*10^3)/(temp+273.15)^2) # Majoube (1971)
  alpha2_l_v <<- exp(52.612*10^-3 - 76.248/(temp+273.15) + (24.844*10^3)/(temp+273.15)^2) # Majoube (1971)
}