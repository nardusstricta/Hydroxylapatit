(1/3 * (-31.33)) - (5/3 * log10(88^-4)) - (4/3 * 7.57)
#Berechung des Gleichgewichtes nach: Ksp = (Ca2+)^5 * (PO43-) * OH-
#Umstellung der Gleichung: log(P) = 10 log(Ca) * 2log(OH-)/Ksp)/6
##
P_oh1 <- function(oh, ca, temp_oh){
  erg_oh1 <- ((10 * log10(ca)) + (2*log10(oh)) + temp_oh)/6 
  return(erg_oh1)
}
#Temperaturabhängigkeit nach McDowell et.al
#Bezogen auf folgende Gleichung: Ksp = (Ca2+)^5 * (PO43-) * OH-
#logKs = -8219.41/Temperatur - 1.6657 - 0.098215 * Temperatur
Ks_temp <- function(temp_P){
  temp_P <- 273.15 + temp_P #Grad Celsius in Kelvin umrechenen
  temp_erg <- -8219.41 / temp_P - 1.6657 - 0.098215 * temp_P
  return(temp_erg)
}
temp_beisp <- seq(0, 40, 1)
plot(temp_beisp, Ks_temp(temp_beisp), ylab="pKs Wert", xlab="Temperatur in C")

#Rotationsmatrix
rotate <- function(x) t(apply(x, 2, rev))

#Lösungsgleichgewicht errechenen von Hydrogenphosphat:
#Funktion der Spezies in Abhängigkeit von der Aktivität, Temperatur und Anfangskonzentration von CaHPO4 
