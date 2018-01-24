####Berechnung des Dissoziationsdiagramm auf Basis der Konzentrationen
#Gleichgewichtskonstanten:

source("rotation_fun.R")
pk1 <- 2.16
pk2 <- 7.21
pk3 <- 12.32
#pks-Werte (-log10 von Gleichgewichtskonstante)
pH <- seq(3, 14, 0.2)
h <- 10 ^ (-pH) 
K1 <- 10 ^ (-pk1)

K2 <- 10 ^ (-pk2)
K3 <- 10 ^ (-pk3)

#Gesamtkonzentration H3 ermittelt über die Konzentrationssumme 
H3 <- h^3 + (h^2 * K1) + (h * K1 * K2) + (K1 * K2 * K3)

# Gleichungen der vier Gleichgewichtskurven im Dissoziationsdiagramm
# Teilen durch die Gesamtkonzentration H3
RAH3 <- (h^3)/H3
RAH2 <- h^2 * K1 / H3
RAH <- h * K1 * K2/H3
RA <- K1 * K2 * K3 / H3
plot(pH, RAH3, type="l", ylim=c(0,1), ylab = "Konzentrationsanteile", main = "Dissoziationsdiagramm der Phosphorsäure")
points(pH, RAH2, type = "l", col="red")
points(pH, RAH, type = "l", col ="blue", lty=2)
points(pH, RA, type = "l", col = "green", lty=2)
legend("topright", c(expression("H"["3"]*"PO"["4"]), expression("H"["2"]*"PO"["4"]^"-"), expression("HPO"["4"]^"2-"), expression("PO"["4"]^"3-")), col = c("black", "red", "blue", "green"), lty = c(1,1,2,2), bty = "n")

###
#Funktion des gelösten Phosphatgehaltes im Gelichgewicht mit Hydroxylapatit###########
#bei gegeben pH und ca+ aktivität in mol/l
###
#Temperaturabhängigkeit nach McDowell et.al
#log Ks = -8219.41/Temperatur - 1.6657 - 0.098215 * Temperatur
Ks_temp <- function(temp_P){
  temp_P <- 273.15 + temp_P #Grad Celsius in Kelvin umrechenen
  temp_erg <- -8219.41 / temp_P - 1.6657 - 0.098215 * temp_P
  return(temp_erg)
}
temp_beisp <- seq(0, 40, 1)
plot(temp_beisp, Ks_temp(temp_beisp))
#Idee: Löslichkeitsprodukt (Ksp) formolieren aus der Gleichgewichtsreaktion von Hydroxylapatit (Ca5(PO4)3OH + 7 H+ <-> 5 Ca2+ + 3H2PO4- + H2O); Ksp = a5 Ca2+ * a3 H2PO4- / a7 H+ (a = aktivität)
#Nach umformen und einsetzen von Ksp = 14.46 ->
#log(aH2PO4-) = 1/3 (14.46 - 7 *pH - 5 log aCa2+)
co2pa_ca(pco02=0.0303975082, temp=25)

ca2(0.0303975082, 25, 8.31753683, 0.84540208)
erg1 <- HPO2(8.31753683,0.0005018758)

p <- co2pa_ca(pco02=0.3039750821, temp=25)
p1 <- ca2(0.3039750821, 25, 7.6659114, p[1])
erg2 <- HPO2(7.6659114,p1)
erg1-erg2
log10(erg2-erg1)
5/3

HPO2 <- function(ph, gammaca = 0.0025, temp_P = 20){
  erg <- 1/3 * (14.46 - (5 * log10(gammaca)) - (7 * ph))
  erg <- 10 ^ erg  * 30.973762 #Umrechung von Mol/l zu Gramm/l
  return(erg *1000) #In Miligramm *1000
}
##
#Beispiel#####
##
HPO2(5)
HPO2(6)
ca_beispiel <- seq(0.0001, 0.009, 0.0001)
pH_beispiel <- seq(3, 10, 0.01)
erg1 <- sapply(pH_beispiel, HPO2, gammaca = 0.0025) #Bei gleichen Ca2+
plot(pH_beispiel, log10(erg1), type = "l", ylab = "Gelöster Phosphor log(mg/l)", xlab = " pH")
#Mit Veränderung der Ca2+ Aktivität
mat_erg <- as.data.frame(matrix(0,length(ca_beispiel), length(pH_beispiel)))
rownames(mat_erg) <- ca_beispiel
colnames(mat_erg) <- pH_beispiel
for (i in 1:length(ca_beispiel)){
  for (j in 1:length(pH_beispiel)){
    mat_erg[i,j] <- HPO2(ph = pH_beispiel[j],  gammaca = ca_beispiel[i])
  }
}
#mat_erg$`5.01`[which(ca_beispiel==0.0025)]
mat_erg <- rotate(mat_erg) # rotation der Matrix 
filled.contour(pH_beispiel, ca_beispiel,   log10(mat_erg), ylab=expression("Ca"^"2+"*"Aktivität in mol l"^"-1"), xlab=("pH-Wert"), main  = "Phosphor in Lösung log(mg/l)", color = terrain.colors)

#Phosphorgehalt Abhängig von der CO2 Konzentration############
CO2_PD <- seq(0.1, 5, 0.1) #CO2 Patialdrücke  in kPA
ph_co2 <- as.data.frame(matrix(0,length(CO2_PD),8)) # Matrix für die Aktivitäten von Ca2+ und den dazugehörigen pH-Wert
colnames(ph_co2) <- c("ca2", "CaHCO", "H", "HCO3", "CO2", "OH", "PH", "paco2")

#Berechung der Aktivitäten und pH-Werten anhand des Calciumsmodelles mit verschienen  CO2 Patialdrücken  
for (i in 1:length(CO2_PD)){
  ph_co2[i,] <- co2pa_ca(pco02=CO2_PD[i], temp1=20)
  print(CO2_PD[i])
}
ph_co2
cacons<-numeric(length(CO2_PD)) #Berechnung der Ca2+ Konzentrationen aus den Aktivitäten 
for (i in 1:length(CO2_PD)){
  cacons[i] <- ca2(pco02=ph_co2$paco2[i], temp1=20, ph=ph_co2$PH[i], ga=ph_co2$ca2[i])
}

plot(ph_co2$paco2, cacons * 40 *1000, type = "l", ylab = expression("Ca" ^ "2+" * " Konzentration in mg" ^ "-1"), xlab = expression("CO" ["2"] * " Patialdruck kPa"))

#Berechung der Phosphatgehalte mit den ca2+ Konzentrationen und dem pH-Wert im Gleichgewicht
plot(ph_co2$PH, cacons, type = "l", xlab = "pH", ylab = expression("ca"^"2+"*"  mol l"^"-1"))
phos <- HPO2(ph = ph_co2$PH, gammaca = cacons)
plot(ph_co2$PH, log10(phos), type = "l")

phos <- as.data.frame(cbind(phos, ph_co2$PH, cacons, CO2_PD))
names(phos) <- c("Phosph", "pH", "Ca2+", "CO2_Patial")

library(ggplot2)
p <- ggplot(phos, aes(pH, Phosph, color = CO2_Patial))
  p + geom_point()



plot(CO2_PD, log10(cacons), type = "l", xlab = "Patialdruck CO2", ylab = "Ca2+ Konzentration mol/l")
#Plot aus Lindsy
phos_03 <- HPO2(ph_co2$PH[1], cacons[1])
##
##Veränderung des pH-Wertes#########
##
#Phosphorgehalt Abhängig von der CO2 Konzentration############
PH_PD <- seq(5, 9, 0.2) #CO2 Patialdrücke  in kPA 
#Berechung der Aktivitäten und pH-Werten anhand des Calciumsmodelles mit verschienen  CO2 Patialdrücken  
ph_co2_3 <- sapply(PH_PD, pH_ca_fun, temp1 = 20)
ph_co2_3 <- as.data.frame(t(ph_co2_3))
colnames(ph_co2_3) <- c("ca2", "CaHCO", "H", "HCO3", "CO2", "OH", "PH", "paco2")
cacons_2 <- numeric(length(PH_PD))
for (i in 1:length(PH_PD)){
  cacons_2[i] <- ca2(pco02=ph_co2_3$paco2[i], temp1=20, ph=ph_co2_3$PH[i], ga=ph_co2_3$ca2[i])
}
phos_2_pH <- HPO2(ph = ph_co2_3$PH, gammaca = cacons_2)
plot(PH_PD, log10(phos_2_pH), type = "l")
plot(ph_co2$PH, log10(phos), type = "l")

