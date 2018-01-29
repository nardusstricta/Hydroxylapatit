####Berechnung des Dissoziationsdiagramm auf Basis der Konzentrationen
#Gleichgewichtskonstanten:

source("rotation_fun.R") # Laden des Calciums-Modell
#Gleichgewichtskonstanten 
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

parts <- function(K){
  K <- 10 ^ (-K)
  H3 <- K^3 + (K^2 * K1) + (K * K1 * K2) + (K1 * K2 * K3)
  P1 <- K^3 * K1/H3
  P2 <- K^2 * K1 / H3
  P3 <- K * K1 * K2/H3
  P4 <- K1 * K2 * K3 / H3
  PF <- c(P1, P2, P3, P4)
  names(PF) <- c("H3PO4", "H2PO4", "HPO24", "PO34")
  return(PF)
}

###
#Funktion des gelösten Phosphatgehaltes im Gelichgewicht mit Hydroxylapatit###########
#bei gegeben pH und ca+ aktivität in mol/l
###
#Beispiel für ca2 Aktivität  und pH-Werte
ca_beispiel <- seq(0.0001, 0.009, 0.0001)
pH_beispiel <- seq(3, 10, 0.01)
temp_beispiel <- seq(0,25, 1)

#Idee: Löslichkeitsprodukt (Ksp) formolieren aus der Gleichgewichtsreaktion von Hydroxylapatit (Ca5(PO4)3OH + 7 H+ <-> 5 Ca2+ + 3H2PO4- + H2O); Ksp = a5 Ca2+ * a3 H2PO4- / a7 H+ (a = aktivität)
#Nach umformen und einsetzen von Ksp = 14.46 ->
#log(aH2PO4-) = 1/3 (14.46 - 7 *pH - 5 log aCa2+)
HPO2 <- function(ph, gammaca = 0.0025, temp_P = 20, g=T){
  erg <- 1/3 * (14.46 - (5 * log10(gammaca)) - (7 * ph))
  erg <- 10 ^ erg   #Umrechung von Mol/l zu Gramm/l
  ifelse(g == T, return(erg* 30.973762 * 1000), return(erg)) #In Miligramm *1000
}
##Funktion für Hpo2-4
HPO24 <- function(ph, gammaca = 0.0025, g= T){
  erg3 <- (1/3 * -7.21) - (5/3 * log10(gammaca)) - (4/3 * ph)
  erg3 <- 10 ^ erg3
  ifelse(g ==T, return(erg3 * 30.973762 *1000), return(erg3))
}
HPO24(5)
HPO24(6)
##
#Beispiel 1#####
##
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
mat_erg <- rotate(mat_erg) # rotation der Matrix 
filled.contour(pH_beispiel, ca_beispiel,   log10(mat_erg), ylab=expression("Ca"^"2+"*"Aktivität in mol l"^"-1"), xlab=("pH-Wert"), main  = "Phosphor in Lösung log(mg/l)", color = terrain.colors)

##
#Beispiel 2#####
##
erg1 <- sapply(pH_beispiel, HPO24, gammaca = 0.0025) #Bei gleichen Ca2+
plot(pH_beispiel, log10(erg1), type = "l", ylab = "Gelöster Phosphor log(mg/l)", xlab = " pH")
#Mit Veränderung der Ca2+ Aktivität
mat_erg <- as.data.frame(matrix(0,length(ca_beispiel), length(pH_beispiel)))
rownames(mat_erg) <- ca_beispiel
colnames(mat_erg) <- pH_beispiel
for (i in 1:length(ca_beispiel)){
  for (j in 1:length(pH_beispiel)){
    mat_erg[i,j] <- HPO24(ph = pH_beispiel[j],  gammaca = ca_beispiel[i])
  }
}
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

cacons<-numeric(length(CO2_PD)) #Berechnung der Ca2+ Konzentrationen aus den Aktivitäten 
for (i in 1:length(CO2_PD)){
  cacons[i] <- ca2(pco02=ph_co2$paco2[i], temp1=20, ph=ph_co2$PH[i], ga=ph_co2$ca2[i])
}

plot(ph_co2$paco2, cacons * 40 *1000, type = "l", ylab = expression("Ca" ^ "2+" * " Konzentration in mg" ^ "-1"), xlab = expression("CO" ["2"] * " Patialdruck kPa"))

#Berechung der Phosphatgehalte mit den ca2+ Konzentrationen und dem pH-Wert im Gleichgewicht
plot(ph_co2$PH, cacons, type = "l", xlab = "pH", ylab = expression("ca"^"2+"*"  mol l"^"-1"))

phos <- HPO2(ph = ph_co2$PH, gammaca = cacons, g = F)
xaxis <- (2*ph_co2$PH) + log10(cacons)
yaxis <- log10(phos) - ph_co2$PH
#Gleichung für die Kalkausfällung:
y <- log10(cacons) + 2 * ph_co2$PH
x <- 9.74 - log10(kpa_fun(CO2_PD))
plot(xaxis, yaxis , type = "l")
#Vergleichen mit den Angaben aus Lindsay:
ba <- 0.0303975082 #niedriger Patialdruck 
bb <- 1.0132502738 #hoher Patialdruck 
abline(v=x[c(1,10)]) #Linie mit den Angaben aus Lindsay
fm <- lm(yaxis ~ xaxis)
summary(fm) # slope ist -5/3 
gdata <- as.data.frame(cbind(xaxis, yaxis, x, CO2_PD, x))
names(gdata) <- c("xaxis", "yaxis", "x", "CO2_PD", "x_lin")
gdata$x_lin[-c(1,10)] <- NA 
library(ggplot2)
p <- ggplot(gdata, aes(xaxis, yaxis)) + geom_point(aes(color= CO2_PD)) 
p + geom_vline(aes(xintercept=x, color = CO2_PD)) +
  scale_colour_gradientn(colours = terrain.colors(10))+
  geom_vline(aes(xintercept=x_lin))

##Veränderung des pH-Wertes und der Temperatur
#pH Bereich ermitteln, durch die Eingabe von dem Patialdruck in das Calcit-Modell
carb_pH <- ph_co2$PH
#Berechung der Aktivitäten und pH-Werten anhand des Calciumsmodelles mit verschienen  CO2 Patialdrücken  
mat_res <- matrix(0,length(carb_pH), length(temp_beispiel))
mat_calcit <- matrix(0,length(carb_pH), length(temp_beispiel))
for(j in 1:length(temp_beispiel)){
      ph_co2_3 <- sapply(carb_pH, pH_ca_fun, temp1 = temp_beispiel[j])
      ph_co2_3 <- as.data.frame(t(ph_co2_3))
      colnames(ph_co2_3) <- c("ca2", "CaHCO", "H", "HCO3", "CO2", "OH", "PH", "paco2")
        for (i in 1:length(carb_pH)){
          mat_calcit[,j]<- ca2(pco02=ph_co2_3$paco2[i], temp1=temp_beispiel[j], ph=ph_co2_3$PH[i], ga=ph_co2_3$ca2[i])
        }
      mat_res[,j] <- HPO2(ph = carb_pH, gammaca = mat_calcit[,j])
}
plot(temp_beispiel, mat_calcit[1,], xlab = "Temperatur", ylab = "Ca Aktität")
matplot(carb_pH, log10(mat_res), type ="l")
phos_2_pH <- HPO2(ph = carb_pH, gammaca = cacons_2)
phos_23_pH <- HPO2(ph = carb_pH, gammaca = mat_res[,which(temp_beispiel==20)])
phos_co2_patial <- HPO2(ph = ph_co2$PH, gammaca = cacons)
plot(carb_pH, log10(phos_23_pH), type = "l")
plot(ph_co2$PH, log10(phos_co2_patial), type = "l") #Ist Gleich auch wenn der pH-Wert unbekannt ist. 
sumcol <- apply(mat_res,2,sum)
plot(temp_beispiel, sumcol, ylab = "Summe des gelösten Phosphat (CO2 variabel)")

