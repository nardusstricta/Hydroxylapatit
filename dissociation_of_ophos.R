#1. Beispiele für die Abhängigkeit des Phosporgehaltes von pH und Ca2-Aktivität in der Bodenlösung
#2. Berechnung der Ca+ Aktivität in dem pH-Bereich, welcher vom CO2 Patialdruck bestimmt wird, mit unterschiedlchen Temperaturen. 

source("functions.R") #Laden des Calciums-Modell und der benötigten Funktionen für die P-Konzentration 


##
#Beispiel Werte für Ca pH 
##

ca_beispiel <- seq(1e-6, 0.01, length.out = 60)
pH_beispiel <- seq(5, 10, length.out = 100)
pH_beispiel2 <- seq(1, 14, length.out = 100)

##
#Abbildung der Dissoziationen der P-Säure:
##

pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/spezi.pdf", width= 6, height=6)
plot(pH_beispiel2, parts(pH_beispiel2)$H3PO4, type="l", ylim=c(0,1), ylab = "Konzentrationsanteile", main = "Dissoziationsdiagramm der Phosphorsäure")
rect(7.004215, 1, 8.248673, 0, border = NA, col= rgb(0.220,0.220,0.220,alpha=0.06))
points(pH_beispiel2, parts(pH_beispiel2)$H2PO4 , type = "l", col="red")
points(pH_beispiel2, parts(pH_beispiel2)$HPO24 , type = "l", col ="blue", lty=2)
points(pH_beispiel2, parts(pH_beispiel2)$PO34, type = "l", col = "green", lty=2)
legend("topleft", c(expression("H"["3"]*"PO"["4"]), expression("H"["2"]*"PO"["4"]^"-"), expression("HPO"["4"]^"2-"), expression("PO"["4"]^"3-"), expression("CaCo"[3]*"-H"[2]*"O-CO"[2]), "System"), col = c("black", "red", "blue", "green", "grey", NA), lty = c(1,1,2,2,1, NA), bty = "n", lwd=c(1,1,1,1,5, NA))
dev.off()

##
#Beispiel Phosphat in Abhängigkeit von pH und Ca Konzentration#####
##

mat_erg1 <- outer(pH_beispiel2, ca_beispiel, HPO2, g=F) #Kombination aus allen Beispiel pH und ca-Werten 

##
#Abbildung pH-Wert gegen Ca-Konzentration:
##

pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/ca_imag.pdf", width= 7, height=5)
filled.contour(pH_beispiel2, ca_beispiel,  log10(mat_erg1), ylab=expression("Ca"^"2+"*"Aktivität in mol l"^"-1"), xlab=("pH-Wert"), main  = "Phosphor in Lösung log(mol/l)", color = terrain.colors)
dev.off()


##
#2. Phosphorgehalt Abhängig vom CO2-Patialdruck####
##

#CO2 Patialdrücke  in kPA Welche für den Boden realistisch sind (atm ist minimum):
CO2_PD <- rep(seq(0.039, 1.5, 0.01), each = 25) #Beispiel Co2-Patialdrücke
temp_beispiel <- rep(seq(1, 25, 1), times = length(seq(0.039, 1.5, 0.01))) #Beispiel Temperaturen


#Anwenden der Funktion (ca_con_s: die mit gegebener Temperatur und CO2-Patialdruck die
#Ca2-Konzentration über das CaCo3-H2O- Modell errechnen kann) auf den Beispiel
# Datensatz (PCo2 = 0.039 - 1.5 kPA; Temperatur = 0-25 Grad Celsius):

mat_calcit <- apply(data.frame(CO2_PD, temp_beispiel), 1, FUN = ca_con_s)
mat_calcit <- as.data.frame(t(mat_calcit))

#Neuer Datensatz "mat_calcit" mit der Ca2-Konzentration, 
#Patialdruck, pH, und der Temperatur

names(mat_calcit) <- c("ca", "pco2", "pH", "temp_b")


#Abbildungen des CaCo3-H2o-Co2 Systemes in Abhängigkeit der Temperatur:

ggplot(mat_calcit, aes(x=temp_b, y=pH))+
  geom_point(aes(color=ca))


#Speichern der Abbildung:

pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/pco2.pdf", width= 7, height=5)
ggplot(mat_calcit, aes(x=pco2, y=ca*40000))+
  geom_point(aes(color=temp_b, size=pH, alpha = pH))+
  labs(y = "Ca2 [mg/l]", x = "CO2 [kPA]", color = "Temperatur\n [C]")+
  scale_colour_gradient(low = "darkblue", high = "darkred")
dev.off()

##
#Modelle####
###

##
#Modell für den pH-Wert
##

fm1 <- lm(pH ~ (I(temp_b^2) + I(temp_b^3) + I(ca^2) +I(ca^3)+ I(ca^4)+ temp_b + ca)^2, data = mat_calcit)


mat_calcit$mod_pH <- predict(fm1, newdata = data.frame(temp_b = mat_calcit$temp_b, ca = mat_calcit$ca))

#Darstellung mit ausgewählten Temperaturen und einem pH-Bereich von 7.25 - 7.3
new <- subset(mat_calcit, (temp_b %in% c(1, 5, 15, 25)) & pH >7.25 & pH < 7.3)

ggplot(new, aes(x=ca, y=pH ))+
  geom_point(aes(color=as.factor(temp_b)))+
  geom_line(aes(ca, mod_pH, color = as.factor(temp_b)))



##
#Modell für den Co2-Patialdruck
##

fm2 <- lm(pco2 ~ (I(temp_b^2) +  I(ca^2)  + I(pH^2) + temp_b + ca + pH)^2, data = mat_calcit)

mat_calcit$mod_pco2 <- predict(fm2, newdata = data.frame(temp_b = mat_calcit$temp_b, ca = mat_calcit$ca, pH = mat_calcit$pH))



##
#Funktion um die pH-Werte bzw. den co2_Patieldruck zu predicten###
##

#pH
pred_ph <- function(ca_I, Temp){
  preds_pH <- predict(fm1, newdata = data.frame(temp_b = Temp, ca = ca_I))
  return(preds_pH)
}

#Test von pH-Modell:

#Max
test <- pred_ph(mat_calcit$ca[which(mat_calcit$ca==min(mat_calcit$ca))], mat_calcit$temp_b[which(mat_calcit$ca==min(mat_calcit$ca))])
test 
mat_calcit$pH[which(mat_calcit$ca==min(mat_calcit$ca))]

#Min
test <- pred_ph(mat_calcit$ca[which(mat_calcit$ca==max(mat_calcit$ca))], mat_calcit$temp_b[which(mat_calcit$ca==max(mat_calcit$ca))])
test 
mat_calcit$pH[which(mat_calcit$ca==max(mat_calcit$ca))]


#CO2
pred_co2 <- function(ca_I, Temp, pH){
  preds_co2 <- predict(fm2, newdata = data.frame(temp_b = Temp, ca = ca_I, pH = pH))
  return(preds_co2)
}

#Test:
test <- pred_co2(mat_calcit$ca[10], mat_calcit$temp_b[10], mat_calcit$pH[10]) 
test 
mat_calcit$pco2[10]


##
#Berechnung der Phosphatwerte in Abhängigkeit der Temperatur####
##

mat_calcit$P_ml  <- HPO2(mat_calcit$mod_pH, mat_calcit$ca)
mat_calcit$P2_ml <- HPO24(mat_calcit$mod_pH, mat_calcit$ca)
mat_calcit$brush <- brush(mat_calcit$mod_pH, mat_calcit$ca)

#Anteile berechnen (in Abhägigkeit des pH-Wertes (Dissoziationsdiagramm))

#Ausschnit aus dem Dissoziationsdiagramm im relevanten pH-Bereich:
plot(mat_calcit$pH, parts(pH=mat_calcit$pH)$H2PO4, type="l", col="red", ylim = c(0,1))
points(mat_calcit$pH, parts(pH=mat_calcit$pH)$HPO24, type="l")
legend("topright", c( expression("H"["2"]*"PO"["4"]^"-"), expression("HPO"["4"]^"2-")), col = c("red", "black"), lty = c(1,1), bty = "n")

#Anteile berechnen:
mat_calcit$P_comp <- P_comp(pH = mat_calcit$pH, P_4 = mat_calcit$P_ml, P_24 = mat_calcit$P2_ml)

#Abbildung der Temperatur und P-Konzentrationen mit Ca als Hilfsvariablen 
ggplot(mat_calcit, aes(x=temp_b, y=log10(P_comp)))+
  geom_point(aes(color=ca))

##
#Abbildung zur überprüfung der Modelle:
## 

#Ca gegen pco2 mit Temperatur als Hilfsvariablen:

new <- subset(mat_calcit, (temp_b %in% c(1, 5, 15, 25)))
ggplot(new, aes(x=ca, y=pco2))+
  geom_point(aes(color=as.factor(temp_b)))+
  geom_line(aes(ca, mod_pco2, color = as.factor(temp_b)))

new_pH <- subset(mat_calcit, (pco2 %in% 1.009)) #Fixierung von pCO2
pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/neu.pdf", width= 7, height=5)
ggplot(new_pH, aes(x=ca, y=log10(P_comp)))+
  geom_point(aes(color=temp_b, size=pH))+
  labs(y = "P (HAP)[mol/l]", x = "Ca2+ [mol/l]", color = "Temperatur\n [C]")
dev.off()

new_pH2 <- subset(mat_calcit, (pH > 7.239 & pH < 7.24)) #Fixierung von pCO2
pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/neu2.pdf", width= 7, height=5)
ggplot(new_pH2, aes(x=ca, y=log10(P_comp)))+
  geom_point(aes(color=temp_b, size=pco2))+
  labs(y = "P (HAP)[mol/l]", x = "Ca2+ [mol/l]", color = "Temperatur\n [C]")
dev.off()

##
#Abbildung nach Lindsay####
##

CO2_PD <- seq(0.1, 5, 0.1) #CO2 Patialdrücke  in kPA Welche für den Boden realistisch sind
ph_co2 <-  ca_con(CO2_PD, temp_b = 25)
phos <- HPO2(ph = ph_co2$pH, gammaca = ph_co2$ca, g = F)

xaxis <- (2*ph_co2$pH) + log10(ph_co2$ca) 
yaxis <- log10(phos) - ph_co2$pH

#Gleichung für die Kalkausfällung nach Lindsay:
y <- log10(ph_co2$ca) + 2 * ph_co2$pH
x <- 9.74 - log10(kpa_fun(ph_co2$pco2))

#Vergleichen mit den Angaben aus Lindsay:
ba <- 0.0303975082 #niedriger Patialdruck 
bb <- 1.0132502738 #hoher Patialdruck 

abline(v=x[c(1,10)]) #Linie mit den Angaben aus Lindsay
fm <- lm(yaxis ~ xaxis)
summary(fm) # slope is -5/3, so wie in Gleichung nach Lindsay angegeben ist

gdata <- as.data.frame(cbind(xaxis, yaxis, x, CO2_PD, x))
names(gdata) <- c("xaxis", "yaxis", "x", "CO2_PD", "x_lin")
gdata$x_lin[-c(1,10)] <- NA 

#Abbildung
pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/gglindsay.pdf", width= 7, height=5)
p <- ggplot(gdata, aes(xaxis, yaxis)) + geom_point(aes(color= CO2_PD)) 
p + geom_vline(aes(xintercept=x, color = CO2_PD, alpha=0.1)) +
  scale_colour_gradientn(colours = terrain.colors(10))+
  geom_vline(aes(xintercept=x_lin))+
  labs(title = "Abbildung nach Lindsay\n", x = "log Ca2 + 2pH", y = "log H2PO4 - pH", color = "Patialdruck\n CO2\n", alpha = "Ausfällung \n von CaCo3\n") 
dev.off()


