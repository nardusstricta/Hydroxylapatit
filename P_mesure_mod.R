##
#gemessene Phosphatwerte Modellieren####
##

source("rotation_fun.R") #Laden des Calciums-Modell und der benötigten Funktionen für die P-Konzentration 
source("dissociation_of_ophos.R") #Laden der Funktion die bei gegebendem Ca-Wert den pH- und Co2 Patialdruck vorhersagt



#Laden der Messwerte:
ICP_Calcium <- read_csv("ICP Calcium.csv")
ICP_Calcium$camoll <- ICP_Calcium$Cammoll/1e+6

##
#Berechnung des pH-Wertes und P-Konzentration:####
##

ICP_Calcium$pH <-  pred_ph(ICP_Calcium$camoll, ICP_Calcium$TemperaturC) #pH 


#Phosphor-Konzentrationen:
ICP_Calcium$P_mod <- HPO2(ph = ICP_Calcium$pH, gammaca = ICP_Calcium$camoll)
ICP_Calcium$P_mod2 <- HPO24(ph = ICP_Calcium$pH, gammaca = ICP_Calcium$camoll)

#Anteile mit einander verrechnen pH-abbhängig:
ICP_Calcium$P_comp <- P_comp(pH=ICP_Calcium$pH, P_4 = ICP_Calcium$P_mod, P_24 = ICP_Calcium$P_mod2)
ICP_Calcium$brushite_mod <- brush(ph = ICP_Calcium$pH, gammaca = ICP_Calcium$camoll)

#Berechung des Co2-Patieldruckes aus dem pH-Wert und ca-Aktivität und der Temperatur:
ICP_Calcium$pco2  <- pred_co2(ca_I = ICP_Calcium$camoll, Temp = ICP_Calcium$TemperaturC, pH = ICP_Calcium$pH)

##
#Abbildung der Modelierten und der gemessen Phosphorkonzentrationen:####
##

melt_mess <- as.data.frame(cbind(ICP_Calcium$TemperaturC, ICP_Calcium$camoll, ICP_Calcium$Mittelwerteh2oExtraktePO4mmoll/1e+6, ICP_Calcium$P_mod, ICP_Calcium$P_mod2, ICP_Calcium$brushite_mod, ICP_Calcium$P_comp)) 

names(melt_mess) <- c("temp_b", "ca", "P_mesure","HAP_mod", "HAP2_mod", "Brushite_mod", "HAP_comp")
melt_mess <- melt(melt_mess, id=c("temp_b", "ca"))

pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/mess_erg.pdf", width= 7, height=5)
pp <- ggplot(melt_mess, aes(temp_b, log10(value)))
pp  + geom_point(aes(color=variable, size=ca)) +
  labs(x = "Temperatur C", y = expression("log HPO"[4]*"- bzw. H2PO"[4]^"2-"), color = "Mineral")
dev.off()

#Abbildung der Ca-Konzentrationen in Abhängigkeit der Temperatur im CaCo3-System:
help_temp <- as.data.frame(cbind(ICP_Calcium$TemperaturC, ICP_Calcium$pco2, ICP_Calcium$camoll))
names(help_temp) <- c ("temp_b", "pco2", "ca")


pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/ca.pdf", width= 7, height=5)
ggplot(mat_calcit, aes(x=temp_b, y=pco2, z=ca))+
  geom_raster(aes(fill = ca))+
  geom_contour(color="white")+
  geom_point(data=help_temp, aes(color="red"))+
  labs(x = "Temperatur [C]", y = expression("CO"[2]*"  Patialdruck  [kPA]"), color ="Messwerte", fill = expression("Ca"^"2+"*"\n[mol l" ^"-1"*"]"))
dev.off()

#Modellierung der HAP-Phosphatkonzentrationen in Abhängigkeit der Temperatur im pH-Bereich vom CaCo3-System 

fm3 <- lm(P_comp  ~ (pH + temp_b + ca + I(ca^2)  +I(temp_b ^2) + I(pH^2))^2, data=mat_calcit)

#HAP Als Funktion der Ca-Konzentration, Temperatur und pH-Wertes
pred_Po <- function(ca_I, Temp, pH){
  preds_Po <- predict(fm3, newdata = data.frame(temp_b = Temp, ca = ca_I, pH = pH))
  return(preds_Po)
}
#Zufügen zu Dataframe:
mat_calcit$P_mod <- pred_Po(mat_calcit$ca, mat_calcit$temp_b, mat_calcit$pH)
#Abbildung im Bereich zwischen pH 7.25 und pH 7.3 mit Temperaturen von 1 5 15 und 25
new1 <- subset(mat_calcit, (temp_b %in%  c(1, 5, 15, 25) & pH >7.25 &pH< 7.3))
pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/mineral_imag2.pdf", width= 7, height=5)
ggplot(new1, aes(x=pH, y= log10(P_comp)))+
  geom_point(aes(color= as.factor(temp_b), size = ca)) +
  geom_line(aes(x=pH, y=log10(P_mod), color=as.factor(temp_b)))+
  labs(x = "pH", y = expression("log HPO"[4]*"- bzw. H2PO"[4]^"2-"* "   [mol l "^"-1"*"]"), color = "Temp\n[C]", size = expression("Ca"^"2+"*"[mol"^"-1"*"]"))
dev.off()



###
#Abbildung der P-Konzentration von HAP und Brushite über  den gesammten pH-Bereich un
#unter Berücksichtigung des CaCO3-H20-CO2 Systemes:
###

#Beispiel für ca2 Aktivität  und pH-Werte erstellen:
ca_beispiel <- seq(1e-6, 0.01, length.out = 60)
pH_beispiel <- seq(5, 10, length.out = 100)

var <- expand.grid(ca_beispiel, pH_beispiel) #Vektor der für jeden ca_Wert einen pH_Wert hat

help_pH <- var$Var2[which(var$Var2 > max(mat_calcit$pH))] #extrahieren von den pH-Werten, welche höher sind als die welche in Abhängigkeit vom Patialdruck  im CaCo3-System herrschen. 
help_pco2 <- seq(min(mat_calcit$pco2), max(mat_calcit$pco2), length.out = 60) #Erstellen von Beispiel Co2-Patialdrücken um die ca-Werte zu modifizieren
help_mat <- as.data.frame(cbind(ca2(pco02 = help_pco2, temp1 =25, ph= help_pH, ga=1), help_pH)) #Berechnung von neuen Ca-Konzentrationen mit gegebenen Co2-Partialdruck und Beispiel pH-Wert 
names(help_mat) <- names(var)

#1.Löschen von allen Werten mit hohen pH-Werten   (> als CaCo3- System)
#2.Ersetzen durch die neu errechenten Ca2-Konzentratien
var[which(var$Var2 > max(mat_calcit$pH)),] <- NA 
var<- rbind(var, help_mat) 

#Berechnung der P-Konzentratioen:
var$HAP <- log10(HPO2(var$Var2, var$Var1)) 
var$HAP2 <- log10(HPO24(var$Var2, var$Var1))
var$HAP_comp <- P_comp(var$Var2, var$HAP, var$HAP2) #Anteile von HAP-Spezies
var$Brushite  <- log10(brush(var$Var2, var$Var1)) #Brushite ohne Berücksichtigung von den Ionenspezies der Phosphatsäure

#Erstellung einer longtable:
melt_var <- melt(var, id=c("Var1", "Var2", "HAP", "HAP2"))
names(melt_var) <- c("ca2", "pH", "HAP", "HAP2",  "Mineral", "P_comp")

#Hilfsmatrix mit den Werten aus der Berechung des CaCo3-Systems
help_calit <- as.data.frame(cbind(mat_calcit$ca, mat_calcit$pH, log10(mat_calcit$P_comp), mat_calcit$temp_b)) 

names(help_calit) <- c("ca2", "pH","P_comp", "temp_b")

#Hilfsmatrix für die Phosphorwerte mit höheren pH-Werten als das CaCo3-System unter der Annahme des kleinsten Co2-Patialdruckes 0.3-kPa (Da dieser den höchsten pH-Wert bestimmt im CaCo3-System )

help_calit_long  <- matrix(NA, length(seq(max(help_calit$pH), 10, by=0.1)), 3)
help_calit_long  <- as.data.frame(help_calit_long )
names(help_calit_long) <- c("ca2", "pH","P_comp")
help_calit_long$pH <- seq(max(help_calit$pH), 10, by=0.1)


#Berechung der Ca-Konzentration mit unterschiedlichen pH-Werten:
help_calit_long$ca2 <- ca2(pco02 = mat_calcit$pco2[which(mat_calcit$pH==max(mat_calcit$pH))], temp1 =25, ph= help_calit_long$pH, ga=1)


#Berechnung der Phosphatwerte
HP1 <- HPO2(ph = help_calit_long$pH, help_calit_long$ca2)
HP2 <- HPO24(ph = help_calit_long$pH, help_calit_long$ca2)
help_calit_long$P_comp <- log10(P_comp(help_calit_long$pH, HP1, HP2))

#Das Gleiche für den Brushite
help_calit_long_b <- help_calit_long
help_calit_long_b$P_comp <- log10(brush(help_calit_long_b$pH, help_calit_long_b$ca2))

#Löschen der Modelierten Phosphorkonzentrationen im Bereich des pH, welcher durch das CaCo3-System bestimmt wird. 
melt_var$P_comp[which(melt_var$pH > min(mat_calcit$pH) & melt_var$pH < max(mat_calcit$pH))] <- NA

#Ergebnisabbildung:
pdf("/home/gabriel/Dokumente/Uni_Master/Aktuelles_Thema/latex2/Bilder/mineral_imag.pdf", width= 7, height=5)
pp1 <- ggplot(melt_var, aes(x=pH, y=P_comp))+
  geom_point(aes(color = Mineral, alpha= ca2))+
  geom_line(data = help_calit, color ="red")+
  geom_line(data = help_calit_long_b, linetype = "dashed", color = "blue")+
  geom_line(data= mat_calcit, aes(x=pH, y=log10(brush)), color = "blue")+
  geom_line(data = help_calit_long, linetype = "dashed", color = "red")+
  labs(x = "pH", y = expression("log HPO"[4]*"- bzw. H2PO"[4]^"2-"* "   [mol l "^"-1"*"]"))
pp1
dev.off()


