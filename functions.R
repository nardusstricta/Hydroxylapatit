#Laden von Paketen:
list.of.packages <- c("reshape2", "ggplot2", "readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(reshape2)
library(ggplot2)
library(readr)

#Funktion für die Berechnung der Phosphorkonzentration in Abhängigkeit von der ca-Aktivität und dem pH-Wert (Gleichgewicht mit HAP und Brushit)

#Idee: Löslichkeitsprodukt (Ksp) formolieren aus der Gleichgewichtsreaktion von Hydroxylapatit (Ca5(PO4)3OH + 7 H+ <-> 5 Ca2+ + 3H2PO4- + H2O); Ksp = a5 Ca2+ * a3 H2PO4- / a7 H+ (a = aktivität)
#Nach umformen und einsetzen von Ksp = 14.46 ->
#log(aH2PO4-) = 1/3 (14.46 - 7 *pH - 5 log aCa2+)

HPO2 <- function(ph, gammaca = 0.0025, temp_P = 20, g=F){
  erg <- 1/3 * (14.46 - (5 * log10(gammaca)) - (7 * ph))
  erg <- 10 ^ erg   
  ifelse(g == T, return(erg* 30.973762 * 1000), return(erg)) #In Miligramm *1000
}



#Hydroxylapatit HpO24:
HPO24 <- function(ph, gammaca = 0.0025, g= F){
  erg3 <- (1/3 * -7.21) - (5/3 * log10(gammaca)) - (4/3 * ph)
  erg3 <- 10 ^ erg3
  ifelse(g ==T, return(erg3 * 30.973762 *1000), return(erg3))
}

#Dicalciumphosphat Dihydrat

brush <- function(ph, gammaca, temp_P,  g=F){
  erg4 <- 0.63 + ph - (log10(gammaca) + 2*ph)
  erg4 <- 10 ^erg4
  ifelse(g ==T, return(erg4 * 30.973762 *1000), return(erg4))
}




#Berechnung der Anteile der Phosphatspezies in Abhängigkeit des pH-Wertes
parts <- function(pH){
  #Gleichgewichtskonstanten 
  pk1 <- 2.16
  pk2 <- 7.21
  pk3 <- 12.32
  #pks-Werte (-log10 von Gleichgewichtskonstante)
  h <- 10 ^ (-pH)
  K1 <- 10 ^ (-pk1)
  K2 <- 10 ^ (-pk2)
  K3 <- 10 ^ (-pk3)
  
  #Gesamtkonzentration H3 ermittelt über die Konzentrationssumme 
 
  H3 <- h^3 + (h^2 * K1) + (h * K1 * K2) + (K1 * K2 * K3)
  
  P1 <- h^3 /H3
  P2 <- h^2 * K1 / H3
  P3 <- h * K1 * K2/H3
  P4 <- K1 * K2 * K3 / H3
  PF <- data.frame(P1, P2, P3, P4)
  names(PF) <- c("H3PO4", "H2PO4", "HPO24", "PO34")
  return(PF)
}


###
#Funktion die bei gegebenen pH-Wert und Phosphatwerte für
# H2PO4 und HPO24 gesammte P-Konzentration ausgibt:
###

P_comp <- function(pH, P_4, P_24){
  part4 <- parts(pH)$H2PO4
  part24 <- parts(pH)$HPO24
  erg_comp <- (part4 * P_4) + (part24 * P_24)
  return(erg_comp)
}


#Umrechnen von kPA zu atm
kpa_fun <- function(kpa){
  kpa/101.3
}

#Temperatur Abhängigkeit der Lösung von CO2 im Wasser
#Für die Temperatur wird der Wert temp1 eingegeben, welcher dann über die Funktion verändert werden kann:

##
#Carroll Funktiton
##

carl <- function(temp=temp1){
  tempK <- 273.15 + temp # in Grad Celsius
  newH <- exp((-6.8346) + ((1.2817*10^4)/tempK)-
                ((3.7668*10^6)/tempK^2) + ((2.997 * 10^8)/(tempK^3))) 
  newkpa <- 1/newH   #Umrechen von Mg Pa auf Kg Pa
  newkpa <- newkpa/1000 
  newkpa <- newkpa*(998/18)
  return(newkpa)
}

##
#Dissoziationstufen der 6 Carbonat Spezies
##
#Ca2
ca2 <- function(pco02, temp1, ph, ga){
  patial_co24 <- (10^(8.28 - log10(carl(temp=temp1)))*((10^-ph)^2))/(pco02*ga)
  return(patial_co24)
}

#cahco3
cahco3 <- function(pco02, temp1, ph, ga){
  patial_co25 <- ((10^3.02)*(10^-ph))/ga
  return(patial_co25)
}

#H+
H <- function(pco02, temp1, ph, ga){
  patial_co21<-(10^-ph)/ga
  return(patial_co21)
}

#HCo3
HCo3 <- function(pco02, temp1, ph, ga){
  patial_co22 <- (10^(-6.36+log10(carl(temp=temp1)))*pco02)/(ga*(10^-ph))
  return(patial_co22)
}

#co23
co23 <- function(pco02, temp1, ph, ga){
  patial_co23 <- (10^(-16.69+log10(carl(temp=temp1)))*pco02)/(ga*((10^-ph)^2))
  return(patial_co23)
}

#oh
oh <- function(pco02, temp1, ph, ga){
  patial_co26 <- (10^-14)/((10^-ph)*ga)
  return(patial_co26)
}

#Elektroneutralitätsformel als Funktion (Dividiert durch die jeweiligen Aktivit?tskoeffizienten usw.)
sumCA <- function(pco02, temp1, ph, ga){
  sumCAer <- (2*ca2(pco02=pco02, temp1=temp1, ph=ph, ga=ga[1]))+
    cahco3(pco02=pco02, temp1=temp1, ph=ph, ga=ga[2])+
    H(pco02=pco02, temp1=temp1, ph=ph, ga=ga[3])-
    HCo3(pco02=pco02, temp1=temp1, ph=ph, ga=ga[4])-
    (2*co23(pco02=pco02, temp1=temp1, ph=ph, ga=ga[5]))-
    oh(pco02=pco02, temp1=temp1, ph=ph, ga=ga[6])
  return(sumCAer)
}


#Ionenst?rke errechen um neuen  Aktivit?tskoeffizienten Gamma zu bekommen--------
#1. Tempereaturabhängigkeit der Parameter A und B errechen mit folgender Funktion:
#A
temp_A <- function(Ax){
  par_a <- c(0.4883, 0.5002, 0.5046, 0.5092, 0.5141, 0.5241, 0.5351, 0.5471, 0.5739)
  par_temp <- c(0,15,20,25,30,40,50,60,80)
  fml <- lm(par_a~par_temp)
  A_ret <- fml$coefficients[2]* Ax + fml$coefficients[1]
  return(A_ret)
}

#B
temp_B <- function(Bx){
  par_b <- c(3.241^10, 3.267^10, 3.276^10, 3.286 ^10, 3.297^10, 3.318^10,3.341^10,3.366^10, 3.420^10)
par_temp <- c(0,15,20,25,30,40,50,60,80)
fmb <- lm(par_b ~ par_temp)
B_ret <- fmb$coefficients[2]*Bx + fmb$coefficients[1]
return(B_ret)
}

#2.Funktion der Ionenstärke 
phmu<-function(pco02=pco02, temp1=temp1, ph=ph, ga=ga){
  #Ionenstärke mu 
  mu <- 0.5*(ca2(pco02=pco02, temp1=temp1, ph=ph,ga=ga[1])*
               4+cahco3(pco02=pco02, temp1 =temp1, ph=ph,ga=ga[2])+
               H(pco02=pco02, temp1=temp1, ph=ph,ga=ga[3])+
               HCo3(pco02=pco02, temp1=temp1, ph=ph,ga=ga[4])+
               co23(pco02=pco02, temp1=temp1, ph=ph,ga=ga[5])*
               4+oh(pco02=pco02, temp1=temp1, ph=ph,ga=ga[6]))
  #temperaturabhängige Koeffizienten usw.
  A<-temp_A(temp1)
  B<-temp_B(temp1)
  d<-c(6e-8, 3e-8, 9e-8, 4e-8, 4.5e-8, 3.5e-8)
  z<-c(2,1,1,1,2,1)
  #Vektor für Gamma erstellen und durch Debye-H?ckel-Formel neu errechen
  ga<-numeric(6) 
  for(i in 1:6){
    ga[i] <- 10^(((-A*z[i]^2) * sqrt(mu))/(1+B * d[i] * (sqrt(mu))))
  }
  return(ga) 
}

#Gleichgewicht lösen mit der Temperatur als Variable und dem CO2 Patialdruck (co2pa_ca)
#Wir haben 6 unbekannte Spezies, wir kennen PCO2 und kennen 7 Zusammenhänge 
#-> Iteration

co2pa_ca <- function(pco02, temp1){#Funktion braucht Patialdurck und Temperatur
  ga <- rep(1,6) #Gamma anfangs immer auf 1 und PH auf 7
  ph <- 7
  repeat{
    ph1 <- ph
    Int <- 7
    repeat{
      er <- sumCA(pco02=pco02, temp1=temp1, ph=ph, ga=ga)
      if(abs(er) < 1e-8|Int==0){
        break
      }else{
        Int <- Int/2
        if(er > 0){
          ph<-ph+Int
        }else{
          ph <- ph - Int}
      }
    }
    if(abs(ph1-ph)<1/10000){#für schnellere Durchläufe reichen auch 1000
      return(c(ga,ph, pco02)) #Bevor auch die zweite repeat Schleife
      #abricht werden 6*Gammas, PH, und der CO2-Patialdruck ausgeben
      break
    }else{
      ga <- phmu(pco02=pco02, temp1=temp1, ph=ph, ga=ga)#neue Gammas abrufen
    }
  }
}


#Funktion die mit gegebener Temperatur und CO2-Patialdruck die Ca2-Konzentration über das CaCo3-H2O- Modell errechnen kann:

#X-Variante:
ca_con_s <- function(x){
  CO2_imp <- x[1]
  temp_b <- x[2]
  zwerg <- co2pa_ca(pco02=CO2_imp, temp1=temp_b)
  erg <- ca2(pco02=zwerg[8], temp1=temp_b, ph=zwerg[7], ga=zwerg[1])
  return(c(erg, zwerg[8], zwerg[7], temp_b))
}

#X, Y Variante:
ca_con <- function(CO2_PD, temp_b){
  mat_com1 <- as.data.frame(t(sapply(CO2_PD, co2pa_ca, temp = temp_b)))
  colnames(mat_com1) <- c("ca2", "CaHCO", "H", "HCO3", "CO2", "OH", "PH", "paco2")
  #Berechnung der Ca2+ Konzentrationen aus den Aktivitäten 
  cacons <- mapply(ca2, mat_com1$paco2, temp_b, mat_com1$PH, mat_com1$ca2)
  erg <- data.frame(cacons, mat_com1$paco2, mat_com1$PH)
  names(erg) <- c("ca", "pco2", "pH")
  return(erg)
}