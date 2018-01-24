matx <- matrix(1:12, 3, 4)
matx[3,2] <- 19
matx
image(matx)
rotate <- function(x) t(apply(x, 2, rev))
imag <- function(m){
  imag
}
image(rotate(matx))

#Carroll Funktiton------------------
carl <- function(temp=20){
  tempK <- 273.15+temp # in Grad Celsius
  newH <- exp((-6.8346)+((1.2817*10^4)/tempK)-
                ((3.7668*10^6)/tempK^2)+((2.997*10^8)/(tempK^3))) 
  newkpa <- 1/newH   #Umrechen von Mg Pa auf Kg Pa
  newkpa <- newkpa/1000 
  newkpa <- newkpa*(998/18)
  return(newkpa)
}
#Dissoziationstufen der 6 Carbonat Spezies-------------------
#1
ca2 <- function(pco02, temp1, ph, ga){
  patial_co24 <- (10^(8.28 - log10(carl(temp=temp1)))*((10^-ph)^2))/(pco02*ga)
  return(patial_co24)
}
#2
cahco3 <- function(pco02, temp1, ph, ga){
  patial_co25 <- ((10^3.02)*(10^-ph))/ga
  return(patial_co25)
}
#3
H <- function(pco02, temp1, ph, ga){
  patial_co21<-(10^-ph)/ga
  return(patial_co21)
}
#4
HCo3 <- function(pco02, temp1, ph, ga){
  patial_co22 <- (10^(-6.36+log10(carl(temp=temp1)))*pco02)/(ga*(10^-ph))
  return(patial_co22)
}
#5
co23 <- function(pco02, temp1, ph, ga){
  patial_co23 <- (10^(-16.69+log10(carl(temp=temp1)))*pco02)/(ga*((10^-ph)^2))
  return(patial_co23)
}
#6
oh <- function(pco02, temp1, ph, ga){
  patial_co26 <- (10^-14)/((10^-ph)*ga)
  return(patial_co26)
}
#Elektroneutralit?tsformel als Funktion (Dividiert durch die jeweiligen Aktivit?tskoeffizienten usw.)------------------
sumCA <- function(pco02, temp1, ph, ga){
  sumCAer <- (2*ca2(pco02=pco02, temp1=temp1, ph=ph, ga=ga[1]))+
    cahco3(pco02=pco02, temp1=temp1, ph=ph, ga=ga[2])+
    H(pco02=pco02, temp1=temp1, ph=ph, ga=ga[3])-
    HCo3(pco02=pco02, temp1=temp1, ph=ph, ga=ga[4])-
    (2*co23(pco02=pco02, temp1=temp1, ph=ph, ga=ga[5]))-
    oh(pco02=pco02, temp1=temp1, ph=ph, ga=ga[6])
  return(sumCAer)
}
#Beispiel
ga <- rep(1, 6)#Gamma wird als Vektor ?bergeben
sumCA(pco02=1, temp=20, ph=7, ga=ga) #0.008107546

#Ionenst?rke errechen um neuen  Aktivit?tskoeffizienten Gamma zu bekommen--------
phmu<-function(pco02=pco02, temp1=temp1, ph=ph, ga=ga){
  #Ionenstärke mu 
  mu <- 0.5*(ca2(pco02=pco02, temp1=temp1, ph=ph,ga=ga[1])*
               4+cahco3(pco02=pco02, temp1 =temp1, ph=ph,ga=ga[2])+
               H(pco02=pco02, temp1=temp1, ph=ph,ga=ga[3])+
               HCo3(pco02=pco02, temp1=temp1, ph=ph,ga=ga[4])+
               co23(pco02=pco02, temp1=temp1, ph=ph,ga=ga[5])*
               4+oh(pco02=pco02, temp1=temp1, ph=ph,ga=ga[6]))
  #temperaturabhängige Koeffizienten usw
  A<-0.504
  B<-0.327e+8
  d<-c(6e-8, 3e-8, 9e-8, 4e-8, 4.5e-8, 3.5e-8)
  z<-c(2,1,1,1,2,1)
  #Vektor f?r Gamma erstellen und durch Debye-H?ckel-Formel neu errechen
  ga<-numeric(6) 
  for(i in 1:6){
    ga[i] <- 10^(((-A*z[i]^2) * sqrt(mu))/(1+B * d[i] * (sqrt(mu))))
  }
  return(ga) 
}

#Wir haben 6 unbekannte Spezies, wir kennen PCO2 und kennen 7 Zusammenh?nge-----
#Gleichungen -> versuchen zu l?sen
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
    if(abs(ph1-ph)<1/1000){#f?r schnellere Durchl?ufe reichen auch 1000
      return(c(ga,ph, pco02)) #Bevor auch die zweite repeat Schleife
      #abricht werden 6*Gammas, PH, und der CO2-Patialdruck ausgeben
      break
    }else{
      ga <- phmu(pco02=pco02, temp1=temp1, ph=ph, ga=ga)#neue Gammas abrufen
    }
  }
}
#Beispiel
co2pa_ca(pco02=0.3, temp=10)#0.7729347 0.9341701 0.9407609 0.9353685
#0.7673852 0.9347748 7.5415645 0.3000000

#Für einen Bekannten pH-Wert und einen unbekannten CO2 Patialdruck#####
#Wir haben 6 unbekannte Spezies, wir kennen PCO2 und kennen 7 Zusammenh?nge-----
#Gleichungen -> versuchen zu l?sen
pH_ca_fun <- function(ph, temp1){#Funktion braucht Patialdurck und Temperatur
  ga <- rep(1,6) #Gamma anfangs immer auf 1 und Patialdruck auf 2.5
  pco02 <- 2.5
  repeat{
    pco02_1 <- pco02
    Int <- 2.5
    repeat{
      er <- sumCA(pco02=pco02, temp1=temp1, ph=ph, ga=ga)
      if(abs(er) < 1e-8|Int==0){
        break
      }else{
        Int <- Int/2
        if(er > 0){
          pco02 <- pco02 + Int
        }else{
          pco02 <- pco02 - Int}
      }
    }
    if(abs(pco02_1 - pco02)<1/10000){#f?r schnellere Durchl?ufe reichen auch 1000
      return(c(ga,ph, pco02)) #Bevor auch die zweite repeat Schleife
      #abricht werden 6*Gammas, PH, und der CO2-Patialdruck ausgeben
      break
    }else{
      ga <- phmu(pco02=pco02, temp1=temp1, ph=ph, ga=ga)#neue Gammas abrufen
    }
  }
}
#Beispiel
pH_ca_fun(ph = 5, temp=10)#0.7729347 0.9341701 0.9407609 0.9353685
#0.7673852 0.9347748 7.5415645 0.3000000
##