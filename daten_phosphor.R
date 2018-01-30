source("rotation_fun.R") 
source("dissociation_of_ophos.R") 
library(readr)
ICP_Calcium <- read_csv("ICP Calcium.csv")
calcit_ph_temp <- read_csv("calcit_ph_temp")
ICP_Calcium$camoll <- ICP_Calcium$Cammoll/1e+6
plot(ICP_Calcium$camoll)
plot(log10(calcit_ph_temp$V1), type ="l")
points(log10(calcit_ph_temp$V26), col="red", type = "l")

#Beispiel Modell für Temp von 25
fmv <- glm(calcit_ph_temp$X1 ~ I(calcit_ph_temp$V25^-2) * I(calcit_ph_temp$V25^-3) * I(calcit_ph_temp$V25^-4)* I(calcit_ph_temp$V25^-5))
summary(fmv)
plot(calcit_ph_temp$V25, 10^ -(calcit_ph_temp$X1), type="l")
lines(10^- (fitted(fmv)) ~ calcit_ph_temp$V25, col = "red")


newdat <- seq(min(calcit_ph_temp$V1), max(calcit_ph_temp$V1), length=200)
preds <- predict.glm(fmv, newdata = data.frame("V1" = newdat))
points(preds, 10^-(calcit_ph_temp$X1), type= "l")
res <- preds-calcit_ph_temp$X1 #Residuen 
plot(res)
sum(abs(res))

calcit_ph_temp1 <- as.matrix(calcit_ph_temp)
pred_ph <- function(ca, Temp){
  Temp <- calcit_ph_temp1[,Temp+1]
  fmv <- glm(calcit_ph_temp1[,1] ~ I(Temp^-2) * I(Temp^-3) * I(Temp^-4)* I(Temp^-5))
  newdat <- ca
  preds <- predict.glm(fmv, newdata = data.frame("Temp" = newdat))
  return(preds)
}
pred_ph(0.0014735733, 2) #Beispiel Daten
erg_ver <- numeric(length(ICP_Calcium$Horizont)) #Berechnen für alle gemessen Temperaturen
erg_ph <- numeric(length(ICP_Calcium$Horizont))
for(i in 1:length(ICP_Calcium$Horizont)){
  erg_ph[i]  <- pred_ph(ca= ICP_Calcium$camoll[i], Temp = ICP_Calcium$TemperaturC[i])
  erg_ver[i] <- HPO2(ph= erg_ph[i], gammaca = ICP_Calcium$camoll[i], g=F)
}
#Unterschied zwischen den Gemessenen und den modellierten Werten:
diff_res <- erg_ver - (ICP_Calcium$Mittelwerteh2oExtraktePO4mmoll/1e+6)
plot(ICP_Calcium$TemperaturC, log10(abs(diff_res))) 
plot(log10(erg_ver), ylim= c(-9,-4))
points(log10(ICP_Calcium$Mittelwerteh2oExtraktePO4mmoll/1e+6), col="red")

plot(log10(ICP_Calcium$camoll),  log10(erg_ver), ylim= c(-9,-4))
points(log10(ICP_Calcium$camoll), log10(ICP_Calcium$Mittelwerteh2oExtraktePO4mmoll/1e+6), col="red")

plot(erg_ph,  log10(erg_ver), ylim= c(-9,-4))
points(erg_ph, log10(ICP_Calcium$Mittelwerteh2oExtraktePO4mmoll/1e+6), col="red")
#Spezies der 
mat_spec <- as.data.frame(sapply(erg_ver, parts))


