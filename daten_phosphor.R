library(readr)
ICP_Calcium <- read_csv("ICP Calcium.csv")
calcit_ph_temp <- read_csv("calcit_ph_temp")
ICP_Calcium$camoll <- ICP_Calcium$Cammoll/1e+6
plot(ICP_Calcium$camoll)
plot(log10(calcit_ph_temp$V1), type ="l")
points(log10(calcit_ph_temp$V26), col="red", type = "l")

fmv <- glm(calcit_ph_temp$X1 ~ I(calcit_ph_temp$V1^-2) * I(calcit_ph_temp$V1^-3) * I(calcit_ph_temp$V1^-4)* I(calcit_ph_temp$V1^-5))
summary(fmv)
plot(calcit_ph_temp$V1, 10^-(calcit_ph_temp$X1), type="l")
lines(10^-(fitted(fmv))~ calcit_ph_temp$V1, col = "red")

newdat <- seq(min(calcit_ph_temp$V1), max(calcit_ph_temp$V1), length=191)
preds <- predict.glm(fmv, newdata = data.frame("V1" = newdat))
points(preds, 10^-(calcit_ph_temp$X1), type= "l")
res <- preds-calcit_ph_temp$X1
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
pred_ph(0.0014735733, 2)
erg_ver <- numeric(length(ICP_Calcium$Horizont))
for(i in 1:length(ICP_Calcium$Horizont)){
  erg_ver[i] <- HPO2(ph=pred_ph(ca= ICP_Calcium$camoll, Temp = ICP_Calcium$TemperaturC), gammaca = ICP_Calcium$camoll)
      }


