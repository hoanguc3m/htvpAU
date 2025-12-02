library(fattvpVAR)
numCores = 15
# setwd("/Backup/Ongoing/WP11/")
setwd("/home/hoanguc3m/Downloads/htvpAU")
##########################################################################

load("G000.RData")
G000_ML_Ent <- MLEnt_TVPGSV(Chain = G000_obj, numCores = numCores)
save(G000_ML_Ent, file = "G000_ML_Ent.RData")
rm(G000_ML_Ent, G000_obj)

load("G001.RData")
G001_ML_Ent <- MLEnt_TVPGSV(Chain = G001_obj, numCores = numCores)
save(G001_ML_Ent, file = "G001_ML_Ent.RData")
rm(G001_ML_Ent, G001_obj)

load("G010.RData")
G010_ML_Ent <- MLEnt_TVPGSV(Chain = G010_obj, numCores = numCores)
save(G010_ML_Ent, file = "G010_ML_Ent.RData")
rm(G010_ML_Ent, G010_obj)

load("G100.RData")
G100_ML_Ent <- MLEnt_TVPGSV(Chain = G100_obj, numCores = numCores)
save(G100_ML_Ent, file = "G100_ML_Ent.RData")
rm(G100_ML_Ent, G100_obj)

load("G011.RData")
G011_ML_Ent <- MLEnt_TVPGSV(Chain = G011_obj, numCores = numCores)
save(G011_ML_Ent, file = "G011_ML_Ent.RData")

load("G101.RData")
G101_ML_Ent <- MLEnt_TVPGSV(Chain = G101_obj, numCores = numCores)
save(G101_ML_Ent, file = "G101_ML_Ent.RData")

load("G110.RData")
G110_ML_Ent <- MLEnt_TVPGSV(Chain = G110_obj, numCores = numCores)
save(G110_ML_Ent, file = "G110_ML_Ent.RData")

load("G111.RData")
G111_ML_Ent <- MLEnt_TVPGSV(Chain = G111_obj, numCores = numCores)
save(G111_ML_Ent, file = "G111_ML_Ent.RData")

##########################################################################

load("T000.RData")
T000_ML_Ent <- MLEnt_TVPTSV(Chain = T000_obj, numCores = numCores)
save(T000_ML_Ent, file = "T000_ML_Ent.RData")
rm(T000_ML_Ent, T000_obj)


load("T001.RData")
T001_ML_Ent <- MLEnt_TVPTSV(Chain = T001_obj, numCores = numCores)
save(T001_ML_Ent, file = "T001_ML_Ent.RData")
rm(T001_ML_Ent, T001_obj)

load("T010.RData")
T010_ML_Ent <- MLEnt_TVPTSV(Chain = T010_obj, numCores = numCores)
save(T010_ML_Ent, file = "T010_ML_Ent.RData")
rm(T010_ML_Ent, T010_obj)

load("T100.RData")
T100_ML_Ent <- MLEnt_TVPTSV(Chain = T100_obj, numCores = numCores)
save(T100_ML_Ent, file = "T100_ML_Ent.RData")
rm(T100_ML_Ent, T100_obj)

load("T011.RData")
T011_ML_Ent <- MLEnt_TVPTSV(Chain = T011_obj, numCores = numCores)
save(T011_ML_Ent, file = "T011_ML_Ent.RData")

load("T101.RData")
T101_ML_Ent <- MLEnt_TVPTSV(Chain = T101_obj, numCores = numCores)
save(T101_ML_Ent, file = "T101_ML_Ent.RData")

load("T110.RData")
T110_ML_Ent <- MLEnt_TVPTSV(Chain = T110_obj, numCores = numCores)
save(T110_ML_Ent, file = "T110_ML_Ent.RData")

load("T111.RData")
T111_ML_Ent <- MLEnt_TVPTSV(Chain = T111_obj, numCores = numCores)
save(T111_ML_Ent, file = "T111_ML_Ent.RData")


setwd("/home/hoanguc3m/Downloads/htvpAU/")
load("G000_ML_Ent.RData"); load("G001_ML_Ent.RData"); load("G010_ML_Ent.RData"); load("G100_ML_Ent.RData")
load("G011_ML_Ent.RData"); load("G101_ML_Ent.RData"); load("G110_ML_Ent.RData"); load("G111_ML_Ent.RData")

load("T000_ML_Ent.RData"); load("T001_ML_Ent.RData"); load("T010_ML_Ent.RData"); load("T100_ML_Ent.RData")
load("T011_ML_Ent.RData"); load("T101_ML_Ent.RData"); load("T110_ML_Ent.RData"); load("T111_ML_Ent.RData")

tab1 <- round(cbind(c(G000_ML_Ent$LL, G100_ML_Ent$LL, G010_ML_Ent$LL, G001_ML_Ent$LL, G110_ML_Ent$LL, G101_ML_Ent$LL, G011_ML_Ent$LL, G111_ML_Ent$LL),
                    c(T000_ML_Ent$LL, T100_ML_Ent$LL, T010_ML_Ent$LL, T001_ML_Ent$LL, T110_ML_Ent$LL, T101_ML_Ent$LL, T011_ML_Ent$LL, T111_ML_Ent$LL)), digits = 2)

rownames(tab1) <- c("000", "100","010","001","110","101","011","111")

load("T000.RData");load("T100.RData");load("T010.RData");load("T001.RData");
load("T110.RData");load("T101.RData");load("T011.RData");load("T111.RData");
T000_nu <- apply(T000_obj$store_nu, MARGIN = 2, FUN = mean)
T100_nu <- apply(T100_obj$store_nu, MARGIN = 2, FUN = mean)
T010_nu <- apply(T010_obj$store_nu, MARGIN = 2, FUN = mean)
T001_nu <- apply(T001_obj$store_nu, MARGIN = 2, FUN = mean)
T110_nu <- apply(T110_obj$store_nu, MARGIN = 2, FUN = mean)
T101_nu <- apply(T101_obj$store_nu, MARGIN = 2, FUN = mean)
T011_nu <- apply(T011_obj$store_nu, MARGIN = 2, FUN = mean)
T111_nu <- apply(T111_obj$store_nu, MARGIN = 2, FUN = mean)
nu_vec <- c(sprintf("(%.2f; %.2f; %.2f)",T000_nu[1],T000_nu[2],T000_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T100_nu[1],T100_nu[2],T100_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T010_nu[1],T010_nu[2],T010_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T001_nu[1],T001_nu[2],T001_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T110_nu[1],T110_nu[2],T110_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T101_nu[1],T101_nu[2],T101_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T011_nu[1],T011_nu[2],T011_nu[3]),
  sprintf("(%.2f; %.2f; %.2f)",T111_nu[1],T111_nu[2],T111_nu[3]))
tab1a <- cbind( matrix(sprintf("%.2f",tab1),nrow=8) , nu_vec)
rownames(tab1a) <- c("000", "100","010","001","110","101","011","111")

xtable::xtable(tab1a)
###################################
setwd("/home/hoanguc3m/Downloads/htvpAU")
load("G000_ML_Ent.RData"); load("G001_ML_Ent.RData"); load("G010_ML_Ent.RData"); load("G100_ML_Ent.RData")
load("G011_ML_Ent.RData"); load("G101_ML_Ent.RData"); load("G110_ML_Ent.RData"); load("G111_ML_Ent.RData")

load("T000_ML_Ent.RData"); load("T001_ML_Ent.RData"); load("T010_ML_Ent.RData"); load("T100_ML_Ent.RData")
load("T011_ML_Ent.RData"); load("T101_ML_Ent.RData"); load("T110_ML_Ent.RData"); load("T111_ML_Ent.RData")

tab1 <- round(cbind(c(G000_ML_Ent$LL, G100_ML_Ent$LL, G010_ML_Ent$LL, G001_ML_Ent$LL, G110_ML_Ent$LL, G101_ML_Ent$LL, G011_ML_Ent$LL, G111_ML_Ent$LL),
                    c(T000_ML_Ent$LL, T100_ML_Ent$LL, T010_ML_Ent$LL, T001_ML_Ent$LL, T110_ML_Ent$LL, T101_ML_Ent$LL, T011_ML_Ent$LL, T111_ML_Ent$LL)), digits = 2)

setwd("/home/hoanguc3m/Downloads/WP11/NonSV")
load("G000_nonSV_ML_Ent_nonSV.RData"); load("G001_nonSV_ML_Ent_nonSV.RData"); load("G010_nonSV_ML_Ent_nonSV.RData"); load("G100_nonSV_ML_Ent_nonSV.RData")
load("G011_nonSV_ML_Ent_nonSV.RData"); load("G101_nonSV_ML_Ent_nonSV.RData"); load("G110_nonSV_ML_Ent_nonSV.RData"); load("G111_nonSV_ML_Ent_nonSV.RData")

load("T000_nonSV_ML_Ent_nonSV.RData"); load("T001_nonSV_ML_Ent_nonSV.RData"); load("T010_nonSV_ML_Ent_nonSV.RData"); load("T100_nonSV_ML_Ent_nonSV.RData")
load("T011_nonSV_ML_Ent_nonSV.RData"); load("T101_nonSV_ML_Ent_nonSV.RData"); load("T110_nonSV_ML_Ent_nonSV.RData"); load("T111_nonSV_ML_Ent_nonSV.RData")
tab2 <- round(cbind(c(G000_nonSV_ML_Ent$LL, G100_nonSV_ML_Ent$LL, G010_nonSV_ML_Ent$LL, G001_nonSV_ML_Ent$LL, G110_nonSV_ML_Ent$LL, G101_nonSV_ML_Ent$LL, G011_nonSV_ML_Ent$LL, G111_nonSV_ML_Ent$LL),
                    c(T000_nonSV_ML_Ent$LL, T100_nonSV_ML_Ent$LL, T010_nonSV_ML_Ent$LL, T001_nonSV_ML_Ent$LL, T110_nonSV_ML_Ent$LL, T101_nonSV_ML_Ent$LL, T011_nonSV_ML_Ent$LL, T111_nonSV_ML_Ent$LL)), digits = 2)


tab3 <- cbind(tab2, tab1)
tab3 <- matrix(sprintf("%.2f",tab3),nrow=8)
rownames(tab3) <- c("000", "100","010","001","110","101","011","111")

xtable::xtable(tab3)
