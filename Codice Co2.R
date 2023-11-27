# Librerie----------------------------------------------------------------
library(car)
library(caret)
library(coefplot)
library(dplyr)
library(factorMerger)
library(forestmodel)
library(gam)
library(ggcorrplot)
library(ggplot2)
library(gvlma)
library(lmtest)
library(MASS)
library(mctest)
library(plyr)
library(sandwich)
library(stringr)
library(tidyverse)
library(VIM)
library(xtable)

# Campionamento, non eseguire -----------------------------------------------------------

# Importazione dataset intero
ds_tot<-read.csv("CENED.csv", header=TRUE, sep=',')

# Calcolo l'intervallo di campionamento
intervallo <- floor(nrow(ds_tot) / 20000)  

# Eseguo il campionamento sistematico
campione_sistematico <- ds_tot[seq(from = sample(1:intervallo, 1), to = nrow(ds_tot), by = intervallo), ]

# Esportazione dataset campionario
write.csv(ds_tot, file = "campionamento_sistematico.csv", row.names = FALSE)

# Importazione dati ---------------------------------------

# Importazione dataset campionamento_sistematico
 setwd("C:/Users/keita/OneDrive/Documenti/universita/data mining e machine learning/progetto DM")
 campione<-read.csv("campionamento_sistematico.csv", header=TRUE, sep=';')

# Rendere fattoriale le variabili
campione$PROVINCIA <- as.factor(campione$PROVINCIA)
campione$DESTINAZIONE_DI_USO <- as.factor(campione$DESTINAZIONE_DI_USO)
campione$MOTIVAZIONE_APE <- as.factor(campione$MOTIVAZIONE_APE)
campione$CLASSE_ENERGETICA <- as.factor(campione$CLASSE_ENERGETICA)
campione$TIPOLOGIA_VENTILAZIONE <- as.factor(campione$TIPOLOGIA_VENTILAZIONE)
campione$TIPOLOGIA_PANNELLO_ST <- as.factor(campione$TIPOLOGIA_PANNELLO_ST)
campione$TIPOLOGIA_PANNELLO_FV <- as.factor(campione$TIPOLOGIA_PANNELLO_FV)
campione$TIPOLOGIA_COMBUSTIBILE <- as.factor(campione$TIPOLOGIA_COMBUSTIBILE)
campione$TIPOLOGIA_GENERATORE <- as.factor(campione$TIPOLOGIA_GENERATORE)

# Selezione edifici abitativi
ds <- campione %>%
  filter(DESTINAZIONE_DI_USO == "E.1(1)" & EDIFICIO_PUBBLICO == "NO")

# Prime analisi -----------------------------------------------------------

par(mfrow = c(3,3))
plot(ds$"PROVINCIA", ds$"EMISSIONI_DI_CO2", xlab="PROVINCIA", ylab="EMISSIONI_CO2")
plot(ds$"MOTIVAZIONE_APE", ds$"EMISSIONI_DI_CO2", xlab="MOTIVAZIONE_APE", ylab="EMISSIONI_CO2")
plot(ds$"VOLUME_LORDO", ds$"EMISSIONI_DI_CO2", xlab="VOLUME_LORDO", ylab="EMISSIONI_CO2")
plot(ds$"SUPERFICIE_VETRATA_OPACA", ds$"EMISSIONI_DI_CO2", xlab="SUP_VET_OPACA", ylab="EMISSIONI_CO2")
plot(ds$"CLASSE_ENERGETICA", ds$"EMISSIONI_DI_CO2", xlab="CLASSE_ENERGETICA", ylab="EMISSIONI_CO2")
plot(ds$"EPH", ds$"EMISSIONI_DI_CO2", xlab="EPH", ylab="EMISSIONI_CO2")
plot(ds$"ETH", ds$"EMISSIONI_DI_CO2", xlab="ETH", ylab="EMISSIONI_CO2")
plot(ds$"ETC", ds$"EMISSIONI_DI_CO2", xlab="ETC", ylab="EMISSIONI_CO2")
plot(ds$"EFER", ds$"EMISSIONI_DI_CO2", xlab="EFER", ylab="EMISSIONI_CO2")
par(mfrow = c(1,1))
par(mfrow = c(2,3))
plot(ds$"EPW", ds$"EMISSIONI_DI_CO2", xlab="EPW", ylab="EMISSIONI_CO2")
plot(ds$"EPT", ds$"EMISSIONI_DI_CO2", xlab="EPT", ylab="EMISSIONI_CO2")
plot(ds$"EF_GLOB_MEDIA_RISCALDAMENTO", ds$"EMISSIONI_DI_CO2", xlab="EF_GLOB_MEDIA_RISCALDAMENTO", ylab="EMISSIONI_CO2")
plot(ds$"EF_GLOB_MEDIA_ACQUA_CALDA_SAN", ds$"EMISSIONI_DI_CO2", xlab="EF_GLOB_MEDIA_ACQUA_CALDA_SAN", ylab="EMISSIONI_CO2")
plot(ds$"EGHW", ds$"EMISSIONI_DI_CO2", xlab="EGHW", ylab="EMISSIONI_CO2")
par(mfrow = c(1,1))

summary(ds)

# Missing data ------------------------------------------------------------

# Creazione dataset covariate
ds_cov <- ds[, c("PROVINCIA",
                 "MOTIVAZIONE_APE",
                 "VOLUME_LORDO",
                 "SUPERFICIE_VETRATA_OPACA",
                 "CLASSE_ENERGETICA",
                 "EPH",
                 "ETH",
                 "ETC",
                 "EFER",
                 "EPW",
                 "EPT",
                 "EF_GLOB_MEDIA_RISCALDAMENTO",
                 "EF_GLOB_MEDIA_ACQUA_CALDA_SAN",
                 "EGHW")]

# Conto valori mancanti per variabile per fare grafico 
missing_data <- ds_cov %>% 
  summarise_all(function(x) sum(is.na(x) | x == "")) %>% 
  gather(variable, missing_count)

# Spazi vuoti di motivazione ape convertiti in na per il grafico
ds_cov$MOTIVAZIONE_APE[ds_cov$MOTIVAZIONE_APE == ""] <- NA

# Grafici
ggplot(missing_data, aes(x = reorder(variable, missing_count), y = missing_count)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(x="Variabile",
       y="Numero di dati mancanti") +
  theme_minimal()

missingness = aggr(ds_cov[, c("MOTIVAZIONE_APE", "SUPERFICIE_VETRATA_OPACA")],
                   col=c('navyblue','yellow'),numbers=TRUE,sortVars=TRUE,
                   labels=c("MOT_APE", "SUP_VET_OP"),cex.axis=.7,gap=2) 

# Visualizzazione oservazioni che hanno valori mancati per le variabili SUPERFICIE_VETRATA_OPACA MOTIVAZIONE_APE 
righe_mancanti <- ds_cov[is.na(ds_cov$MOTIVAZIONE_APE) | ds_cov$MOTIVAZIONE_APE == "" |
                           is.na(ds_cov$SUPERFICIE_VETRATA_OPACA) | ds_cov$SUPERFICIE_VETRATA_OPACA == "", ]
print(righe_mancanti)

# Dataset dopo missing data
ds <- ds[, c("EMISSIONI_DI_CO2",
    "CODICE_IDENTIFICATIVO_PRATICA",
                       "PROVINCIA",
                       "MOTIVAZIONE_APE",
                       "VOLUME_LORDO",
                       "SUPERFICIE_VETRATA_OPACA",
                       "CLASSE_ENERGETICA",
                       "EPH",
                       "ETH",
                       "ETC",
                       "EFER",
                       "EPW",
                       "EPT",
                       "EF_GLOB_MEDIA_RISCALDAMENTO",
                       "EF_GLOB_MEDIA_ACQUA_CALDA_SAN",
                       "EGHW")]
ds <- ds %>%
  filter(
    !(MOTIVAZIONE_APE == "" | is.na(MOTIVAZIONE_APE)),
    !(SUPERFICIE_VETRATA_OPACA == "" | is.na(SUPERFICIE_VETRATA_OPACA))
  )

# Eliminare 8 osservazioni problematiche, errori di inserimento
ds <- ds[-c(3229, 7935, 10104, 16509), ]
ds <- ds[ds$EFER >= 0 & ds$EFER < 4000000, ]

# Controllo missing
missingness = aggr(ds[, c("MOTIVAZIONE_APE", "SUPERFICIE_VETRATA_OPACA")],
                   col=c('navyblue','yellow'),numbers=TRUE,sortVars=TRUE,
                   labels=c("MOT_APE", "SUP_VET_OP"),cex.axis=.7,gap=2) 

# Optimal grouping --------------------------------------------------------

# Numero di livelli per variabile
num_livelli <- sapply(ds, function(x) if(is.factor(x)) length(levels(x)) else NA)
print(num_livelli)

# Grafici 
par(mfrow = c(1,3))
plot(ds$MOTIVAZIONE_APE, ds$EMISSIONI_DI_CO2)
plot(ds$PROVINCIA, ds$EMISSIONI_DI_CO2)
plot(ds$CLASSE_ENERGETICA, ds$EMISSIONI_DI_CO2)
par(mfrow = c(1,1))

# Classe_energetica
ds <- ds %>%
  mutate(CLASSE_ENERGETICA = case_when(
    CLASSE_ENERGETICA %in% c("A", "A+") ~ "Alta",
    CLASSE_ENERGETICA %in% c("B", "C", "D") ~ "Media",
    CLASSE_ENERGETICA %in% c("E", "F", "G") ~ "Bassa",
    TRUE ~ CLASSE_ENERGETICA  
  ))

# Provincia
ds <- ds %>%
  mutate(PROVINCIA = case_when(
    PROVINCIA %in% c("MI", "MB", "VA") ~ "Province densità elevata",
    PROVINCIA %in% c("LC", "BG", "CO", "LO", "BS") ~ "Province densità media",
    PROVINCIA %in% c("MN","PV", "CR", "SO") ~ "Province densità bassa",
    TRUE ~ PROVINCIA  
  ))

ds$CLASSE_ENERGETICA <- as.factor(ds$CLASSE_ENERGETICA)
ds$PROVINCIA <- as.factor(ds$PROVINCIA)

# MOTIVAZIONE_APE 
reduce_levels_MOT <- mergeFactors(response = ds$EMISSIONI_DI_CO2, factor = ds$MOTIVAZIONE_APE)
mergingHistory(reduce_levels_MOT, showStats = TRUE ) %>%head(5)

# plot aggregations based on mimimum AIC
plot(reduce_levels_MOT, panel = "GIC",title = "", panelGrid = FALSE )

# plot aggregations based on mimimum SBC
plot(reduce_levels_MOT, panel = "GIC", penalty = log(NROW(ds$MOTIVAZIONE_APE)), title = "", panelGrid = FALSE )

# plot aggregations using means
plot(reduce_levels_MOT, palette = "Reds")
plot(reduce_levels_MOT, responsePanel = "boxplot", colorCluster = TRUE)

# save reduced levels
MOTIVAZIONE_APE_REDUCED=cutTree(reduce_levels_MOT)
length(MOTIVAZIONE_APE_REDUCED)
class(MOTIVAZIONE_APE_REDUCED)
table(MOTIVAZIONE_APE_REDUCED)

# add to dataset
ds$MOTIVAZIONE_APE_REDUCED=as.numeric(MOTIVAZIONE_APE_REDUCED)
ds$MOTIVAZIONE_APE_REDUCED=as.factor(ds$MOTIVAZIONE_APE_REDUCED)
table(MOTIVAZIONE_APE_REDUCED,ds$MOTIVAZIONE_APE_REDUCED)

# try if the gruoped var has been grouped well 
par(mfrow = c(1,3))
plot(ds$MOTIVAZIONE_APE_REDUCED, ds$EMISSIONI_DI_CO2)
plot(ds$PROVINCIA, ds$EMISSIONI_DI_CO2)
plot(ds$CLASSE_ENERGETICA, ds$EMISSIONI_DI_CO2)
par(mfrow = c(1,1))

# ggplot
ggplot(ds, aes(x=MOTIVAZIONE_APE_REDUCED, y=EMISSIONI_DI_CO2, color=MOTIVAZIONE_APE_REDUCED)) +
  geom_boxplot(notch = TRUE, 
               outlier.colour="darkorange") +        
  labs(x = "MOTIVAZIONE_APE_REDUCED") +
  theme(legend.title = element_blank()) +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Barplot
ds_count <- ds %>%
  group_by(MOTIVAZIONE_APE_REDUCED) %>%
  tally(name = "Count")

ggplot(ds_count, aes(x = MOTIVAZIONE_APE_REDUCED, y = Count)) + 
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  coord_flip() +  # Rendi il grafico orizzontale (rimuovi questa riga per un barplot verticale)
  labs(x = "MOTIVAZIONE_APE_REDUCED",
       y = "Numero di Occorrenze") +
  theme_minimal()

# Primo fit ---------------------------------------------------------------

# FIT del primo modello
fit1=lm(EMISSIONI_DI_CO2 ~
      PROVINCIA+MOTIVAZIONE_APE_REDUCED+VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+CLASSE_ENERGETICA+EPH+ETH+ETC+
      EFER+EPW+EPT+EGHW+EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)
summary(fit1)
length(fit1$residuals)
drop1(fit1, test="F")

# Xtable
xtable(fit1)

# plots importanti
par(mfrow=c(2,2))
plot(fit1)
par(mfrow=c(1,1))

# Grafici coefficienti
coefplot(fit1, decreasing = TRUE, sort = "magnitude",intercept=FALSE)
print(forest_model(fit1))

# Multicollinearità -------------------------------------------------------

# Variabili fattoriali
cov=attr(terms(fit1), "term.labels") 
cov

ds_fac <- ds[,cov] %>% dplyr::select_if(is.factor)
colnames(ds_fac)

combos <- combn(ncol(ds_fac),2)
adply(combos, 2, function(x) {
  test <- chisq.test(ds_fac[, x[1]], ds_fac[, x[2]])
  tab  <- table(ds_fac[, x[1]], ds_fac[, x[2]])
  out <- data.frame("Row" = colnames(ds_fac)[x[1]]
                    , "Column" = colnames(ds_fac[x[2]])
                    , "Chi.Square" = round(test$statistic,3)
                    , "df"= test$parameter
                    , "p.value" = round(test$p.value, 3)
                    , "n" = sum(table(ds_fac[,x[1]], ds_fac[,x[2]]))
                    , "u1" =length(unique(ds_fac[,x[1]]))-1
                    , "u2" =length(unique(ds_fac[,x[2]]))-1
                    , "nMinu1u2" =sum(table(ds_fac[,x[1]], ds_fac[,x[2]]))* min(length(unique(ds_fac[,x[1]]))-1 , length(unique(ds_fac[,x[2]]))-1) 
                    , "Chi.Square norm"  =test$statistic/(sum(table(ds_fac[,x[1]], ds_fac[,x[2]]))* min(length(unique(ds_fac[,x[1]]))-1 , length(unique(ds_fac[,x[2]]))-1)) 
  )
  
  
  return(out)
  
})

# Variabili quantitative
fit_tol <- lm(EMISSIONI_DI_CO2 ~
                VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+EPH+ETH+ETC+EFER+EPW+EPT+EGHW
              +EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)
imcdiag(fit_tol)

fit_tol_1 <- lm(EMISSIONI_DI_CO2 ~
                  VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+EPH+ETH+ETC+EFER+EPW+EGHW
                +EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)
imcdiag(fit_tol_1)

# Fit del modello senza variabili problematiche
fit2=lm(EMISSIONI_DI_CO2 ~
          PROVINCIA+MOTIVAZIONE_APE_REDUCED+VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+CLASSE_ENERGETICA+EPH+ETH+ETC+
          EFER+EPW+EGHW+EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)
summary(fit2)
length(fit2$residuals)
anova(fit2)

par(mfrow=c(2,2))
plot(fit2)
par(mfrow=c(1,1))

# Test assunzioni
gvlma(fit2)

# Linearità ---------------------------------------------------------------

# BOX-COX
ds$y_pos <- ds$EMISSIONI_DI_CO2+1
summary(ds$EMISSIONI_DI_CO2)
summary(log(ds$y_pos))
fit3 <- lm(y_pos~
             PROVINCIA+MOTIVAZIONE_APE_REDUCED+VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+CLASSE_ENERGETICA+EPH+ETH+ETC+
             EFER+EPW+EGHW+EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)
boxcox_reg1<-boxcox(fit3)

title("Lambda")
lambda=round(boxcox_reg1$x[which.max(boxcox_reg1$y)], 2)
lambda

lambda_fit3 <- lm(log(y_pos)~
                    PROVINCIA+MOTIVAZIONE_APE_REDUCED+VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+CLASSE_ENERGETICA+EPH+ETH+ETC+
                    EFER+EPW+EGHW+EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)
summary(lambda_fit3)

par(mfrow = c(2, 2)) 
plot(lambda_fit3)
par(mfrow = c(1, 1))

# spline
lambda_fit3_gam_1 <- gam(log(y_pos)~
                         PROVINCIA+MOTIVAZIONE_APE_REDUCED+VOLUME_LORDO+SUPERFICIE_VETRATA_OPACA+CLASSE_ENERGETICA+EPH+ETH+ETC+
                         EFER+EPW+EGHW+EF_GLOB_MEDIA_RISCALDAMENTO+EF_GLOB_MEDIA_ACQUA_CALDA_SAN, data=ds)

lambda_fit3_gam <- gam(log(y_pos)~
                         PROVINCIA+
                         MOTIVAZIONE_APE_REDUCED+
                         s(VOLUME_LORDO)+
                         s(SUPERFICIE_VETRATA_OPACA)+
                         CLASSE_ENERGETICA+
                         s(EPH)+
                         s(ETH)+
                         s(ETC)+
                         s(EFER)+
                         s(EPW)+
                         s(EGHW)+
                         s(EF_GLOB_MEDIA_RISCALDAMENTO)+
                         s(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data=ds)
summary(lambda_fit3_gam)

plot(lambda_fit3_gam)
par(mfrow=c(1,1)) 

summary(ds)
ds$EF_GLOB_MEDIA_ACQUA_CALDA_SAN = ds$EF_GLOB_MEDIA_ACQUA_CALDA_SAN+1
ds$EF_GLOB_MEDIA_RISCALDAMENTO = ds$EF_GLOB_MEDIA_RISCALDAMENTO+1
ds$ETH=ds$ETH+1

fit4 <- lm(log(y_pos)~    PROVINCIA+
             MOTIVAZIONE_APE_REDUCED+
             log(VOLUME_LORDO)+
             SUPERFICIE_VETRATA_OPACA+
             CLASSE_ENERGETICA+
             I((EPH)^(0.5))+
             log(ETH)+
             I(ETC^2)+
             ETC+
             EFER+
             EPW+
             I(EGHW^2)+
             EGHW+
             I((1/EF_GLOB_MEDIA_RISCALDAMENTO))+
             log(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data=ds)
summary(fit4)

par(mfrow=c(2,2))
plot(fit4)
par(mfrow=c(1,1))

plot(fit4$fitted.values, log(ds$EMISSIONI_DI_CO2))

resettest(fit4, power = 2, type = "fitted", data = ds)

# Linearità Loess ---------------------------------------------------------

ggplot(data = ds, aes(x = VOLUME_LORDO, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw()

ggplot(data = ds, aes(x = EPH, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw()

ggplot(data = ds, aes(x = ETH, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw()

ggplot(data = ds, aes(x = ETC, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw()

ggplot(data = ds, aes(x = EF_GLOB_MEDIA_ACQUA_CALDA_SAN, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 30))

ggplot(data = ds, aes(x = EF_GLOB_MEDIA_RISCALDAMENTO, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw()

ggplot(data = ds, aes(x = EGHW, y = log(y_pos))) + geom_point(color = "gray") +
  geom_smooth(method = "loess", span = 0.2, col = "red") +
  geom_smooth(method = "loess", span = 0.5, col = "purple3") +
  theme_bw()

ds_loess <- ds 
ds_loess$EGHW = ds_loess$EGHW+1
ds_loess$ETC=ds_loess$ETC+1

lambda_fit3_gam_3 <- lm(log(y_pos)~
                          PROVINCIA+
                          MOTIVAZIONE_APE_REDUCED+
                          I((1/VOLUME_LORDO))+
                          SUPERFICIE_VETRATA_OPACA+
                          CLASSE_ENERGETICA+
                          I((EPH)^(0.5))+
                          log(ETH)+
                          (I(1/ETC))+
                          EFER+
                          EPW+
                          I((1/EGHW))+
                          I((1/EF_GLOB_MEDIA_RISCALDAMENTO))+
                          I((1/EF_GLOB_MEDIA_ACQUA_CALDA_SAN)), data=ds_loess)
summary(lambda_fit3_gam_3)

par(mfrow=c(2,2))
plot(lambda_fit3_gam_3)
par(mfrow=c(1,1))

resettest(lambda_fit3_gam_3, power = 2, type = "fitted", data = ds_loess)

# Punti influenti ---------------------------------------------------------

par(mfrow=c(1,1))
influencePlot(fit4,  main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

# Calcola i valori di leverage utilizzando hatvalues(fit1)
leverage_values <- hatvalues(fit4)

# Calcola la media dei valori di leverage moltiplicata per 2
mean_leverage <- 2 * mean(leverage_values)

# Crea un dataframe con i dati
data <- data.frame(index = 1:length(leverage_values), leverage = leverage_values)

#Conto punti di leva
count_leverage <- which(hatvalues(fit4)>2*mean(hatvalues(fit4))) 

# Crea il grafico utilizzando ggplot2
ggplot(data, aes(x = index, y = leverage)) +
  geom_line(stat = "identity", color = "blue") +  
  geom_hline(yintercept = mean_leverage, linetype = "dashed", color = "red") +
  labs(ylab = "leverage", xlab = "index")

# Eliminazione punti influenti
dffits = dffits(fit4)

thresh3 = 2*sqrt(length(fit4$coefficients)/length(dffits))
dffits[dffits > thresh3]

ds$dffits = dffits
Noinflu = ds[ds$dffits <= thresh3, ]

fit5 = lm(log(y_pos)~ PROVINCIA+
            MOTIVAZIONE_APE_REDUCED+
            log(VOLUME_LORDO)+
            SUPERFICIE_VETRATA_OPACA+
            CLASSE_ENERGETICA+
            I((EPH)^(0.5))+
            log(ETH)+
            I(ETC^2)+
            ETC+
            EFER+
            EPW+
            I(EGHW^2)+
            EGHW+
            I((1/EF_GLOB_MEDIA_RISCALDAMENTO))+
            log(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data=Noinflu)

summary(fit5)
drop1(fit5, test = "F")

par(mfrow = c(2,2))
plot(fit5)
par(mfrow = c(1,1))

plot(fit5$fitted.values, log(Noinflu$y_pos))
abline(a=0, b=1, col="red", lty=2)
resettest(fit5, power = 2, type = "fitted", data = Noinflu)

# Eteroschedasticità ------------------------------------------------------

plot(fit5$fitted.values, fit5$residuals)

bptest(fit5)
ncvTest(fit5)

coeftest(fit5, vcov = vcovHC(fit5))

std_correct = fit5 %>% 
  vcovHC() %>% 
  diag() %>% 
  sqrt()
std_correct

# proviamo a togliere SUPERFICIE_VETRATA_OPACA EPW ETC            

fit6 = lm(log(y_pos)~ PROVINCIA+
                 MOTIVAZIONE_APE_REDUCED+
                 log(VOLUME_LORDO)+
                 CLASSE_ENERGETICA+
                 I((EPH)^(0.5))+
                 log(ETH)+
                 EFER+
                 I(EGHW^2)+
                 EGHW+
                 I((1/EF_GLOB_MEDIA_RISCALDAMENTO))+
                 log(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data=Noinflu)

plot(fit6$fitted.values, log(Noinflu$y_pos))
abline(a=0, b=1, col="red", lty=2)
bptest(fit6)
ncvTest(fit6)

coeftest(fit6, vcov = vcovHC(fit6))

gvlma(fit6)

anova(fit5, fit6, test="F")

# Model Selection ---------------------------------------------------------

step <- stepAIC(fit6, direction="both")
stepSBC <- stepAIC(fit6, direction = "both", k = log(nrow(Noinflu)))

fit7 <- lm(log(y_pos) ~ PROVINCIA + MOTIVAZIONE_APE_REDUCED + log(VOLUME_LORDO) + 
                           CLASSE_ENERGETICA + I((EPH)^(0.5)) + log(ETH) + I(EGHW^2) + 
                           EGHW + I((1/EF_GLOB_MEDIA_RISCALDAMENTO)) + log(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data=Noinflu)
summary(fit7)
drop1(fit7, test="F")

par(mfrow = c(2,2))
plot(fit7)
par(mfrow = c(1,1))
gvlma(fit7)

# Robust inference --------------------------------------------------------

BOOT.MOD=Boot(fit7, R=1999)
summary(BOOT.MOD, high.moments=TRUE)

# confint boot
Confint(BOOT.MOD, level=c(.95), type="perc")
hist(BOOT.MOD, legend="separate")

# Confronto modello iniziale e finale -------------------------------------

plot(fit7$fitted.values, log(Noinflu$y_pos))
abline(a=0, b=1, col="red", lty=2)
plot(fit1$fitted.values, ds$EMISSIONI_DI_CO2)
abline(a=0, b=1, col="red", lty=2)

gvlma(fit7)
print(forest_model(fit7))

# Interpretazione coefficienti variabili categoriche
(exp(fit7$coefficients[c(1:8, 10, 11)])-1)*100

# Modello logistico -------------------------------------------------------

Binario=ds
# Crea la nuova variabile 'co2'
Binario <- Binario %>%
  mutate(co2 = ifelse(EMISSIONI_DI_CO2 < 40, 1, 0))
table(Binario$co2)
Binario$co2 <- factor(Binario$co2)

# Grafico variabile dipendente binaria
ggplot(Binario, aes(x=co2, fill=co2)) +
  geom_bar(position='dodge') +
  geom_text(stat='count', aes(label=scales::percent(..count../sum(..count..))), vjust=-0.5) +
  labs(x="status", y="Percentuale") +
  theme_minimal()

#Modello
logit = glm(co2 ~ PROVINCIA + MOTIVAZIONE_APE_REDUCED + CLASSE_ENERGETICA + log(VOLUME_LORDO)  + I((EPH)^(0.5)) + log(ETH) + I(EGHW^2) + 
              EGHW + I((1/EF_GLOB_MEDIA_RISCALDAMENTO)) + log(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data = Binario, family = "binomial")
summary(logit)
drop1(logit, test="LRT")
exp(logit$coefficients)

# Quasi separation
ggplot(Binario, aes(x = CLASSE_ENERGETICA, y = co2)) + 
  geom_jitter(width = 0.2, height = 0.2, size = 3, alpha = 0.6) +
  labs(x = "Classe Energetica", y = "CO2 Emissions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(Binario$CLASSE_ENERGETICA, Binario$co2)

# Modello senza classe energetica
logit1 = glm(co2 ~ PROVINCIA + MOTIVAZIONE_APE_REDUCED + log(VOLUME_LORDO)  + I((EPH)^(0.5)) + log(ETH) + I(EGHW^2) + 
              EGHW + I((1/EF_GLOB_MEDIA_RISCALDAMENTO)) + log(EF_GLOB_MEDIA_ACQUA_CALDA_SAN), data = Binario, family = "binomial")
summary(logit1)
drop1(logit1, test="LRT")
print(forest_model(logit1))
exp(logit1$coefficients)

#ACCURANCY CONFUSION MATRIX 

# Probabilità predette
probs <- predict(logit1, type = "response")

# Conversioe probabilità predette in classi predette usando una soglia di 0.5
predizioni <- ifelse(probs > 0.5, 1, 0)

# Fattorizzazione variabili
risposta_reale <- as.factor(Binario$co2)
predizioni <- as.factor(predizioni)

# Matrice di confusione
cm <- confusionMatrix(predizioni, risposta_reale)
print(cm)

#R2 mac fadden
null = glm(co2 ~ 1, data=Binario,family = binomial)
R2=1-(logit1$deviance/null$deviance)
R2

# Test lrt
anova(null, logit, test = "LRT")