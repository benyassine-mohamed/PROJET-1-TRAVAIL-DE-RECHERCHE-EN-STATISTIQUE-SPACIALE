
##############################
####Importation de données####
##############################

#Installation des packages utiles :

install.packages("knitr")
install.packages("rgdal")
install.packages("raster")
install.packages("gstat")
install.packages("fields")
install.packages("maps")
install.packages("sp")
install.packages("rspatial")

library(knitr)
opts_chunk$set(fig.width = 5, fig.height = 5, fig.cap='',  collapse = TRUE)
library(rgdal)
library(raster)
library(gstat)
library(fields)
library(rspatial)


if (!require("rspatial")) devtools::install_github('rspatial/rspatial')

library(rspatial)
#Pollution de l'air des données Californie
#La lecture des données

x <- sp_data("airqual")
head(x)

#Pour obtenir des chiffres plus faciles à lire, je multiplie OZDLYAV par 1000 
x$OZDLYAV <- x$OZDLYAV * 1000

#OZDLYAV : Concentration moyenne


#Creation d'un SpatialPointsDataFrame

library(sp)
coordinates(x) <- ~LONGITUDE + LATITUDE
proj4string(x) <- CRS('+proj=longlat +datum=NAD83')

TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
library(rgdal)

aq =spTransform(x,TA)


#Polygone, modèle de trame 

cageo <- sp_data('counties.rds')
ca <- spTransform(cageo, TA)
r <- raster(ca)
res(r) <- 10  # 10 km  ( les unités de CRS sont en km )
g <- as(r, 'SpatialGrid')



###############################
#### Ajuster un variogramme####
###############################

install.packages("gstat")
library(gstat)
require(gstat)

gs <- gstat(formula=OZDLYAV~1, locations=aq)
v <- variogram(gs, width=20)
head(v)
plot(v)
head(v)


####################################
#### Modelisation du variogramme####
####################################


fve <- fit.variogram(v, vgm(85, "Exp", 75, 20))
fve

plot(variogramLine(fve, 400), type='l', ylim=c(0,120))
points(v[,2:3], pch=20, col='red')

#Sohérique au lieu d'exponentielle

fvs <- fit.variogram(v, vgm(85, "Sph", 75, 20))
fvs

plot(variogramLine(fvs, 400), type='l', ylim=c(0,120) ,col='blue', lwd=2)
points(v[,2:3], pch=20, col='red')


####################################
#### Interpolation du krigeage #####
####################################

k <- gstat(formula=OZDLYAV~1, locations=aq, model=fve)
# Les valeurs prédites
kp <- predict(k, g)

#spplot(kp)

# Variance

ok = brick(kp)
ok = mask(ok,ca)

names(ok) <- c('prediction', 'variance')
plot(ok)

#Comparer avec d'autres méthodes
#Utilisons à nouveau le gstat pour faire l'interpolation des IDW. L'approche de base d'abord.

library(gstat)
idm <- gstat(formula=OZDLYAV~1, locations=aq)
idp <- interpolate(r, idm)

idp <- mask(idp, ca)
plot(idp)

# Optimisation 
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
f1 <- function(x, test, train) {
  nmx <- x[1]
  idp <- x[2]
  if (nmx < 1) return(Inf)
  if (idp < .001) return(Inf)
  m <- gstat(formula=OZDLYAV~1, locations=train, nmax=nmx, set=list(idp=idp))
  p <- predict(m, newdata=test, debug.level=0)$var1.pred
  RMSE(test$OZDLYAV, p)
}
set.seed(20150518)
i <- sample(nrow(aq), 0.2 * nrow(aq))
tst <- aq[i,]
trn <- aq[-i,]
opt <- optim(c(8, .5), f1, test=tst, train=trn)
opt

#our optimal model

m <- gstat(formula=OZDLYAV~1, locations=aq, nmax=opt$par[1], set=list(idp=opt$par[2]))
idw <- interpolate(r, m)
## [inverse distance weighted interpolation]
idw <- mask(idw, ca)
plot(idw)

# Validation croisée

library(dismo)
nfolds <- 5
k <- kfold(aq, nfolds)
ensrmse <- tpsrmse <- krigrmse <- idwrmse <- rep(NA, 5)

for (i in 1:nfolds) {
  test <- aq[k!=i,]
  train <- aq[k==i,]
  m <- gstat(formula=OZDLYAV~1, locations=train, nmax=opt$par[1], set=list(idp=opt$par[2]))
  p1 <- predict(m, newdata=test, debug.level=0)$var1.pred
  idwrmse[i] <-  RMSE(test$OZDLYAV, p1)
  m <- gstat(formula=OZDLYAV~1, locations=train, model=fve)
  p2 <- predict(m, newdata=test, debug.level=0)$var1.pred
  krigrmse[i] <-  RMSE(test$OZDLYAV, p2)
  m <- Tps(coordinates(train), train$OZDLYAV)
  p3 <- predict(m, coordinates(test))
  tpsrmse[i] <-  RMSE(test$OZDLYAV, p3)
  w <- c(idwrmse[i], krigrmse[i], tpsrmse[i])
  weights <- w / sum(w)
  ensemble <- p1 * weights[1] + p2 * weights[2] + p3 * weights[3]
  ensrmse[i] <-  RMSE(test$OZDLYAV, ensemble)
}
rmi <- mean(idwrmse)
rmk <- mean(krigrmse)
rmt <- mean(tpsrmse)
rms <- c(rmi, rmt, rmk)
rms

rme <- mean(ensrmse)
rme