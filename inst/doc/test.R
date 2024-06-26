## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE, echo = FALSE-------------------------------------
library(multipleOutcomes)
library(dplyr)
library(survival)
options(digits = 2)

printObject <- function(obj){
  message(gsub('_', '::', deparse(substitute(obj))))
  message(paste0('Relative Efficiency: ', format(attr(obj, 'Rel. Eff.'), digits = 3)))
  print(obj)
}

## ----echo=FALSE---------------------------------------------------------------
data(pharmacoSmoking, package = 'asaur')

asaur_pharmacoSmoking <- 
  pated(
    Surv(time = ttr, event = relapse) ~ grp, 
    age ~ grp, 
    yearsSmoking ~ grp, 
    priorAttempts ~ grp, 
    longestNoSmoke ~ grp,
    gender ~ grp, 
    I(race == 'black') ~ grp,
    I(race == 'hispanic') ~ grp,
    I(race == 'white') ~ grp,
    I(employment == 'ft') ~ grp, 
    I(employment == 'pt') ~ grp,
    I(levelSmoking == 'heavy') ~ grp,
    family = c('coxph', rep('gaussian', 4), rep('binomial', 7)),
    data = pharmacoSmoking %>% mutate(grp = ifelse(grp == 'combination', 1, 0))
  )

printObject(asaur_pharmacoSmoking)


## ----echo=FALSE---------------------------------------------------------------
data(glioma, package = 'coin')

coin_glioma <- 
  pated(
    Surv(time = time, event = event) ~ group,
    age ~ group, 
    sex ~ group,
    I(histology == 'GBM') ~ group,
    family = c('coxph', 'gaussian', rep('binomial', 2)), 
    data = glioma %>% mutate(group = ifelse(group == 'Control', 0, 1), event = 1 * event)
  )

printObject(coin_glioma)

## ----echo=FALSE---------------------------------------------------------------
data(burn, package = 'iBST')

iBST_burn <- 
  pated(
    Surv(time = T3, event = D3) ~ Z1,
    Z2 ~ Z1, 
    Z3 ~ Z1,
    Z5 ~ Z1, 
    Z6 ~ Z1, 
    Z7 ~ Z1,
    Z8 ~ Z1,
    Z9 ~ Z1,
    Z10 ~ Z1,
    I(Z11 == 1) ~ Z1, 
    I(Z11 == 2) ~ Z1,
    I(Z11 == 3) ~ Z1,
    
    Z4 ~ Z1,
    family = c('coxph', rep('binomial', 11), 'gaussian'), 
    data = burn
  )

printObject(iBST_burn)

## ----echo=FALSE---------------------------------------------------------------
data(d.oropha.rec, package = 'invGauss')

invGauss_d.oropha.rec <- 
  pated(
    Surv(time = time, event = status) ~ treatm, 
    I(sex == 1) ~ treatm,
    
    age ~ treatm,
    tstage ~ treatm,
    nstage ~ treatm,
    family = c('coxph', rep('gaussian', 1), rep('gaussian', 3)),
    data = d.oropha.rec %>% mutate(treatm = ifelse(treatm == 2, 1, 0))
  )

printObject(invGauss_d.oropha.rec)

## ----echo=FALSE---------------------------------------------------------------
data(aids.id, package = 'JM')

JM_aids.id <- 
  pated(
    Surv(time = Time, event = death) ~ drug,
    CD4 ~ drug, 
    gender ~ drug,
    I(prevOI == 'AIDS') ~ drug,
    I(AZT == 'intolerance') ~ drug,
    family = c('coxph', 'gaussian', rep('binomial', 3)),
    data = aids.id %>% mutate(drug = ifelse(drug == 'ddC', 1, 0))
  )

printObject(JM_aids.id)

## ----echo=FALSE---------------------------------------------------------------
data(actg, package = 'multipleOutcomes')

mlr3proba_actg <- 
  pated(
    Surv(time = time, event = event) ~ tx,
    strat2 ~ tx,
    sex ~ tx,
    I(ivdrug == 1) ~ tx,
    I(raceth == 1) ~ tx, 
    I(raceth == 2) ~ tx,
    I(raceth == 3) ~ tx,
    hemophil ~ tx,
    I(karnof == 100) ~ tx,
    I(karnof == 90) ~ tx,
    I(karnof == 80) ~ tx,
    I(karnof == 70) ~ tx,
    cd4 ~ tx,
    priorzdv ~ tx,
    age ~ tx,
    family = c('coxph', rep('binomial', 11), rep('gaussian', 3)),
    data = actg %>% mutate(event = 1 * (censor + censor_d > 0))
  )

printObject(mlr3proba_actg)

## ----echo=FALSE---------------------------------------------------------------
data(dataOvarian1, package = 'joint.Cox')

set.seed(123)
dat <- dataOvarian1 %>% dplyr::select(t.event, event, CXCL12, NCOA3, PDPN, TEAD1, TIMP2, YWHAB)
n <- 500
ctrl <- dat[sample(nrow(dat), n, TRUE), ] %>% mutate(grp = 0)
trt <- dat[sample(nrow(dat), n, TRUE), ] %>% mutate(t.event = t.event / .8, grp = 1)

joint.Cox_dataOvarian1 <- 
  pated(
    Surv(time = t.event, event = event) ~ grp, 
    CXCL12 ~ grp, 
    NCOA3 ~ grp,
    PDPN ~ grp,
    TEAD1 ~ grp,
    TIMP2 ~ grp,
    YWHAB ~ grp,
    family = c('coxph', rep('gaussian', 6)),
    data = rbind(ctrl, trt)
  )

printObject(joint.Cox_dataOvarian1)

## ----echo=FALSE---------------------------------------------------------------
data(Pbc3, package = 'pec')

pec_Pbc3 <- 
  pated(
    Surv(time = days, event = event) ~ tment, 
    sex ~ tment, 
    I(stage == 1) ~ tment, 
    I(stage == 2) ~ tment, 
    I(stage == 3) ~ tment, 
    I(stage == 4) ~ tment, 
    gibleed ~ tment,
    
    age ~ tment,
    crea ~ tment,
    #alb ~ tment,
    bili ~ tment,
    alkph ~ tment,
    asptr ~ tment,
    weight ~ tment,
    family = c('coxph', rep('binomial', 6), rep('gaussian', 6)),
    data = Pbc3 %>% mutate(event = ifelse(status == 0, 0, 1))
  )

printObject(pec_Pbc3)

## ----echo=FALSE---------------------------------------------------------------
data(cost, package = 'pec')
set.seed(10)
n <- 300
ctrl <- cost[sample(nrow(cost), n, TRUE), ] %>% mutate(trt = 0)
trt <- cost[sample(nrow(cost), n, TRUE), ] %>% mutate(time = time / .8, trt = 1)
dat <- rbind(ctrl, trt)

set.seed(1)
pec_cost <- 
  pated(
    Surv(time = time, event = status) ~ trt,
    age ~ trt,
    strokeScore ~ trt,
    cholest ~ trt,
    
    sex ~ trt, 
    hypTen ~ trt,
    ihd ~ trt,
    prevStroke ~ trt,
    othDisease ~ trt,
    alcohol ~ trt,
    diabetes ~ trt,
    smoke ~ trt,
    atrialFib ~ trt,
    hemor ~ trt,
    family = c('coxph', rep('gaussian', 3), rep('binomial', 10)),
    data = dat
  )

printObject(pec_cost)

## ----echo=FALSE---------------------------------------------------------------
data(GBSG2, package = 'pec')

pec_GBSG2 <- 
  pated(
    Surv(time = time, event = cens) ~ horTh,
    #age ~ horTh,
    tsize ~ horTh,
    pnodes ~ horTh,
    progrec ~ horTh,
    #estrec ~ horTh,
    
    #menostat ~ horTh,
    I(tgrade == 'I') ~ horTh,
    I(tgrade == 'II') ~ horTh,
    I(tgrade == 'III') ~ horTh,
    family = c('coxph', rep('gaussian', 3), rep('binomial', 3)),
    data = GBSG2 %>% mutate(horTh = ifelse(horTh == 'yes', 1, 0))
  )

printObject(pec_GBSG2)

## ----echo=FALSE---------------------------------------------------------------
data(follic, package = 'randomForestSRC')

randomForestSRC_follic <- 
  pated(
    Surv(time = time, event = status) ~ ch, 
    age ~ ch, 
    hgb ~ ch, 
    #clinstg ~ ch,
    family = c('coxph', rep('gaussian', 2)),#, 'binomial'),
    data = follic %>% mutate(clinstg = ifelse(clinstg == 1, 1, 0), ch = ifelse(ch == 'Y', 1, 0), status = ifelse(status == 0, 0, 1))
  )

printObject(randomForestSRC_follic)

