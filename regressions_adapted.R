
# use the pvs calculated in script "cognitiva_codes.R"
if (exists("use_new_pvs")==F) {
  use_new_pvs <- T
}

# run additional regressions for robustness/curiosity
explore <- F
# use probability-scaled plausible values (instead of those scaled in the traditional way)
use_probsums <- T

#install.packages("survey")
#install.packages("missMDA")
#install.packages("ltm")
# OECD regressions
library(survey)
library(missMDA)
library(ltm)

pvreg = function(formula, design, ...) {
  formula = as.character(formula)
  lspl = unlist(strsplit(formula[2], "\\+"))
  npv = length(lspl)
  m   = 1 + 1 /npv
  regs = lapply (lspl, function(x){
    rm = formula(paste(x, "~", formula[3]))  
    summary(svyglm(rm, design=design))
  })
  degf = regs[[1]]$df.residual
  betas = lapply(regs, function(x)x$coef)
  betas = Map(function(x){x[,2]=x[,2]^2;x}, betas)
  b = sapply(betas, function(x)x[,1])
  impuvar = apply(b, 1, var)
  br = Reduce('+',betas) / length(betas)
  br[,2] = sqrt(br[,2] + m*impuvar)
  br[,3] = br[,1] / br[,2]
  br[,4] = 2 * pt(-abs(br[,3]), degf)
  list(summary=br, regressions=regs)
}

regfun = function(x, test){
  ind = paste0(x, collapse='+')
  dep = paste0(paste0(paste(paste0('PV',1:5),test,sep='_'),collapse='+'),'~')
  formula = as.formula(paste0(dep,ind))
  pvreg(formula, frd, na.rm=T)
}

#de-activated: setwd('~/WD/oecd')
s = readRDS('msStudentsWithPCA.rds')
s = subset(s, country_id != "UK-ENG" & main_study_sample=="MS" & targetLanguage_id=="EN")
s = subset(s, !is.na(JKzone))
dim(s)

table(s$PL1_LIST, s$PL1_READ)
sum(table(s$PL1_LIST, s$PL1_READ))

#############################################################################
### replace old with new pvs

if (use_new_pvs == T & use_probsums==F) {
  pvdf <- readRDS("pvdf.rds")
  s <- merge(s, pvdf, all = T)
  colnames(s)
  s$PV1_LIST <- s$listpv1
  s$PV2_LIST <- s$listpv2
  s$PV3_LIST <- s$listpv3
  s$PV4_LIST <- s$listpv4
  s$PV5_LIST <- s$listpv5
  s$PV1_READ <- s$readpv1
  s$PV2_READ <- s$readpv2
  s$PV3_READ <- s$readpv3
  s$PV4_READ <- s$readpv4
  s$PV5_READ <- s$readpv5
  s <- s[,1:524]
}

if (use_new_pvs == T & use_probsums==T) {
  pvdf <- readRDS("pvdf.rds")
  s <- merge(s, pvdf, all = T)
  colnames(s)
  s$PV1_LIST <- s$listpsum1
  s$PV2_LIST <- s$listpsum2
  s$PV3_LIST <- s$listpsum3
  s$PV4_LIST <- s$listpsum4
  s$PV5_LIST <- s$listpsum5
  s$PV1_READ <- s$readpsum1
  s$PV2_READ <- s$readpsum2
  s$PV3_READ <- s$readpsum3
  s$PV4_READ <- s$readpsum4
  s$PV5_READ <- s$readpsum5
  s <- s[,1:524]
}


str(s$PV5_READ)
#View(s[,395:420])

### until here: replace old with new pvs
#############################################################################







# gen avg pv for depvar
s$read <- (s$PV1_READ + s$PV2_READ + s$PV3_READ + s$PV4_READ + s$PV5_READ) / 5
s$list <- (s$PV1_LIST + s$PV2_LIST + s$PV3_LIST + s$PV4_LIST + s$PV5_LIST) / 5
s$writ <- (s$PV1_WRIT_C + s$PV2_WRIT_C + s$PV3_WRIT_C + s$PV4_WRIT_C + s$PV5_WRIT_C) / 5
s$listsocio <- s$list * s$I08_ST_A_S19B
# summarise
table(s$PL1_LIST) / sum(as.numeric(table(s$PL1_LIST)))
table(s$PL1_READ) / sum(as.numeric(table(s$PL1_READ)))
table(s$PL1_WRIT_C) / sum(as.numeric(table(s$PL1_WRIT_C)))
c(mean(s$PV1_LIST[is.na(s$PV1_LIST)==FALSE]), sd(s$PV1_LIST[is.na(s$PV1_LIST)==FALSE]))
c(mean(s$PV1_READ[is.na(s$PV1_READ)==FALSE]), sd(s$PV1_READ[is.na(s$PV1_READ)==FALSE]))
c(mean(s$PV1_WRIT_C[is.na(s$PV1_WRIT_C)==FALSE]), sd(s$PV1_WRIT_C[is.na(s$PV1_WRIT_C)==FALSE]))
c(mean(s$PV1_WRIT2[is.na(s$PV1_WRIT2)==FALSE]), sd(s$PV1_WRIT2[is.na(s$PV1_WRIT2)==FALSE]))
c(mean(s$I08_ST_A_S19B[is.na(s$I08_ST_A_S19B)==FALSE]), sd(s$I08_ST_A_S19B[is.na(s$I08_ST_A_S19B)==FALSE]))
c(mean(s$PV1_READ[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]), sd(s$PV1_READ[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]))
c(mean(s$PV1_LIST[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]), sd(s$PV1_LIST[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]))
c(mean(s$read[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]), sd(s$read[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]))
c(mean(s$list[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]), sd(s$list[is.na(s$PV1_READ)==FALSE&is.na(s$PV1_LIST)==FALSE]))


#View(s[,350:400])
#View(s[,500:dim(s)[2]])
# students
table(s$PL1_LIST,s$PL1_READ)



# new pv diff
s$PV1_dreadlist <- s$PV1_READ - s$PV1_LIST
s$PV2_dreadlist <- s$PV2_READ - s$PV2_LIST
s$PV3_dreadlist <- s$PV3_READ - s$PV3_LIST
s$PV4_dreadlist <- s$PV4_READ - s$PV4_LIST
s$PV5_dreadlist <- s$PV5_READ - s$PV5_LIST
s$dreadlist <- (s$PV1_dreadlist + s$PV2_dreadlist + s$PV3_dreadlist + s$PV4_dreadlist + s$PV5_dreadlist) / 5 
s$FSW_dreadlist <- 0
s$FSW_dreadlist[is.na(s$dreadlist)==FALSE] <- 1
summary(s$dreadlist)


### generate numeric PL variables, in case they turn convenient

# run the following line if you need to delete the old numeric PL variables (if you already have them in the dataset s)
#s <- s[,1:535]

for (myskill in c("READ", "LIST", "WRIT_C")) {
  for (i in 1:5) {
    s$temp <- NA
    s$temp[eval(str2expression(paste0("s$PL",i,"_",myskill)))=="-A1"] <- 0
    s$temp[eval(str2expression(paste0("s$PL",i,"_",myskill)))=="A1"] <- 1
    s$temp[eval(str2expression(paste0("s$PL",i,"_",myskill)))=="A2"] <- 2
    s$temp[eval(str2expression(paste0("s$PL",i,"_",myskill)))=="B1"] <- 3
    s$temp[eval(str2expression(paste0("s$PL",i,"_",myskill)))=="B2"] <- 4
    colnames(s)[ncol(s)] <- paste0("numPL",i,"_",myskill)
  }
}
head(s[,536:dim(s)[2]])
table(s$numPL1_LIST[is.na(s$numPL1_READ)==F])/sum(table(s$numPL1_LIST[is.na(s$numPL1_READ)==F]))
table(s$numPL1_READ[is.na(s$numPL1_LIST)==F])/sum(table(s$numPL1_READ[is.na(s$numPL1_LIST)==F]))

#############################################################################
### replace old with new pls

if (use_new_pvs == T) {
  pvdf <- readRDS("pvdf.rds")
  s <- merge(s, pvdf, all = T)
  colnames(s)
  s$numPL1_LIST <- s$listlevpv1
  s$numPL2_LIST <- s$listlevpv2
  s$numPL3_LIST <- s$listlevpv3
  s$numPL4_LIST <- s$listlevpv4
  s$numPL5_LIST <- s$listlevpv5
  s$numPL1_READ <- s$readlevpv1
  s$numPL2_READ <- s$readlevpv2
  s$numPL3_READ <- s$readlevpv3
  s$numPL4_READ <- s$readlevpv4
  s$numPL5_READ <- s$readlevpv5
  s <- s[,1:550]
}

str(s$PV5_READ)
#View(s[,395:420])

### until here: replace old with new pls
#############################################################################





### calculate thresholds used by SurveyLang for transformation from PV to PL

myv <- as.data.frame(cbind(0,0))
rownames(myv)[nrow(myv)] <- "row0"
colnames(myv) <- c("min","max")
for (myskill in c("READ", "LIST", "WRIT_C")) {
  for (mylevel in 0:4) {
    for (i in 1:5) {
      print(paste0("s$PV",i,"_",myskill, "level: ", mylevel))
      temp <- cbind( min(eval(str2expression(paste0("s$PV",i,"_",myskill)))[eval(str2expression(paste0("s$numPL",i,"_",myskill)))==as.numeric(mylevel)&is.na(eval(str2expression(paste0("s$numPL",i,"_",myskill))))==FALSE&is.na(eval(str2expression(paste0("s$PV",i,"_",myskill))))==FALSE]) , max(eval(str2expression(paste0("s$PV",i,"_",myskill)))[eval(str2expression(paste0("s$numPL",i,"_",myskill)))==as.numeric(mylevel)&is.na(eval(str2expression(paste0("s$numPL",i,"_",myskill))))==FALSE&is.na(eval(str2expression(paste0("s$PV",i,"_",myskill))))==FALSE]) )
      colnames(temp) <- c("min","max")
      rownames(temp) <- paste0("s$PV",i,"_",myskill,"_numPL_",mylevel)
      myv <- rbind(myv,temp)
    }
  }
}

myv <- myv[-1,]
head(myv)
write.csv(myv,"PLthresholds.csv")
max(s$PV1_LIST[s$numPL1_LIST==0&is.na(s$numPL1_LIST)==FALSE&is.na(s$PV1_LIST)==FALSE])


### generate z-standardised plausible values

for (myskill in c("READ","LIST","WRIT_C")) {
  for (myval in 1:5) {
    tempv <- eval(str2expression(paste0("s$PV",myval,"_",myskill)))
    s$temp <- (tempv - mean(tempv, na.rm=T)) / sd(tempv, na.rm=T)
    colnames(s)[ncol(s)] <- paste0("PV",myval,"z_",myskill)
  }
}

#View(s[,550:dim(s)[2]])


### generate CEFR-anchored plausible values
scol <- dim(s)[2]


# calculate thresholds for CEFR-level-based standardisation

thrb_READ <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="1"]), max(myv$max[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="0"]))   )
thrc_READ <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="2"]), max(myv$max[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="1"]))   )
thrd_READ <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="3"]), max(myv$max[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="2"]))   )
thre_READ <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="4"]), max(myv$max[substr(rownames(myv),7,10)=="READ" & substr(rownames(myv),18,18)=="3"]))   )
thrb_LIST <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="1"]), max(myv$max[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="0"]))   )
thrc_LIST <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="2"]), max(myv$max[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="1"]))   )
thrd_LIST <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="3"]), max(myv$max[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="2"]))   )
thre_LIST <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="4"]), max(myv$max[substr(rownames(myv),7,10)=="LIST" & substr(rownames(myv),18,18)=="3"]))   )
thrb_WRIT_C <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="1"]), max(myv$max[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="0"]))   )
thrc_WRIT_C <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="2"]), max(myv$max[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="1"]))   )
thrd_WRIT_C <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="3"]), max(myv$max[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="2"]))   )
thre_WRIT_C <- mean(   c(min(myv$min[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="4"]), max(myv$max[substr(rownames(myv),7,10)=="WRIT" & substr(rownames(myv),20,20)=="3"]))   )

# extreme ("entry" and "exit") thresholds (a and f)
thra_READ <- mean(quantile(s$PV1_READ[s$numPL1_READ==0], .01, na.rm = T), quantile(s$PV2_READ[s$numPL2_READ==0], .01, na.rm = T), quantile(s$PV3_READ[s$numPL3_READ==0], .01, na.rm = T), quantile(s$PV4_READ[s$numPL4_READ==0], .01, na.rm = T), quantile(s$PV5_READ[s$numPL5_READ==0], .01, na.rm = T)) 
thrf_READ <- mean(quantile(s$PV1_READ[s$numPL1_READ==4], .99, na.rm = T), quantile(s$PV2_READ[s$numPL2_READ==4], .99, na.rm = T), quantile(s$PV3_READ[s$numPL3_READ==4], .99, na.rm = T), quantile(s$PV4_READ[s$numPL4_READ==4], .99, na.rm = T), quantile(s$PV5_READ[s$numPL5_READ==4], .99, na.rm = T)) 
thra_WRIT_C <- mean(quantile(s$PV1_WRIT_C[s$numPL1_WRIT_C==0], .01, na.rm = T), quantile(s$PV2_WRIT_C[s$numPL2_WRIT_C==0], .01, na.rm = T), quantile(s$PV3_WRIT_C[s$numPL3_WRIT_C==0], .01, na.rm = T), quantile(s$PV4_WRIT_C[s$numPL4_WRIT_C==0], .01, na.rm = T), quantile(s$PV5_WRIT_C[s$numPL5_WRIT_C==0], .01, na.rm = T)) 
thrf_WRIT_C <- mean(quantile(s$PV1_WRIT_C[s$numPL1_WRIT_C==4], .99, na.rm = T), quantile(s$PV2_WRIT_C[s$numPL2_WRIT_C==4], .99, na.rm = T), quantile(s$PV3_WRIT_C[s$numPL3_WRIT_C==4], .99, na.rm = T), quantile(s$PV4_WRIT_C[s$numPL4_WRIT_C==4], .99, na.rm = T), quantile(s$PV5_WRIT_C[s$numPL5_WRIT_C==4], .99, na.rm = T)) 
thra_LIST <- mean(quantile(s$PV1_LIST[s$numPL1_LIST==0], .01, na.rm = T), quantile(s$PV2_LIST[s$numPL2_LIST==0], .01, na.rm = T), quantile(s$PV3_LIST[s$numPL3_LIST==0], .01, na.rm = T), quantile(s$PV4_LIST[s$numPL4_LIST==0], .01, na.rm = T), quantile(s$PV5_LIST[s$numPL5_LIST==0], .01, na.rm = T)) 
thrf_LIST <- mean(quantile(s$PV1_LIST[s$numPL1_LIST==4], .99, na.rm = T), quantile(s$PV2_LIST[s$numPL2_LIST==4], .99, na.rm = T), quantile(s$PV3_LIST[s$numPL3_LIST==4], .99, na.rm = T), quantile(s$PV4_LIST[s$numPL4_LIST==4], .99, na.rm = T), quantile(s$PV5_LIST[s$numPL5_LIST==4], .99, na.rm = T)) 
thr_READ <- c(thra_READ, thrb_READ, thrc_READ, thrd_READ, thre_READ, thrf_READ)
thr_LIST <- c(thra_LIST, thrb_LIST, thrc_LIST, thrd_LIST, thre_LIST, thrf_LIST)
thr_WRIT_C <- c(thra_WRIT_C, thrb_WRIT_C, thrc_WRIT_C, thrd_WRIT_C, thre_WRIT_C, thrf_WRIT_C)



# to choose the correlation-with-zPV-maximising thresholds (a and f) for for extreme levels ("entry" and "exit" thresholds), you need the results from the code "choosepars.R" within this R project
#source("choosepars.R")
# for listening, a=.34 and f=.38
# for reading, a=.38 and f=.46



# gen cefr-anchored PV levels

s <- s[,1:scol]
if (use_new_pvs==F) {
  for (myskill in c("READ","LIST","WRIT_C")) {
    for (myval in 1:5) {
      tempv <- eval(str2expression(paste0("s$PV",myval,"_",myskill)))
      templ <- eval(str2expression(paste0("s$numPL",myval,"_",myskill)))
      thr <- eval(str2expression(paste0("thr_",myskill)))
      for (mylev in 0:4) {
        #myskill <- "READ"
        #myval <- 5
        #mylev <- 2
        c1factor <- ifelse(mylev==4, 2, 1)
        s$temp[templ==mylev&is.na(templ)==F&is.na(tempv)==F] <- mylev + c1factor * (tempv[templ==mylev&is.na(templ)==F&is.na(tempv)==F]-thr[mylev+1]) / (thr[mylev+2] - thr[mylev+1])
        print(c( paste0(myskill, myval, mylev),c1factor, quantile(s$temp[templ==mylev&is.na(templ)==F&is.na(tempv)==F], probs=.01), quantile(s$temp[templ==mylev&is.na(templ)==F&is.na(tempv)==F], probs=.99) ))
      }
      colnames(s)[ncol(s)] <- paste0("PV",myval,"cefr_",myskill)
    }
  }
}

s$PV1cefr_READ[1:10]

if (use_new_pvs==T & use_probsums==F) {
  
  for (myskill in c("READ","LIST","WRIT_C")) {
    for (myval in 1:5) {
      tempv <- eval(str2expression(paste0("s$PV",myval,"_",myskill)))
      templ <- eval(str2expression(paste0("s$numPL",myval,"_",myskill)))
      thr <- eval(str2expression(paste0("thr_",myskill)))
      for (mylev in 0:4) {
        #myskill <- "READ"
        #myval <- 5
        #mylev <- 2
        c1factor <- ifelse(mylev==4, 2, 1)
        s$temp[templ==mylev&is.na(templ)==F&is.na(tempv)==F] <- mylev + c1factor * (tempv[templ==mylev&is.na(templ)==F&is.na(tempv)==F] - min(tempv[templ==mylev&is.na(templ)==F&is.na(tempv)==F], na.rm=T) ) / ( max(tempv[templ==mylev&is.na(templ)==F&is.na(tempv)==F], na.rm=T) - min(tempv[templ==mylev&is.na(templ)==F&is.na(tempv)==F], na.rm=T) )
        print(c( paste0(myskill, myval, mylev),c1factor, quantile(s$temp[templ==mylev&is.na(templ)==F&is.na(tempv)==F], probs=.01), quantile(s$temp[templ==mylev&is.na(templ)==F&is.na(tempv)==F], probs=.99) ))
      }
      colnames(s)[ncol(s)] <- paste0("PV",myval,"cefr_",myskill)
    }
  }
  
}


if (use_new_pvs==T & use_probsums==T) {
  s$PV1cefr_READ <- s$PV1_READ
  s$PV2cefr_READ <- s$PV2_READ
  s$PV3cefr_READ <- s$PV3_READ
  s$PV4cefr_READ <- s$PV4_READ
  s$PV5cefr_READ <- s$PV5_READ
  s$PV1cefr_LIST <- s$PV1_LIST
  s$PV2cefr_LIST <- s$PV2_LIST
  s$PV3cefr_LIST <- s$PV3_LIST
  s$PV4cefr_LIST <- s$PV4_LIST
  s$PV5cefr_LIST <- s$PV5_LIST
}


# compare plots and cor
with(s, plot(PV1cefr_READ, PV1cefr_LIST))
c(mean(s$PV1cefr_READ, na.rm=T), sd(s$PV1cefr_READ, na.rm=T))
c(mean(s$PV1cefr_LIST, na.rm=T), sd(s$PV1cefr_LIST, na.rm=T))
cor(s$PV1_READ, s$PV1_LIST, use="complete.obs")
with(s, plot(PV1_READ, PV1_LIST))
with(s, plot(PV1z_READ, PV1z_LIST))
cor(s$PV1_READ, s$PV1_LIST, use="complete.obs")
cor(s$PV1cefr_READ, s$PV1cefr_LIST, use="complete.obs")
table(s$numPL1_LIST, s$numPL1_READ)
table(floor(s$PV1_READ*10)/10)
table(floor(s$PV1_READ), floor(s$PV1_LIST))


### gen diff of cefr-anchored and z-standardised PV levels


for (myskill in c("cefr_READ-cefr_LIST","z_READ-z_LIST")) {
  for (myval in 1:5) {
    myskills <- unlist(strsplit(myskill, split="-"))
    temp1 <- eval(str2expression(paste0("s$PV",myval,myskills[1])))
    temp2 <- eval(str2expression(paste0("s$PV",myval,myskills[2])))
    s$temp <- temp1 - temp2
    colnames(s)[ncol(s)] <- paste0("PV",myval,"_",gsub("_","",myskills[1]),gsub("_","",myskills[2]))
  }
}


# gen diff of numeric plausible levels (PLs)


for (myskill in c("_READ-_LIST","_WRIT_C-_READ")) {
  for (myval in 1:5) {
    myskills <- unlist(strsplit(myskill, split="-"))
    temp1 <- eval(str2expression(paste0("s$numPL",myval,myskills[1])))
    temp2 <- eval(str2expression(paste0("s$numPL",myval,myskills[2])))
    s$temp <- temp1 - temp2
    colnames(s)[ncol(s)] <- paste0("dPL",myval,"_",gsub("_","",myskills[1]),gsub("_","",myskills[2]))
  }
}



#  gen avg diffs for cefr-anchored and z-standardised PV levels
s$PV_zREADzLIST = (s$PV1_zREADzLIST + s$PV2_zREADzLIST + s$PV3_zREADzLIST + s$PV4_zREADzLIST + s$PV5_zREADzLIST) / 5
s$PV_cefrREADcefrLIST = (s$PV1_cefrREADcefrLIST + s$PV2_cefrREADcefrLIST + s$PV3_cefrREADcefrLIST + s$PV4_cefrREADcefrLIST + s$PV5_cefrREADcefrLIST) / 5
s$dPL_READLIST = (s$dPL1_READLIST + s$dPL2_READLIST + s$dPL3_READLIST + s$dPL4_READLIST + s$dPL5_READLIST) / 5
with(s, plot(PV_zREADzLIST, dPL_READLIST))
with(s, plot(PV_zREADzLIST, s$PV_cefrREADcefrLIST))
with(s, plot(s$PV_cefrREADcefrLIST, dPL_READLIST))


# couple of histograms and sumstats to check 
hist(s$PV1_dreadlist, breaks = 50)
dim(s[is.na(s$PV1_dreadlist)==F,])
c(mean(s$PV1_dreadlist, na.rm=T), sd(s$PV1_dreadlist, na.rm=T))
grand_mean <- mean(c(s$PV1_dreadlist,s$PV2_dreadlist,s$PV3_dreadlist,s$PV4_dreadlist,s$PV5_dreadlist), na.rm=T)
sd_within <- mean(c(sd(s$PV1_dreadlist, na.rm=T),sd(s$PV2_dreadlist, na.rm=T),sd(s$PV3_dreadlist, na.rm=T),sd(s$PV4_dreadlist, na.rm=T),sd(s$PV5_dreadlist, na.rm=T)), na.rm=T)
sd_between <- sqrt(   .2*(mean(s$PV1_dreadlist, na.rm=T)^2 + mean(s$PV2_dreadlist, na.rm=T)^2 + mean(s$PV3_dreadlist, na.rm=T)^2 + mean(s$PV4_dreadlist, na.rm=T)^2 + mean(s$PV5_dreadlist, na.rm=T)^2) - grand_mean^2   )
grand_sd <- sd_within + sd_between
grand_se <- grand_sd / sqrt(dim(s[is.na(s$PV1_dreadlist)==F,])[1])
print(paste0("grand_mean: ", round(grand_mean,5), "; grand_sd: ", round(grand_sd,5), "; grand_se: ", round(grand_se,5)))

summary(s$PV1_dreadlist)
summary(s$PV2_dreadlist)
sum(s$PV1_dreadlist<0, na.rm = T) / sum(s$PV1_dreadlist<9999999 , na.rm = T)

hist(s$PV1_cefrREADcefrLIST, breaks = 50)
c(mean(s$PV1_cefrREADcefrLIST, na.rm=T), sd(s$PV1_cefrREADcefrLIST, na.rm=T))
summary(s$PV1_cefrREADcefrLIST)
sum(s$PV1_cefrREADcefrLIST<0, na.rm = T) / sum(s$PV1_cefrREADcefrLIST<9999999 , na.rm = T)
c(mean(s$PV1_zREADzLIST, na.rm=T), sd(s$PV1_zREADzLIST, na.rm=T))
meanw_read <- mean(s$FSW_READ[is.na(s$PV1_cefrREADcefrLIST)==F])
c(mean(s$PV1_cefrREADcefrLIST*s$FSW_READ/meanw_read, na.rm=T), sd(s$PV1_cefrREADcefrLIST*s$FSW_READ/meanw_read, na.rm=T))

# sumstats by country
host_table <- data.frame(mean=0, stdev=0, stderr=0)
for (nat in c('BE de', 'BE fr', 'BE nl', 'BG','EE','EL','ES','FR','HR','MT','NL','PL','PT','SE','SI')) {
  print(nat)
  natresults <- as.data.frame(cbind(mean(s$PV1_dreadlist[s$cnt==nat], na.rm=T), sd(s$PV1_dreadlist[s$cnt==nat], na.rm=T), sd(s$PV1_dreadlist[s$cnt==nat], na.rm=T)/sqrt(length(s$PV1_dreadlist[s$cnt==nat&is.na(s$PV1_dreadlist)==F])) ))
  colnames(natresults) <- c("mean", "stdev", "stderr")
  rownames(natresults) <- nat
  host_table <- rbind(host_table,natresults)
}
write.csv(host_table,"national_means.csv")

hist(s$PV_zREADzLIST)
c(mean(s$PV_zREADzLIST, na.rm=T), sd(s$PV_zREADzLIST, na.rm=T))
summary(s$PV_zREADzLIST)
hist(s$dPL_READLIST)
c(mean(s$dPL_READLIST, na.rm=T), sd(s$dPL_READLIST, na.rm=T))



#############################################################################
### doing some work on control variables
#############################################################################

colnames(s)
# motivation items, Q33
alpha_motivation <- cronbach.alpha(s[,129:138], na.rm=T)
alpha_academic <- cronbach.alpha(s[,c(142,146,147)], na.rm=T)
alpha_sport <- cronbach.alpha(s[,c(147,156,165)], na.rm=T)
alpha_vocat <- cronbach.alpha(s[,c(146,155,164)], na.rm=T)
alpha_cultu <- cronbach.alpha(s[,c(142,151,160)], na.rm=T)
alpha_fl <- cronbach.alpha(s[,c(145,154,163)], na.rm=T)
alpha_tl <- cronbach.alpha(s[,c(144,153,162)], na.rm=T)
alpha_allfl <- cronbach.alpha(s[,c(144,153,162,145,154,163)], na.rm=T)
alpha_spvoc <- cronbach.alpha(s[,c(146,155,164,147,156,165)], na.rm=T)
alpha_homexposure <- cronbach.alpha(s[,c(364, 362, 366, 367)], na.rm=T)
alpha_oosread <- cronbach.alpha(s[,c(103,122, 123, 124)], na.rm=T)
alpha_classlistening <- cronbach.alpha(s[,c(290, 300, 289, 299, 293, 303)], na.rm=T)
nsb_classlistening <- 10/length(c(292, 302, 288, 298))
sbalpha_classlistening <- as.numeric(alpha_classlistening[1])*nsb_classlistening / (1+as.numeric(alpha_classlistening[1])*(nsb_classlistening-1))
alpha_classreading <- cronbach.alpha(s[,c(292, 302, 288, 298)], na.rm=T)
nsb_classreading <- 10/length(c(292, 302, 288, 298))
sbalpha_classreading <- as.numeric(alpha_classreading[1])*nsb_classreading / (1+as.numeric(alpha_classreading[1])*(nsb_classreading-1))
alpha_ictclass <- cronbach.alpha(s[,c(248, 249, 250)], na.rm=T) # you could add 305, 306 as needed
nsb_ictclass <- 10/length(c(248, 249, 250))
sbalpha_ictclass <- as.numeric(alpha_ictclass[1])*nsb_ictclass / (1+as.numeric(alpha_ictclass[1])*(nsb_ictclass-1))
alpha_ictoos <- cronbach.alpha(s[,c(103, 108, 110, 115)], na.rm=T)
alpha_tluse <- cronbach.alpha(s[,c(242, 243, 244)], na.rm=T) # you could add 240, 241 (but if you don't, you can use the ready-made index) 
alpha_goodclass <- cronbach.alpha(s[,c(267, 268, 269, 270, 271, 273, 274, 275)], na.rm=T) 
nsb_goodclass <- 10/length(c(267, 268, 269, 270, 271, 273, 274, 275))
sbalpha_goodclass <- as.numeric(alpha_goodclass[1])*nsb_goodclass / (1+as.numeric(alpha_goodclass[1])*(nsb_goodclass-1))
alpha_face2face <- cronbach.alpha(s[,c(106, 107, 113, 114)], na.rm=T) # you could add 305, 306 as needed
nsb_face2face <- 10/length(c(106, 107, 113, 114))
sbalpha_face2face <- as.numeric(alpha_face2face[1])*nsb_face2face / (1+as.numeric(alpha_face2face[1])*(nsb_face2face-1))
alpha_txtbclass <- cronbach.alpha(s[,c(251, 257, 260)], na.rm=T) # you could add 305, 306 as needed
nsb_txtbclass <- 10/length(c(251, 257, 260))
sbalpha_txtbclass <- as.numeric(alpha_txtbclass[1])*nsb_txtbclass / (1+as.numeric(alpha_txtbclass[1])*(nsb_txtbclass-1))
alpha_readoos <- cronbach.alpha(s[,c(103, 110, 122, 123, 124)], na.rm=T)
nsb_readoos <- 10/length(c(103, 110, 122, 123, 124))
sbalpha_readoos <- as.numeric(alpha_readoos[1])*nsb_readoos / (1+as.numeric(alpha_readoos[1])*(nsb_readoos-1))
alpha_listoos <- cronbach.alpha(s[,c(117, 118, 119, 120, 121)], na.rm=T)
nsb_listoos <- 10/length(c(117, 118, 119, 120, 121))
sbalpha_listoos <- as.numeric(alpha_listoos[1])*nsb_listoos / (1+as.numeric(alpha_listoos[1])*(nsb_listoos-1))




s$d57 <- s$SQt57i05 - s$SQt57i03  
s$d58 <- s$SQt58i05 - s$SQt58i03   
s$d61 <- s$SQt61i05 - s$SQt61i03  
alpha_classdiff <- cronbach.alpha(as.data.frame(cbind(s$d57, s$d58, s$d61)), na.rm=T)

# multiple imputation function, run for motivation and other indexes

multipleimpute <- function(indexcols) {
  colrelavg <- as.data.frame(rbind(rep(NA,length(indexcols))))
  colnames(colrelavg) <- indexcols
  s$tempsum <- 0
  s$tempcount <- 0
  s$temp <- NA
  for (r in indexcols) {
    colrelavg[1,match(r, indexcols)] <- mean(s[,r], na.rm=T)
    s$tempsum[is.na(s[,r])==F] <- s$tempsum[is.na(s[,r])==F] + s[is.na(s[,r])==F,r]
    s$tempcount[is.na(s[,r])==F] <- s$tempcount[is.na(s[,r])==F] + 1
  }
  colrelavg <- colrelavg / mean(as.numeric(colrelavg))
  s$temp[s$tempcount>=1] <- s$tempsum[s$tempcount>=1] / s$tempcount[s$tempcount>=1]
  for (r in indexcols) {
    s[is.na(s[,r])==T,r] <- s$temp[is.na(s[,r])==T] * colrelavg[1,match(r, indexcols)]
  }
  return(s)
}

# motivation index Q33
s <- multipleimpute( c(129,130,131,132,133,134,135,136,137, 138) )

# interest in fl subjects
s <- multipleimpute(c(144,153,162,145,154,163))
#View(s[,c(144,153,162,145,154,163)])

# interests in sports and vocational
s <- multipleimpute(c(146,155,164,147,156,165))
#View(s[,c(146,155,164,147,156,165)])

# oos reading
s <- multipleimpute(c(103,122, 123, 124))
#View(s[,c(103,122, 123, 124)])

# class listening
s <- multipleimpute(c(283, 290, 300, 282, 289, 299, 286, 293, 303))

# class reading
s <- multipleimpute(c(285, 292, 302, 281, 288, 298))

# travel to tl countries
s <- multipleimpute(c(219, 221))

# ict in class
s <- multipleimpute(c(248, 249, 250))

# the good class
s <- multipleimpute(c(240, 241, 242, 243, 244))
s <- multipleimpute(c(267, 268, 269, 270, 271, 273, 274, 275))

# face2face
s <- multipleimpute(c(106, 107, 113, 114))

# txtb in class
s <- multipleimpute(c(251, 257, 260))


# read out of class
s <- multipleimpute(c(103, 110, 122, 123, 124))

# list out of class
s <- multipleimpute(c(117, 118, 119, 120, 121))



### pca and other indices

# extrinsic motivation - works not

s$natest <- is.na(s[,129]) + is.na(s[,130]) + is.na(s[,131]) + is.na(s[,132]) + is.na(s[,133]) + is.na(s[,134]) + is.na(s[,135]) + is.na(s[,136]) + is.na(s[,137])
pca <- prcomp(s[s$natest==0,129:137], center = TRUE,scale. = TRUE)
summary(pca)
str(pca)
s$motivationQ33 <- NA
s$motivationQ33[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# interest in fl subjects
s$natest <- is.na(s[,144])
pca <- prcomp(s[s$natest==0,c(144,153,162,145,154,163)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$interestflQ34 <- NA
s$interestflQ34[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# interest in vocational/sport subjects
s$natest <- is.na(s[,146])
pca <- prcomp(s[s$natest==0,c(146,155,164,147,156,165)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$interestvsQ34 <- NA
s$interestvsQ34[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# oos reading - works not
s$natest <- is.na(s[,103])
pca <- prcomp(s[s$natest==0,c(103, 122, 123, 124)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$oosreadQ29 <- NA
s$oosreadQ29[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)
s$oosreadQ29 <- NA
s$oosreadQ29max <- apply(cbind(s[,103], s[,122], s[,123], s[,124]),FUN=max, 1)
summary(s$oosreadQ29max)

s$homeexposureQ4 <- s$I03_ST_A_S25B + s$I03_ST_A_S04B + s$I03_ST_A_S26B + s$I03_ST_A_S27B

# class listening
s$natest <- is.na(s[,283])
pca <- prcomp(s[s$natest==0,c(290, 300, 289, 299, 293, 303)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$classoralQ58 <- NA
s$classoralQ58[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)


# class reading
s$natest <- is.na(s[,285])
pca <- prcomp(s[s$natest==0,c(292, 302, 288, 298)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$classwrittenQ58 <- NA
s$classwrittenQ58[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)
cor(s$classwrittenQ58, s$classoralQ58, use="complete.obs")

# class listening & reading combined
s$natest <- is.na(s[,283]) | is.na(s[,285])
pca <- prcomp(s[s$natest==0,c(283, 290, 300, 282, 289, 299, 286, 293, 303,    285, 292, 302, 281, 288, 298)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$classall <- NA
s$classall[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# ict use in class
s$natest <- is.na(s[,248])
pca <- prcomp(s[s$natest==0,c(248, 249, 250)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$classictQ51 <- NA
s$classictQ51[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)


# the good class with use of tl included -> compare with: I09_ST_M_S33B and I09_IN_M_S50A
s$natest <- is.na(s[,240]) | is.na(s[,267]) 
pca <- prcomp(s[s$natest==0,c(240, 241, 242, 243, 244,    267, 268, 269, 270, 271, 273, 274, 275)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$goodclassQ54 <- NA
s$goodclassQ54[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)


# the good class without use of tl -> compare with: I09_ST_M_S33B
s$natest <- is.na(s[,240]) | is.na(s[,267]) 
pca <- prcomp(s[s$natest==0,c(267, 268, 269, 270, 271, 273, 274, 275)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$goodclass2Q54 <- NA
s$goodclass2Q54[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# use of tl -> compare with: I09_IN_M_S50A
s$natest <- is.na(s[,240]) | is.na(s[,267]) 
pca <- prcomp(s[s$natest==0,c(240, 241, 242, 243, 244)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$usetl <- NA
s$usetl[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# face to face English use
s$natest <- is.na(s[,106])
pca <- prcomp(s[s$natest==0,c(106, 107, 113, 114)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$face2faceQ29 <- NA
s$face2faceQ29[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

# txtb use in class
s$natest <- is.na(s[,251])
pca <- prcomp(s[s$natest==0,c(251, 257, 260)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$txtbQ52 <- NA
s$txtbQ52[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

s$natest <- is.na(s[,103])
pca <- prcomp(s[s$natest==0,c(103, 110, 122, 123, 124)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$readoosQ31 <- NA
s$readoosQ31[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)



s$natest <- is.na(s[,117])
pca <- prcomp(s[s$natest==0,c(117, 118, 119, 120, 121)], center = TRUE,scale. = TRUE)
summary(pca)
pca$rotation
s$listoosQ31 <- NA
s$listoosQ31[s$natest==0] <- pca$x[,1] / sd(pca$x[,1], na.rm = T)

s$SQt43i01_cleaned <- NA
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01<40.5] <- 40
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=40.5&s$SQt43i01<45.5] <- 45
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=45.5&s$SQt43i01<50.5] <- 50
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=50.5&s$SQt43i01<55.5] <- 55
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=55.5&s$SQt43i01<60.5] <- 60
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=60.5&s$SQt43i01<80.5] <- 80
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=80.5&s$SQt43i01<90.5] <- 90
s$SQt43i01_cleaned[is.na(s$SQt43i01)==F&s$SQt43i01>=90.5&s$SQt43i01<120.5] <- 120
s$I01_ST_M_S44A_cont <- s$SQt43i01_cleaned * s$SQt44i01 / 60
s$I01_ST_M_S44A_cont[is.na(s$I01_ST_M_S44A_cont)==F&s$I01_ST_M_S44A_cont>10] <- 10
s$I01_ST_M_S44A_cont[is.na(s$I01_ST_M_S44A_cont)==F&s$I01_ST_M_S44A_cont==0] <- NA
summary(s$I01_ST_M_S44A_cont)
summary(s$I01_ST_M_S44A)
s$I01_ST_M_S44A_cont <- log(s$I01_ST_M_S44A_cont)
cor(s$I01_ST_M_S44A, s$PV_cefrREADcefrLIST, use="complete.obs")
#with(s, plot(I01_ST_M_S44A_cont, PV_cefrREADcefrLIST))
#s$I01_ST_M_S40B_log <- log( rowSums(s[,194:208], na.rm = T) )


# travel to tl countries
s$travelQ45 <- s$SQt45i01 + s$SQt45i03

# difficult to learn, read - speak
s$diffQ48 <- s$SQt48i05 - s$SQt48i03



############################################################
### set up weights
############################################################


s$unweighted <- 1
#s <- s[,c(1:626)]
average_weight <- function(mycountry, fsw="READ") {
  if (fsw=="READ") {avw <- mean(s$FSW_READ[s$country_id==mycountry&is.na(s$PV1_dreadlist)==F], na.rm=T)}
  if (fsw=="QUES") {avw <- mean(s$FSW_QUES[s$country_id==mycountry&is.na(s$PV1_dreadlist)==F], na.rm=T)}
  if (fsw=="LIST") {avw <- mean(s$FSW_LIST[s$country_id==mycountry&is.na(s$PV1_dreadlist)==F], na.rm=T)}
  mydf <- data.frame(country_id=mycountry, avw=as.numeric(avw))
  colnames(mydf) <- c("country_id", paste0("avw",fsw))
  return(mydf)
} 
avw_read <- do.call(rbind, lapply(unique(s$country_id), average_weight, fsw="READ"))
avw_list <- do.call(rbind, lapply(unique(s$country_id), average_weight, fsw="LIST"))
avw_ques <- do.call(rbind, lapply(unique(s$country_id), average_weight, fsw="QUES"))
dim(s)
s <- merge(s, avw_read, by="country_id", all.x = T)
s <- merge(s, avw_list, by="country_id", all.x = T)
s <- merge(s, avw_ques, by="country_id", all.x = T)
dim(s)
s$countryavg_wread <- s$FSW_READ / s$avwREAD
s$countryavg_wlist <- s$FSW_LIST / s$avwLIST
s$countryavg_wques <- s$FSW_QUES / s$avwQUES
#lapply(list(s$countryavg_wread,s$countryavg_wlist,s$countryavg_wques), summary)
#lapply(list(s$countryavg_wread,s$countryavg_wlist,s$countryavg_wques), sd, na.rm=T)
#with(s[is.na(s$PV1_dreadlist)==F,],plot(countryavg_wread,countryavg_wlist, xlim=c(0,6.1), ylim=c(0,6.1)))

s$countryavg_w <- s$countryavg_wread
mean(s$countryavg_w[is.na(s$PV1_dreadlist)==F])
sd(s$countryavg_w[is.na(s$PV1_dreadlist)==F])
#plot(s$countryavg_wques[is.na(s$PV1_dreadlist)==F],s$countryavg_w[is.na(s$PV1_dreadlist)==F], xlim=c(0,6.1), ylim=c(0,6.1))

temp_weighted_diff <- s$PV1_dreadlist*s$countryavg_w
mean(s$PV1_dreadlist, na.rm=T)
mean(s$countryavg_w[is.na(s$PV1_dreadlist)==F])
mean(temp_weighted_diff, na.rm=T)



#############################################
# save smaller dataset for cognitive work
#############################################

sergio <- s[is.na(s$PV1_LIST)==F&is.na(s$PV1_READ)==F,c("respondent_id","I08_ST_A_S19B","I03_ST_A_S03A","I08_ST_A_S15A","I03_ST_A_S26B","I08_ST_A_S01A","I08_ST_A_S02A","cnt","I01_ST_M_S40A","I01_ST_M_S44A_cont","listoosQ31","I03_ST_A_S31A","classoralQ58","classwrittenQ58","travelQ45","classictQ51","I09_IN_M_S50A","goodclass2Q54","countryavg_w", "face2faceQ29", "txtbQ52")]
#&s$I03_ST_A_S04B==0
dim(sergio)
sergio <- sergio[complete.cases(sergio)==T,]
dim(sergio)
saveRDS(sergio, "sergio.rds")


######################################################
### start to try out REGRESSIONS
######################################################


#ses - home location - immigrant backgr - use of tl at home - gender - age - duration FL learning - FL lessons/week
#I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I01_ST_M_S39A + I01_ST_M_S44B
# simple linear models for experimentation
myest <- lm(dreadlist ~ I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A_cont + travelQ45 + face2faceQ29 + classoralQ58 + classwrittenQ58 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54, data = s)
summary(myest)
myest <- lm(dreadlist ~ I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + cnt:I01_ST_M_S40A + cnt:I01_ST_M_S44A_cont + travelQ45 + face2faceQ29 + classoralQ58 + classwrittenQ58 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54, data = s)
summary(myest)


#########################################################################################
### conditional model
#########################################################################################

myest0r <- lm(read ~ list + I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A + travelQ45 + face2faceQ29 + classoralQ58 + classwrittenQ58 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54, data = s)
summary(myest0r)
myest0l <- lm(list ~ read + I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A + travelQ45 + face2faceQ29 + classoralQ58 + classwrittenQ58 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54, data = s)
summary(myest0l)

run_conditional_models <- function(pv) {
  myr <- paste0("PV",pv,"_READ")
  myl <- paste0("PV",pv,"_LIST")
  myformula_r <- paste0(myr, " ~ ", myl, " + I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A + travelQ45 + face2faceQ29 + classoralQ58 + classwrittenQ58 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54")
  myformula_l <- paste0(myl, " ~ ", myr, " + I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A + travelQ45 + face2faceQ29 + classoralQ58 + classwrittenQ58 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54")
  myest_r <- lm(myformula_r, data = s, weights = countryavg_wread)
  myest_l <- lm(myformula_l, data = s, weights = countryavg_wread)
  output <- data.frame(b_r=myest_r$coefficients, sd_r=sqrt(diag(vcov(myest_r))), b_l=myest_l$coefficients, sd_l=sqrt(diag(vcov(myest_l))))
  colnames(output) <- paste0(c("b_r", "sd_r", "b_l", "sd_l"),pv)
  return(output)
}
cmod <- do.call(cbind,lapply(1:5, run_conditional_models))
str(cmod)

cmod$br <- .2 * (cmod$b_r1 + cmod$b_r2 + cmod$b_r3 + cmod$b_r4 + cmod$b_r5) 
cmod$bl <- .2 * (cmod$b_l1 + cmod$b_l2 + cmod$b_l3 + cmod$b_l4 + cmod$b_l5) 
cmod$sdr <- sqrt(.2 * (cmod$b_r1^2 + cmod$b_r2^2 + cmod$b_r3^2 + cmod$b_r4^2 + cmod$b_r5^2) - (.2 * (cmod$b_r1^2 + cmod$b_r2^2 + cmod$b_r3^2 + cmod$b_r4^2 + cmod$b_r5^2))^2   + (.2 * (cmod$sd_r1^2 + cmod$sd_r2^2 + cmod$sd_r3^2 + cmod$sd_r4^2 + cmod$sd_r5^2) ) )
cmod$sdl <- sqrt(.2 * (cmod$b_l1^2 + cmod$b_l2^2 + cmod$b_l3^2 + cmod$b_l4^2 + cmod$b_l5^2) - (.2 * (cmod$b_l1^2 + cmod$b_l2^2 + cmod$b_l3^2 + cmod$b_l4^2 + cmod$b_l5^2))^2   + (.2 * (cmod$sd_l1^2 + cmod$sd_l2^2 + cmod$sd_l3^2 + cmod$sd_l4^2 + cmod$sd_l5^2) ) )
cmod[,c("br","bl","sdr","sdl")]

#########################################################################################
### mega-loop for exploratory regressions
#########################################################################################


exploreg <- function (myx, myvar) {
  rm(myest)
  print(myx)
  s$temp <- (s[,myx] - mean(s[,myx], na.rm=T)) / sd(s[,myx], na.rm=T)
  try(myest <- lm(as.formula(paste0(myvar, " ~ temp + I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I01_ST_M_S39A + I01_ST_M_S44B + cnt")), data = s))
  if (exists("myest")==T) {coeffs <- as.data.frame(coef(summary(myest)))[2,c(1,3)]}
  if (exists("myest")==T) {coeffs$df <- myest$df.residual} else {coeffs <- as.data.frame(cbind(NA,NA,NA))}
  rownames(coeffs) <- colnames(s)[myx]
  colnames(coeffs) <- c("Estimate", "t value", "df")
  return(coeffs)
}

if (explore == T) {
  mycoefflist <- lapply(c(12:397), exploreg_cefr, myvar = "PV_cefrREADcefrLIST")
  myexplocoeffs_cefr <- do.call("rbind", mycoefflist)
  head(myexplocoeffs_cefr)
  write.csv(myexplocoeffs_cefr,"preliminary_myexplocoeffs_cefrplus.csv")
  
  mycoefflist <- lapply(c(12:397), exploreg, myvar = "PV_zREADzLIST")
  myexplocoeffs_z <- do.call("rbind", mycoefflist)
  head(myexplocoeffs_z)
  write.csv(myexplocoeffs_z,"preliminary_myexplocoeffs_zplus.csv")
  
  mycoefflist <- lapply(c(12:397), exploreg, myvar = "dPL_READLIST")
  myexplocoeffs_l <- do.call("rbind", mycoefflist)
  head(myexplocoeffs_l)
  write.csv(myexplocoeffs_l,"preliminary_myexplocoeffs_lplus.csv")
}




cor((s$PV1_LIST-s$list), (s$PV1_READ-s$read))
cor((s$PV3_LIST-s$list), (s$PV3_READ-s$read))


