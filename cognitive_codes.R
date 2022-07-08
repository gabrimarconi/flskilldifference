


# in the item bootstrapping, sample the items separately for each individual (instead of once for everyone)
bootindividually <- F
# conventional calculation of proficiency estimation variance (instead of bootstrapping)
noboot <- T
# plausible values for reading conditional on listening scores
secondscore <- F
# proportion of correct answers to calculate cognitive scores (instead of count)
propscore <- F
# run latent regression
twosteps <- T
# use probability to be at level B1 (instead of expected value of CEFR level)
proB1<- F
# when this coefficient is 0, a student is classified at level X if she has a 50% chance to answer correctly a question with the average difficulty for that level
levcoeff <- as.numeric(.0)
# combinations/correlations:
#TTF .936/938
#TFF .928/.927
#TFT .903/.904
#TTT .907/902
#FTF .949/928
#FFF .940/.900 -> fast combination to get something
#FFT .898/.886
#FTT .918/.905
# levcoeff = .7 -> 2/3 of right answers to qualify at a certain level

############################################### 
### generate file "k" with the info needed for loglinear rasch model
############################################### 

sergio <- readRDS("sergio.rds")
cogno_sco <- readRDS("cogno_sco.rds")
#View(cogno_sco[,660:665])
cogno_sco <- cogno_sco[,1:138]
k <- merge(cogno_sco,sergio)
a <- rowMeans(k[,9:66],na.rm=T)
k <- k[is.na(a)==F,]
#k$a <- rowMeans(k[,67:138],na.rm=T)
#k <- k[is.na(k$a)==F,]
k$bookl <- substr(k$booklet_L,1,3)
k$bookr <- substr(k$booklet_R,1,3)
k$sumL <- rowSums(k[,9:66], na.rm=T)
k$sumR <- rowSums(k[,67:138], na.rm=T)
k$totL <- apply(as.matrix(1:dim(k)[1]), 1, function(i) {length(na.omit(as.numeric(k[i,9:66])))})
k$totR <- apply(as.matrix(1:dim(k)[1]), 1, function(i) {length(na.omit(as.numeric(k[i,67:138])))})
if (propscore==T) {
  k$sumL <- k$sumL/k$totL
  k$sumR <- k$sumR/k$totR
}
cor(k$sumL, k$sumR, use="complete.obs")
cor(k$sumL[k$bookl=="EL3"], k$sumR[k$bookr=="ER3"], use="complete.obs")

dim(k)
colnames(k)
k$respondent_id <- as.character(k$respondent_id)
#View(k[,1:50])
#View(k[,(dim(k)[2]-49):dim(k)[2]])

###############################################
### bootstrap sampling
###############################################

if (bootindividually==F) {
  
  ### bootstrapping function: creates vector of question weights for bootstrap samples and sums of correct items per bootstrap sample
  
  bootstrappamelo <- function(mycols) {
    #mycols <- 9:66
    # generate response data matrix and question vector
    mymat <- k[,mycols]
    questions <- as.data.frame(as.character(mycols))
    colnames(questions) <- c("q")
    # generate bootstrap sample of questions from question vector
    questionw <- as.data.frame(table(sample(questions$q,length(questions$q),replace = T)))
    colnames(questionw) <- c("q", "w")
    questionw <- merge(questions,questionw, all.x = T)
    questionw <- questionw[order(as.numeric(questionw$q)),]
    # calculate total number of items answered in bootstrap sample
    if (propscore==T) {
      totItems <- apply(   as.matrix(1:dim(mymat)[1]), 1, function(i) {sum(na.omit(as.numeric(as.numeric(questionw$w) + mymat[i,] - mymat[i,])))}   )
    }
    # calculate sum of correct responses and divide by total items
    questionw$w[is.na(questionw$w)==T] <- 0
    mymat[is.na(mymat)] <- 0
    bootsum <- as.matrix(mymat)%*%as.matrix(as.numeric(questionw$w))
    if (propscore==T) {bootsum <- bootsum / totItems}
    # compile and return output
    myoutput <- list(questionw, bootsum)
    return(myoutput)
  }
  
  ### listening
  
  boot1L <- bootstrappamelo(9:66)
  boot2L <- bootstrappamelo(9:66)
  boot3L <- bootstrappamelo(9:66)
  boot4L <- bootstrappamelo(9:66)
  boot5L <- bootstrappamelo(9:66)
  w1L <- boot1L[[1]]
  w2L <- boot2L[[1]]
  w3L <- boot3L[[1]]
  w4L <- boot4L[[1]]
  w5L <- boot5L[[1]]
  colnames(w1L) <- c("q","w1L")
  colnames(w2L) <- c("q","w2L")
  colnames(w3L) <- c("q","w3L")
  colnames(w4L) <- c("q","w4L")
  colnames(w5L) <- c("q","w5L")
  k$sumL1 <- boot1L[[2]]
  k$sumL2 <- boot2L[[2]]
  k$sumL3 <- boot3L[[2]]
  k$sumL4 <- boot4L[[2]]
  k$sumL5 <- boot5L[[2]]
  summary(k$sumL3)
  cor(k$sumL1, k$sumL2, use="complete.obs")
  
  ### reading
  
  boot1R <- bootstrappamelo(67:138)
  boot2R <- bootstrappamelo(67:138)
  boot3R <- bootstrappamelo(67:138)
  boot4R <- bootstrappamelo(67:138)
  boot5R <- bootstrappamelo(67:138)
  w1R <- boot1R[[1]]
  w2R <- boot2R[[1]]
  w3R <- boot3R[[1]]
  w4R <- boot4R[[1]]
  w5R <- boot5R[[1]]
  colnames(w1R) <- c("q","w1R")
  colnames(w2R) <- c("q","w2R")
  colnames(w3R) <- c("q","w3R")
  colnames(w4R) <- c("q","w4R")
  colnames(w5R) <- c("q","w5R")
  k$sumR1 <- boot1R[[2]]
  k$sumR2 <- boot2R[[2]]
  k$sumR3 <- boot3R[[2]]
  k$sumR4 <- boot4R[[2]]
  k$sumR5 <- boot5R[[2]]
  
}


if (bootindividually==T) {
  
  ### bootstrapping 2: over the items of each single respondent
  
  # function to calculate bootstrap sample for single respondent, with relative weights and sum of correct answers
  rowweights <- function(myrow, mycols) {
    #myrow <- 9
    #mycols <- 9:66
    # generate response data matrix and question vector
    print(myrow)
    mymat <- k[myrow,mycols]
    questionlist <- as.data.frame(as.character(mycols))
    colnames(questionlist) <- "q"
    questions <- as.data.frame(   na.omit(cbind(questionlist$q, t(mymat)))   )
    colnames(questions) <- c("q","y")
    # generate bootstrap sample of questions from question vector
    nrq <- length(questions$q)
    questionw <- as.data.frame(table(sample(questions$q,nrq,replace = T)))
    colnames(questionw) <- c("q", "w")
    # merge question vector with weights and calculate sample sum of correct answers
    questionw <- merge(questions,questionw, all.x = T)
    sumLboot <- sum(as.numeric(questionw$y)*as.numeric(questionw$w), na.rm = T)
    if (propscore==T) {sumLboot <- sumLboot / nrq}
    questionlist <- merge(questionlist, questionw, all.x = T)
    questionlist <- questionlist[order(as.numeric(questionlist$q)),]
    #compile and return output
    myoutput <- t(as.matrix(c(questionlist$w,sumLboot)))
    return(myoutput)
  }
  # function to extend previous function to all respondents, get weights for each question, and get sum of correct answers for each respondent
  calcwmat <- function (mycolumns) {
    #mycolumns <- 9:66
    wmatlist <- lapply(1:dim(k)[1], rowweights, mycols=mycolumns)
    wmat <- as.data.frame(do.call("rbind", wmatlist))
    colnames(wmat) <- c(as.character(mycolumns), "sumLboot")
    wmat$respondent_id <- k$respondent_id
    logwmat <- lapply(1:(dim(wmat)[2]-1), function(i) {as.matrix(cbind(wmat[,i],wmat$sumLboot,wmat$respondent_id,i+(mycolumns[1]-1)))})
    logwmat <- as.data.frame(do.call("rbind", logwmat))
    colnames(logwmat) <- c("w", "sumLboot", "respondent_id", "q")
    logwmat$w[is.na(logwmat$w)==T] <- 0
    return(logwmat)
  }
  
  ### listening
  listwmat1 <- calcwmat(9:66)
  listwmat2 <- calcwmat(9:66)
  listwmat3 <- calcwmat(9:66)
  listwmat4 <- calcwmat(9:66)
  listwmat5 <- calcwmat(9:66)
  #View(listwmat1)
  # add provisionary empty sum variables
  k$sumL1 <- NA
  k$sumL2 <- NA
  k$sumL3 <- NA
  k$sumL4 <- NA
  k$sumL5 <- NA
  ### read
  readwmat1 <- calcwmat(67:138)
  readwmat2 <- calcwmat(67:138)
  readwmat3 <- calcwmat(67:138)
  readwmat4 <- calcwmat(67:138)
  readwmat5 <- calcwmat(67:138)
  #View(readwmat1)
  # add provisionary empty sum variables
  k$sumR1 <- NA
  k$sumR2 <- NA
  k$sumR3 <- NA
  k$sumR4 <- NA
  k$sumR5 <- NA
  
}




### find out which listening items belong to which booklet

whichlistlets <- apply(as.matrix(9:66),1,function (i)  {ifelse(table(k$bookl,ifelse(is.na(k[,i])==F,1,NA))>0,1,0)}   )
colnames(whichlistlets) <- paste0("c",9:66)
# View(whichbooklets)
# results: 
 # booklet 1 -> 9:26, 37:42 (9 own items + 15 link items)
 # booklet 2 -> 14:18, 23:31, 37:48, 55:58 (30 link items)
 # booklet 3 -> 27:36, 43:58 (10 own items + 15 link items)

# count response pattern number -> result: they're so many that the aggregated loglinear model does not save much matrix space compared to estimating but by but
#prova <- apply(as.matrix(1:dim(k)[1]),1,function (i) {paste0(k[i,9:66], collapse = "_")})
#length(unique(prova))




#############################################
### listening all booklets
#############################################

### compile dataset for listening rasch - booklet 1
logdf <- lapply(c(9:66), function(i) {as.matrix(cbind(k[,i],k$respondent_id,i,k$booklet_L,k$booklet_R,k[139:174]))})
logdf <- as.data.frame(do.call("rbind", logdf))
colnames(logdf) <- c("y", "respondent_id", "q", "booklet_L", "booklet_R", colnames(logdf[(length(colnames(logdf))-35):length(colnames(logdf))]))
logdf <- logdf[is.na(logdf$y)==F,]

### set appropriate variable type
logdf$q <- as.character(logdf$q)
logdf$y <- as.numeric(logdf$y)
for (i in c(6:11,13:23)) {
  print(colnames(logdf)[i])
  logdf[,i] <- as.numeric(logdf[,i])
  # binarize continuous indexes (if needed to save computing power)
  #logdf$temp <- ifelse(logdf[,i] > quantile(logdf[,i], probs = .5, na.rm = T), 1, 0)
  #colnames(logdf) <- c(   colnames(logdf)[1:(dim(logdf)[2]-1)], paste0(colnames(logdf)[i],"_b")   )
}


### add bootstrap weights and (in the case of individual bootstrap) sums

if (bootindividually==F) {
  ### merge with bootstrap weights
  logdf <- merge(logdf, w1L, all.x = T)
  logdf <- merge(logdf, w2L, all.x = T)
  logdf <- merge(logdf, w3L, all.x = T)
  logdf <- merge(logdf, w4L, all.x = T)
  logdf <- merge(logdf, w5L, all.x = T)
}

if (bootindividually==T) {
  ### merge with bootstrap weights
  logdf$w1L <- NA
  logdf$w2L <- NA
  logdf$w3L <- NA
  logdf$w4L <- NA
  logdf$w5L <- NA
  # list1
  logdf <- merge(logdf, listwmat1, all.x = T)
  logdf$w1L <- as.numeric(logdf$w)
  logdf$sumL1 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list2
  logdf <- merge(logdf, listwmat2, all.x = T)
  logdf$w2L <- as.numeric(logdf$w)
  logdf$sumL2 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list3
  logdf <- merge(logdf, listwmat3, all.x = T)
  logdf$w3L <- as.numeric(logdf$w)
  logdf$sumL3 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list4
  logdf <- merge(logdf, listwmat4, all.x = T)
  logdf$w4L <- as.numeric(logdf$w)
  logdf$sumL4 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list5
  logdf <- merge(logdf, listwmat5, all.x = T)
  logdf$w5L <- as.numeric(logdf$w)
  logdf$sumL5 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read1
  logdf <- merge(logdf, readwmat1[duplicated(readwmat1$respondent_id)==F,1:3], all.x = T)
  logdf$sumR1 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read2
  logdf <- merge(logdf, readwmat2[duplicated(readwmat2$respondent_id)==F,1:3], all.x = T)
  logdf$sumR2 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read3
  logdf <- merge(logdf, readwmat3[duplicated(readwmat3$respondent_id)==F,1:3], all.x = T)
  logdf$sumR3 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read4
  logdf <- merge(logdf, readwmat4[duplicated(readwmat4$respondent_id)==F,1:3], all.x = T)
  logdf$sumR4 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read5
  logdf <- merge(logdf, readwmat5[duplicated(readwmat5$respondent_id)==F,1:3], all.x = T)
  logdf$sumR5 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
}


# generate weight for base estimates
logdf$w0L <- 1

### activate your favourite sumL:
logdf$sumL <- as.numeric(logdf$sumL)
if (propscore==F) {
  logdf$sumL0 <- as.numeric( as.numeric(logdf$sumL) - logdf$y )
  logdf$sumL1 <- as.numeric( as.numeric(logdf$sumL1) - logdf$y )
  logdf$sumL2 <- as.numeric( as.numeric(logdf$sumL2) - logdf$y )
  logdf$sumL3 <- as.numeric( as.numeric(logdf$sumL3) - logdf$y )
  logdf$sumL4 <- as.numeric( as.numeric(logdf$sumL4) - logdf$y )
  logdf$sumL5 <- as.numeric( as.numeric(logdf$sumL5) - logdf$y )
}
if (propscore==T) {
  logdf$sumL0 <- as.numeric( as.numeric(logdf$sumL) )
  logdf$sumL1 <- as.numeric( as.numeric(logdf$sumL1) )
  logdf$sumL2 <- as.numeric( as.numeric(logdf$sumL2) )
  logdf$sumL3 <- as.numeric( as.numeric(logdf$sumL3) )
  logdf$sumL4 <- as.numeric( as.numeric(logdf$sumL4) )
  logdf$sumL5 <- as.numeric( as.numeric(logdf$sumL5) )
}
logdf$sumR <- as.numeric(logdf$sumR)
logdf$sumR1 <- as.numeric(logdf$sumR1)
logdf$sumR2 <- as.numeric(logdf$sumR2)
logdf$sumR3 <- as.numeric(logdf$sumR3)
logdf$sumR4 <- as.numeric(logdf$sumR4)
logdf$sumR5 <- as.numeric(logdf$sumR5)

colnames(logdf)


### generate PVs

generatepv <- function(myscore, mysecondscore=NULL, myw, savediffic=F) {
  #myw <- logdf$w0L
  #myscore <- logdf$sumL
  #mysecondscore <- logdf$sumR1
  #savediffic <- T
  logdf$bootw <- myw
  logdf$score <- myscore
  logdf$scoresq <- logdf$score^2  
  if (secondscore==T) {logdf$secondscore <- mysecondscore}
  
  # rasch logit model
  print("input variables for rasch logit")
  if (twosteps==F & secondscore==F) {
    mylogit <- glm(y ~ score:bookl + scoresq + q +  booklet_L + cnt +  I08_ST_A_S19B + I03_ST_A_S03A +  I08_ST_A_S15A +  I03_ST_A_S26B +  I08_ST_A_S01A +  I08_ST_A_S02A + I01_ST_M_S40A +  I01_ST_M_S44A_cont +  interestflQ34 +  interestvsQ34 +  classoralQ58 +  classwrittenQ58 + travelQ45 +  classictQ51 +  I09_IN_M_S50A +  goodclass2Q54 +  diffQ48, data = logdf, family = binomial, weights = bootw)
    #mean(predict(mylogit,logdf,type="response"))
    }
  if (twosteps==F & secondscore==T) {
    mylogit <- glm(y ~ score:bookl + scoresq + q +  booklet_L + cnt +  I08_ST_A_S19B + I03_ST_A_S03A +  I08_ST_A_S15A +  I03_ST_A_S26B +  I08_ST_A_S01A +  I08_ST_A_S02A + I01_ST_M_S40A +  I01_ST_M_S44A_cont +  interestflQ34 +  interestvsQ34 +  classoralQ58 +  classwrittenQ58 + travelQ45 +  classictQ51 +  I09_IN_M_S50A +  goodclass2Q54 +  diffQ48 + secondscore, data = logdf, family = binomial, weights = bootw)
  }
  if (twosteps==T & secondscore==T) {
    mylogit <- glm(y ~ score:bookl + q +  booklet_L, data = logdf, family = binomial, weights = bootw)
    logdf$list_step1 <- ifelse(logdf$bookl=="EL1", myscore*mylogit$coefficients[65], 0) + ifelse(logdf$bookl=="EL2", myscore*mylogit$coefficients[66], 0) + ifelse(logdf$bookl=="EL3", myscore*mylogit$coefficients[67], 0) +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99))]   )   # fitting the booklet
    myols <- lm(list_step1 ~ I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A_cont + classoralQ58 + classwrittenQ58 + travelQ45 + face2faceQ29 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54 , data = logdf, weights = bootw)
    logdf$list <- as.numeric(myols$fitted.values)
    logdf$list_err <- as.numeric(myols$residuals)
  }
  if (twosteps==T & secondscore==F) {
    mylogit <- glm(y ~ score:bookl + q +  booklet_L, data = logdf, family = binomial, weights = bootw)
    logdf$list_step1 <- ifelse(logdf$bookl=="EL1", myscore*mylogit$coefficients[65], 0) + ifelse(logdf$bookl=="EL2", myscore*mylogit$coefficients[66], 0) + ifelse(logdf$bookl=="EL3", myscore*mylogit$coefficients[67], 0) +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99))]   )   # fitting the booklet
    myols <- lm(list_step1 ~ I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + cnt:I01_ST_M_S40A + cnt:I01_ST_M_S44A_cont + classoralQ58 + classwrittenQ58 + travelQ45 + face2faceQ29 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54 , data = logdf, weights = bootw)
    logdf$list <- as.numeric(myols$fitted.values)
    logdf$list_err <- as.numeric(myols$residuals)
  }
  print("logit estimation finished")
  summary(mylogit)
  length(t(as.matrix(mylogit$coefficients)))
  
  # fit ability 
  # (only needed if twosteps==F; for twosteps==T, ability is fitted in the previous part of the process)
  if (twosteps==F & secondscore==F) {
    logdf$list <- ifelse(logdf$bookl=="EL1", myscore*mylogit$coefficients[97], 0) + ifelse(logdf$bookl=="EL2", myscore*mylogit$coefficients[98], 0) + ifelse(logdf$bookl=="EL3", myscore*mylogit$coefficients[99], 0) + myscore^2*mylogit$coefficients[2] +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99))]   ) +   # fitting the booklet
      ifelse(   is.na(match(logdf$cnt,substr(names(mylogit$coefficients),4,99)))==T,0,   mylogit$coefficients[match(logdf$cnt,substr(names(mylogit$coefficients),4,99))]   ) +   # fitting the country
      as.matrix(logdf[,c(6:11,13:23)]) %*% t(as.matrix(mylogit$coefficients))[80:96]  # fitting the continuous indexes
  }
  if (twosteps==F & secondscore==T) {
    logdf$list <- ifelse(logdf$bookl=="EL1", myscore*mylogit$coefficients[98], 0) + ifelse(logdf$bookl=="EL2", myscore*mylogit$coefficients[99], 0) + ifelse(logdf$bookl=="EL3", myscore*mylogit$coefficients[100], 0) + myscore^2*mylogit$coefficients[2] + mysecondscore*mylogit$coefficients[97] +  # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99))]   ) +   # fitting the booklet
      ifelse(   is.na(match(logdf$cnt,substr(names(mylogit$coefficients),4,99)))==T,0,   mylogit$coefficients[match(logdf$cnt,substr(names(mylogit$coefficients),4,99))]   ) +   # fitting the country
      as.matrix(logdf[,c(6:11,13:23)]) %*% t(as.matrix(mylogit$coefficients))[80:96]  # fitting the continuous indexes
  }
  print(summary(logdf$list))
  
  # get difficulties
  # levels: A1 (111/213) -> 9:13, 19:22; A2 (123/221/321) -> 14:18, 23:26, 37:42; B1 (231/432/531) -> 27:31, 43:48, 55:60; B2 (242/442/543) -> 32:36, 49:54, 61:66
  diffic <- c(1:8, mylogit$coefficients[58] + as.numeric(mylogit$coefficients[1]),mylogit$coefficients[1] , mylogit$coefficients[2:57] + as.numeric(mylogit$coefficients[1]))
  A1 <- mean(diffic[c(9:13, 19:22)], na.rm=T)
  A2 <- mean(diffic[c(14:18, 23:26, 37:42)], na.rm=T)
  B1 <- mean(diffic[c(27:31, 43:48, 55:60)], na.rm=T)
  B2 <- mean(diffic[c(32:36, 49:54, 61:66)], na.rm=T)
  C1 <- min(diffic[c(32:36, 49:54, 61:66)], na.rm=T)
  print(paste0(c("A1: ", "A2: ", "B1: ", "B2: ", "C1: "), c(A1, A2, B1, B2, C1)))
  if (savediffic==T) {write.csv( (-as.numeric(c(A1,A2,B1,B2,C1))+levcoeff), "difficlist.csv" )}
  logdf$listlev[logdf$list>as.numeric(-C1)+levcoeff] <- 5
  logdf$listlev[logdf$list>as.numeric(-B2)+levcoeff] <- 4
  logdf$listlev[logdf$list<as.numeric(-B2)+levcoeff] <- 3
  logdf$listlev[logdf$list<as.numeric(-B1)+levcoeff] <- 2
  logdf$listlev[logdf$list<as.numeric(-A2)+levcoeff] <- 1
  logdf$listlev[logdf$list<as.numeric(-A1)+levcoeff] <- 0
  
  # generate proficiency score
  logdf$listpsum <- exp(logdf$list+as.numeric(A1))/(1 + exp(logdf$list+as.numeric(A1))) + exp(logdf$list+as.numeric(A2))/(1 + exp(logdf$list+as.numeric(A2))) + exp(logdf$list+as.numeric(B1))/(1 + exp(logdf$list+as.numeric(B1))) + exp(logdf$list+as.numeric(B2))/(1 + exp(logdf$list+as.numeric(B2))) + exp(logdf$list+as.numeric(C1))/(1 + exp(logdf$list+as.numeric(C1)))
  if(proB1==T) {logdf$listpsum <- exp(logdf$list+as.numeric(B1))/(1 + exp(logdf$list+as.numeric(B1)))}
  print(summary(logdf$listlev))
  print(summary(logdf$listpsum))
  
  # compile and return output
  if (noboot==F) {
    measure_list <- logdf[duplicated(logdf$respondent_id)==F,c("respondent_id","list","listlev","listpsum")]
  } else {
    measure_list <- logdf[duplicated(logdf$respondent_id)==F,c("respondent_id","list","listlev","listpsum", "list_err")]
  }
  
  print("all calculations finished")
  return(measure_list)
}
list000 <- generatepv(logdf$sumL0, logdf$sumR, logdf$w0L, savediffic = T)
if (noboot == T) {
  colnames(list000) <- c("respondent_id", "list000", "listlev000", "listpsum0", "list_err")
}
if (noboot == F) {
  colnames(list000) <- c("respondent_id", "list000", "listlev000", "listpsum0")
  listpv1 <- generatepv(logdf$sumL1, logdf$sumR1, logdf$w1L)
  listpv2 <- generatepv(logdf$sumL2, logdf$sumR2, logdf$w2L)
  listpv3 <- generatepv(logdf$sumL3, logdf$sumR3, logdf$w3L)
  listpv4 <- generatepv(logdf$sumL4, logdf$sumR4, logdf$w4L)
  listpv5 <- generatepv(logdf$sumL5, logdf$sumR5, logdf$w5L)
  colnames(listpv1) <- c("respondent_id", "listpv1", "listlevpv1", "listpsum1")
  colnames(listpv2) <- c("respondent_id", "listpv2", "listlevpv2", "listpsum2")
  colnames(listpv3) <- c("respondent_id", "listpv3", "listlevpv3", "listpsum3")
  colnames(listpv4) <- c("respondent_id", "listpv4", "listlevpv4", "listpsum4")
  colnames(listpv5) <- c("respondent_id", "listpv5", "listlevpv5", "listpsum5")
}


# smaller logit
#mylogit <- glm(y ~ sumLn +  booklet_L + q, data = logdf, family = binomial)
#summary(mylogit)
#length(t(as.matrix(mylogit$coefficients)))
#
#logdf$list_smaller <- logdf$sumL*mylogit$coefficients[2] +   # fitting the number of correct answers
#  ifelse(   is.na(match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_L,substr(names(mylogit$coefficients),10,99))]   )   # fitting the booklet
#  summary(logdf$list_smaller)
#
#measure_list_smaller <- logdf[duplicated(logdf$respondent_id)==F,c("respondent_id","list_smaller")]


#############################################
### reading all booklets
#############################################

# compile dataset for reading rasch
logdf <- lapply(c(67:138), function(i) {as.matrix(cbind(k[,i],k$respondent_id,i,k$booklet_L,k$booklet_R,k[139:174]))})
logdf <- as.data.frame(do.call("rbind", logdf))
colnames(logdf) <- c("y", "respondent_id", "q", "booklet_L", "booklet_R", colnames(logdf[(length(colnames(logdf))-35):length(colnames(logdf))]))
logdf <- logdf[is.na(logdf$y)==F,]

# set appropriate variable type
logdf$q <- as.character(logdf$q)
logdf$y <- as.numeric(logdf$y)
for (i in c(6:11,13:23)) {
  print(colnames(logdf)[i])
  logdf[,i] <- as.numeric(logdf[,i])
  # binarize continuous indexes (if needed to save computing power)
  #logdf$temp <- ifelse(logdf[,i] > quantile(logdf[,i], probs = .5, na.rm = T), 1, 0)
  #colnames(logdf) <- c(   colnames(logdf)[1:(dim(logdf)[2]-1)], paste0(colnames(logdf)[i],"_b")   )
}


### add bootstrap weights and (in the case of individual bootstrap) sums

if (bootindividually==F) {
  ### merge with bootstrap weights
  logdf <- merge(logdf, w1R, all.x = T)
  logdf <- merge(logdf, w2R, all.x = T)
  logdf <- merge(logdf, w3R, all.x = T)
  logdf <- merge(logdf, w4R, all.x = T)
  logdf <- merge(logdf, w5R, all.x = T)
}

if (bootindividually==T) {
  ### merge with bootstrap weights
  logdf$w1R <- NA
  logdf$w2R <- NA
  logdf$w3R <- NA
  logdf$w4R <- NA
  logdf$w5R <- NA
  # read1
  logdf <- merge(logdf, readwmat1, all.x = T)
  logdf$w1R <- as.numeric(logdf$w)
  logdf$sumR1 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read2
  logdf <- merge(logdf, readwmat2, all.x = T)
  logdf$w2R <- as.numeric(logdf$w)
  logdf$sumR2 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read3
  logdf <- merge(logdf, readwmat3, all.x = T)
  logdf$w3R <- as.numeric(logdf$w)
  logdf$sumR3 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read4
  logdf <- merge(logdf, readwmat4, all.x = T)
  logdf$w4R <- as.numeric(logdf$w)
  logdf$sumR4 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # read5
  logdf <- merge(logdf, readwmat5, all.x = T)
  logdf$w5R <- as.numeric(logdf$w)
  logdf$sumR5 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list1
  logdf <- merge(logdf, listwmat1[duplicated(listwmat1$respondent_id)==F,1:3], all.x = T)
  logdf$sumL1 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list2
  logdf <- merge(logdf, listwmat2[duplicated(listwmat2$respondent_id)==F,1:3], all.x = T)
  logdf$sumL2 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list3
  logdf <- merge(logdf, listwmat3[duplicated(listwmat3$respondent_id)==F,1:3], all.x = T)
  logdf$sumL3 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list4
  logdf <- merge(logdf, listwmat4[duplicated(listwmat4$respondent_id)==F,1:3], all.x = T)
  logdf$sumL4 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
  # list5
  logdf <- merge(logdf, listwmat5[duplicated(listwmat5$respondent_id)==F,1:3], all.x = T)
  logdf$sumL5 <- as.numeric(logdf$sumLboot)
  logdf <- logdf[,1:(dim(logdf)[2]-2)]
}


# generate weight for base estimates
logdf$w0R <- 1


### activate your favourite sumR:
logdf$sumR <- as.numeric(logdf$sumR)
if (propscore==F) {
  logdf$sumR0 <- as.numeric( as.numeric(logdf$sumR) - logdf$y )
  logdf$sumR1 <- as.numeric( as.numeric(logdf$sumR1) - logdf$y )
  logdf$sumR2 <- as.numeric( as.numeric(logdf$sumR2) - logdf$y )
  logdf$sumR3 <- as.numeric( as.numeric(logdf$sumR3) - logdf$y )
  logdf$sumR4 <- as.numeric( as.numeric(logdf$sumR4) - logdf$y )
  logdf$sumR5 <- as.numeric( as.numeric(logdf$sumR5) - logdf$y )
}

if (propscore==T) {
  logdf$sumR0 <- as.numeric( as.numeric(logdf$sumR) )
  logdf$sumR1 <- as.numeric( as.numeric(logdf$sumR1) )
  logdf$sumR2 <- as.numeric( as.numeric(logdf$sumR2) )
  logdf$sumR3 <- as.numeric( as.numeric(logdf$sumR3) )
  logdf$sumR4 <- as.numeric( as.numeric(logdf$sumR4) )
  logdf$sumR5 <- as.numeric( as.numeric(logdf$sumR5) )
}
logdf$sumL <- as.numeric(logdf$sumL)
logdf$sumL1 <- as.numeric(logdf$sumL1)
logdf$sumL2 <- as.numeric(logdf$sumL2)
logdf$sumL3 <- as.numeric(logdf$sumL3)
logdf$sumL4 <- as.numeric(logdf$sumL4)
logdf$sumL5 <- as.numeric(logdf$sumL5)


logdf <- merge(logdf, listpv1, all.x = T)
logdf <- merge(logdf, listpv2, all.x = T)
logdf <- merge(logdf, listpv3, all.x = T)
logdf <- merge(logdf, listpv4, all.x = T)
logdf <- merge(logdf, listpv5, all.x = T)
#View(logdf)
dim(listpv1)
#length(unique(listpv1$respondent_id))

colnames(logdf)

### generate PVs for reading

generatepv <- function(myscore, mysecondscore=NULL, myw, savediffic=F) {
  #myw <- logdf$w0R
  #myscore <- logdf$sumR0
  #mysecondscore <- logdf$sumL1
  #savediffic <- T
  logdf$bootw <- myw
  logdf$score <- myscore
  logdf$scoresq <- logdf$score^2  
  logdf$secondscore <- mysecondscore
  
  # rasch logit model
  print("input variables for rasch logit")
  if (twosteps==F & secondscore==F) {
    mylogit <- glm(y ~ score:bookr + scoresq + q +  booklet_R + cnt +  I08_ST_A_S19B + I03_ST_A_S03A +  I08_ST_A_S15A +  I03_ST_A_S26B +  I08_ST_A_S01A +  I08_ST_A_S02A + I01_ST_M_S40A +  I01_ST_M_S44A_cont +  interestflQ34 +  interestvsQ34 +  classoralQ58 +  classwrittenQ58 + travelQ45 +  classictQ51 +  I09_IN_M_S50A +  goodclass2Q54 +  diffQ48, data = logdf, family = binomial, weights = myw)
  }
  if (twosteps==F & secondscore==T) {
    mylogit <- glm(y ~ score:bookr + scoresq + q +  booklet_R + cnt +  I08_ST_A_S19B + I03_ST_A_S03A +  I08_ST_A_S15A +  I03_ST_A_S26B +  I08_ST_A_S01A +  I08_ST_A_S02A + I01_ST_M_S40A +  I01_ST_M_S44A_cont +  interestflQ34 +  interestvsQ34 +  classoralQ58 +  classwrittenQ58 + travelQ45 +  classictQ51 +  I09_IN_M_S50A +  goodclass2Q54 +  diffQ48 + secondscore, data = logdf, family = binomial, weights = myw)
  }
  if (twosteps==T & secondscore==T) {
    mylogit <- glm(y ~ score:bookr + q +  booklet_R, data = logdf, family = binomial, weights = bootw)
    logdf$read_step1 <- ifelse(logdf$bookr=="ER1", myscore*mylogit$coefficients[90], 0) + ifelse(logdf$bookr=="ER2", myscore*mylogit$coefficients[91], 0) + ifelse(logdf$bookr=="ER3", myscore*mylogit$coefficients[92], 0) +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99))]   )   # fitting the booklet
    myols <- lm(read_step1 ~ I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + I01_ST_M_S40A + I01_ST_M_S44A_cont + classoralQ58 + classwrittenQ58 + travelQ45 + face2faceQ29 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54 + secondscore, data = logdf, weights = bootw)
    logdf$read <- as.numeric(myols$fitted.values)
    logdf$read_err <- as.numeric(myols$residuals)
  }
  if (twosteps==T & secondscore==F) {
    mylogit <- glm(y ~ score:bookr + q +  booklet_R, data = logdf, family = binomial, weights = bootw)
    logdf$read_step1 <- ifelse(logdf$bookr=="ER1", myscore*mylogit$coefficients[90], 0) + ifelse(logdf$bookr=="ER2", myscore*mylogit$coefficients[91], 0) + ifelse(logdf$bookr=="ER3", myscore*mylogit$coefficients[92], 0) +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99))]   )   # fitting the booklet
    myols <- lm(read_step1 ~ I08_ST_A_S19B + I03_ST_A_S03A + I08_ST_A_S15A + I03_ST_A_S26B + I08_ST_A_S01A + I08_ST_A_S02A + I03_ST_A_S31A + cnt + cnt:I01_ST_M_S40A + cnt:I01_ST_M_S44A_cont + classoralQ58 + classwrittenQ58 + travelQ45 + face2faceQ29 + classictQ51 + txtbQ52 + I09_IN_M_S50A + goodclass2Q54, data = logdf, weights = bootw)
    logdf$read <- as.numeric(myols$fitted.values)
    logdf$read_err <- as.numeric(myols$residuals)
  }
  print("logit estimation finished")
  summary(mylogit)
  length(t(as.matrix(mylogit$coefficients)))

  # fit ability (only needed if twosteps==F; for twosteps==T, ability is fitted in the previous part of the process)
  if (twosteps==F & secondscore==F) {
    logdf$read <- ifelse(logdf$bookr=="ER1", myscore*mylogit$coefficients[122], 0) + ifelse(logdf$bookr=="ER2", myscore*mylogit$coefficients[123], 0) + ifelse(logdf$bookr=="ER3", myscore*mylogit$coefficients[124], 0) + myscore^2*mylogit$coefficients[2] +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99))]   ) +   # fitting the booklet
      ifelse(   is.na(match(logdf$cnt,substr(names(mylogit$coefficients),4,99)))==T,0,   mylogit$coefficients[match(logdf$cnt,substr(names(mylogit$coefficients),4,99))]   ) +   # fitting the country
      as.matrix(logdf[,c(6:11,13:23)]) %*% t(as.matrix(mylogit$coefficients))[105:121]  # fitting the continuous indexes
  }
  if (twosteps==F & secondscore==T) {
    logdf$read <- ifelse(logdf$bookr=="ER1", myscore*mylogit$coefficients[123], 0) + ifelse(logdf$bookr=="ER2", myscore*mylogit$coefficients[124], 0) + ifelse(logdf$bookr=="ER3", myscore*mylogit$coefficients[125], 0) + myscore^2*mylogit$coefficients[2] + mysecondscore*mylogit$coefficients[122] +   # fitting the number of correct answers
      ifelse(   is.na(match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99))]   ) +   # fitting the booklet
      ifelse(   is.na(match(logdf$cnt,substr(names(mylogit$coefficients),4,99)))==T,0,   mylogit$coefficients[match(logdf$cnt,substr(names(mylogit$coefficients),4,99))]   ) +   # fitting the country
      as.matrix(logdf[,c(6:11,13:23)]) %*% t(as.matrix(mylogit$coefficients))[105:121]  # fitting the continuous indexes
  }
  print(summary(logdf$read))
  
  # get difficulties
  # levels: A1 (112/211/312) -> 67:74, 80:84; A2 (223/321/423/523) -> 75:79, 85:99; B1 (532/631/731) -> 100:110, 123:126; B2 (642/741/841) -> 111:122, 127:138
  diffic <- c(1:66, mylogit$coefficients[41:73] + as.numeric(mylogit$coefficients[1]), mylogit$coefficients[1] , mylogit$coefficients[3:40] + as.numeric(mylogit$coefficients[1]))
  A1 <- mean(diffic[c(67:74, 80:84)], na.rm=T)
  A2 <- mean(diffic[c(75:79, 85:99)], na.rm=T)
  B1 <- mean(diffic[c(100:110, 123:126)], na.rm=T)
  B2 <- mean(diffic[c(111:122, 127:138)], na.rm=T)
  C1 <- min(diffic[c(111:122, 127:138)], na.rm=T)
  print(paste0(c("A1: ", "A2: ", "B1: ", "B2: ", "C1: "), c(A1, A2, B1, B2, C1)))
  if (savediffic==T) {write.csv( (-as.numeric(c(A1,A2,B1,B2,C1))+levcoeff), "difficread.csv" )}
  logdf$readlev[logdf$read>as.numeric(-C1)+levcoeff] <- 5
  logdf$readlev[logdf$read>as.numeric(-B2)+levcoeff] <- 4
  logdf$readlev[logdf$read<as.numeric(-B2)+levcoeff] <- 3
  logdf$readlev[logdf$read<as.numeric(-B1)+levcoeff] <- 2
  logdf$readlev[logdf$read<as.numeric(-A2)+levcoeff] <- 1
  logdf$readlev[logdf$read<as.numeric(-A1)+levcoeff] <- 0
  
  # calculate proficiency scores
  logdf$readpsum <- exp(logdf$read+as.numeric(A1))/(1 + exp(logdf$read+as.numeric(A1))) + exp(logdf$read+as.numeric(A2))/(1 + exp(logdf$read+as.numeric(A2))) + exp(logdf$read+as.numeric(B1))/(1 + exp(logdf$read+as.numeric(B1))) + exp(logdf$read+as.numeric(B2))/(1 + exp(logdf$read+as.numeric(B2))) + exp(logdf$read+as.numeric(C1))/(1 + exp(logdf$read+as.numeric(C1)))
  if(proB1==T) {logdf$readpsum <- exp(logdf$read+as.numeric(B1))/(1 + exp(logdf$read+as.numeric(B1)))}
  print(summary(logdf$readlev))
  print(summary(logdf$readpsum))

  # compile and return output
  if (noboot==F) {
    measure_read <- logdf[duplicated(logdf$respondent_id)==F,c("respondent_id","read", "readlev", "readpsum")]
  } else {
    measure_read <- logdf[duplicated(logdf$respondent_id)==F,c("respondent_id","read", "readlev", "readpsum","read_err")]
  }
  
  print("all calculations finished")
  return(measure_read)
}
read000 <- generatepv(logdf$sumR0, logdf$sumL, logdf$w0R, savediffic = T)
if (noboot == T) {
  colnames(read000) <- c("respondent_id", "read000", "readlev000", "readpsum0", "read_err")
}
if (noboot == F) {
  colnames(read000) <- c("respondent_id", "read000", "readlev000", "readpsum0")
  readpv1 <- generatepv(logdf$sumR1, logdf$sumL1, logdf$w1R)
  readpv2 <- generatepv(logdf$sumR2, logdf$sumL2, logdf$w2R)
  readpv3 <- generatepv(logdf$sumR3, logdf$sumL3, logdf$w3R)
  readpv4 <- generatepv(logdf$sumR4, logdf$sumL4, logdf$w4R)
  readpv5 <- generatepv(logdf$sumR5, logdf$sumL5, logdf$w5R)
  colnames(readpv1) <- c("respondent_id", "readpv1", "readlevpv1", "readpsum1")
  colnames(readpv2) <- c("respondent_id", "readpv2", "readlevpv2", "readpsum2")
  colnames(readpv3) <- c("respondent_id", "readpv3", "readlevpv3", "readpsum3")
  colnames(readpv4) <- c("respondent_id", "readpv4", "readlevpv4", "readpsum4")
  colnames(readpv5) <- c("respondent_id", "readpv5", "readlevpv5", "readpsum5")
}



# smaller logit
#mylogit <- glm(y ~ sumRn +  booklet_R + q, data = logdf, family = binomial)
#summary(mylogit)
#length(t(as.matrix(mylogit$coefficients)))
#
#logdf$read_smaller <- logdf$sumR*mylogit$coefficients[2] +   # fitting the number of correct answers
#  ifelse(   is.na(match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99)))==T,0,   mylogit$coefficients[match(logdf$booklet_R,substr(names(mylogit$coefficients),10,99))]   )   # fitting the booklet
#summary(logdf$read_smaller)
#
#measure_read_smaller <- logdf[duplicated(logdf$respondent_id)==F,c("respondent_id","read_smaller")]


###########################################
### tryout (correlate the various measures) and store
###########################################


if (noboot==F) {
  pvdf <- merge(list000, listpv1, all.x = T)
  pvdf <- merge(pvdf, listpv2, all.x = T)
  pvdf <- merge(pvdf, listpv3, all.x = T)
  pvdf <- merge(pvdf, listpv4, all.x = T)
  pvdf <- merge(pvdf, listpv5, all.x = T)
  pvdf <- merge(pvdf, read000, all.x = T)
  pvdf <- merge(pvdf, readpv1, all.x = T)
  pvdf <- merge(pvdf, readpv2, all.x = T)
  pvdf <- merge(pvdf, readpv3, all.x = T)
  pvdf <- merge(pvdf, readpv4, all.x = T)
  pvdf <- merge(pvdf, readpv5, all.x = T)
  str(pvdf)
  
  cor(pvdf$list000, pvdf$read000, use="complete.obs")
  with(pvdf, plot(read000, list000))
  with(pvdf, plot(readpsum0, listpsum0))
  cor(pvdf$list_err, pvdf$read_err, use="complete.obs")
  cov((pvdf$list_err+pvdf$list000), (pvdf$read_err+pvdf$read000), use="complete.obs")
  cov(pvdf$listpv1, pvdf$readpv1, use="complete.obs")
  cor(pvdf$listpv2, pvdf$readpv2, use="complete.obs")
  cor(pvdf$listpv3, pvdf$readpv3, use="complete.obs")
  cor(pvdf$listpv4, pvdf$readpv4, use="complete.obs")
  cor(pvdf$listpv5, pvdf$readpv5, use="complete.obs")
  cor(pvdf$list000, pvdf$listpv3, use="complete.obs")
  cor(pvdf$listpv1, pvdf$listpv2, use="complete.obs")
  cor(pvdf$read000, pvdf$readpv3, use="complete.obs")
  cor(pvdf$readpv1, pvdf$readpv4, use="complete.obs")
  
  
  mean (cor(pvdf$readpv1, pvdf$readpv2, use="complete.obs"), cor(pvdf$readpv1, pvdf$readpv3, use="complete.obs"), cor(pvdf$readpv1, pvdf$readpv4, use="complete.obs"), cor(pvdf$readpv1, pvdf$readpv5, use="complete.obs"), cor(pvdf$readpv2, pvdf$readpv3, use="complete.obs"), cor(pvdf$readpv2, pvdf$readpv4, use="complete.obs"), cor(pvdf$readpv2, pvdf$readpv5, use="complete.obs"), cor(pvdf$readpv3, pvdf$readpv4, use="complete.obs"), cor(pvdf$readpv3, pvdf$readpv5, use="complete.obs"), cor(pvdf$readpv4, pvdf$readpv5, use="complete.obs"))
  mean (cor(pvdf$listpv1, pvdf$listpv2, use="complete.obs"), cor(pvdf$listpv1, pvdf$listpv3, use="complete.obs"), cor(pvdf$listpv1, pvdf$listpv4, use="complete.obs"), cor(pvdf$listpv1, pvdf$listpv5, use="complete.obs"), cor(pvdf$listpv2, pvdf$listpv3, use="complete.obs"), cor(pvdf$listpv2, pvdf$listpv4, use="complete.obs"), cor(pvdf$listpv2, pvdf$listpv5, use="complete.obs"), cor(pvdf$listpv3, pvdf$listpv4, use="complete.obs"), cor(pvdf$listpv3, pvdf$listpv5, use="complete.obs"), cor(pvdf$listpv4, pvdf$listpv5, use="complete.obs"))
  
  mean (cor(pvdf$readpsum1, pvdf$readpsum2, use="complete.obs"), cor(pvdf$readpsum1, pvdf$readpsum3, use="complete.obs"), cor(pvdf$readpsum1, pvdf$readpsum4, use="complete.obs"), cor(pvdf$readpsum1, pvdf$readpsum5, use="complete.obs"), cor(pvdf$readpsum2, pvdf$readpsum3, use="complete.obs"), cor(pvdf$readpsum2, pvdf$readpsum4, use="complete.obs"), cor(pvdf$readpsum2, pvdf$readpsum5, use="complete.obs"), cor(pvdf$readpsum3, pvdf$readpsum4, use="complete.obs"), cor(pvdf$readpsum3, pvdf$readpsum5, use="complete.obs"), cor(pvdf$readpsum4, pvdf$readpsum5, use="complete.obs"))
  mean (cor(pvdf$listpsum1, pvdf$listpsum2, use="complete.obs"), cor(pvdf$listpsum1, pvdf$listpsum3, use="complete.obs"), cor(pvdf$listpsum1, pvdf$listpsum4, use="complete.obs"), cor(pvdf$listpsum1, pvdf$listpsum5, use="complete.obs"), cor(pvdf$listpsum2, pvdf$listpsum3, use="complete.obs"), cor(pvdf$listpsum2, pvdf$listpsum4, use="complete.obs"), cor(pvdf$listpsum2, pvdf$listpsum5, use="complete.obs"), cor(pvdf$listpsum3, pvdf$listpsum4, use="complete.obs"), cor(pvdf$listpsum3, pvdf$listpsum5, use="complete.obs"), cor(pvdf$listpsum4, pvdf$listpsum5, use="complete.obs"))
  
  mean(pvdf$read000 - pvdf$list000)
  mean(pvdf$readlev000 - pvdf$listlev000)
  table(pvdf$readlev000)/length(pvdf$readlev000)
  table(pvdf$listlev000)/length(pvdf$listlev000)
  mean(pvdf$readlevpv5 - pvdf$listlevpv5)
  mean(pvdf$readpsum0 - pvdf$listpsum0)
  sd(pvdf$readpsum5 - pvdf$listpsum5)
}

if (noboot==T) {
  pvdf <- merge(list000, read000, all = T)
  s2l <- var(pvdf$list_err, na.rm = T)
  s2r <- var(pvdf$read_err, na.rm = T)
  rho <- cov(pvdf$list_err, pvdf$read_err, use = "complete.obs")
  generate_pvs_noboot <- function(i) {
    err_read <- rnorm(dim(pvdf)[1], mean = 0, sd = sqrt(s2r))
    err_list <- rnorm(dim(pvdf)[1], mean = rho * err_read, sqrt(s2l - rho^2))
    pvread <- pvdf$read000 + err_read
    pvlist <- pvdf$list000 + err_list
    output <- list(pvread, pvlist)
    return(output)
  }
  determine_pl <- function(x) {
    newvector <- rep(NA,length(x))
    newvector[x>as.numeric(-dd[5])+levcoeff] <- 5
    newvector[x<as.numeric(-dd[5])+levcoeff] <- 4
    newvector[x<as.numeric(-dd[4])+levcoeff] <- 3
    newvector[x<as.numeric(-dd[3])+levcoeff] <- 2
    newvector[x<as.numeric(-dd[2])+levcoeff] <- 1
    newvector[x<as.numeric(-dd[1])+levcoeff] <- 0
    return(newvector)
  }
  
  pvs <- lapply(1:5, generate_pvs_noboot)
  pvdf$readpv1 <- pvs[[1]][[1]]
  pvdf$readpv2 <- pvs[[2]][[1]]
  pvdf$readpv3 <- pvs[[3]][[1]]
  pvdf$readpv4 <- pvs[[4]][[1]]
  pvdf$readpv5 <- pvs[[5]][[1]]
  pvdf$listpv1 <- pvs[[1]][[2]]
  pvdf$listpv2 <- pvs[[2]][[2]]
  pvdf$listpv3 <- pvs[[3]][[2]]
  pvdf$listpv4 <- pvs[[4]][[2]]
  pvdf$listpv5 <- pvs[[5]][[2]]
  dd <- -read.csv("difficread.csv")[,2]
  pvdf$readlevpv1 <- determine_pl(pvdf$readpv1)
  pvdf$readlevpv2 <- determine_pl(pvdf$readpv2)
  pvdf$readlevpv3 <- determine_pl(pvdf$readpv3)
  pvdf$readlevpv4 <- determine_pl(pvdf$readpv4)
  pvdf$readlevpv5 <- determine_pl(pvdf$readpv5)
  pvdf$readpsum1 <- exp(pvdf$readpv1+as.numeric(dd[1]))/(1 + exp(pvdf$readpv1+as.numeric(dd[1]))) + exp(pvdf$readpv1+as.numeric(dd[2]))/(1 + exp(pvdf$readpv1+as.numeric(dd[2]))) + exp(pvdf$readpv1+as.numeric(dd[3]))/(1 + exp(pvdf$readpv1+as.numeric(dd[3]))) + exp(pvdf$readpv1+as.numeric(dd[4]))/(1 + exp(pvdf$readpv1+as.numeric(dd[4]))) + exp(pvdf$readpv1+as.numeric(dd[5]))/(1 + exp(pvdf$readpv1+as.numeric(dd[5])))
  pvdf$readpsum2 <- exp(pvdf$readpv2+as.numeric(dd[1]))/(1 + exp(pvdf$readpv2+as.numeric(dd[1]))) + exp(pvdf$readpv2+as.numeric(dd[2]))/(1 + exp(pvdf$readpv2+as.numeric(dd[2]))) + exp(pvdf$readpv2+as.numeric(dd[3]))/(1 + exp(pvdf$readpv2+as.numeric(dd[3]))) + exp(pvdf$readpv2+as.numeric(dd[4]))/(1 + exp(pvdf$readpv2+as.numeric(dd[4]))) + exp(pvdf$readpv2+as.numeric(dd[5]))/(1 + exp(pvdf$readpv2+as.numeric(dd[5])))
  pvdf$readpsum3 <- exp(pvdf$readpv3+as.numeric(dd[1]))/(1 + exp(pvdf$readpv3+as.numeric(dd[1]))) + exp(pvdf$readpv3+as.numeric(dd[2]))/(1 + exp(pvdf$readpv3+as.numeric(dd[2]))) + exp(pvdf$readpv3+as.numeric(dd[3]))/(1 + exp(pvdf$readpv3+as.numeric(dd[3]))) + exp(pvdf$readpv3+as.numeric(dd[4]))/(1 + exp(pvdf$readpv3+as.numeric(dd[4]))) + exp(pvdf$readpv3+as.numeric(dd[5]))/(1 + exp(pvdf$readpv3+as.numeric(dd[5])))
  pvdf$readpsum4 <- exp(pvdf$readpv4+as.numeric(dd[1]))/(1 + exp(pvdf$readpv4+as.numeric(dd[1]))) + exp(pvdf$readpv4+as.numeric(dd[2]))/(1 + exp(pvdf$readpv4+as.numeric(dd[2]))) + exp(pvdf$readpv4+as.numeric(dd[3]))/(1 + exp(pvdf$readpv4+as.numeric(dd[3]))) + exp(pvdf$readpv4+as.numeric(dd[4]))/(1 + exp(pvdf$readpv4+as.numeric(dd[4]))) + exp(pvdf$readpv4+as.numeric(dd[5]))/(1 + exp(pvdf$readpv4+as.numeric(dd[5])))
  pvdf$readpsum5 <- exp(pvdf$readpv5+as.numeric(dd[1]))/(1 + exp(pvdf$readpv5+as.numeric(dd[1]))) + exp(pvdf$readpv5+as.numeric(dd[2]))/(1 + exp(pvdf$readpv5+as.numeric(dd[2]))) + exp(pvdf$readpv5+as.numeric(dd[3]))/(1 + exp(pvdf$readpv5+as.numeric(dd[3]))) + exp(pvdf$readpv5+as.numeric(dd[4]))/(1 + exp(pvdf$readpv5+as.numeric(dd[4]))) + exp(pvdf$readpv5+as.numeric(dd[5]))/(1 + exp(pvdf$readpv5+as.numeric(dd[5])))
  if(proB1==T) {
    pvdf$readpsum1 <- exp(pvdf$readpv1+as.numeric(dd[3]))/(1 + exp(pvdf$readpv1+as.numeric(dd[3])))
    pvdf$readpsum2 <- exp(pvdf$readpv2+as.numeric(dd[3]))/(1 + exp(pvdf$readpv2+as.numeric(dd[3])))
    pvdf$readpsum3 <- exp(pvdf$readpv3+as.numeric(dd[3]))/(1 + exp(pvdf$readpv3+as.numeric(dd[3])))
    pvdf$readpsum4 <- exp(pvdf$readpv4+as.numeric(dd[3]))/(1 + exp(pvdf$readpv4+as.numeric(dd[3])))
    pvdf$readpsum5 <- exp(pvdf$readpv5+as.numeric(dd[3]))/(1 + exp(pvdf$readpv5+as.numeric(dd[3])))
  }
  
  dd <- -read.csv("difficlist.csv")[,2]
  pvdf$listlevpv1 <- determine_pl(pvdf$listpv1)
  pvdf$listlevpv2 <- determine_pl(pvdf$listpv2)
  pvdf$listlevpv3 <- determine_pl(pvdf$listpv3)
  pvdf$listlevpv4 <- determine_pl(pvdf$listpv4)
  pvdf$listlevpv5 <- determine_pl(pvdf$listpv5)
  pvdf$listpsum1 <- exp(pvdf$listpv1+as.numeric(dd[1]))/(1 + exp(pvdf$listpv1+as.numeric(dd[1]))) + exp(pvdf$listpv1+as.numeric(dd[2]))/(1 + exp(pvdf$listpv1+as.numeric(dd[2]))) + exp(pvdf$listpv1+as.numeric(dd[3]))/(1 + exp(pvdf$listpv1+as.numeric(dd[3]))) + exp(pvdf$listpv1+as.numeric(dd[4]))/(1 + exp(pvdf$listpv1+as.numeric(dd[4]))) + exp(pvdf$listpv1+as.numeric(dd[5]))/(1 + exp(pvdf$listpv1+as.numeric(dd[5])))
  pvdf$listpsum2 <- exp(pvdf$listpv2+as.numeric(dd[1]))/(1 + exp(pvdf$listpv2+as.numeric(dd[1]))) + exp(pvdf$listpv2+as.numeric(dd[2]))/(1 + exp(pvdf$listpv2+as.numeric(dd[2]))) + exp(pvdf$listpv2+as.numeric(dd[3]))/(1 + exp(pvdf$listpv2+as.numeric(dd[3]))) + exp(pvdf$listpv2+as.numeric(dd[4]))/(1 + exp(pvdf$listpv2+as.numeric(dd[4]))) + exp(pvdf$listpv2+as.numeric(dd[5]))/(1 + exp(pvdf$listpv2+as.numeric(dd[5])))
  pvdf$listpsum3 <- exp(pvdf$listpv3+as.numeric(dd[1]))/(1 + exp(pvdf$listpv3+as.numeric(dd[1]))) + exp(pvdf$listpv3+as.numeric(dd[2]))/(1 + exp(pvdf$listpv3+as.numeric(dd[2]))) + exp(pvdf$listpv3+as.numeric(dd[3]))/(1 + exp(pvdf$listpv3+as.numeric(dd[3]))) + exp(pvdf$listpv3+as.numeric(dd[4]))/(1 + exp(pvdf$listpv3+as.numeric(dd[4]))) + exp(pvdf$listpv3+as.numeric(dd[5]))/(1 + exp(pvdf$listpv3+as.numeric(dd[5])))
  pvdf$listpsum4 <- exp(pvdf$listpv4+as.numeric(dd[1]))/(1 + exp(pvdf$listpv4+as.numeric(dd[1]))) + exp(pvdf$listpv4+as.numeric(dd[2]))/(1 + exp(pvdf$listpv4+as.numeric(dd[2]))) + exp(pvdf$listpv4+as.numeric(dd[3]))/(1 + exp(pvdf$listpv4+as.numeric(dd[3]))) + exp(pvdf$listpv4+as.numeric(dd[4]))/(1 + exp(pvdf$listpv4+as.numeric(dd[4]))) + exp(pvdf$listpv4+as.numeric(dd[5]))/(1 + exp(pvdf$listpv4+as.numeric(dd[5])))
  pvdf$listpsum5 <- exp(pvdf$listpv5+as.numeric(dd[1]))/(1 + exp(pvdf$listpv5+as.numeric(dd[1]))) + exp(pvdf$listpv5+as.numeric(dd[2]))/(1 + exp(pvdf$listpv5+as.numeric(dd[2]))) + exp(pvdf$listpv5+as.numeric(dd[3]))/(1 + exp(pvdf$listpv5+as.numeric(dd[3]))) + exp(pvdf$listpv5+as.numeric(dd[4]))/(1 + exp(pvdf$listpv5+as.numeric(dd[4]))) + exp(pvdf$listpv5+as.numeric(dd[5]))/(1 + exp(pvdf$listpv5+as.numeric(dd[5])))
  if(proB1==T) {
    pvdf$listpsum1 <- exp(pvdf$listpv1+as.numeric(dd[3]))/(1 + exp(pvdf$listpv1+as.numeric(dd[3])))
    pvdf$listpsum2 <- exp(pvdf$listpv2+as.numeric(dd[3]))/(1 + exp(pvdf$listpv2+as.numeric(dd[3])))
    pvdf$listpsum3 <- exp(pvdf$listpv3+as.numeric(dd[3]))/(1 + exp(pvdf$listpv3+as.numeric(dd[3])))
    pvdf$listpsum4 <- exp(pvdf$listpv4+as.numeric(dd[3]))/(1 + exp(pvdf$listpv4+as.numeric(dd[3])))
    pvdf$listpsum5 <- exp(pvdf$listpv5+as.numeric(dd[3]))/(1 + exp(pvdf$listpv5+as.numeric(dd[3])))
  }
  str(pvdf)
  
  cor(pvdf$list000, pvdf$read000, use="complete.obs")
  with(pvdf, plot(read000, list000))
  with(pvdf, plot(readpsum0, listpsum0))
  cor(pvdf$listpv1, pvdf$readpv1, use="complete.obs")
  cor(pvdf$listpv2, pvdf$readpv2, use="complete.obs")
  cor(pvdf$listpv3, pvdf$readpv3, use="complete.obs")
  cor(pvdf$listpv4, pvdf$readpv4, use="complete.obs")
  cor(pvdf$listpv5, pvdf$readpv5, use="complete.obs")
  cor(pvdf$list000, pvdf$listpv3, use="complete.obs")
  cor(pvdf$listpv1, pvdf$listpv2, use="complete.obs")
  cor(pvdf$read000, pvdf$readpv3, use="complete.obs")
  cor(pvdf$readpv1, pvdf$readpv4, use="complete.obs")
  
  
  mean (cor(pvdf$readpv1, pvdf$readpv2, use="complete.obs"), cor(pvdf$readpv1, pvdf$readpv3, use="complete.obs"), cor(pvdf$readpv1, pvdf$readpv4, use="complete.obs"), cor(pvdf$readpv1, pvdf$readpv5, use="complete.obs"), cor(pvdf$readpv2, pvdf$readpv3, use="complete.obs"), cor(pvdf$readpv2, pvdf$readpv4, use="complete.obs"), cor(pvdf$readpv2, pvdf$readpv5, use="complete.obs"), cor(pvdf$readpv3, pvdf$readpv4, use="complete.obs"), cor(pvdf$readpv3, pvdf$readpv5, use="complete.obs"), cor(pvdf$readpv4, pvdf$readpv5, use="complete.obs"))
  mean (cor(pvdf$listpv1, pvdf$listpv2, use="complete.obs"), cor(pvdf$listpv1, pvdf$listpv3, use="complete.obs"), cor(pvdf$listpv1, pvdf$listpv4, use="complete.obs"), cor(pvdf$listpv1, pvdf$listpv5, use="complete.obs"), cor(pvdf$listpv2, pvdf$listpv3, use="complete.obs"), cor(pvdf$listpv2, pvdf$listpv4, use="complete.obs"), cor(pvdf$listpv2, pvdf$listpv5, use="complete.obs"), cor(pvdf$listpv3, pvdf$listpv4, use="complete.obs"), cor(pvdf$listpv3, pvdf$listpv5, use="complete.obs"), cor(pvdf$listpv4, pvdf$listpv5, use="complete.obs"))
  
  mean (cor(pvdf$readpsum1, pvdf$readpsum2, use="complete.obs"), cor(pvdf$readpsum1, pvdf$readpsum3, use="complete.obs"), cor(pvdf$readpsum1, pvdf$readpsum4, use="complete.obs"), cor(pvdf$readpsum1, pvdf$readpsum5, use="complete.obs"), cor(pvdf$readpsum2, pvdf$readpsum3, use="complete.obs"), cor(pvdf$readpsum2, pvdf$readpsum4, use="complete.obs"), cor(pvdf$readpsum2, pvdf$readpsum5, use="complete.obs"), cor(pvdf$readpsum3, pvdf$readpsum4, use="complete.obs"), cor(pvdf$readpsum3, pvdf$readpsum5, use="complete.obs"), cor(pvdf$readpsum4, pvdf$readpsum5, use="complete.obs"))
  mean (cor(pvdf$listpsum1, pvdf$listpsum2, use="complete.obs"), cor(pvdf$listpsum1, pvdf$listpsum3, use="complete.obs"), cor(pvdf$listpsum1, pvdf$listpsum4, use="complete.obs"), cor(pvdf$listpsum1, pvdf$listpsum5, use="complete.obs"), cor(pvdf$listpsum2, pvdf$listpsum3, use="complete.obs"), cor(pvdf$listpsum2, pvdf$listpsum4, use="complete.obs"), cor(pvdf$listpsum2, pvdf$listpsum5, use="complete.obs"), cor(pvdf$listpsum3, pvdf$listpsum4, use="complete.obs"), cor(pvdf$listpsum3, pvdf$listpsum5, use="complete.obs"), cor(pvdf$listpsum4, pvdf$listpsum5, use="complete.obs"))
  
  mean(pvdf$read000 - pvdf$list000)
  mean(pvdf$readlev000 - pvdf$listlev000)
  table(pvdf$readlev000)/length(pvdf$readlev000)
  table(pvdf$listlev000)/length(pvdf$listlev000)
  mean(pvdf$readlevpv5 - pvdf$listlevpv5)
  mean(pvdf$readpsum0 - pvdf$listpsum0)
  sd(pvdf$readpsum5 - pvdf$listpsum5)
  
}


saveRDS(pvdf, "pvdf.rds")

#prova_smaller <- merge(measure_list_smaller,measure_read_smaller)
#str(prova_smaller)
#cor(prova_smaller$list_smaller, prova_smaller$read_smaller, use="complete.obs")
#with(prova_smaller, plot(read_smaller, list_smaller))

  





