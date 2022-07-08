
dim(s)[1]

##################
### proper estimates
##################



### regressions

rgrss = list(
  model1 = c(
    'I08_ST_A_S19B',
    'I03_ST_A_S03A',
    'I08_ST_A_S15A',
    'I03_ST_A_S26B',
    'I08_ST_A_S01A',
    'I08_ST_A_S02A',
    'I03_ST_A_S31A',
    'cnt',
    'I01_ST_M_S40A',
    'I01_ST_M_S44A_cont',
    'classoralQ58',
    'classwrittenQ58',
    'travelQ45',
    'face2faceQ29',
    'classictQ51',
    'txtbQ52',
    'I09_IN_M_S50A',
    'goodclass2Q54'
  ),
  model2 = c('I08_ST_A_S19B',
             'I03_ST_A_S03A',
             'I08_ST_A_S15A',
             'I03_ST_A_S26B',
             'I08_ST_A_S01A',
             'I08_ST_A_S02A',
             'I03_ST_A_S31A',
             'cnt',
             'I01_ST_M_S40A',
             'I01_ST_M_S44A_cont',
             'classoralQ58',
             'classwrittenQ58',
             'travelQ45',
             'face2faceQ29'
  ), 
  model3= c(
    'I08_ST_A_S19B',
    'I03_ST_A_S03A',
    'I08_ST_A_S15A',
    'I03_ST_A_S26B',
    'I08_ST_A_S01A',
    'I08_ST_A_S02A',
    'I03_ST_A_S31A',
    'cnt',
    'I01_ST_M_S40A',
    'I01_ST_M_S44A_cont',
    'travelQ45',
    'face2faceQ29'
  ), 
  model4 = c(
    'I08_ST_A_S19B',
    'I03_ST_A_S03A',
    'I08_ST_A_S15A',
    'I03_ST_A_S26B',
    'I08_ST_A_S01A',
    'I08_ST_A_S02A',
    'I03_ST_A_S31A',
    'cnt',
    'I01_ST_M_S40A',
    'I01_ST_M_S44A_cont'
  )
)


#unweighted
fd = svydesign(ids = ~JKrep, strata = ~JKzone, nest = TRUE, data = s, weights = ~unweighted)
frd = as.svrepdesign(design = fd, type='Fay', fay.rho=.5)
cefrREADcefrLIST_SQ = lapply(rgrss, regfun, test='dreadlist')
saveRDS(cefrREADcefrLIST_SQ, file='regr_cefrREADcefrLIST_SQ.rds')
cefrREADcefrLIST_SQ_summ = lapply(cefrREADcefrLIST_SQ, function(x)x$summary)
cefrREADcefrLIST_SQ_df = lapply(cefrREADcefrLIST_SQ, function(x)x$regressions[[1]]$df.null)
cefrREADcefrLIST_SQ_vcov = lapply(cefrREADcefrLIST_SQ, function(x)x$regressions[[1]]$cov.unscaled)
rm(cefrREADcefrLIST_SQ)

#weighted
fd = svydesign(ids = ~JKrep, strata = ~JKzone, nest = TRUE, data = s, weights = ~countryavg_w)
frd = as.svrepdesign(design = fd, type='Fay', fay.rho=.5)
cefrREADcefrLIST_SQW = lapply(rgrss, regfun, test='dreadlist')
saveRDS(cefrREADcefrLIST_SQW, file='regr_cefrREADcefrLIST_SQW.rds')
cefrREADcefrLIST_SQW_summ = lapply(cefrREADcefrLIST_SQW, function(x)x$summary)
cefrREADcefrLIST_SQW_df = lapply(cefrREADcefrLIST_SQW, function(x)x$regressions[[1]]$df.null)
cefrREADcefrLIST_SQW_vcov = lapply(cefrREADcefrLIST_SQW, function(x)x$regressions[[1]]$cov.unscaled)
rm(cefrREADcefrLIST_SQW)



csvtable <- rbind(cefrREADcefrLIST_SQ_summ[[1]],c(cefrREADcefrLIST_SQ_df[[1]],NA,NA,NA),cefrREADcefrLIST_SQ_summ[[2]],c(cefrREADcefrLIST_SQ_df[[2]],NA,NA,NA),cefrREADcefrLIST_SQ_summ[[3]],c(cefrREADcefrLIST_SQ_df[[3]],NA,NA,NA),cefrREADcefrLIST_SQ_summ[[4]],c(cefrREADcefrLIST_SQ_df[[4]],NA,NA,NA))
write.csv(csvtable, "cefrREADcefrLIST_SQ_summ.csv")
write.csv(cefrREADcefrLIST_SQ_vcov[[4]], "cefrREADcefrLIST_SQ_cov.csv")
csvtableW <- rbind(cefrREADcefrLIST_SQW_summ[[1]],c(cefrREADcefrLIST_SQW_df[[1]],NA,NA,NA),cefrREADcefrLIST_SQW_summ[[2]],c(cefrREADcefrLIST_SQW_df[[2]],NA,NA,NA),cefrREADcefrLIST_SQW_summ[[3]],c(cefrREADcefrLIST_SQW_df[[3]],NA,NA,NA),cefrREADcefrLIST_SQW_summ[[4]],c(cefrREADcefrLIST_SQW_df[[4]],NA,NA,NA))
write.csv(csvtableW, "cefrREADcefrLIST_SQW_summ.csv")
write.csv(cefrREADcefrLIST_SQW_vcov[[4]], "cefrREADcefrLIST_SQW_cov.csv")



z_estimate <- F
if (z_estimate==T) {
  # with z stdsation, unweighted
  fd = svydesign(ids = ~JKrep, strata = ~JKzone, nest = TRUE, data = s, weights = ~unweighted)
  frd = as.svrepdesign(design = fd, type='Fay', fay.rho=.5)
  zREADzLIST_SQ = lapply(rgrss, regfun, test='zREADzLIST')
  saveRDS(zREADzLIST_SQ, file='regr_zREADzLIST_SQ.rds')
  zREADzLIST_SQ_summ = lapply(zREADzLIST_SQ, function(x)x$summary)
  zREADzLIST_SQ_df = lapply(zREADzLIST_SQ, function(x)x$regressions[[1]]$df.null)
  zREADzLIST_SQ_vcov = lapply(zREADzLIST_SQ, function(x)x$regressions[[1]]$cov.unscaled)
  rm(zREADzLIST_SQ)
  
  csvtablez <- rbind(zREADzLIST_SQ_summ[[1]],c(zREADzLIST_SQ_df[[1]],NA,NA,NA),zREADzLIST_SQ_summ[[2]],c(zREADzLIST_SQ_df[[2]],NA,NA,NA),zREADzLIST_SQ_summ[[3]],c(zREADzLIST_SQ_df[[3]],NA,NA,NA),zREADzLIST_SQ_summ[[4]],c(zREADzLIST_SQ_df[[4]],NA,NA,NA))
  write.csv(csvtablez, "zREADzLIST_SQ_summ.csv")
  write.csv(zREADzLIST_SQ_vcov[[4]], "zREADzLIST_SQ_cov.csv")
}
  




rgrssnat = list(
  model1 = c(
    'I08_ST_A_S19B',
    'I03_ST_A_S03A',
    'I08_ST_A_S15A',
    'I03_ST_A_S26B',
    'I08_ST_A_S01A',
    'I08_ST_A_S02A',
    'I03_ST_A_S31A',
    'travelQ45',
    'face2faceQ29',
    'I01_ST_M_S40A',
    'I01_ST_M_S44A_cont'
  )
)

countryestimate <- F
if (countryestimate==T) {
  host_table <- data.frame(coeff=0, stderr=0,t=0,p=0)
  for (nat in c('BE de', 'BE fr', 'BE nl', 'BG','EE','EL','ES','FR','HR','MT','NL','PL','PT','SE','SI')) {
    print(nat)
    fd <- svydesign(ids = ~JKrep, strata = ~JKzone, nest = TRUE, data = s[news$cnt==nat,], weights = ~countryavg_w)
    frd <- as.svrepdesign(design = fd, type='Fay', fay.rho=.5)
    estnat_SQ <- lapply(rgrssnat, regfun, test='dreadlist')
    natresults <- as.data.frame(lapply(estnat_SQ, function(x)x$summary))[length(rgrssnat$model1):(length(rgrssnat$model1)+1),]
    colnames(natresults) <- c("coeff", "stderr", "t", "p")
    rownames(natresults) <- paste0(nat, rownames(natresults))
    host_table <- rbind(host_table,natresults)
  }
  write.csv(host_table,"national_results.csv")
}





#########################################################################################
### conditional model
#########################################################################################


rgrss_r = list(
  model1 = c(
    'I08_ST_A_S19B',
    'I03_ST_A_S03A',
    'I08_ST_A_S15A',
    'I03_ST_A_S26B',
    'I08_ST_A_S01A',
    'I08_ST_A_S02A',
    'I03_ST_A_S31A',
    'cnt',
    'myl',
    'I01_ST_M_S40A',
    'I01_ST_M_S44A_cont'
  )
)

rgrss_l = list(
  model1 = c(
    'I08_ST_A_S19B',
    'I03_ST_A_S03A',
    'I08_ST_A_S15A',
    'I03_ST_A_S26B',
    'I08_ST_A_S01A',
    'I08_ST_A_S02A',
    'I03_ST_A_S31A',
    'cnt',
    'myr',
    'I01_ST_M_S40A',
    'I01_ST_M_S44A_cont'
  )
)



pv <- as.numeric(1)

news <- s

# conditioning proficiencies to be used as independent vars
news$myr <- eval(str2expression(paste0("news$PV",as.character(pv),"_READ")))
news$myl <- eval(str2expression(paste0("news$PV",as.character(pv),"_LIST")))
print(pv)
print(summary(news$myr))
print(summary(news$myl))

# conditional models, weighted
fd = svydesign(ids = ~JKrep, strata = ~JKzone, nest = TRUE, data = news, weights = ~countryavg_w)
frd = as.svrepdesign(design = fd, type='Fay', fay.rho=.5)
newest_r <- lapply(rgrss_r, regfun, test='READ')
print(newest_r$model1$summary[1][1])
newest_l <- lapply(rgrss_l, regfun, test='LIST')
print(newest_l$model1$summary[1][1])

# compile and return output
cmod <- data.frame(b_r=newest_r$model1$summary[,1], sd_r=newest_r$model1$summary[,2], b_l=newest_l$model1$summary[,1], sd_l=newest_l$model1$summary[,2] )
colnames(cmod) <- paste0(c("b_r", "sd_r", "b_l", "sd_l"),pv)


for (pv in 2:5) {
  #pv <- as.numeric(2)
  
  news <- s
  
  # conditioning proficiencies to be used as independent vars
  news$myr <- eval(str2expression(paste0("news$PV",as.character(pv),"_READ")))
  news$myl <- eval(str2expression(paste0("news$PV",as.character(pv),"_LIST")))
  print(pv)
  print(summary(news$myr))
  print(summary(news$myl))
  
  # conditional models, weighted
  fd = svydesign(ids = ~JKrep, strata = ~JKzone, nest = TRUE, data = news, weights = ~countryavg_w)
  frd = as.svrepdesign(design = fd, type='Fay', fay.rho=.5)
  newest_r <- lapply(rgrss_r, regfun, test='READ')
  print(newest_r$model1$summary[1][1])
  newest_l <- lapply(rgrss_l, regfun, test='LIST')
  print(newest_l$model1$summary[1][1])
  
  # compile and return output
  output <- data.frame(b_r=newest_r$model1$summary[,1], sd_r=newest_r$model1$summary[,2], b_l=newest_l$model1$summary[,1], sd_l=newest_l$model1$summary[,2] )
  colnames(output) <- paste0(c("b_r", "sd_r", "b_l", "sd_l"),pv)
  cmod <- cbind(cmod,output)
}



cmod$br <- .2 * (cmod$b_r1 + cmod$b_r2 + cmod$b_r3 + cmod$b_r4 + cmod$b_r5) 
cmod$bl <- .2 * (cmod$b_l1 + cmod$b_l2 + cmod$b_l3 + cmod$b_l4 + cmod$b_l5) 
cmod$sdr_within <- sqrt(round(.2 * (cmod$sd_r1^2 + cmod$sd_r2^2 + cmod$sd_r3^2 + cmod$sd_r4^2 + cmod$sd_r5^2), 5) )
cmod$sdl_within <- sqrt(round(.2 * (cmod$sd_l1^2 + cmod$sd_l2^2 + cmod$sd_l3^2 + cmod$sd_l4^2 + cmod$sd_l5^2), 5) )
cmod$sdr_between <- sqrt(round(.2 * (cmod$b_r1^2 + cmod$b_r2^2 + cmod$b_r3^2 + cmod$b_r4^2 + cmod$b_r5^2) - cmod$br^2 , 5 ))
cmod$sdl_between <- sqrt(round(.2 * (cmod$b_l1^2 + cmod$b_l2^2 + cmod$b_l3^2 + cmod$b_l4^2 + cmod$b_l5^2) - cmod$bl^2 , 5 ))
cmod$sdr <- cmod$sdr_within + cmod$sdr_between
cmod$sdl <- cmod$sdl_within + cmod$sdl_between

cmod_bsd <- cmod[,c("br","bl","sdr","sdl")]
write.csv(cmod_bsd,"cmod.csv")

