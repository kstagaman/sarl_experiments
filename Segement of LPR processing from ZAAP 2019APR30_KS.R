##code to run LPR and LSR independently
## L. Truong, April 30, 2019
rm(list=ls())

##all the packages you'll  need
library(MASS)
library(reshape2)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(dplyr)
library(NMF) #
library(tidyr)
library(stringr)
library(grid)
require(car)
require(pracma) #
library(splines)
library(gridExtra)
library(cowplot)
library(rlang)
library(tidyverse)
library(gplots) #
library(drc) #
library(purrrlyr) #
select <- dplyr::select
rename <- dplyr::rename


Args <- commandArgs(TRUE)
work.drive <- Args[1]
work.path <- Args[2]
user_Path <- Args[3]
Result_Path <- Args[4]
format_files <- list.files(paste(work.drive,work.path,user_Path,Result_Path,sep="/"),recursive=T,pattern="*.xls")
fnam_csv <- grep("Single_Tab", format_files, value=T)
key_csv<- grep("Table", format_files, value=T)

data.dir <- paste(work.drive,work.path,user_Path,Result_Path,sep="/") # dir(data.dir)
data.file <- substr(basename(fnam_csv), 1, nchar(basename(fnam_csv)) - 4)#"Collapsed_Tanguay_ZF_Toxcast_1078_2013MAY27.csv"#"Collapsed_Tanguay_ZF_Toxcast_1078_2013MAY1.csv"#"Collapsed_Tanguay_ZF_Toxcast_1078_FINAL.csv" # "Tanguay_ZF_Toxcast_data_1073chem_2013FEB07.csv" #"Tanguay_ZF_Toxcast_data_1073chem_2013JAN31.csv"
TXid.file <- key_csv # {filename, NULL}
remove_dead <- T

output.prefix <- toupper(format(Sys.Date(),"%Y%b%d"))

#########################################################################
####
###manipulate and create the dev tox file needed for LPR
####
#########################################################################
print("Starting Analysis")
org.dat<- read.table(file=paste0(data.dir,data.file, ".xls"), sep="\t", header=T, skipNul=T,
                     as.is= T, na.strings="NA", strip.white = T)

###change the newest ZAAP output into the format everyone is use to
org.dat$well <- paste0(org.dat$row_id, org.dat$col_id)
##make sure the well column has "0" in front of everything
org.dat$well <- factor(mapply(gsub, list('\\d+'), sprintf('%.2d', as.numeric(regmatches(org.dat$well, gregexpr('\\d+', org.dat$well)))), org.dat$well))
corresponding_barcodes <- unique(org.dat$plate_id)

dat <- data.frame()
print(corresponding_barcodes)
###Fill in control data based on solvent
for(pi in corresponding_barcodes){
  dat.temp <- org.dat %>% filter(plate_id == pi)


  ##in the case that people put in treatment, but no conc value
  index.ctrl <- which(is.na(dat.temp$value))
  dat.temp[index.ctrl,'value'] <- 0


  ##fill in conc for the rest.
  idex.in <- which(is.na(dat.temp$chemical_id) & !is.na(dat.temp$solvent_map))

  if(length(unique(dat.temp$chemical_id))>2){
    dat.rep <- data.frame()
    ci <- unique(dat.temp$chemical_id)[!is.na(unique(dat.temp$chemical_id))]
    for(i in ci){
      chem.name <- unique(dat.temp[which(dat.temp$chemical_id == i),'chemical_name'])
      chem.bottle.id <- unique(dat.temp[which(dat.temp$chemical_id == i),'bottle_id'])
      chem.cas <- unique(dat.temp[which(dat.temp$chemical_id == i),'casrn'])
      blank <- dat.temp[idex.in,]
      blank$chemical_id <- i
      blank$chemical_name <- chem.name
      blank$bottle_id <- chem.bottle.id
      blank$casrn <- chem.cas
      blank$value <- 0
      dat.rep <- rbind(dat.rep, blank)
    }
    dat.temp <- dat.temp[-idex.in,]
    dat.temp <- rbind(dat.temp, dat.rep)
  } else {
    dat.temp[idex.in,'chemical_id']<- unique(dat.temp$chemical_id)[!is.na(unique(dat.temp$chemical_id))]
    dat.temp[idex.in,'chemical_name']<- unique(dat.temp$chemical_name)[unique(dat.temp$chemical_name) !=""]
    dat.temp[idex.in,'bottle_id']<- unique(dat.temp$bottle_id)[unique(dat.temp$bottle_id) !=""]
    if(length(unique(dat.temp$casrn)) > 1){
      dat.temp[idex.in,'casrn']<- unique(dat.temp$casrn)[-which(unique(dat.temp$casrn)=="")]
    } else{
      dat.temp[idex.in,'casrn']<- unique(dat.temp$casrn)
    }

    dat.temp[idex.in,'value']<- 0
  }
  dat <- rbind(dat, dat.temp)

}

conc_units <- unique(dat$units)

print(table(dat$value, dat$bottle_id))
print(table(dat$value, dat$chemical_name))

##########################
### mapping of TX IDs
##########################
TX.id.dat2 <- unique(dat[c('bottle_id','casrn', 'chemical_id','chemical_name')])
TX.id.dat2$casrn <- as.numeric(as.character(TX.id.dat2$casrn))

TX.id.dat.org <- read.table(file=paste0(data.dir,TXid.file), sep="\t", header=T, skipNul=T, as.is=T)
TX.id.dat.org$casrn <- as.numeric(as.character(TX.id.dat.org$casrn))


##only if there is more unique ids
if(length(unique(TX.id.dat2$bottle_id)) < length(unique(TX.id.dat.org$bottle_id))){
  TX.id.dat <- TX.id.dat.org %>% inner_join(TX.id.dat2)
}else{
  TX.id.dat <- TX.id.dat2
}
head(TX.id.dat)

TX.id.dat <- TX.id.dat %>% dplyr::rename(bottle.id = bottle_id,
                                         chemical.name = chemical_name,
                                         chemical.id = chemical_id,
                                         casrn = casrn)

###clean up the mapping so all the comma's are removed
TX.id.dat$chemical.name <- gsub(",", "_", TX.id.dat$chemical.name)
TX.id.dat <- TX.id.dat[!duplicated(TX.id.dat),]
print(TX.id.dat)

dev.tox <- dat %>% select(plate_id, bottle_id, well,value, units, date_of_plate,
                          MO24, DP24, SM24, NC24, MORT, YSE_, AXIS, EYE_, SNOU, JAW_,
                          OTIC, PE__, BRAI, SOMI, PFIN, CFIN, PIG_, CIRC, TRUN, SWIM,
                          NC__, TR__, DNC_, comment) %>%
  rename(plate.id = plate_id, bottle.id = bottle_id,
         conc = value, date = date_of_plate)

##create folder structure
dir.create(file.path(paste(data.dir, "Devtox_by_bottleid_per_day/",sep="")))
dir.create(file.path(paste(data.dir, "Devtox_by_bottleid/",sep="")))
dir.create(file.path(paste(data.dir, "Devtox_by_chemicalid/",sep="")))
dir.create(file.path(paste(data.dir, "GUI files/",sep="")))
dir.create(file.path(paste(data.dir, "GUI files/NoRemoval/",sep="")))
dir.create(file.path(paste(data.dir, "GUI files/RemovedAffected/",sep="")))
dir.create(file.path(paste(data.dir, "Reporting_Material/",sep="")))


#####################
###for LPR
#####################
format_files <- list.files(paste(work.drive,work.path,user_Path,Result_Path,sep="/"),recursive=T,pattern="*.xls")

if(length(grep("LPR", format_files)) >0 ){

  dir.create(file.path(paste(data.dir, "LPR/", sep="")))
  conclists <- dev.tox[-which(names(dev.tox) == "comment"),] %>% join(TX.id.dat)
  conclists <- conclists %>% select(plate.id, bottle.id, chemical.id, chemical.name, casrn,
                                    well, conc, units, date, MO24, DP24, SM24, NC24, MORT,
                                    YSE_, AXIS, EYE_, SNOU, JAW_, OTIC, PE__, BRAI, SOMI, PFIN, CFIN,
                                    PIG_, CIRC, TRUN, SWIM, NC__, TR__, DNC_, comment)


  ##change the conclists so there is two digits always
  Any_End<-vector(length=dim(conclists)[1])
  for(i in 1:dim(conclists)[1]){
    if(sum(conclists[i,-c(1:9,32,33)], na.rm=T)>0){
      Any_End[i]<-1
    }else{
      Any_End[i]<-0
    }
    if(i%%10000==0){
      cat(i,"\n")
    }
  }

  conclists<-cbind(conclists,Any_End)

  conclists$well <- factor(mapply(gsub, list('\\d+'),
                                  sprintf('%.02d', as.numeric(regmatches(conclists$well, gregexpr('\\d+', conclists$well)))), conclists$well))

  ###summary counts
  morphsum <- conclists %>% group_by(chemical.id, chemical.name, conc) %>%
    summarise(morph.normal = length(which(Any_End == 0)),
              org_ss = length(conc))
  morphsum$perc.norm <- round(morphsum$morph.normal/morphsum$org_ss*100,3)
  morphsum$behavior.include <- 1
  morphsum[which(morphsum$perc.norm <=60 ),'behavior.include'] <- 0


  ###read in LPR files
  read.data <- function(file){
    dat <- read.table(file=paste(data.dir, "/", file,sep=""),sep="\t",  header=T, skipNul=T)
    print(file)
    return(dat)
  }

  print("Done reading behavior files")
  light_dark <- do.call(rbind, lapply(grep("LPR", list.files(data.dir, recursive=T, ".xls"),value=T), read.data))
  light_dark$well <- unlist(lapply(strsplit(as.character(light_dark$vd_aname), "\\-"), "[", 1))

  if(length(grep("STR", format_files)) >0){
    startle <- read.table(file=paste0(data.dir,grep("STR", format_files, value=T)), sep="\t", header=T, skipNul = T)
    startle$well <- unlist(lapply(strsplit(as.character(startle$vd_aname), "\\-"), "[", 1))
  }

  assays <- c("LPR", "LSR")

  for(assay in assays){
    print(paste0("Starting analysis of ", assay))
    if(assay == "LPR"){set2 <- light_dark}
    if(assay == "LSR"){set2 <- startle}

    set2$vd_samplecode <- as.integer(as.character(set2$vd_samplecode))
    set2 <- set2 %>% rename(plate.id = vd_samplecode) %>%
      inner_join(conclists, by=c("plate.id", "well")) %>%
      unite(plate_well, plate.id, well, remove = F)

    set2$total <- set2$vd_smldist + set2$vd_lardist
    set2<-subset(set2,set2$DNC_==0) #######remove DNC

    ###write out barcodes + box
    boxset <- set2 %>% distinct(vp_locationcode, plate.id)
    write.csv(boxset,file=paste(data.dir,"behavior_5day_plates_Zebraboxes_key_",Sys.Date(),"_",data.file, ".csv",sep=""), row.names=F)

    ###############################################
    ####do some actual analysis now
    ###############################################
    all <- set2 %>% select(well, plate.id, chemical.id, chemical.name, date, conc,vd_start, vd_inadist,
                           vd_smldist, vd_lardist, plate_well, total) %>%
      rename(starttime = vd_start)
    raw_files2 <- all %>% arrange(plate_well, conc, starttime)
    chemical_name <- all %>% distinct(chemical.name)

    ##formatting all raw data
    raw_files <- set2
    raw_files[which(raw_files$Any_End== 1), 'total'] <- NA
    raw_files <- raw_files %>% select(well, plate.id, chemical.id, chemical.name, vp_locationcode,
                                      conc,vd_start, vd_inadist, vd_smldist, vd_lardist, plate_well, total) %>%
      rename(starttime = vd_start) %>%
      arrange(plate_well, conc, starttime)

    behav5d<- dcast(raw_files,chemical.id+conc+plate.id+well ~ starttime, value.var="total")
    colnames(behav5d) <- c("chemical.id", "conc", "plate", "well", paste0("t", seq(0,ncol(behav5d)-5)))
    behav5d <- behav5d[order(behav5d[,1], behav5d$conc, behav5d$plate),]

    ##No removal for GUI
    # raw_files2$conc <- paste0(raw_files2$conc, " ", raw_files2$Unit)
    behav5d.unfiltered<- dcast(raw_files2,chemical.id+conc+plate.id+well ~ starttime, value.var="total")
    colnames(behav5d.unfiltered) <- c("chemical.id", "conc", "plate", "well", paste0("t", seq(0,ncol(behav5d.unfiltered)-5)))
    behav5d.unfiltered <- behav5d.unfiltered[order(behav5d.unfiltered[,1], behav5d.unfiltered$conc, behav5d.unfiltered$plate),]

    ###No removal for GUI per day
    behav5d.date <- dcast(raw_files2,chemical.id+conc+date+plate.id+well ~ starttime, value.var="total")
    colnames(behav5d.date) <- c("chemical.id", "conc", "date", "plate", "well", paste0("t", seq(0,ncol(behav5d.date)-6)))
    behav5d.date <- behav5d.date[order(behav5d.date[,1], behav5d.date$conc, behav5d.date$plate),]

    ##write out all the files
    write.csv(raw_files,file=paste(data.dir,"GUI files/RemovedAffected/behavior_5day_",assay, "_",Sys.Date(),"_",data.file, ".csv",sep=""), row.names=F)
    write.csv(raw_files2,file=paste(data.dir,"GUI files/NoRemoval/behavior_5day_",assay, "_NoRemovals_",Sys.Date(),"_",data.file, ".csv",sep=""), row.names=F)
    write.csv(behav5d,file=paste(data.dir,"GUI files/RemovedAffected/behavior_5d_FOR_GUI_", assay, "_",Sys.Date(),"_",data.file, ".csv", sep=""), row.names=F)
    write.csv(behav5d.unfiltered,file=paste(data.dir,"GUI files/NoRemoval/behavior_5d_FOR_GUI_", assay, "_noRemoval_",Sys.Date(),"_",data.file, ".csv", sep=""), row.names=F)
    write.csv(behav5d.date,file=paste(data.dir,"GUI files/NoRemoval/behavior_5d_FOR_GUI_", assay, "_Per_Day_noRemoval_",Sys.Date(),"_",data.file, ".csv", sep=""), row.names=F)

    ##summary by Chemical and conc
    temp2 <- raw_files %>% group_by(chemical.id, chemical.name, conc, starttime) %>%
      summarise(mean_tot = mean(total, na.rm=T),
                sd=sd(total, na.rm=T),
                se = sd(total,na.rm=T)/sqrt(length(total)),
                affected = length(which(is.na(total))),
                count=length(which(!is.na(total) & total != -1)),
                org_ss = length(which(!is.na(total) & total != -1 | is.na(total))))

    temp2$conc <- as.numeric(as.character(temp2$conc))
    temp2$conc <- factor(temp2$conc, levels = sort(unique(as.numeric(as.character(temp2$conc)))))

    ##summary by plate
    tempx <- raw_files %>% group_by(chemical.id, vp_locationcode, chemical.name, plate.id, conc, starttime) %>%
      summarise(mean_tot = mean(total, na.rm=T),
                sd=sd(total, na.rm=T),
                se = sd(total,na.rm=T)/sqrt(length(total)),
                affected = length(which(is.na(total))),
                count=length(which(!is.na(total) & total != -1)),
                org_ss = length(which(!is.na(total) & total != -1 | is.na(total))))

    # tempx$conc <- as.numeric(as.character(tempx$conc))
    tempx$conc <- factor(tempx$conc, levels = sort(unique(as.numeric(as.character(tempx$conc)))))
    tempx$plate_Location <- paste0(tempx$vp_locationcode, "-", tempx$plate.id)


    ##processing for only LPR
    if(assay == "LPR"){
      temp2 <- temp2[which(temp2$starttime %in% c(366:1434)),]
      tempx <- tempx[which(tempx$starttime %in% c(366:1434)),]

      temp_qc <- ddply(temp2, c('chemical.id', "chemical.name", 'conc'), summarise,
                       count=max(count),
                       org_ss=max(org_ss))
      temp_qc$perc <- round(temp_qc$count/temp_qc$org_ss,2)
      temp_qc$ss_more_70perc <- as.numeric(as.character(temp_qc$perc)) >=.3
      temp_qc$chemical.id <- as.character(temp_qc$chemical.id)
      # temp_qc$conc <- as.character(as.numeric(temp_qc$conc))

      ##now remove concentrations that don't have enough samples
      temp_qc_fail <- temp_qc[which(temp_qc$ss_more_70perc == FALSE),c('chemical.id', 'conc')]

      if(nrow(temp_qc_fail) >= 1 ){
        for(i in 1:nrow(temp_qc_fail)){
          temp2 <- temp2[-which(temp2$chemical.id == temp_qc_fail[i,1] &
                                  temp2$conc == temp_qc_fail[i,2]),]
          tempx <- tempx[-which(tempx$chemical.id == temp_qc_fail[i,1] &
                                  tempx$conc == temp_qc_fail[i,2]),]
        }
      } ##removes concs with not enough samples

      ##collapse down the results with affect counts
      temp2_collapse <- temp2 %>%
        group_by(chemical.id, chemical.name, conc) %>%
        summarise(affected = max(affected),
                  count= max(count),
                  initalcount = max(org_ss))


      #############################################
      ##stats according to Zhang et al
      #############################################
      ##look at controls only
      behav5d <- read.csv("gf-vs-cvz-wChorions-10day-01and02-CVZ/GUI files/RemovedAffected/behavior_5d_FOR_GUI_LPR_2019-04-18_Single_Tab_File  04-18-2019 141044.csv")
      behav5d$conc <- as.numeric(as.character(behav5d$conc))
      behav5d.c <- behav5d %>% filter(conc == 0)

      quant <- function(x){quantile(x, c(0.05,0.95), na.rm=T)}

      ###lets analyze only the 2nd epoch
      l1 <- paste0("t", seq(61,89))
      d1 <- paste0("t", seq(90,119))
      lstart <- 61
      lend <- 89
      dstart <- 90
      dend <- 119
      ep2 <- paste0("t", seq(lstart,dend))

      qdiff <- function(x, probs = c(0.95), na.rm = TRUE) {
        the_quantiles <- quantile(x = x, probs = probs, na.rm = na.rm)
        return(log(the_quantiles) - min(0))
      }

      ###use guozhu's code
      upper = 0.95
      lower = 0.05
      # test5d <- behav5d[which(behav5d$chemical.id == "C09375"),]
      DE <- function(x) log(quantile(x,upper,na.rm=TRUE)-quantile(x,lower,na.rm=TRUE)+1)


      ##calculate h vals per plate
      hcal <- melt(behav5d, id=c("chemical.id", "well", "conc", "plate"))

      hcal <- hcal  %>% group_by(chemical.id, conc, plate, variable) %>%
        summarise(h = DE(value))

      # hcal[which(hcal$h <= 0),'h'] <- 0

      hcal.wide <- dcast(hcal,chemical.id + conc + plate ~ variable, value.var="h")

      phases <- c("light", "dark")
      ci <- unique(hcal$chemical.id)[!is.na(unique(hcal$chemical.id))]
      cc <- sort(unique(as.numeric(as.character(hcal$conc)))) ##concs


      ##########LETS DO SOME STATS
      AUC.df <- data.frame()

      headers <- c("chemical.id", "conc",
                   "Light_Cycle_Pval", "Dark_Cycle_Pval",
                   "Light_Cycle_AUC","Dark_Cycle_AUC",
                   # "Light_AUC_Ratio", "Dark_AUC_Ratio",
                   "Light_AUC_RelativeRatio",
                   "Dark_AUC_RelativeRatio")
      tab2 <- data.frame()

      for(chi in ci){
        hdata <- hcal %>% filter(chemical.id == chi)
        ctrl <- hdata %>% filter(conc == 0)
        # fails <- temp_qc_fail[which(temp_qc_fail$Chemical.ID == chi),'conc']
        cc <- sort(unique(as.numeric(as.character(hdata$conc)))) ##concs
        some1 <- matrix(nrow=length(cc), ncol=length(headers), dimnames=list(cc, headers))

        for(c in cc){
          hdata1 <- hdata %>% filter(conc==c)
          # plates <- unique(hdata1$plate)

          ##setup things
          some1[,1] <- chi
          some1[paste(c),2] <- paste(c)

          ##check to see if there is enough good fish to do this snalysis
          if(morphsum[which(morphsum$chemical.id == chi & morphsum$conc == c),'behavior.include'] == 1){
            ##dark AUC
            ksdark <-suppressWarnings(ks.test(t(ctrl[which(ctrl$variable %in% d1),'h']),
                                              t(hdata1[which(hdata1$variable %in% d1),'h']), alternative="two.sided"))#Kolmogorov-Smirnov Goodness-of-Fit
            some1[paste(c),'Dark_Cycle_Pval'] <- round(ksdark$p.value, 3)

            ##light AUC
            kslight <-suppressWarnings(ks.test(t(ctrl[which(ctrl$variable %in% l1),'h']),
                                               t(hdata1[which(hdata1$variable %in% l1),'h']), alternative="two.sided"))#Kolmogorov-Smirnov Goodness-of-Fit
            some1[paste(c),'Light_Cycle_Pval'] <- round(kslight$p.value, 3)

            for(phase in phases){
              if(phase == "light"){
                h1 = paste0("t", lstart)
                hend = paste0("t", lend)
                between = paste0("t", seq(lstart+1, lend-1))
              } else{
                h1 = paste0("t", dstart)
                hend = paste0("t", dend)
                between = paste0("t", seq(dstart+1, dend-1))
              }
              hdat <- hdata1 %>% filter(variable %in% c(h1,hend,between))


              # for(pl in plates){
              #   hdat1 <- hdat %>% filter(plate == pl)

              hdat$variable <- gsub("[a-zA-Z]", "\\2", as.character(hdat$variable))
              auc = trapz(as.numeric(hdat$variable), hdat$h)

              hcalc <- data.frame(chi, c,  auc, phase)
              AUC.df <- rbind(AUC.df, hcalc)
              if(phase == "dark"){
                some1[paste(c),'Dark_Cycle_AUC'] <- round(mean(AUC.df$auc),3)
              }
              if(phase == "light"){
                some1[paste(c),'Light_Cycle_AUC'] <- round(mean(AUC.df$auc),3)
              }
              # } ##plates
            }##light phases

          }else {
            some1[paste(c),2] <- paste(c)
            some1[paste(c),c("Light_Cycle_Pval", "Dark_Cycle_Pval",
                             "Light_Cycle_AUC", "Dark_Cycle_AUC",
                             # "Light_AUC_Ratio", "Dark_AUC_Ratio",
                             "Light_AUC_RelativeRatio", "Dark_AUC_RelativeRatio")] <- NA
            auc <- NA
            for(phase in phases){
              hcalc <- data.frame(chi, c, auc, phase)
              AUC.df <- rbind(AUC.df, hcalc)
            }
          } ##if there isnt enough data points
        }## concs

        AUC.ch <- AUC.df[which(AUC.df$chi == chi),]
        c.c <- unique(AUC.ch$c)

        # AUC.ch1 <- AUC.ch %>% filter(c == co)
        AUC.ch1.wide <- dcast(AUC.ch, chi + c  ~ phase, value.var="auc")

        for(co in c.c){
          ##conc
          sub1 <- AUC.ch1.wide[which(AUC.ch1.wide$c == co),]
          Lratio <- sub1[1,'light']
          Dratio <- sub1[1,'dark']

          ##control
          AUC.ctrl <- AUC.ch1.wide[which(AUC.ch1.wide$c == 0),]
          Lratio.c <- AUC.ctrl[1,'light']
          Dratio.c <- AUC.ctrl[1,'dark']

          ##relative ratio
          LRratio <- (Lratio- Lratio.c) / Lratio.c
          DRratio<- (Dratio-Dratio.c) / Dratio.c

          # some1[paste(co), 'Light_AUC_Ratio'] <- round(Lratio,3)
          # some1[paste(co), 'Dark_AUC_Ratio'] <- round(Dratio,3)

          some1[paste(co), 'Light_AUC_RelativeRatio'] <- round(LRratio,3)
          some1[paste(co), 'Dark_AUC_RelativeRatio'] <- round(DRratio,3)
        }
        tab2<- rbind(tab2, some1)
      }


      tab_melt2 <- melt(tab2, id=c("chemical.id", "conc"))
      tab_melt2 <- tab_melt2 %>% separate(variable, into = c("interval", "cycle", "readout"), sep="_")
      tab_melt2 <- tab_melt2 %>% select("chemical.id", "conc", "interval", "readout", "value")
      tab_melt2 <- tab_melt2 %>% spread(readout, value)

      stats_table2 <- merge(temp_qc, tab_melt2, by=c("chemical.id", "conc"))
      stats_table2[,c('conc', 'AUC', 'Pval', 'RelativeRatio')] = apply(stats_table2[,c('conc','AUC', 'Pval', 'RelativeRatio')],2, function(x) as.numeric(as.character(x)))
      stats_table2 <- stats_table2[order(stats_table2$chemical.id, stats_table2$interval, stats_table2$conc ),]
      stats_table2[which(stats_table2$Pval <= 0.01 & stats_table2$ss_more_70perc == TRUE),'significance p<0.01'] <- "YES"
      stats_table2[which(stats_table2$`significance p<0.01` == "YES" & stats_table2$RelativeRatio >=.1 ),'activity'] <- "HYPER"
      stats_table2[which(stats_table2$`significance p<0.01` == "YES" & stats_table2$RelativeRatio <= -.1),'activity'] <- "HYPO"
      stats_table2[which(!is.na(stats_table2$activity)),'Sig'] <- "YES"
      write.csv(stats_table2, file=paste0(data.dir,  "/", assay, "/5dpf_Behavior_Light_Dark_Stats_Summary_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)


      ###GRAPHS
      hcal2 <- hcal %>% filter(variable %in% ep2)
      hcal2 <- ddply(hcal2, c("chemical.id", "conc", "variable"), summarise, ave.h = mean(h))
      hcal2$time <- gsub("[a-zA-Z]", "\\2", hcal2$variable)
      hcal2$time <- as.numeric(as.character(hcal2$time))
      hcal2$conc <- as.factor(hcal2$conc)


      ###calculate average per min
      hcal3 <- hcal2
      hcal3[which(hcal3$time %in% seq(61,65)),'sec'] <- 0
      hcal3[which(hcal3$time %in% seq(66,70)),'sec'] <- 30
      hcal3[which(hcal3$time %in% seq(71,75)),'sec'] <- 60
      hcal3[which(hcal3$time %in% seq(76,80)),'sec'] <- 90
      hcal3[which(hcal3$time %in% seq(81,85)),'sec'] <- 120
      hcal3[which(hcal3$time %in% seq(86,90)),'sec'] <- 150
      hcal3[which(hcal3$time %in% seq(91,95)),'sec'] <- 180
      hcal3[which(hcal3$time %in% seq(96,100)),'sec'] <- 210
      hcal3[which(hcal3$time %in% seq(101,105)),'sec'] <- 240
      hcal3[which(hcal3$time %in% seq(106,110)),'sec'] <- 270
      hcal3[which(hcal3$time %in% seq(111,115)),'sec'] <- 300
      hcal3[which(hcal3$time %in% seq(116,119)),'sec'] <- 330

      hcal4 <- ddply(hcal3, c("chemical.id", "conc", "sec"), summarise, h.sec = sum(ave.h))
      # hcal4$conc <- as.numeric(as.character(hcal4$conc))


      # print(ggplot(hcal4[which(hcal4$chemical.id == "C09375"),], aes(sec, h.sec, colour=conc))+
      #         geom_path(lwd=1.2) +
      #         # facet_grid(.~ chemical.id)+
      #         scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
      #                                     "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
      #                                     "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
      #                                     "darkblue", "darkslateblue", "darkseagreen1"))+
      #         labs(x="Time", y= "Differential Entropy")
      # )

      #####################################################################
      #############calculate AUC
      sub_chemids <- unique(temp2$chemical.id)[!is.na(unique(temp2$chemical.id))]
      control_conc <- 0


      dark_t <- c(paste0("t", seq(90,119)),paste0("t", seq(150,179)),paste0("t", seq(210,239)))
      # d1 <- paste0("t", seq(90,119))
      # d2 <- paste0("t", seq(150,179))
      # d3 <- paste0("t", seq(210,239))

      light_t <- c(paste0("t", seq(61,89)),paste0("t", seq(120,149)),paste0("t", seq(180,209)))
      # l1 <- paste0("t", seq(61,89))
      # l2 <-paste0("t", seq(120,149))
      # l3 <- paste0("t", seq(180,209))

      all_t <- paste0("t", seq(61,239))
      t_trans <- "t90"

      ## Stats
      tab <- data.frame()
      for(vals in sub_chemids){
        headers <- c("chemical.id", "conc", "Light_Cycle_Pval", "Dark_Cycle_Pval",
                     "Light_Cycle_AUC", "Dark_Cycle_AUC",
                     "Light_Cycle_PercDiff", "Dark_Cycle_PercDiff",
                     "trans_dist", "trans_PercDiff", "trans_AOV_pval")
        dat <- behav5d[which(behav5d$chemical.id == vals),]
        dat$conc <- as.numeric(as.character(dat$conc))

        c1 <- dat[which(dat$conc == control_conc),]
        c1.clean <- c1[complete.cases(c1),]
        c1_all_aucs <- apply(c1.clean[,match(c(dark_t,light_t), names(c1.clean))], 2, trapz)
        c1_dark_aucs <- apply(c1.clean[,match(dark_t, names(c1.clean))], 2, trapz)

        # d1_cauc  <- apply(c1.clean[,match(d1, names(c1.clean))], 2, trapz) ##columns
        # d2_cauc <- apply(c1.clean[,match(d2, names(c1.clean))], 2, trapz)
        # d3_cauc <- apply(c1.clean[,match(d3, names(c1.clean))], 2, trapz)
        # d1_caucr  <- apply(c1.clean[,match(d1, names(c1.clean))], 1, trapz) ##every fish
        # d2_caucr  <- apply(c1.clean[,match(d2, names(c1.clean))], 1, trapz) ##every fish
        # d3_caucr  <- apply(c1.clean[,match(d3, names(c1.clean))], 1, trapz) ##every fish

        # l1_cauc  <- apply(c1.clean[,match(l1, names(c1.clean))], 2, trapz) ##columns
        # l2_cauc  <- apply(c1.clean[,match(l2, names(c1.clean))], 2, trapz) ##columns
        # l3_cauc  <- apply(c1.clean[,match(l3, names(c1.clean))], 2, trapz) ##columns
        # l1_caucr  <- apply(c1.clean[,match(l1, names(c1.clean))], 1, trapz) ##every fish
        # l2_caucr  <- apply(c1.clean[,match(l2, names(c1.clean))], 1, trapz) ##every fish
        # l3_caucr  <- apply(c1.clean[,match(l3, names(c1.clean))], 1, trapz) ##every fish


        c1_light_aucs <- apply(c1.clean[,match(light_t, names(c1.clean))], 2, trapz) ##columnes

        c_transition <- mean(c1.clean[,t_trans])

        ##ANoVA for transition time
        trans_dat <- dat %>% select(chemical.id, conc,well, plate, t_trans)
        concs_av <- trans_dat %>% distinct( conc)
        fails <- temp_qc_fail[which(temp_qc_fail$chemical.id == vals),'conc']

        if(length(fails)>0){
          trans_dat$conc <- as.numeric(as.character(trans_dat$conc))
          trans_dat <- trans_dat[-which(trans_dat$conc %in% fails),]
        }

        if(length(unique(trans_dat$conc)) >1){
          trans_dat$conc <- as.factor(trans_dat$conc)
          trans_stats <- data.frame(TukeyHSD(aov(t90 ~ conc, trans_dat))$conc)
          trans_stats$Comparison <- row.names(trans_stats)
          trans_stats <- trans_stats[grep("0", trans_stats$Comparison),]
          trans_stats$conc <- unlist(lapply(strsplit(as.character(trans_stats$Comparison), "\\-"), "[", 1))
          trans_stats$conc <- as.numeric(as.character(trans_stats$conc))
          trans_stats <- rbind(rep(0,6), trans_stats)
          trans_stats <- trans_stats[1:length(unique(dat$conc)),]
          trans_stats$conc <- unique(dat$conc)
          trans_stats[is.na(trans_stats)] <- 0
          trans_stats <- data.frame(trans_stats)
        } else{
          trans_stats <- matrix(nrow=length(unique(dat$conc)),ncol=6,
                                dimnames=list(c(),c("diff", "lwr", "upr", "p.adj", "Comparison", "conc")), 0)
          trans_stats[,6] <-  unique(dat$conc)
          trans_stats <- data.frame(trans_stats)
        }



        ###look at the epochs
        cc <- sort(unique(as.numeric(as.character(dat$conc)))) ##concs
        some <- matrix(nrow=length(cc), ncol=length(headers), dimnames=list(cc, headers))
        for(c in cc){
          t1 <- dat[which(dat$conc == c),]
          t <- t1[complete.cases(t1),]
          if(nrow(t) > 1){

            # ###calculate AUC per minute per light or dark period
            all_aucs <- apply(t[,c(61:239)], 2, trapz)
            dark_aucs <- apply(t[,match(dark_t, names(t))], 2, trapz)
            light_aucs <- apply(t[,match(light_t, names(t))], 2, trapz)

            ##setup things
            some[,1] <- vals
            some[paste(c),2] <- paste(c)

            ##dark AUC
            kdark <-suppressWarnings(ks.test(c1_dark_aucs,  dark_aucs, alternative="two.sided"))#Kolmogorov-Smirnov Goodness-of-Fit
            some[paste(c),'Dark_Cycle_Pval'] <- round(kdark$p.value, 3)
            some[paste(c),'Dark_Cycle_AUC'] <- round(sum(dark_aucs),3)

            ##light AUC
            klight <-ks.test(c1_light_aucs, light_aucs, alternative="two.sided")#Kolmogorov-Smirnov Goodness-of-Fit
            some[paste(c),'Light_Cycle_Pval'] <- round(klight$p.value, 3)
            some[paste(c),'Light_Cycle_AUC'] <- round(sum(light_aucs),3)

            ##calculate % difference
            diff_dark <- (sum(dark_aucs) - sum(c1_dark_aucs)) / (sum(c1_dark_aucs))*100
            diff_light <- (sum(light_aucs) - sum(c1_light_aucs)) / (sum(c1_light_aucs))*100
            some[paste(c),'Dark_Cycle_PercDiff'] <- round(diff_dark,3)
            some[paste(c),'Light_Cycle_PercDiff'] <- round(diff_light,3)

            ##calculate light to dark transition - 6 sec after light on
            treat_transition <- mean(t[,t_trans])
            some[paste(c),'trans_dist'] <- round(treat_transition,3)
            some[paste(c),'trans_PercDiff'] <- round((treat_transition - c_transition)/c_transition *100,3)
            some[paste(c), 'trans_AOV_pval'] <- round(trans_stats[which(trans_stats$conc == c),'p.adj'],3)

          } else{
            some[,1] <- vals
            some[paste(c),2] <- paste(c)
            some[paste(c),c("Light_Cycle_Pval", "Dark_Cycle_Pval",
                            "Light_Cycle_AUC", "Dark_Cycle_AUC", "Light_Cycle_PercDiff", "Dark_Cycle_PercDiff",
                            "trans_dist", "trans_PercDiff","trans_dist", "trans_PercDiff", "trans_AOV_pval")] <- NA
          }
        } ##for all the concentrations
        tab <- rbind(tab, some)
      } ##runs through each chemical

      #####################################
      ####new stats 9/20/17
      #####################################
      tab$chemical.id <- as.character(tab$chemical.id)
      tab_melt <- melt(tab[-c(9:11)], id=c("chemical.id", "conc"))
      tab_melt <- tab_melt %>% separate(variable, into = c("interval", "cycle", "readout"), sep="_")
      tab_melt <- tab_melt %>% spread(readout, value)

      stats_table <- merge(temp_qc, tab_melt, by=c("chemical.id", "conc"))
      stats_table[,c('conc', 'AUC', 'Pval','PercDiff')] = apply(stats_table[,c('conc','AUC', 'Pval','PercDiff')],2, function(x) as.numeric(as.character(x)))
      stats_table <- stats_table[order(stats_table$chemical.id, stats_table$interval, stats_table$conc ),]
      stats_table[which(stats_table$Pval <= 0.01 & stats_table$ss_more_70perc == TRUE),'significance p<0.01'] <- "YES"
      stats_table[which(stats_table$`significance p<0.01` == "YES" & stats_table$PercDiff >= 50 ),'significance p<0.01 and 50% cutoff'] <- "YES"
      stats_table[which(stats_table$`significance p<0.01` == "YES" & stats_table$PercDiff <= -50),'significance p<0.01 and 50% cutoff'] <- "YES"
      stats_table[which(stats_table$Pval <= 0.01 & stats_table$ss_more_70perc == TRUE & stats_table$PercDiff >50 ),'activity'] <- "HYPER"
      stats_table[which(stats_table$Pval <= 0.01 & stats_table$ss_more_70perc == TRUE &  stats_table$PercDiff < -50),'activity'] <- "HYPO"

      ##add in conc units:
      stats_table <- merge(stats_table, conc_table[c(1,3,5,6)])
      stats_table <- stats_table[c("chemical.id", "chemical.name", "conc", "units", "interval",
                                   "count", "org_ss", "perc", "ss_more_70perc", "AUC", "PercDiff",
                                   "Pval", "significance p<0.01", "significance p<0.01 and 50% cutoff","activity")]

      ##add comments why there is NA
      stats_table[is.na(stats_table$interval), 'Comments'] <- "Not enough animals for stats test"
      write.csv(stats_table, file=paste0(data.dir,  "/", assay, "/5dpf_Behavior_Light_Dark_AUCmethod_Stats_Summary_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)

      ###########stats table for the transition
      stats_table_trans <- tab %>% select(chemical.id, conc, trans_dist, trans_PercDiff, trans_AOV_pval) %>%
        inner_join(temp_qc, by=c("chemical.id", "conc")) %>%
        select(chemical.id, chemical.name, conc, count, org_ss, ss_more_70perc, trans_dist,trans_PercDiff, trans_AOV_pval)
      stats_table_trans[,c('trans_dist', 'trans_PercDiff', 'trans_AOV_pval')] = apply(stats_table_trans[,c('trans_dist', 'trans_PercDiff', 'trans_AOV_pval')],2, function(x) as.numeric(as.character(x)))
      stats_table_trans[which(stats_table_trans$trans_AOV_pval <= 0.01 & stats_table_trans$ss_more_70perc == TRUE & stats_table_trans$conc != 0),'significance p<0.01'] <- "YES"

      #####################################################################################
      ##### Create report quality table for VP
      #####################################################################################
      if(length(which(stats_table2$`significance p<0.01` == "YES"))>0){
        vp_tab <- stats_table2 %>% filter(`significance p<0.01` == "YES") %>%
          select(chemical.name, conc, interval, count, AUC, RelativeRatio, Pval,activity)
        # rename(chemical.name = Chemical, interval = Interval, count = N, activity = Activity)
        # vp_tab <- stats_table2[which(stats_table2$activity %in% c("HYPER", "HYPO")),]
        vp_tab[vp_tab$Pval < 0.01 & vp_tab$Pval>0.01,'Pval'] <- "<0.01"
        vp_tab[vp_tab$Pval <= 0.01 & vp_tab$Pval>0.001,'Pval'] <- "<0.001"
        vp_tab[vp_tab$Pval < 0.001 & vp_tab$Pval>0, 'Pval'] <- "<0.0001"
        colnames(vp_tab) <- c("Chemical", "conc", "Interval", "N","AUC", "RelativeRatio", "P value", "Activity")
        vp_tab <- vp_tab[order(vp_tab$Chemical, vp_tab$conc),]
        write.csv(vp_tab,paste0(data.dir,"Reporting_Material/VP_", assay, "_Significance_Table_Summary_New_LPR_System_", data.file, ".csv"),quote=T,row.names=F)
      }

      if(length(which(stats_table_trans$`significance p<0.01` == "YES"))>0){
        vp_tab_trans <- stats_table_trans[which(stats_table_trans$`significance p<0.01` == "YES"),-c(5,6)]
        vp_tab_trans <- vp_tab_trans[-c(1,8)]
        colnames(vp_tab_trans) <- c("Chemical", "conc", "N", "Average Distance (mm)", "Percent change", "P value")
        vp_tab_trans <- vp_tab_trans[order(vp_tab_trans$Chemical, vp_tab_trans$conc),]
        write.csv(vp_tab_trans,paste0(data.dir,"Reporting_Material/VP_", assay, "_Transition_Significance_Table_Summary_New_LPR_System_", data.file, ".csv"),quote=T,row.names=F)
      }

      #####################################################################################
      ##### Compute LEL tables
      #####################################################################################
      vp_AUC <- data.frame(c(sapply(unique(stats_table$chemical.id),function(x) rep(x,length(unique(stats_table$interval))+1))),
                           rep(NA,length(unique(stats_table$chemical.id))*length(c(unique(stats_table$interval), "Transition"))),
                           rep(c(unique(stats_table$interval), "Transition"),length(unique(stats_table$chemical.id))),
                           stringsAsFactors=F)

      names(vp_AUC)<-c("chemical.id","conc","interval")
      vp_AUC<-merge(vp_AUC,unique(chem.key[,c("chemical.id","chemical.name")]))


      for (i in 1:nrow(vp_AUC)) { ##this adds in the LEL
        ##light dark
        vp_AUC.grp<-stats_table2$conc[stats_table2$chemical.id==vp_AUC$chemical.id[i] &
                                        stats_table2$interval==vp_AUC$interval[i] & stats_table2$Sig == "YES"]
        vp_AUC.grp<-vp_AUC.grp[!is.na(vp_AUC.grp)]
        if (length(vp_AUC.grp)>0) {vp_AUC$conc[i]<-vp_AUC.grp[1]}
      }

      ##fill in transition time
      tab$trans_AOV_pval <- as.numeric(as.character(tab$trans_AOV_pval))
      trans_df <- tab[c(1,2,11)]
      trans_df <- trans_df[which(trans_df$trans_AOV_pval < 0.01 & trans_df$conc != 0),]
      trans_df$conc <- as.numeric(as.character(trans_df$conc))
      trans_df <- trans_df %>% group_by(chemical.id) %>% summarise(conc = min(conc))
      for(id in unique(trans_df$chemical.id)){
        vp_AUC[which(vp_AUC$chemical.id == id & vp_AUC$interval == "Transition"),'conc'] <- trans_df[which(trans_df$chemical.id == id), 'conc']
      }


      ##summarize # of "good" concs
      good.conc <- ddply(stats_table2, c("chemical.id", "chemical.name", "interval"), summarise,
                         total.no.conc = length(unique(conc)),
                         no.conc.used.for.analysis = length(which(ss_more_70perc == "TRUE")),
                         no.conc.w.sig.hits = length(which(Sig == "YES")))

      vp_AUC <- vp_AUC %>% rename(LEL = conc)
      print(head(vp_AUC))

      good.conc <- good.conc %>% join(vp_AUC[-which(vp_AUC$interval == "Transition"),])
      # LEL.5d.AUC <- vp_AUC[-which(vp_AUC$Interval == "Transition"),] %>% rename(interval = Interval) %>%
      # join(good.conc) %>% spread(interval, CONC) %>% rename(LEL.dark = Dark, LEL.light = Light)



      # LEL.5d.AUC <- good.conc %>% mutate(id=1:n()) %>%
      #   spread(interval, LEL) %>% rename(LEL.dark = Dark, LEL.light = Light)
      #
      LEL.5d.AUC <- vp_AUC[-which(vp_AUC$interval == "Transition"),] %>%
        join(good.conc) %>%
        spread(interval, LEL) %>% rename(LEL.dark = Dark, LEL.light = Light)
      #
      LEL.5d.AUC[is.na(LEL.5d.AUC)]<-0

      ###determine if a chem is a hit
      LEL.5d.AUC[which(LEL.5d.AUC$LEL.dark > 0 | LEL.5d.AUC$LEL.light > 0),'LPR.hit'] <- "YES"


      #####write out the csv files
      # write.csv(stats_table, file=paste0(data.dir,  "/", assay, "/5dpf_Behavior_Light_Dark_Stats_Summary_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
      write.csv(tab, file=paste0(data.dir,  "/", assay, "/5dpf_Behavior_raw_stats_table_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
      write.csv(stats_table_trans, file=paste0(data.dir,  "/", assay, "/5dpf_Behavior_Transition_Stats_Summary_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
      write.csv(LEL.5d.AUC, file=paste0(data.dir, "/", assay, "/5dpf_Behavior_LEL_Table_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
      write.csv(good.conc, file=paste0(data.dir, "/", assay, "/5dpf_Behavior_LEL_and_metadata_Table_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)

    }##close off LPR


    ################################
    ## DO SOME WORK FOR LSR
    ################################
    if(assay == "LSR"){
      dir.create(file.path(paste(data.dir, "LSR/", sep="")))
      SR <- paste0("t", seq(0,9))

      ## Stats
      tab <- data.frame()
      for(vals in sub_chemids){
        headers <- c("chemical.id", "conc", "AUC.Pval", "AUC", "PercDiff", "Peak.no", "Peak_dist", "Peak.pval", "PeakDiff")

        dat <- behav5d %>% filter(chemical.id == vals) %>% select(chemical.id, conc, plate, well,SR)
        dat$conc <- as.numeric(as.character(dat$conc))
        cc <- as.numeric(as.character(unique(dat$conc)))
        c1 <- dat[which(dat$conc == control_conc),]
        c1 <- c1[complete.cases(c1),]
        c1_aucs <- apply(c1[,match(c(SR), names(c1))], 2, trapz)
        c_mean <- apply(c1[,match(c(SR), names(c1))], 2, mean)

        SR.peak <- names(grep(max(c_mean), c_mean, value=T)) ##what time point is the peak
        c1_peak <- get(SR.peak, c1)


        some <- matrix(nrow=length(cc), ncol=length(headers), dimnames=list(cc, headers))
        for(c in cc){
          t <- dat[which(dat$conc == c),]
          t <- t[complete.cases(t),]

          if(nrow(t) > 1){
            t_aucs <- apply(t[,match(SR, names(t))], 2, trapz)
            t_peak <- get(SR.peak, t)

            ##setup things
            some[,1] <- vals
            some[paste(c),2] <- paste(c)

            ##AUC
            SR.auc <-ks.test(c1_aucs, t_aucs, alternative="two.sided")#Kolmogorov-Smirnov Goodness-of-Fit
            some[paste(c),'AUC.Pval'] <- round(SR.auc$p.value, 3)
            some[paste(c),'AUC'] <- round(mean(t_aucs),3)
            some[paste(c),'PercDiff'] <- round((mean(t_aucs) - mean(c1_aucs)) / (mean(c1_aucs))*100,3)

            ##looking at the peak - 3 secs in
            SR.Peak.compare <- ks.test(c1_peak, t_peak, alternative = "two.sided")
            some[paste(c),'Peak.no'] <- as.numeric(sapply(strsplit(SR.peak, "t"), "[[", 2))+1
            some[paste(c),'Peak_dist'] <- round(mean(t_peak),3)
            some[paste(c),'Peak.pval'] <- round(SR.Peak.compare$p.value, 3)
            some[paste(c),'PeakDiff'] <- round((mean(t_peak) - mean(c1_peak)) / (mean(c1_peak))*100,3)
          } else{
            some[,1] <- vals
            some[paste(c),2] <- paste(c)
            some[paste(c),'Peak.no'] <- as.numeric(sapply(strsplit(SR.peak, "t"), "[[", 2))
            some[paste(c),c("AUC.Pval", "AUC", "PercDiff", "Peak_dist", "Peak.pval", "PeakDiff")] <- NA
          }
        } ## for each conc
        tab <- rbind(tab, some)
      }## for each chemical ID

      temp2 <- temp2[which(temp2$starttime %in% c(1467:1476)),]
      tempx <- tempx[which(tempx$starttime %in% c(1467:1476)),]

      #####################################
      ####new stats 9/20/17
      #####################################
      stats_table_LSR <- merge(temp_qc, tab, by=c("chemical.id", "conc"),all.x=T)

      ##convert the 3 columns into numeric
      stats_table_LSR[,c('conc', 'AUC', 'AUC.Pval','PercDiff')] = apply(stats_table_LSR[,c('conc','AUC', 'AUC.Pval','PercDiff')],2, function(x) as.numeric(as.character(x)))
      stats_table_LSR[which(stats_table_LSR$AUC.Pval <= 0.01 & stats_table_LSR$ss_more_70perc == TRUE),'significance p<0.01'] <- "YES"
      stats_table_LSR[which(stats_table_LSR$AUC.Pval <= 0.01 & stats_table_LSR$ss_more_70perc == TRUE & stats_table_LSR$PercDiff >50),'significance p<0.01 and 50% cutoff'] <- "YES"
      stats_table_LSR[which(stats_table_LSR$AUC.Pval <= 0.01 & stats_table_LSR$ss_more_70perc == TRUE & stats_table_LSR$PercDiff < -50),'significance p<0.01 and 50% cutoff'] <- "YES"
      stats_table_LSR[which(stats_table_LSR$AUC.Pval <= 0.01 & stats_table_LSR$ss_more_70perc == TRUE & stats_table_LSR$PercDiff >50 ),'activity'] <- "HYPER"
      stats_table_LSR[which(stats_table_LSR$PAUC.Pvalval <= 0.01 & stats_table_LSR$ss_more_70perc == TRUE &  stats_table_LSR$PercDiff < -50),'activity'] <- "HYPO"

      ##add in conc units:
      stats_table_LSR <- merge(stats_table_LSR, conc_table[c(1,3,5,6)])
      stats_table_LSR <- stats_table_LSR[c("chemical.id", "chemical.name", "conc", "units",
                                           "count", "org_ss", "perc", "ss_more_70perc", "AUC", "PercDiff",
                                           "AUC.Pval", "significance p<0.01", "significance p<0.01 and 50% cutoff","activity")]

      ##add comments why there is NA
      stats_table_LSR[is.na(stats_table_LSR$AUC), 'Comments'] <- "Not enough animals for stats test"

      ###########stats table for the transition
      stats_table_LSR_peak <- tab %>% select(chemical.id, conc, Peak.no, Peak_dist, Peak.pval, PeakDiff) %>%
        inner_join(temp_qc, by=c("chemical.id", "conc")) %>%
        select(chemical.id, chemical.name, conc, count, org_ss, ss_more_70perc, Peak.no, Peak_dist,Peak.pval, PeakDiff)
      stats_table_LSR_peak$Peak.pval <- as.numeric(as.character(stats_table_LSR_peak$Peak.pval))
      stats_table_LSR_peak[,c('Peak_dist', 'Peak.pval', 'PeakDiff')] = apply(stats_table_LSR_peak[,c('Peak_dist', 'Peak.pval', 'PeakDiff')],2, function(x) as.numeric(as.character(x)))
      stats_table_LSR_peak[which(stats_table_LSR_peak$Peak.pval < 0.05 & stats_table_LSR_peak$ss_more_70perc == TRUE & stats_table_LSR_peak$conc != 0),'significance p<0.05'] <- "YES"
      stats_table_LSR_peak$conc <- as.numeric(as.character(stats_table_LSR_peak$conc))

      #####################################################################################
      ##### Create report quality table for VP
      #####################################################################################
      if(length(which(stats_table_LSR$`significance p<0.01` == "YES"))>0){
        vp_LSR_tab <- stats_table_LSR[which(stats_table_LSR$`significance p<0.01 and 50% cutoff` == "YES"),c(2:5, 9:11, 14)]
        vp_LSR_tab[vp_LSR_tab$AUC.Pval < 0.05 & vp_LSR_tab$AUC.Pval>0.01,'Pval'] <- "<0.01"
        vp_LSR_tab[vp_LSR_tab$AUC.Pval <= 0.01 & vp_LSR_tab$AUC.Pval>0.001,'Pval'] <- "<0.001"
        vp_LSR_tab[vp_LSR_tab$AUC.Pval < 0.001 & vp_LSR_tab$AUC.Pval>0, 'Pval'] <- "<0.0001"
        colnames(vp_LSR_tab) <- c("Chemical", "conc", "conc.units", "N", "AUC", "Percent change", "P value", "Activity")
        write.csv(vp_LSR_tab,paste0(data.dir,"Reporting_Material/VP_", assay, "_Significance_Table_Summary_New_LPR_System_", data.file, ".csv"),quote=T,row.names=F)
      }

      if(length(which(stats_table_LSR_peak$`significance p<0.01` == "YES"))>0){
        vp_LSR_peak_tab <- stats_table_LSR_peak[which(stats_table_LSR_peak$`significance p<0.01` == "YES"),c(2:4,8,10,9)]
        vp_LSR_peak_tab[vp_LSR_peak_tab$Peak.pval < 0.01 & vp_LSR_peak_tab$Peak.pval>0.01,'Peak.pval'] <- "<0.01"
        vp_LSR_peak_tab[vp_LSR_peak_tab$Peak.pval <= 0.01 & vp_LSR_peak_tab$Peak.pval>0.001,'Peak.pval'] <- "<0.001"
        vp_LSR_peak_tab[vp_LSR_peak_tab$Peak.pval < 0.001 & vp_LSR_peak_tab$Peak.pval>0, 'Peak.pval'] <- "<0.0001"
        colnames(vp_LSR_peak_tab) <- c("Chemical", "conc", "N", "Peak Distance", "Peak Diff", "P value")
        write.csv(vp_LSR_peak_tab,paste0(data.dir,"Reporting_Material/VP_", assay, "_Peak_Significance_Table_Summary_", data.file, ".csv"),quote=T,row.names=F)
      }

      ###make LEL Table
      LSR_LEL <- data.frame(c(sapply(unique(stats_table$chemical.id),function(x) rep(x,2))),
                            rep(NA,length(unique(stats_table$chemical.id))*2),
                            rep(c("LSR.AUC", "LSR.Peak"),length(unique(stats_table$chemical.id))),
                            stringsAsFactors=F)

      names(LSR_LEL)<-c("chemical.id","conc","Method")
      LSR_LEL<-merge(LSR_LEL,unique(chem.key[,c("chemical.id","chemical.name")]))


      ##fill in the table
      stats_table_LSR$AUC.Pval <- as.numeric(as.character(stats_table_LSR$AUC.Pval))
      LSR_df <- stats_table_LSR[c(1,2,3,11)]
      LSR_df <- LSR_df[which(LSR_df$AUC.Pval < 0.05 & LSR_df$conc != 0),]
      LSR_df$conc <- as.numeric(as.character(LSR_df$conc))
      LSR_df <- LSR_df %>% group_by(chemical.id) %>% summarise(conc = min(conc))

      LSR_peak_df <- stats_table_LSR_peak %>% select(chemical.id, chemical.name,conc, Peak.pval) %>%
        filter(Peak.pval < 0.05) %>% filter(conc != 0) %>% group_by(chemical.id) %>% summarise(conc = min(conc))

      for(id in unique(LSR_df$chemical.id)){
        LSR_LEL[which(LSR_LEL$chemical.id == id & LSR_LEL$Method == "LSR.AUC"),'conc'] <- LSR_df[which(LSR_df$chemical.id == id), 'conc']
      }

      for(id in unique(LSR_peak_df$chemical.id)){
        LSR_LEL[which(LSR_LEL$chemical.id == id & LSR_LEL$Method == "LSR.Peak"),'conc'] <- LSR_peak_df[which(LSR_peak_df$chemical.id == id), 'conc']
      }
      LSR.LEL <- LSR_LEL %>% spread(Method, conc)
      LSR.LEL[is.na(LSR.LEL)]<-0

      ##write out csv
      write.csv(stats_table_LSR_peak, file=paste0(data.dir, "/", assay, "/5dpf_Behavior_Peak_Stats_Summary_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
      write.csv(stats_table_LSR, file=paste0(data.dir, "/", assay, "/5dpf_Behavior_AUC_Stats_Summary_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
      write.csv(LSR.LEL, file=paste0(data.dir, "/", assay, "/5dpf_Behavior_LEL_Table_",assay, "_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)
    } ##this ends startle


    ###################################
    ##plotting for everyone
    ###################################
    limits <- aes(ymax = mean_tot + se, ymin=mean_tot - se)
    pdf(file=paste0(data.dir,"/", assay, "/5dpf_Behavior_plate_View_Plots_",data.file,"_", assay, "_",output.prefix, ".pdf",sep=""), width=16, height=8)

    for(low in sub_chemids){
      run1 <- unique(raw_files[which(raw_files$chemical.id == low),'plate.id'])
      chemical_name <- TX.id.dat[which(TX.id.dat$chemical.id == low),'chemical.name']
      raw_files$conc <- as.factor(raw_files$conc)
      # raw_files$total

      print(ggplot(raw_files[which(raw_files$chemical.id == low),],
                   aes(starttime, total, group = conc))+
              geom_point(lwd=2)+
              geom_line(lwd=2)+
              facet_wrap(~ plate.id)+
              ggtitle(paste0(low, ": ",chemical_name))+
              labs(x="Time (min)", y = "Total Distance Moved (cm)")

      )

      print(ggplot(raw_files[which(raw_files$chemical.id == low),],
                   aes(starttime, total, colour = well, group = conc))+
              geom_point(lwd=2, aes(group = well))+
              geom_line(lwd=2, aes(group = well))+
              # facet_wrap(~well)+
              ggtitle(paste0(low, ": ",chemical_name))+
              labs(x="Time (min)", y = "Total Distance Moved (cm)")
      )



      for(i in run1){
        print(ggplot(raw_files[which(raw_files$plate.id == i & raw_files$chemical.id == low),],
                     aes(starttime, total, colour=conc, group=conc)) +
                geom_point()+
                geom_line()+
                # ylim(0,10)+
                # geom_path()+
                ggtitle(paste0(chemical_name, ": Barcode ", i))+
                facet_wrap(~well, ncol=12)+
                theme(legend.position="bottom")+
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1"))

        )

        print(ggplot(raw_files[which(raw_files$plate.id == i & raw_files$chemical.id == low & raw_files$starttime < 1520),], aes(starttime, total, colour=conc, group=conc)) +
                geom_point()+
                geom_line()+
                ylim(0,10)+
                # geom_path()+
                ggtitle(paste0(chemical_name, ": Barcode ", i))+
                facet_wrap(~well, ncol=12)+
                theme(legend.position="bottom")+
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1",
                                            "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                            "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                            "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                            "tan", "slateblue3", "violetred4"))

        )
      }##per run
    } ##for all ids
    dev.off()
    ##################################
    ####graph the behavior plots
    pdf(file=paste0(data.dir,"/", assay, "/5dpf_Behavior_Plots_",data.file, "_", assay, "_", output.prefix, ".pdf",sep=""), width=16, height=8)
    if(assay == "LPR"){
      print(ggplot(hcal4, aes(sec, h.sec, colour=conc))+
              geom_path(lwd=1.2) +
              facet_grid(chemical.id~.)+
              scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                          "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                          "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                          "darkblue", "darkslateblue", "darkseagreen1",
                                          "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                          "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                          "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                          "tan", "slateblue3", "violetred4"))+
              labs(x="Time", y= "Differential Entropy")
      )
    }

    for(low in sub_chemids){
      temp1 <- temp2[which(temp2$chemical.id == low),]
      temp1.1 <- tempx[which(tempx$chemical.id == low),]

      limits <- aes(ymax = mean_tot + se, ymin=mean_tot - se)
      chemical_name <- TX.id.dat[which(TX.id.dat$chemical.id == low),'chemical.name']
      count <- unique(temp1$count)

      ##new plots
      #### AUC barplots
      if(assay == "LPR"){
        hcal_1 <- hcal4 %>% filter(chemical.id == low)

        print(ggplot(hcal_1, aes(sec, h.sec, colour=conc))+
                geom_path(lwd=1.2) +
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1",
                                            "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                            "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                            "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                            "tan", "slateblue3", "violetred4"))+
                labs(x="Time", y= "Differential Entropy", title=chemical_name)
        )



        print(ggplot(hcal_1, aes(sec, h.sec, colour=conc))+
                geom_path(lwd=1.2) +
                facet_grid(.~ conc)+
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1"))+
                labs(x="Time", y= "Differential Entropy", title=chemical_name)
        )
      }

      print(ggplot(temp1, aes(starttime, mean_tot, group=conc,  colour=conc)) +
              geom_point(lwd=2)+
              geom_line(lwd=2)+
              geom_errorbar(limits, width=0.2)+
              ggtitle(paste0(low, ":", chemical_name))+
              ylim(0,8)+
              xlab("StartTime (s)")+
              ylab("total Movment (mm)")+
              # scale_x_discrete(labels=c("1468" = "2", "1470" = "4", "1472" = "6", "1474" = "8", "1476" ="10"))+
              # scale_color_manual(values=c("royal blue","red"))+ ##for vitamin E
              # scale_color_manual(values=c("red","royal blue"))+
              theme(legend.position="bottom")+
              scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                          "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                          "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                          "darkblue", "darkslateblue", "darkseagreen1",
                                          "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                          "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                          "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                          "tan", "slateblue3", "violetred4")))


      ##no error bar
      print(ggplot(temp1, aes(starttime, mean_tot, group=conc,  colour=conc)) +
              geom_point(lwd=2)+
              geom_line(lwd=2)+
              # ylim(0,5)+##for vitamin E
              # scale_color_manual(values=c("royal blue","red"))+ ##for vitamin E
              #   geom_path(lwd=2) +
              #         geom_errorbar(limits, width=0.2)+
              ggtitle(paste0(low, ":", chemical_name, " - n = ", count))+
              xlab("StartTime (min)")+
              ylab("total Movment (mm)")+
              theme(legend.position="bottom")+
              scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                          "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                          "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                          "darkblue", "darkslateblue", "darkseagreen1",
                                          "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                          "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                          "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                          "tan", "slateblue3", "violetred4")))

      print(ggplot(temp1, aes(starttime, mean_tot, group=conc,  colour=conc)) +
              geom_point(lwd=2)+
              geom_line(lwd=2)+
              geom_errorbar(limits, width=0.2)+
              ggtitle(paste0(low, ":",chemical_name))+
              xlab("StartTime (min)")+
              ylab("total Movment (mm)")+
              theme(legend.position="bottom")+
              facet_grid(.~ conc)+
              scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                          "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                          "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                          "darkblue", "darkslateblue", "darkseagreen1",
                                          "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                          "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                          "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                          "tan", "slateblue3", "violetred4")))

      ##plate location
      print(ggplot(temp1.1, aes(starttime, mean_tot, group=conc,  colour=conc)) +
              geom_point(lwd=2)+
              geom_path(lwd=2) +
              geom_errorbar(limits, width=0.2)+
              ggtitle(paste0(low, ":", chemical_name))+
              xlab("StartTime (sec)")+
              ylab("total Movment (mm)")+
              theme(legend.position="bottom")+
              facet_grid(plate_Location ~ conc)+
              scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                          "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                          "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                          "darkblue", "darkslateblue", "darkseagreen1",
                                          "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                          "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                          "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                          "tan", "slateblue3", "violetred4")))

      #### AUC barplots
      if(assay == "LPR"){
        stats_table$conc <- as.factor(stats_table$conc)
        print(ggplot(stats_table[which(stats_table$chemical.id == low),],
                     aes(interval, AUC, group=conc)) +
                geom_bar(aes(fill=conc), position="dodge", stat="identity")+
                theme_bw() +
                geom_text(aes(label=count, colour=conc), position=position_dodge(width=0.9), vjust = -0.25)+
                labs(x="", y="Average Area Under the Curve")+
                theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=1))+
                ggtitle(paste0(low, ":", chem.key$chemical.name[chem.key$chemical.id==low]))+
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1",
                                            "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                            "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                            "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                            "tan", "slateblue3", "violetred4"))+
                scale_fill_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                           "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                           "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                           "darkblue", "darkslateblue", "darkseagreen1",
                                           "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                           "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                           "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                           "tan", "slateblue3", "violetred4")))




        ##cdf plots
        if(length(which(is.nan(temp1$mean_tot)))<= 0){
          temp1 <- temp1
        }else{
          temp1 <- temp1[-which(is.nan(temp1$mean_tot)),]
        }

        maxx <- .5*max(temp1$mean_tot)
        print(x<-ggplot(temp1, aes(mean_tot, colour = conc)) +
                stat_ecdf(size=2) +
                xlab("Total Movement (mm)")+
                ylab("Portion of Time")+
                ggtitle(chemical_name) +
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "green", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1"))+
                annotate("text", x=-1.2, y=.125, label="Cycle 1", color = "black", cex=6)+
                annotate("text", x=-1.2, y=.375, label="Cycle 2", color = "black", cex=6)+
                annotate("text", x=-1.2, y=.625, label="Cycle 3", color = "black", cex=6)+
                annotate("text", x=-1.2, y=.875, label="Cycle 4", color = "black", cex=6)+
                annotate("text", x=maxx-2, y=0, label="Hypo", color = "black", cex=6)+
                annotate("text", x=maxx+2, y=0, label="Hyper", color = "black", cex=6)+
                geom_segment(aes(x=maxx-2.5,y=0, xend=maxx-4, yend=0), colour="black", arrow=arrow(), lwd=2)+
                geom_segment(aes(x=maxx+2.75,y=0, xend=maxx+4, yend=0), colour="black", arrow=arrow(), lwd=2)+
                theme(plot.title = element_text(size=20, face="bold",margin = margin(), debug = FALSE))
        )



      }

      if(assay == "LSR"){
        stats_table_LSR$conc <- as.factor(stats_table_LSR$conc)
        print(ggplot(stats_table_LSR[which(stats_table_LSR$chemical.id == low),],
                     aes(conc, AUC, group=conc)) +
                geom_bar(aes(fill=conc), position="dodge", stat="identity")+
                theme_bw() +
                geom_text(aes(label=count, colour=conc), position=position_dodge(width=0.9), vjust = -0.25)+
                labs(x="concentration", y="Average Area Under the Curve")+
                theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=1))+
                ggtitle(paste0(low, ":", chem.key$chemical.name[chem.key$chemical.id==low]))+
                scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                            "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                            "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                            "darkblue", "darkslateblue", "darkseagreen1",
                                            "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                            "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                            "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                            "tan", "slateblue3", "violetred4"))+
                scale_fill_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
                                           "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
                                           "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
                                           "darkblue", "darkslateblue", "darkseagreen1",
                                           "darkolivegreen", "darkorchid", "forestgreen", "darkred",
                                           "deeppink", "darkslategray3", "lightcyan3", "hotpink",
                                           "khaki", "pink3", "plum", "seagreen1", "tomato", "steelblue",
                                           "tan", "slateblue3", "violetred4"))+
                theme(legend.position="none"))
      } ## for LSR
      print(paste0("Making 5d Plots...", low))
    } ##for all the 5d plots
    dev.off()

    ######################Publication quality figures
    #   pdf(file=paste0(data.dir,"Reporting_Material/5dpf_", assay, "_Behavior_Plots_",data.file, "_", output.prefix, ".pdf",sep=""), width=16, height=8)
    #
    #   for(low in sub_chemids){
    #     chemical_name <- chem.key$chemical.name[chem.key$chemical.id==low]
    #     chemical_name <- str_replace_all(chemical_name,"_"," ")
    #
    #
    #     if(assay == "LPR"){
    #     ggplot(stats_table[which(stats_table$chemical.id == low),],
    #            aes(interval, AUC, group=conc)) +
    #       geom_bar(aes(fill=conc), position="dodge", stat="identity")+
    #       theme_bw() +
    #       geom_text(aes(label=count, colour=conc), position=position_dodge(width=0.9), vjust = -0.25)+
    #       labs(x="", y="Average Area Under the Curve")+
    #       theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=1))+
    #       ggtitle(chemical_name)+
    #        scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
    #                                   "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
    #                                   "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
    #                                   "darkblue", "darkslateblue", "darkseagreen1"))+
    #       scale_fill_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
    #                                   "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
    #                                   "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
    #                                   "darkblue", "darkslateblue", "darkseagreen1"))+
    #       theme(axis.text=element_text(size=12, face="bold"),
    #             axis.title=element_text(size=14,face="bold"))
    #     } ##printing for LPR
    #
    #     if(assay == "LSR"){
    #       print(ggplot(stats_table_LSR[which(stats_table_LSR$chemical.id == low),],
    #                    aes(conc, AUC, group=conc)) +
    #               geom_bar(aes(fill=conc), position="dodge", stat="identity")+
    #               theme_bw() +
    #               geom_text(aes(label=count, colour=conc), position=position_dodge(width=0.9), vjust = -0.25)+
    #               labs(x="concentration", y="Average Area Under the Curve")+
    #               theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=1))+
    #               ggtitle(chemical_name)+
    #               scale_color_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
    #                                           "gray", "purple", "red", "blue", "yellow", "red", "royllel.alblue4",
    #                                           "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
    #                                           "darkblue", "darkslateblue", "darkseagreen1"))+
    #               scale_fill_manual(values=c("black", "#CC6666", "#9999CC", "cornflowerblue",
    #                                          "gray", "purple", "red", "blue", "yellow", "red", "royalblue4",
    #                                          "darkmagenta", "darkolivegreen", "brown2", "darkcyan",
    #                                          "darkblue", "darkslateblue", "darkseagreen1"))+
    #               theme(legend.position="none")+
    #                theme(axis.text=element_text(size=12, face="bold"),
    #               axis.title=element_text(size=14,face="bold")))
    #     } ##printing for LSR
    # } ## ends printing for the report
    # graphics.off()
  } ##end of assays
}# if there are LPR files

##################################################################################################
###
#####Create a hits table
###
##################################################################################################
LEL.morph.stripped$chemical.id  <- as.character(LEL.morph.stripped$chemical.id)
LEL.morph.stripped <- LEL.morph.stripped %>% join(chem.key)

##heatmap of just morph
hit.tab <- LEL.morph.stripped
hit.tab[hit.tab == 0] <- NA
col_order <- c("chemical.id","chemical.name", "any.except.MORT", "any.effect", "MO24", "DP24",
               "SM24", "NC24", "MORT", "YSE_", "AXIS", "EYE_", "SNOU", "JAW_", "OTIC", "PE__",
               "BRAI", "SOMI", "PFIN","CFIN", "PIG_", "CIRC", "TRUN", "SWIM", "NC__",
               "TR__")
hit.tab <- hit.tab[, col_order]


# hit.tab <- hit.tab%>% select(starts_with("chemical.name"))
# hit.tab <- hit.tab[c(1,26,2:25,27:34,36,37)]
write.csv(hit.tab, file=paste0(data.dir, "/ZF_Morph_Hit_Call_Table_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)

##add colors for the columns
dd <- c(rep("royalblue", 2), rep("darkmagenta", 22))
#####################################################################################

#############MAKE A HEATMAP OF ALL THE DATA
heatmapfile_wide <- hit.tab %>% select(chemical.name, any.except.MORT, any.effect, MO24, DP24, SM24, NC24,
                                       MORT, YSE_, AXIS, EYE_, SNOU, JAW_, OTIC, PE__, BRAI, SOMI, PFIN,
                                       CFIN, PIG_, CIRC, TRUN, SWIM, NC__, TR__)

end_points_start <- "2"
heatmapfile<- heatmapfile_wide[,end_points_start:ncol(heatmapfile_wide)]

heatmapfile_matrix <- data.matrix(heatmapfile)

row.names(heatmapfile_matrix) <- hit.tab$chemical.name
colnames(heatmapfile_matrix) <- c("Any.Except.Mort", "Any.Effect", "24 hpf Mortality", "24 hpf Developmental Progression",
                                  "24 hpf Spontaneous Movement", "24 hpf Notochord", "120 hpf Mortality", "Yolk Sac Edema",
                                  "Axis", "Eye", "Snout", "Jaw", "Otolith", "Pericardial Edema", "Brain", "Somites",
                                  "Pectoral Fin", "Caudal Fin", "Pigmentation", "Circulation", "Trunk", "Swim Bladder",
                                  "Notochord", "Touch-Response")

##check if the full data matrix is only NA
if(any(!is.na(heatmapfile_matrix))== FALSE){
  heatmapfile_matrix[is.na(heatmapfile_matrix)] <- 0

  aheatmap(heatmapfile_matrix,
           color="white",
           breaks=1,
           scale="none",
           annRow = LEL.morph.stripped2["chemical.name"], annColors = ann_colors,
           Rowv=NA,Colv=NA,cellheight = 20, cellwidth=20, border=T, annLegend=NA,
           fontsize=5, labRow=row.names(heatmapfile_matrix),
           cexCol=1.5, cexRow=1.5,
           legend=NA,
           filename=paste0(data.dir, "Overview_Morphology_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""),
           width=11, height=8.5)
} else{
  aheatmap(heatmapfile_matrix,
           color=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
           breaks=1,
           scale="none",
           annRow = LEL.morph.stripped2["chemical.name"], annColors = ann_colors,
           Rowv=NA,Colv=NA,cellheight = 20, cellwidth=20, border=T, annLegend=NA,
           fontsize=5, labRow=row.names(heatmapfile_matrix),
           cexCol=1.5, cexRow=1.5,
           filename=paste0(data.dir, "Overview_Morphology_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""),
           width=11, height=20)
}

if(nrow(heatmapfile_matrix)>=2){
  pdf(paste(data.dir, "Overview_Morphology_Heatmap_Pretty_",data.file, "_", output.prefix, ".pdf",sep=""),
      height=11, width=20,bg = "white")
  heatmap <- heatmap.2(heatmapfile_matrix,
                       na.rm=TRUE,
                       cellnote = heatmapfile_matrix,
                       notecol="white",
                       notecex=.75,
                       keysize = .4,
                       key=FALSE,
                       key.title = "LEL",
                       denscol=NULL,
                       margins=c(20,40),
                       cexRow=1.5,
                       cexCol=1.3,
                       trace = "none",
                       colsep=1:ncol(heatmapfile_matrix),
                       rowsep=1:nrow(heatmapfile_matrix),
                       sepcolor="grey70",
                       dendrogram = c("none"),
                       ColSideColors = dd,
                       cex.main = 4,
                       # breaks = c(0.000969, 0.009, 0.02, 0.06, 0.1, 0.5,1,2,18,3),
                       scale="none",
                       Colv=NA, ##no cluster = NA
                       Rowv=NA, ##no cluster = NA
                       col=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
                       main="Morphology Heatmap",
                       xlab=""
  )

  leg_col <- unique(dd)
  legend("topright", legend = c("Any Effect", "Morphology"), title="Assays", fill=leg_col,
         border=FALSE, bty="y", y.intersp=0.75, cex=.9)
  graphics.off()
}



# if(nrow(vp_AUC)>1  & length(grep("Motion", format_files)) > 0){
if(length(grep("LPR", format_files)) >0 ){
  if(nrow(vp_AUC)>1  & exists("LEL.EPR")==TRUE){

    LEL.5d.hits <- vp_AUC[-which(vp_AUC$interval == "Transition"),] %>%
      spread(interval, LEL) %>% rename(LEL.dark = Dark, LEL.light = Light)


    hit.tab <- LEL.morph.stripped %>% inner_join(LEL.EPR) %>% inner_join(LEL.5d.hits) %>% inner_join(LSR.LEL)
    hit.tab[hit.tab == 0] <- NA
    # hit.tab <- hit.tab %>% select(-ks.pval, -delta, -n.fish, -n.MO24, -perc.MO24, -ordinality, -binary)


    col_order <- c("chemical.id","chemical.name", "any.except.MORT", "any.effect", "MO24", "DP24",
                   "SM24", "NC24", "MORT", "YSE_", "AXIS", "EYE_", "SNOU", "JAW_", "OTIC", "PE__",
                   "BRAI", "SOMI", "PFIN","CFIN", "PIG_", "CIRC", "TRUN", "SWIM", "NC__",
                   "TR__", "B", "E", "R", "LEL.dark", "LEL.light", "LSR.AUC", "LSR.Peak")
    hit.tab <- hit.tab[, col_order]


    # hit.tab <- hit.tab%>% select(starts_with("chemical.name"))
    # hit.tab <- hit.tab[c(1,26,2:25,27:34,36,37)]
    write.csv(hit.tab, file=paste0(data.dir, "/ZF_Hit_Call_Table_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)

    ##add colors for the columns
    dd <- c(rep("darkmagenta", 24), rep("darkorange3", 3), rep("gold2", 2), rep("lightgreen", 2))
    #####################################################################################

    #############MAKE A HEATMAP OF ALL THE DATA
    heatmapfile_wide <- hit.tab %>% select(chemical.name, any.except.MORT, any.effect, MO24, DP24, SM24, NC24,
                                           MORT, YSE_, AXIS, EYE_, SNOU, JAW_, OTIC, PE__, BRAI, SOMI, PFIN,
                                           CFIN, PIG_, CIRC, TRUN, SWIM, NC__, TR__, B, E, R, LEL.dark, LEL.light,
                                           LSR.AUC, LSR.Peak)
    # heatmapfile_wide <- hit.tab[c(2,25,26,3:24,27:29, 30:34)]

    end_points_start <- "2"
    heatmapfile<- heatmapfile_wide[,end_points_start:ncol(heatmapfile_wide)]

    ##fill empty space with 100
    # heatmapfile <- sapply(heatmapfile, function (x)
    # as.numeric(gsub("(0)|( +$)", "NA", x)))

    heatmapfile_matrix <- data.matrix(heatmapfile)

    row.names(heatmapfile_matrix) <- hit.tab$chemical.name
    colnames(heatmapfile_matrix) <- c("Any.Except.Mort", "Any.Effect", "24 hpf Mortality", "24 hpf Developmental Progression",
                                      "24 hpf Spontaneous Movement", "24 hpf Notochord", "120 hpf Mortality", "Yolk Sac Edema",
                                      "Axis", "Eye", "Snout", "Jaw", "Otolith", "Pericardial Edema", "Brain", "Somites",
                                      "Pectoral Fin", "Caudal Fin", "Pigmentation", "Circulation", "Trunk", "Swim Bladder",
                                      "Notochord", "Touch-Response", "EPR-Background", "EPR-Excitatory", "EPR-Refractory",
                                      "LPR-Dark Cycle", "LPR-Light Cycle", "LSR-All", "LSR-Peak")

    if(nrow(heatmapfile_matrix)==1){
      aheatmap(heatmapfile_matrix,
               breaks=3,
               color=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
               # color = c('white', brewer.pal(n=9, name="Blues")),
               scale="none",
               annRow = LEL.morph.stripped2["chemical.name"],
               annColors = ann_colors,
               Rowv=NA,Colv=NA,cellheight = 20,
               cellwidth=20, border=T, annLegend=NA,
               # main="Heatmap of Developmental Toxicity LEL",
               fontsize=5, labRow=row.names(heatmapfile_matrix),
               # fontsize=5, labRow=paste0(LEL.morph.stripped1$chemical.id, ":", LEL.morph.stripped1$chemical.name),
               cexCol=1.5, cexRow=1.5,
               filename=paste0(data.dir, "Overview_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""),
               width=11, height=8.5)
    }  else{
      aheatmap(heatmapfile_matrix,
               color=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
               # color = c('white', brewer.pal(n=9, name="Blues")),
               scale="none",
               annRow = LEL.morph.stripped2["chemical.name"], annColors = ann_colors,
               Rowv=NA,Colv=NA,cellheight = 20, cellwidth=20, border=T, annLegend=NA,
               # main="Heatmap of Developmental Toxicity LEL",
               fontsize=5, labRow=row.names(heatmapfile_matrix),
               # fontsize=5, labRow=paste0(LEL.morph.stripped1$chemical.id, ":", LEL.morph.stripped1$chemical.name),
               cexCol=1.5, cexRow=1.5,
               filename=paste0(data.dir, "Overview_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""),
               width=11, height=8.5)
    }

    if(nrow(heatmapfile_matrix)>=2){
      pdf(paste(data.dir, "Overview_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""), height=8.5, width=11,bg = "white")
      heatmap <- heatmap.2(heatmapfile_matrix,
                           na.rm=TRUE,
                           cellnote = heatmapfile_matrix,
                           notecol="white",
                           notecex=.75,
                           keysize = .4,
                           key=FALSE,
                           key.title = "LEL",
                           denscol=NULL,
                           margins=c(40,23),
                           cexRow=1.5,
                           cexCol=1,
                           trace = "none",
                           colsep=1:ncol(heatmapfile_matrix),
                           rowsep=1:nrow(heatmapfile_matrix),
                           sepcolor="grey70",
                           # labRow=TRUE,
                           # labCol=TRUE,
                           # vline="none",
                           # tracecol="gray",
                           #vline="white",
                           dendrogram = c("none"),
                           ColSideColors = dd,
                           cex.main = 4,
                           # breaks = c(0.000969, 0.009, 0.02, 0.06, 0.1, 0.5,1,2,18,3),
                           scale="none",
                           Colv=NA, ##no cluster = NA
                           Rowv=NA, ##no cluster = NA
                           col=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
                           main="Overall Heatmap",
                           xlab=""
      )

      leg_col <- unique(dd)
      legend("bottom", legend = c("Morphology", "EPR", "LPR", "LSR"), title="Assays", fill=leg_col,
             border=FALSE, bty="y", y.intersp=0.75, cex=.9)

      graphics.off()
    } ##end of graphs
  }else{
    ##tthis is if there is no EPR
    LEL.5d.hits <- vp_AUC[-which(vp_AUC$interval == "Transition"),] %>%
      spread(interval, LEL) %>% rename(LEL.dark = Dark, LEL.light = Light)


    hit.tab <- LEL.morph.stripped %>%  inner_join(LEL.5d.hits) %>% inner_join(LSR.LEL)
    hit.tab[hit.tab == 0] <- NA
    # hit.tab <- hit.tab %>% select(-ks.pval, -delta, -n.fish, -n.MO24, -perc.MO24, -ordinality, -binary)


    col_order <- c("chemical.id","chemical.name", "any.except.MORT", "any.effect", "MO24", "DP24",
                   "SM24", "NC24", "MORT", "YSE_", "AXIS", "EYE_", "SNOU", "JAW_", "OTIC", "PE__",
                   "BRAI", "SOMI", "PFIN","CFIN", "PIG_", "CIRC", "TRUN", "SWIM", "NC__",
                   "TR__", "LEL.dark", "LEL.light", "LSR.AUC", "LSR.Peak")
    hit.tab <- hit.tab[, col_order]


    # hit.tab <- hit.tab%>% select(starts_with("chemical.name"))
    # hit.tab <- hit.tab[c(1,26,2:25,27:34,36,37)]
    write.csv(hit.tab, file=paste0(data.dir, "/ZF_Hit_Call_Table_", data.file, "_", output.prefix, ".csv",sep=""), row.names=F)

    ##add colors for the columns
    dd <- c(rep("darkmagenta", 24),
            # rep("darkorange3", 3),
            rep("gold2", 2),
            rep("lightgreen", 2))
    #####################################################################################

    #############MAKE A HEATMAP OF ALL THE DATA
    heatmapfile_wide <- hit.tab %>% select(chemical.name, any.except.MORT, any.effect, MO24, DP24, SM24, NC24,
                                           MORT, YSE_, AXIS, EYE_, SNOU, JAW_, OTIC, PE__, BRAI, SOMI, PFIN,
                                           CFIN, PIG_, CIRC, TRUN, SWIM, NC__, TR__,
                                           # B, E, R,
                                           LEL.dark, LEL.light,
                                           LSR.AUC, LSR.Peak)
    # heatmapfile_wide <- hit.tab[c(2,25,26,3:24,27:29, 30:34)]

    end_points_start <- "2"
    heatmapfile<- heatmapfile_wide[,end_points_start:ncol(heatmapfile_wide)]

    ##fill empty space with 100
    # heatmapfile <- sapply(heatmapfile, function (x)
    # as.numeric(gsub("(0)|( +$)", "NA", x)))

    heatmapfile_matrix <- data.matrix(heatmapfile)

    row.names(heatmapfile_matrix) <- hit.tab$chemical.name
    colnames(heatmapfile_matrix) <- c("Any.Except.Mort", "Any.Effect", "24 hpf Mortality", "24 hpf Developmental Progression",
                                      "24 hpf Spontaneous Movement", "24 hpf Notochord", "120 hpf Mortality", "Yolk Sac Edema",
                                      "Axis", "Eye", "Snout", "Jaw", "Otolith", "Pericardial Edema", "Brain", "Somites",
                                      "Pectoral Fin", "Caudal Fin", "Pigmentation", "Circulation", "Trunk", "Swim Bladder",
                                      "Notochord", "Touch-Response",
                                      # "EPR-Background", "EPR-Excitatory", "EPR-Refractory",
                                      "LPR-Dark Cycle", "LPR-Light Cycle", "LSR-All", "LSR-Peak")

    aheatmap(heatmapfile_matrix,
             color=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
             # color = c('white', brewer.pal(n=9, name="Blues")),
             scale="none",
             annRow = LEL.morph.stripped2["chemical.name"], annColors = ann_colors,
             Rowv=NA,Colv=NA,cellheight = 20, cellwidth=20, border=T, annLegend=NA,
             # main="Heatmap of Developmental Toxicity LEL",
             fontsize=5, labRow=row.names(heatmapfile_matrix),
             # fontsize=5, labRow=paste0(LEL.morph.stripped1$chemical.id, ":", LEL.morph.stripped1$chemical.name),
             cexCol=1.5, cexRow=1.5,
             filename=paste0(data.dir, "Overview_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""),
             width=11, height=8.5)

    if(nrow(heatmapfile_matrix)>=2){
      pdf(paste(data.dir, "Overview_Heatmap_",data.file, "_", output.prefix, ".pdf",sep=""), height=8.5, width=11,bg = "white")
      heatmap <- heatmap.2(heatmapfile_matrix,
                           na.rm=TRUE,
                           cellnote = heatmapfile_matrix,
                           notecol="white",
                           notecex=.75,
                           keysize = .4,
                           key=FALSE,
                           key.title = "LEL",
                           denscol=NULL,
                           margins=c(40,23),
                           cexRow=1.5,
                           cexCol=1,
                           trace = "none",
                           colsep=1:ncol(heatmapfile_matrix),
                           rowsep=1:nrow(heatmapfile_matrix),
                           sepcolor="grey70",
                           # labRow=TRUE,
                           # labCol=TRUE,
                           # vline="none",
                           # tracecol="gray",
                           #vline="white",
                           dendrogram = c("none"),
                           ColSideColors = dd,
                           cex.main = 4,
                           # breaks = c(0.000969, 0.009, 0.02, 0.06, 0.1, 0.5,1,2,18,3),
                           scale="none",
                           Colv=NA, ##no cluster = NA
                           Rowv=NA, ##no cluster = NA
                           col=colorRampPalette(rev(c("lightblue",  "midnightblue")))(50),
                           main="Overall Heatmap",
                           xlab=""
      )

      leg_col <- unique(dd)
      legend("bottom", legend = c("Morphology", "LPR", "LSR"), title="Assays", fill=leg_col,
             border=FALSE, bty="y", y.intersp=0.75, cex=.9)

      graphics.off()
    }
