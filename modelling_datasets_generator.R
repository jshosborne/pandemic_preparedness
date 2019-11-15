#code to translate master extraction sheet for MERS-CoV into separate experimental run datasets

dat_dir <- commandArgs()[4]

# master<-read.csv(paste0(gsub('09_modelling_datasets', '', dat_dir), "/04_extraction_sheet/MERSCoV_extraction_master_copy_3.6.19.csv"),
#                  stringsAsFactors = FALSE)

master <- read.csv(paste0(dat_dir, '/who_lit_merge_03_20.csv'), as.is=T, stringsAsFactors = FALSE)

animal<-subset(master, master$organism_type == "mammal")
animal<-animal[is.na(animal$patient_type),]

master<-master[master$patient_type!="absent",]
master<-master[master$patient_type!="import",]

#if there are any human cases that have patient_type == NA, reclassify as "unspecified"
for (i in 1:nrow(master)){
  if(toString(master$organism_type[i])=="human" && is.na(master$patient_type[i])){
    master$patient_type[i]<-"unspecified"
  }
}
human<-subset(master, master$organism_type == "human")

# Experiment 1 - No secondary human cases. detections: mers-cov & mers-like
exp1<-rbind(animal, subset(human, human$patient_type != 'secondary'))

write.csv(exp1,
          paste0(dat_dir,
                 "/experiment_subsets/experiment1_all_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)

# Experiment 2 - detections: mers-cov only
exp2<-subset(exp1, exp1$pathogen == "MERS-CoV")

write.csv(exp2,
          paste0(dat_dir,
                 "/experiment_subsets/experiment2_mers_cov_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)

#Experiment 3 - all index detections: mers-cov & mers-like
exp3<-subset(human, human$patient_type == 'index')

write.csv(exp3,
          paste0(dat_dir,
                 "/experiment_subsets/experiment3_index_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)
		  
# Experiment 4 - all index detections: mers-cov only
exp4<-subset(human, human$patient_type == 'index' & human$pathogen == "MERS-CoV")

write.csv(exp4,
          paste0(dat_dir,
                 "/experiment_subsets/experiment4_mers_cov_index_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)
		  
# Experiment 5 - only human detections: mers-cov & mers-like
exp5<-subset(human, human$patient_type != 'secondary')

write.csv(exp5,
          paste0(dat_dir,
                 "/experiment_subsets/experiment5_human_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)

# Experiment 6 - only human detections: mers-cov only
exp6<-subset(human, human$patient_type != 'secondary' & human$pathogen == "MERS-CoV")

write.csv(exp6,
          paste0(dat_dir,
                 "/experiment_subsets/experiment6_mers_cov_human_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)

# Experiment 7 - pcr detections: mers-cov & mers-like
exp7<-subset(exp5, exp5$diagnostic == 'PCR')

write.csv(exp7,
          paste0(dat_dir,
                 "/experiment_subsets/experiment7_pcr_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)
		  
# Experiment 8 - pcr detections: mers-cov only 
exp8<-subset(exp6, exp6$diagnostic == 'PCR')

write.csv(exp8,
          paste0(dat_dir,
                 "/experiment_subsets/experiment8_mers_cov_pcr_",
                 Sys.Date(),
                 ".csv"),
          row.names = FALSE)

write.csv(data.frame(), paste0(dat_dir, '/experiment_subsets/finished.csv'))