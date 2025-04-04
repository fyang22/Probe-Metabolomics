#################### Read mzml data fuction #################### 
read_mzml <- function(file_path){
    list_of_files <- list.files(path = file_path,
                            recursive = TRUE,
                            pattern = "\\.mzML$",
                            full.names = TRUE)
    sps <- Spectra(list_of_files)
    return(sps)
}

#Fuction: Normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

#Fuction: Normalize MS/MS, 
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * (5 / 100 ) #remove the Int < 5%
}

#################### filter spectra contain probe fragments ####################
sps_probe <- function(sps,mass_probe,mass_fragment){
    sps <- filterPrecursorMz(sps, c(mass_probe,1200))
    sps_normalized <- addProcessing(sps, norm_int)
    sps_normalized <- filterIntensity(sps_normalized, intensity = low_int)
    has_probe_fragment <- containsMz(sps_normalized,
                                    mz = mass_fragment,
                                    tolerance = 0.005,
                                    ppm = 5,
                                    which = "all")
    sps_sub <- sps_normalized[has_probe_fragment]
    return(sps_sub)
} 

##################### write spectra table #######################
# sps_table_maker <- function(sps_probe, mass_probe, mass_group, polarity){
#     df_sps <- spectraData(sps_probe,
#                             c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
#     df_sps <- as.data.frame(df_sps)
#     df_sps$mass_metabolite <- df_sps$precursorMz - mass_probe + mass_group
#     df_sps$mass_ion <- if(polarity == 1) {
#         round(df_sps$mass_metabolite + proton, 4)
#     } else if(polarity == 0) {
#         round(df_sps$mass_metabolite - proton, 4)
#     } else {
#         message("No polarity")
#         NA  # Return NA or handle missing case
#     }
#     write.csv(df_sps,here("output", paste0(substitute(df_sps), ".csv")))
# }


##################### function for ms1 match ####################
MS1_match <- function(df_probe, df_db,db_param,mz){
    MS1_match <- matchValues(df_probe, 
                         df_db, 
                         db_param, mzColname= mz)
    MS1_matched <- matchedData(MS1_match)[whichQuery(MS1_match),]
    MS1_matched <- MS1_matched[!is.na(MS1_matched$score),]
    MS1_matched <- MS1_matched[order(MS1_matched$ppm_error,decreasing = FALSE),]
    #MS1_matched <- MS1_matched[!duplicated(MS1_matched$ID),]
    MS1_matched$ID <- seq.int(nrow(MS1_matched))
    return(MS1_matched)
}




# Read control mzml files & filter precursor mz range
proton <- 1.00728
light_probe <- 354.2056
delta_mass_13C <- 1.00728*6
mass_group <- 17.0033
fragment_light <- c(105.0335,148.0757)#,160.0756)
fragment_heavy <- c(111.0509,154.0931,160.0756)
heavy_probe <- light_probe + delta_mass_13C
print(heavy_probe)
adduct = "[M+H]+"

# Read light probe
sps_light <- read_mzml(here("data","raw","light"))

# filter spectra by mz range and fragmentation
sps_probe_light <- sps_probe(sps_light,light_probe,fragment_light)
df_probe_light  <- spectraData(sps_probe_light,
                            c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
#df_probe_light$precursorMz <- round(df_probe_light$precursorMz,4)
df_probe_light$mass_metabolite <- df_probe_light$precursorMz - light_probe + mass_group
df_probe_light<-df_probe_light[!duplicated(df_probe_light$mass_metabolite),]
df_probe_light$mass_protonated <- df_probe_light$mass_metabolite + proton
df_probe_light$mass_deprotonated <- df_probe_light$mass_metabolite - proton
write.csv(df_probe_light,here("output", paste0(substitute(df_probe_light), "4.csv")))

df_probe_light <- Mass2MzParam(adducts = adduct,
                     tolerance = 0.002, ppm = 15)

# df_sub_sps_light_mass_molecule mtch with HMDB
df_probe_light<- read.csv(here("data","hmdb_cleanup_v02062023.csv"),header = TRUE, sep = ",")

probe_match <- MS1_match(df_probe_light,df_hmdb,parm2,"mass_protonated")
write.csv(probe_match,here("output", "light_hmdb_matched5.csv"))
# df_mtches <- as.data.frame(mtches_aglycon)
# df_mtches <- df_mtches  %>%
#   select(rtime, precursorMz,
#   target_rtime,target_precursorMz,target_aglycon_theo,target_exactmass,
#   adduct,ppm_error,ID) %>%       # Select columns
#   rename(EA_rtime = rtime,          # Rename
#          EA_mz = precursorMz,
#          rtime = target_rtime,
#          mz = target_precursorMz,
#          ppm_aglycon = ppm_error
#          )
# write.csv(df_mtches, here("output","aglycon_match_cleanup.csv"))

# Read heavy probe
sps_heavy <- read_mzml(here("data","raw","heavy"))
# filter spectra by mz range and fragmentation
sps_heavy <- filterPrecursorMz(sps_heavy, c(heavy_probe,1200))
has_heavy_probe_fp <- containsMz(sps_heavy,
                                mz = c(111.0509,154.0931,160.0756), # ms/ms fragmentation patterns of heavy probe
                                tolerance = 0.005, ppm = 5,
                                which = "all"
                                )
sub_sps_heavy <- sps_heavy[has_heavy_probe_fp]
df_sub_sps_heavy <- spectraData(sub_sps_heavy,
                            c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
df_sub_sps_heavy$light_probe_est <- (df_sub_sps_heavy$precursorMz-delta_mass_13C)
write.csv(df_sub_sps_heavy,here("output","sps_heavy.csv"))