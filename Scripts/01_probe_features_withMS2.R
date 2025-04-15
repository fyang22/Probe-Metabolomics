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
    x > max(x, na.rm = TRUE) * (10 / 100 ) #remove the Int < 5%
}

#################### filter spectra contain probe fragments ####################
sps_probe <- function(sps,mass_probe,mass_fragment){
    sps <- filterPrecursorMzRange(sps, c(mass_probe,1200))
    sps_normalized <- addProcessing(sps, norm_int)
    sps_normalized <- filterIntensity(sps_normalized, intensity = low_int)
    has_probe_fragment <- containsMz(sps_normalized,
                                    mz = mass_fragment,
                                    tolerance = 0.005,
                                    ppm = 10,
                                    which = "all")
    sps_sub <- sps_normalized[has_probe_fragment]
    return(sps_sub)
} 



# Read control mzml files & filter precursor mz range
proton <- 1.00728
light_probe <- 353.19832
delta_mass_13C <- 6.0201
mass_group <- 17.0033
#fragment_light <- c(105.0335,148.0757,160.0756)
fragment_light <-c(148.0757,160.0756)
#fragment_heavy <- c(111.0509,154.0931,160.0756)
fragment_heavy <- c(154.0931,160.0756)
heavy_probe <- light_probe + delta_mass_13C
print(heavy_probe)
adduct = "[M+H]+"

df_hmdb<- read.csv(here("data","hmdb_cleanup_v02062023.csv"))#,header = TRUE, $
df_hmdb$exactmass_light <- df_hmdb$exactmass + light_probe - mass_group
df_hmdb$exactmass_heavy <- df_hmdb$exactmass + heavy_probe - mass_group
write.csv(df_hmdb,here("hmdb_CA_probe.csv"))


parm2 <- Mass2MzParam(adducts = adduct,
                     tolerance = 0.005, ppm = 5)

# Read light probe
sps_light <- read_mzml(here("data","raw","light"))
sps_light_sub <- sps_probe(sps_light,light_probe,fragment_light)

df_sps_light <- spectraData(sps_light_sub,c("rtime","dataOrigin","precursorMz"))
write.csv(df_sps_light,here("output","sps_probe_light.csv"))

MS1_match <- matchValues(df_sps_light, 
                         df_hmdb, 
                         parm2, mzColname= "precursorMz",
                         massColname = "exactmass_light")

MS1_matched <- matchedData(MS1_match)
MS1_matched <- MS1_matched[!is.na(MS1_matched$ppm_error),]
MS1_matched <- MS1_matched[order(MS1_matched$ppm_error,decreasing = FALSE),]
MS1_matched <- MS1_matched[!duplicated(MS1_matched$target_accession),]
MS1_matched <- MS1_matched[,c("rtime","dataOrigin","precursorMz",
                        "target_accession","target_name","target_exactmass",
                        "target_exactmass_light","target_exactmass_heavy","target_chemical_formula",
                        "adduct","ppm_error")]
write.csv(MS1_matched,here("output", "light_hmdb_match.csv"))


# Read heavy probe
sps_heavy <- read_mzml(here("data","raw","heavy"))
sps_heavy_sub <- sps_probe(sps_heavy,heavy_probe,fragment_heavy)

df_sps_heavy <- spectraData(sps_heavy_sub,c("rtime","dataOrigin","precursorMz"))
write.csv(df_sps_heavy,here("output","sps_probe_heavy.csv"))

MS1_match2 <- matchValues(df_sps_heavy, 
                         df_hmdb, 
                         parm2, mzColname= "precursorMz",
                         massColname = "exactmass_heavy")

MS1_matched <- matchedData(MS1_match2)
MS1_matched <- MS1_matched[!is.na(MS1_matched$ppm_error),]
MS1_matched <- MS1_matched[order(MS1_matched$ppm_error,decreasing = FALSE),]
MS1_matched <- MS1_matched[!duplicated(MS1_matched$target_accession),]
MS1_matched <- MS1_matched[,c("rtime","dataOrigin","precursorMz",
                        "target_accession","target_name","target_exactmass",
                        "target_exactmass_light","target_exactmass_heavy","target_chemical_formula",
                        "adduct","ppm_error")]
write.csv(MS1_matched,here("output", "heavy_hmdb_match.csv"))
