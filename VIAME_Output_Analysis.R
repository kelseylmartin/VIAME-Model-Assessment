#### Cleaning VIAME Data Output
## Authors: Kelsey Martin, Jack Prior, Sara Thompson
# Nov. 5, 2025


# installing packages
inst <- c("readr","dplyr","tidyverse","lubridate","data.table","microbenchmark","tidyr","knitr","ggplot2", "ggpubr", "gridExtra", "FSA", "FSAdata")
for (p in inst) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p)
    require(p,character.only = TRUE)}
}

### Important functions
# to select multiple stations for min count @ multiple confidences
#########################################################################
# if all CSVs come from same source web or desktop use this read function
read_plus <- function(flnm) {
  tryCatch({
    read_csv(flnm,skip=2,col_names = FALSE, col_select = c(1,2,3,4,5,6,7,8,9,10,11), show_col_types = FALSE) %>%  #Original run has 2 header rows
      mutate(filename = flnm) %>% 
      mutate(across(everything(), as.character))
  }, error = function(e) {
    print(paste("Error reading file: ", flnm)) # this has been happening when files don't have any observations. Saw it only in NGI files
    return(NULL)  # If there's an error reading the file, return NULL
  })
}



# if mixing online and desktop produced CSVs use the below section which removes time columns that are incompatible for joining. VIAME time format isn't great aanyway
read_plus_mixed <- function(flnm) {
  read_csv(flnm,skip=2,col_names = FALSE, col_select = c(1,3,4,5,6,7,8,9,10,11), show_col_types = FALSE) %>%  #Original run has 2 header rows
    mutate(filename = flnm) %>% 
    mutate(across(everything(), as.character))
}

## Naming your directories. These will be consistent throughout the script and will prevent each user from having to change paths for every file
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
data_root <- file.path(script.dir, "Data")
wrkdir  <- file.path(data_root, "/") # working directory (folder where data files are located)
trkdir  <- file.path(data_root, "Tracks/") # track directory (folder where all track data files are located)
trthdir <- file.path(data_root, "Truth/") # truth directory (folder where all groundtruthed data files are located)

# Create output directories - this will create the appropriate output folders wherever you have this R script located
dir.create(paste0(script.dir, "/Output")) # output directory
dir.create(paste0(script.dir, "/Output/Part I - Counts")) # part I output directory (count data per model run and confidence)
dir.create(paste0(script.dir, "/Output/Part II - Groundtruthing")) # part II output directory (groundtruthed data)
dir.create(paste0(script.dir, "/Output/Part III - Data Analysis")) # part III output directory (data analysis)
dir.create(paste0(script.dir, "/Output/Part III - Data Analysis/Figures")) # part III output directory (data analysis figures)
outdir <- paste0(script.dir, "/Output/")  # (where you want new data output files to go) 


## Bringing in important csv files
Allreadtimes <- read.csv(file.path(wrkdir, "All GFISHER Read Times.csv"))
Readtimekey <- read.csv(file.path(wrkdir, "final_filled_key_v2.csv")) %>% mutate(Videotime = sub("\\..*", "", Timestamp))
Species_List <- read.csv(file.path(wrkdir, "Species List.csv"))
fwri_ref_key <- read.csv(paste0(wrkdir, "env3LABS_93to24.csv")) %>%
  dplyr::mutate(SITE_ID = gsub("_|-", "", SITE_ID),
                FWRI_Ref = ifelse(LAB != "PC", REFERENCE, SITE_ID)) %>%
  dplyr::select(FWRI_Ref, SITE_ID, YEAR, LAB, HAB_STRAT) %>%
  dplyr::rename(Deployment = FWRI_Ref) %>%
  dplyr::filter(SITE_ID != "")


# run everything above this line first
# ----------------------------------------------------------------------------------------------------------------------------------------------
##################################################
## Part I: Importing Tracks and Creating Counts ##
##################################################
# MANUAL INPUT!! Putting all of the different VIAME model runs into a list that we can loop through below. This can grow or shrink depending on the different model runs that you want to look at and it will not affect the code below. I currently have my folders structured where the West files are in a separate folder within my data folder. Adjust below how you have your folders structured.
dts <- list.files(paste0(trkdir, ""), full.names = T)


# looping through each folder you have selected, in our case this will loop through each year_model
for (t in 1:length(dts)) {
  # singling out each of the datasets
  dt <- dts[[t]]
  filename <- basename(dt)
  
  # creating a variable to identify which model run and version we are looking at (to be used later in the code for exporting final datasets)
  model.run <- as.character(sub('.*_', '', filename))
  year <- as.numeric(sub('_.*',"",filename))
  print(paste0("Loading in track files for ", year, " v", model.run))
  
  # loading all CSVs into a table
  tbl.raw <- list.files(dt, pattern = "*.csv", full.names = T) 
  tbl <- tbl.raw %>% 
    map_df(~read_plus(.))
  # note that the table will not be produced if all columns are not formatted correctly 

  # renaming and fixing columns and column headers
  tbl$filename <- basename(tbl$filename)
  tbl$filename <- gsub("\\.csv$", "", tbl$filename)
  tbl$filename <- gsub("_\\d+\\.\\d+_(tracks|detections)$|_cam.*_(tracks|detections)$|_(tracks|detections)_cam.*$| cam.*_(tracks|detections)$|_(tracks|detections)$","",tbl$filename)
  tbl$filename <- gsub("_|-", "", tbl$filename)
  
  # Add the 2019-SPECIFIC fix for prefixes. This will be necessary to join with the groundtruthing and read time files
  if (year == 2019) {
    tbl$filename <- gsub("^19W_", "2019W", tbl$filename)
    tbl$filename <- gsub("^19E_", "2019E", tbl$filename)
  }
  
  # in 2022, PC deployment names don't include the year for some reason so we first look for the cases where the filename starts with a letter and then we add the year 
  if (year == 2022) {
    indices_to_fix <- grepl("^[A-Za-z]", tbl$filename)
    tbl$filename[indices_to_fix] <- paste0("2022-", tbl$filename[indices_to_fix])
  }
  
  col.names <- c("TrackID", "VidIdent", "UniqFrame", "TL_X", "TL_Y", "BR_X", "BR_Y", "DetLen_Conf", "Tar_Len", "SP", "CP", "Deployment")
  names(tbl) <- col.names
  
  # had to import all columns as characters so now will reassign certain columns as numbers
  tbl[,c("TrackID", "UniqFrame", "TL_X", "TL_Y", "BR_X", "BR_Y", "DetLen_Conf", "Tar_Len", "CP")] <- tbl[,c("TrackID", "UniqFrame", "TL_X", "TL_Y", "BR_X", "BR_Y", "DetLen_Conf", "Tar_Len", "CP")] %>%  
    mutate(across(everything(), as.numeric))
    
  # separating out stations and number of distinct occurrences for each station
  stations <- dplyr::distinct(tbl,Deployment)
  
  # Create number of distinct track file
  ndist<-tbl %>% dplyr::group_by(Deployment) %>% dplyr::summarise(ndist = n_distinct(TrackID))
  write.csv(ndist, file = paste0(outdir, model.run, "_ndist.csv"), row.names = F)
  
  ########## Shortening the track files ###########
  Allreadtimes$ReferenceID <- gsub("_|-", "", Allreadtimes$ReferenceID)
  
  # Get the minimum Start frame for each deployment
  start_lookup <- Allreadtimes %>%
    left_join(Readtimekey, by = c("StartTime" = "Videotime")) %>%
    dplyr::group_by(ReferenceID) %>%
    dplyr::summarise(Start = min(Frame, na.rm = TRUE), .groups = "drop")
  
  # Get the maximum End frame for each deployment
  end_lookup <- Allreadtimes %>%
    left_join(Readtimekey, by = c("EndTime" = "Videotime")) %>%
    dplyr::group_by(ReferenceID) %>%
    dplyr::summarise(End = max(Frame, na.rm = TRUE), .groups = "drop")
  
  # Join them into one final lookup table
  frame_lookup <- full_join(start_lookup, end_lookup, by = "ReferenceID")
  
  # Handle -Inf/Inf from min/max on empty sets (replaces them with NA)
  frame_lookup <- frame_lookup %>%
    mutate(Start = ifelse(is.infinite(Start), NA, Start),
           End = ifelse(is.infinite(End), NA, End))
  
  # Filter the tracks dataframe to keep only frames within the start and end range. But this only works is there is data for the start and end frames. In other words, NAs will be filtered out so we have to add them back in.
  tbl.final <- tbl %>%
    left_join(frame_lookup, by = c("Deployment" = "ReferenceID"))
  tbl.final <- tbl.final %>%
    dplyr::filter((UniqFrame >= Start & UniqFrame <= End) | (is.na(Start) & is.na(End)))
  
  ########## Calculating Various Metrics ###########
  cplong <- tbl.final %>%
    dplyr::select(TrackID, UniqFrame, Deployment, CP) %>%
    pivot_longer(
      cols = starts_with("C"),
      names_to = "Groups_p",
      values_to = "Probs",
      values_drop_na = TRUE
    )
  
  splong <- tbl.final %>%
    dplyr::select(TrackID, UniqFrame, Deployment, SP) %>%
    pivot_longer(
      cols = starts_with("S"),
      names_to = "Groups_s",
      values_to = "Species",
      values_drop_na = TRUE
    )
  
  vertrnk <- merge(cplong, splong, by = c("TrackID", "UniqFrame", "Deployment"), all.x = T)
  vertrnk <- vertrnk %>% dplyr::select(TrackID, UniqFrame, Deployment, Species, Probs)
  
  #create multiple outputs by adjusting filter of "Probs" at multiple intervals to desired resolution of confidence
  for (c in c(1:9,95)) { 
    confidence <- as.numeric(paste0("0.", c))
    subtbl.name <- paste0("out", confidence)
    print(paste0(year, " v", model.run, " ", confidence))
    
    out <- vertrnk %>%
    dplyr::group_by(Deployment,TrackID,UniqFrame) %>% 
    dplyr::mutate(grpmax = max(Probs)) %>% 
    ungroup() %>%
    filter(Probs>confidence) %>%
    arrange(Deployment,UniqFrame) %>%
    mutate(count=1)
   
    #Pull first instances by species and calculate the time (S) based on VIAME processing at 5 frame per second
    #adjust the desired confidence by adjusting out0.X and rerun from here to the end for each confidence interval
    out <- out %>% dplyr::group_by(Deployment)
    first <- out %>%
      dplyr::distinct(Species,.keep_all = TRUE) %>%
      dplyr::mutate(First = UniqFrame/5,Spec_Viame_Dash = Species) %>% #5 fps
      dplyr::select(Deployment,Spec_Viame_Dash, UniqFrame, First)
    first <- first %>%
      dplyr::inner_join(Species_List, by = "Spec_Viame_Dash") %>% 
      dplyr::select(Deployment,Species,UniqFrame,First)
    
    #Use ftable version of crosstabs to get species specific frequencies by station and unique camera frame 
    #adjust the desired confidence by adjusting out0.X
    freqtable <-out %>% 
      dplyr::select(Deployment,UniqFrame,Species) %>% 
      ftable() %>% 
      as.data.frame() %>% 
      arrange(Deployment,UniqFrame)
    
    #Summarize mean and min counts
    #MinCount
    sta_trkmin <- suppressMessages(freqtable %>% 
      dplyr::group_by(Deployment,Species) %>%
      dplyr::summarise(max(Freq)))
    names(sta_trkmin)[2] <- "Spec_Viame_Dash"
    names(sta_trkmin)[3] <- "MaxNCount"
    sta_trkminw <- sta_trkmin %>% pivot_wider(names_from = Spec_Viame_Dash, values_from = MaxNCount)
    #Join back to species table with no biocode
    sta_trkmin <- dplyr::inner_join(sta_trkmin, Species_List, by = "Spec_Viame_Dash")
    sta_trkmin <- dplyr::select(sta_trkmin,Deployment,MaxNCount,Species)
    #write.csv(sta_trkmin, file = paste0(outdir, model.run, "_", confidence, "min.csv"), row.names = F) #exporting just min count, adjust file name for model type and confidence
    
    #MeanCount
    sta_trkmean <- suppressMessages(freqtable %>% 
      dplyr::group_by(Deployment,Species) %>%
      dplyr::summarise(sum(Freq)))
    names(sta_trkmean)[2] <- "Spec_Viame_Dash"
    names(sta_trkmean)[3] <- "Sums"
    sta_trkmean <- dplyr::mutate(sta_trkmean, MeanCount = Sums/7500) %>% 
      dplyr::select(Deployment, Spec_Viame_Dash, MeanCount) #25 min video * 60 sec/min * 5 fps evaluation in Viame = 7500 frames evaluated 
    sta_trkmeanw <- sta_trkmean %>% pivot_wider(names_from = Spec_Viame_Dash, values_from = MeanCount)
    #Join back to species table with no biocode
    sta_trkmean <- dplyr::inner_join(sta_trkmean, Species_List, by = "Spec_Viame_Dash")
    sta_trkmean <- dplyr::select(sta_trkmean,Deployment,MeanCount,Species)
    #write.csv(sta_trkmean, file = paste0(outdir, model.run, "_", confidence, "mean.csv"), row.names = F)  # exporting just mean count, adjust file name for model type and confidence
    
    #TotalCount
    sta_trktotal <- suppressMessages(freqtable %>% 
      dplyr::group_by(Deployment,Species) %>%
      dplyr::summarise(sum(Freq)))
    names(sta_trktotal)[2] <- "Spec_Viame_Dash"
    names(sta_trktotal)[3] <- "Sums"
    sta_trktotal <- dplyr::mutate(sta_trktotal, TotalCount = Sums) %>% 
      dplyr::select(Deployment, Spec_Viame_Dash, TotalCount)
    sta_trktotalw <- sta_trktotal %>% pivot_wider(names_from = Spec_Viame_Dash, values_from = TotalCount)
    #Join back to species table with no biocode
    sta_trktotal <- dplyr::inner_join(sta_trktotal, Species_List, by = "Spec_Viame_Dash")
    sta_trktotal <- dplyr::select(sta_trktotal,Deployment,TotalCount,Species)
    #write.csv(sta_trktotal, file = paste0(outdir, model.run, "_", confidence, "total.csv"), row.names = F) # exporting just total count, adjust file name for model type and confidence
    
    #Pull MaxN, MeanCount and Time at first arrival into single file to export
    sta_vars <- dplyr::inner_join(sta_trkmean, sta_trkmin, by = c("Deployment","Species"))
    sta_vars <- dplyr::inner_join(sta_vars, first, c("Deployment","Species"))
    sta_vars <- dplyr::inner_join(sta_vars, ndist, c("Deployment"))
    sta_vars <- dplyr::inner_join(sta_vars, sta_trktotal, by = c("Deployment","Species"))
    
    #Include metadata on model run for version of model being used
    sta_vars <- sta_vars %>% dplyr::mutate(Version = paste0("v", model.run))
    sta_vars$year <- as.numeric(year)
    sta_vars$Deployment <- as.character(gsub('_|-', '', sta_vars$Deployment))
    sta_vars$Confidence <- as.character(confidence)
    
    # exporting
    write.csv(sta_vars, file = paste0(outdir, "Part I - Counts/", year, "_v", model.run, "_", confidence, ".csv"), row.names = F)
 
  }
  if (t == length(dts)) {print("Part I Complete!")}
}




##################################################
############ Part II: Groundtruthing #############
##################################################
# pull in all manual count files
# No longer using the method of having separate groundtruthing files for each site/cruise. From now on, using the full database files from PASC and PC (and eventually FWC)
truth <- read.csv(paste0(trthdir, "/maxn3LABS_93to24.csv")) %>% 
  #dplyr::filter(LAB %in% c("PASC", "PC")) %>%  # can add a filter here to only look at data from a particular lab
  dplyr::select(-LAB)

names(truth) <- toupper(names(truth))
truth$REFERENCE <- gsub("_|-", "", truth$REFERENCE)
truth[-1] <- lapply(truth[-1], function(x) as.numeric(as.character(x)))
truth[is.na(truth)] <- 0
truth[truth == ""] <- 0

# below is the code needed to keep the species that have greater than or equal to 20 occurrences
presence_absence_matrix <- (truth > 0)*1 
occurrences <- colSums(presence_absence_matrix, na.rm = TRUE)
species_to_keep <- names(occurrences[occurrences > 20]) # keeping species that have greater than 20 occurrences
truth_final <- truth[, species_to_keep]


# Make stationkey row name, remove no-fish stations, THEN BRING STATIONKEY BACK
truth_final <- truth_final %>% 
  column_to_rownames(var = "REFERENCE") %>% 
  dplyr::filter(rowSums(.) > 0) %>% 
  rownames_to_column(var = "Deployment") # Bring Deployment back
truth_final <- merge(truth_final, fwri_ref_key, by = "Deployment", all.x = T) %>%
  dplyr::mutate(Deployment = ifelse(LAB %in% c("PC", "FWRI") | (!is.na(YEAR) & YEAR > 2023), SITE_ID, Deployment))

# Get a simple list of all species in the truth file
truth_species_list <- names(truth_final)[-c(1, 299:302)] # Exclude 'Deployment' col and fwri_key_ref cols

#### depracated code to match up species differences between VIAME and truth species lists
# truthspec <- data.frame(truth_species_list)
# names(truthspec) <- c("Truth_Species")
# truthspec <- truthspec %>% dplyr::arrange(Truth_Species)
# viamespec <- data.frame(unique(p2dat$Species))
# names(viamespec) <- c("VIAME")
# viamespec <- viamespec %>% dplyr::arrange(VIAME)
# write.csv(viamespec, paste(outdir, "VIAME Species Check.csv"), row.names = F)
# write.csv(truthspec, paste(outdir, "Truth Species Check.csv"), row.names = F)


# pull in the files created in part 1
p1files <- paste0(outdir, "Part I - Counts/")
p1out <- list.files(p1files, pattern = "*.csv", full.names = T) 
p2dat <- p1out %>% 
  map_df(~suppressMessages(read_csv(., col_names = T, show_col_types = F)))
p2dat$Species <- trimws(p2dat$Species)

# Get AI metadata (the "universe" of deployments we care about)
deploydat <- unique(p2dat[c("Deployment", "year", "Version", "Confidence")])


# PREP SPECIES MATCHING FILE
spec_match <- read.csv(paste0(wrkdir, "Species List Match.csv")) %>% 
  dplyr::mutate(Species = trimws(VIAME), 
                Truth_Species = trimws(Truth)) %>% 
  dplyr::select(Species, Truth_Species) # Only need these two

# Join AI data with matching file 
# Use inner_join: we can only compare AI species that are in the match file.
p2dat.matched <- inner_join(p2dat, spec_match, by = "Species")

# Get a final, definitive list of all unique "Truth_Species" 
# from BOTH the AI and Manual files to loop through.
ai_species_list <- unique(p2dat.matched$Truth_Species)
all_species_to_compare <- unique(c(truth_species_list, ai_species_list))

# LOOP, JOIN, AND FIND FPs/FNs
allspec <- list()
print(paste0("Comparing ", length(all_species_to_compare), " species..."))

for (s in 1:length(all_species_to_compare)){
  spec <- all_species_to_compare[s]
  
  # --- AI Data for this species ---
  # Get all AI counts for this one species
  subdat_ai <- p2dat.matched %>%
    dplyr::filter(Truth_Species == spec) %>%
    dplyr::rename(VIAME_MaxN = MaxNCount) %>%
    # We only need these columns for the join
    dplyr::select(Deployment, Version, Confidence, Species, Truth_Species, VIAME_MaxN)
  
  # --- Manual Data for this species ---
  # Check if this species even exists in the truth file
  if (spec %in% names(truth_final)) {
    subdat_manual <- truth_final %>%
      dplyr::select(Deployment, all_of(spec)) %>%
      dplyr::rename(Manual = all_of(spec))
  } else {
    # If not, create a dummy tibble with 0s for this species
    subdat_manual <- tibble(Deployment = truth_final$Deployment, Manual = 0)
  }
  
  # --- Create the "Master Comparison" ---
  # We join the AI metadata (deploydat) with the manual counts.
  # This creates our "universe": one row for every deployment the AI ran,
  # populated with the correct manual count (which is 0 if not seen).
  # This is the "subtruth" object from your script, but created more safely.
  master_list <- inner_join(deploydat, subdat_manual, by = "Deployment")
  
  # --- Full Outer Join ---
  # This is the key: full_join keeps records from BOTH. This section of code basically ensures that we ALWAYS have columns for model version and confidence regardless of whether VIAME actually recorded any counts for this species
  # It keeps:
  #   1. All deployments from master_list (for FNs)
  #   2. All deployments from subdat_ai (for FPs)
  groundtruth <- full_join(subdat_ai, master_list, by = c("Deployment", "Version", "Confidence")) %>%
    dplyr::mutate(
      # Fill in NAs
      VIAME_MaxN = ifelse(is.na(VIAME_MaxN), 0, VIAME_MaxN),
      Manual = ifelse(is.na(Manual), 0, Manual),
      # Add species name back (it's NA for false negatives)
      Truth_Species = spec) %>%
    # We don't need the original "Species" name anymore
    dplyr::select(Deployment, year, Version, Confidence, Truth_Species, 
                  VIAME_MaxN, Manual)
  
  allspec[[s]] <- data.frame(groundtruth)
}

# 5. FINALIZE AND EXPORT
print("Finalizing master file...")
groundtruth_master <- do.call("rbind", allspec) %>%
  # filtering out all occurrences where both VIAME_MaxN and Manual counts were 0
  dplyr::filter(VIAME_MaxN != 0 | Manual != 0) %>% 
  dplyr::rename(Species = Truth_Species) %>%
  dplyr::mutate(Model_info = paste(year, Version, Confidence, sep = "_"))

# making sure that there is a unique row for each deployment, species, and model_info
groundtruth_master <- groundtruth_master %>%
  dplyr::distinct(Deployment, Species, Model_info, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(year)) # this is simply here in case the truth file is not up to date. If not, then there is no point in comparing to the manual counts so we will not need to proceed with that year.

# adding back in the info from the fwri_ref_key dataset
fwri_ref_key1 <- fwri_ref_key %>%
  dplyr::mutate(Deployment = ifelse(LAB %in% c("PC", "FWRI")| (!is.na(YEAR) & YEAR > 2023), SITE_ID, Deployment))
groundtruth_master <- merge(groundtruth_master, fwri_ref_key1, by = "Deployment", all.x = T)

# exporting master file
write.csv(groundtruth_master, 
          paste0(outdir, "Part II - Groundtruthing/Combined_Species_AllModels.csv"), 
          row.names = F)

print("Part II Complete!")




##################################################
############# Part III: Data Analysis ############
##################################################
# bringing in the above master script
combined_master <- read.csv(paste0(outdir, "Part II - Groundtruthing/Combined_Species_AllModels.csv")) %>% 
  dplyr::filter(!LAB %in% "FWRI") %>%
  dplyr::arrange(Species)

# important functions
# This ONE function calculates all standard metrics (Precision, Accuracy, Recall, FPR, FNR, and Ratios)
calculate_metrics <- function(df, species = "none") {
  
  # Define the core calculation logic.
  summarize_metrics <- function(.data) {
    .data %>% dplyr::summarize(
      # --- 1. FUNDAMENTAL COUNTS ---
      Agree = ifelse(Manual == VIAME_MaxN, 1, 0), # Agreement flag
      Difference = Manual-VIAME_MaxN, # Difference flag
      Relaxed = ifelse(Difference == 1 | Difference == -1 | Difference == 0, 1, 0), # Relaxed agreement flag for counts within 1/-1 of VIAME counts
      TP = sum(Manual > 0 & VIAME_MaxN > 0, na.rm = T), # True Positive
      FP = sum(Manual == 0 & VIAME_MaxN > 0, na.rm = T), # False Positive
      FN = sum(Manual > 0 & VIAME_MaxN == 0, na.rm = T), # False Negative
      TN = sum(Manual == 0 & VIAME_MaxN == 0, na.rm = T), # True Negative
      
      # --- 2. STANDARD METRICS (RATES) ---
      Precision = TP / (TP + FP),
      Recall_TPR = TP / (TP + FN), # True Positive Rate (Sensitivity)
      FPR = FP / (FP + TN),        # False Positive Rate
      FNR = FN / (TP + FN),        # False Negative Rate (1 - Recall)
      Accuracy = (TP + TN) / (TP + FP + FN + TN),
      
      # --- 3. RATIOS ---
      False_P_Ratio = FP / (TP + FN + TN),
      False_N_Ratio = FN / (TP + FP + TN),
      
      # --- 4. RAW COUNTS FOR CONTEXT ---
      Total_Actual_Positives = TP + FN,
      Total_Actual_Negatives = FP + TN,
      
      .groups = 'drop'
    )
  }
  
  # 1. Aggregate across all species
  if (species == "none") {
    df_out <- df %>%    
      dplyr::group_by(year, Version, Confidence) %>%
      summarize_metrics()
    
    # 2. Calculate for each species individually
  } else if (species == "all") {
    df_out <- df %>%    
      dplyr::group_by(year, Version, Confidence, Species) %>%
      summarize_metrics()
    
    # 3. Calculate for one specific species
  } else if (any(species %in% df$Species) == TRUE){
    df_out <- df %>%
      dplyr::filter(Species == species) %>%
      dplyr::group_by(year, Version, Confidence, Species) %>%
      summarize_metrics()
    
  } else {
    print("No species detected with that name. Check spelling and try again.")
    return(NULL) # Return NULL on failure
  }
  
  # Handle NaN (0/0) and Inf (X/0)
  df_out <- df_out %>%
    dplyr::mutate(across(where(is.numeric), ~ifelse(is.nan(.), NA, .))) %>%
    dplyr::mutate(across(where(is.numeric), ~ifelse(is.infinite(.), NA, .)))
  return(df_out)
}

# Create a function that will generate a percent agreement metric across deployments.
percent_metric <- function(df, variable1, variable2, group, metric){
  new.obj <- df %>%
    group_by({{variable1}}, {{variable2}}, {{group}}) %>%
    summarise(
      percentage_1s = mean({{metric}}) * 100,
      percentage_0s = (1 - mean({{metric}})) * 100,
      count = n())
  return(new.obj)
}


########## Processing steps ##########
# options for the function below:
  # species = "none" - gives you overall false positives/false negatives by model run, version, and confidence (not species-specific). DEFAULT - if you put nothing, this is what it will give you
  # species = "all" - gives you the false positives/false negatives for ALL species by model run, version, and confidence
  # species = "SPECIES_NAME" - gives you the false positives/false negatives for the specific species that you have chosen by model run, version, and confidence
metrics <- calculate_metrics(combined_master, species = "all")
percent_agreement <- percent_metric(metrics, year, Species, Confidence, Agree)
relaxed_agreement <- percent_metric(metrics, year, Species, Confidence, Relaxed)
percent_agreement$Agreement <- "Exact"
relaxed_agreement$Agreement <- "Relaxed"
PA.comp <- rbind(percent_agreement, relaxed_agreement)

# exporting
write.csv(metrics, paste0(outdir, "Part III - Data Analysis/Combined_False_Pos_Neg.csv"))



###############################################
########## Creating Analysis Reports ##########
###############################################
# Would you like to cut out large schools (i.e., 999, 299, 399)?
  print("Would you like to remove counts with large schools (i.e., 299, 399, and 999)? (y, n)")
  remove_large_schools <- rstudioapi::showPrompt(
    title = "Manual Input Required",
    message = "Would you like to remove counts with large schools (i.e., 299, 399, and 999)? Please enter y or n."
  )
  # checking to see if the user cancelled the value selection
  if (is.null(remove_large_schools)) {
    stop("Script cancelled by user.", call. = FALSE)
  }
  
  
out.runs <- list()
spec_runs <- list()
for (s in 1:length(unique(combined_master$Species))){
  # creating the analysis report directory for each model 
  analysrepdir <- paste0(script.dir, "/Output/Part III - Data Analysis/Analysis Reports/")
  dir.create(analysrepdir) # part III output directory (data analysis figures)
  
  # subsetting data to this model run
  species <- unique(combined_master$Species)[s]
  pa_spec <- PA.comp %>% dplyr::filter(Species == species)
  spec_pretty <- gsub("_", " ", str_to_sentence(species))
  spec_master <- combined_master %>% dplyr::filter(Species == species)
  sub_false <- metrics %>% dplyr::filter(Species == species)
  anrepmkd <- paste0(script.dir, "/VIAME Output Analysis Report.Rmd")
  rmarkdown::render(anrepmkd,
    output_dir = analysrepdir,
    output_format = "html_document",
    output_file = paste(spec_pretty, "Analysis Report.html"),
    params = list(
      spec_master = spec_master,
      sub_false = sub_false,
      species = species,
      outdir = outdir))
  
  # converting all confidences to a data frame for precision and difference
  spec_prec <- data.frame(do.call("rbind", precision))
  spec_diff <- data.frame(do.call("rbind", difference))
  if (sum(spec_diff$Frequency, na.rm = T) == 0){  
  spec_prec_diff <- merge(spec_prec, spec_diff, by = c("year", "Confidence"), all.x = T)
  spec_prec_diff <- spec_prec_diff %>% 
      dplyr::mutate(year = year, Confidence = Confidence, n = n, validn = validn, R = NA, PercAgree = NA, ASD = NA, ACV = NA, AAD = NA, APE = NA, Percent = NA) %>% dplyr::select(-Count, -Frequency)
  out.prec <- merge(mod_false, spec_prec_diff, by = c("year", "Confidence"), all.x = T) 
  } else {
  spec_diff <- spec_diff %>% 
    dplyr::filter(Count %in% c("-1", "0", "1")) %>%
    dplyr::group_by(year, Confidence) %>% 
    dplyr::summarise(Percent = sum(as.numeric(Frequency), na.rm = T))
  spec_prec_diff <- merge(spec_prec, spec_diff, by = c("year", "Confidence"), all.x = T)
  out.prec <- merge(mod_false, spec_prec_diff, by = c("year", "Confidence"), all.x = T)
  }
  out.runs[[s]] <- out.prec
}

spec_runs <- do.call("rbind", out.runs)
write.csv(spec_runs, paste0(outdir, "Part III - Data Analysis/", "Combined_Species_Prec.csv"))




