#### Annotation Library Tracking for VIAME CSV files
## Authors:Jack Prior, Sara Thomas, Kelsey Martin
# Jan. 21, 2026


# --- 1. SETUP & PACKAGE INSTALLATION ---
inst <- c("readr","dplyr","tidyverse","lubridate","data.table","microbenchmark",
          "tidyr","knitr","ggplot2", "ggpubr", "gridExtra", "FSA", "FSAdata")

for (p in inst) {
  if(!require(p, character.only = TRUE)) {
    install.packages(p)
    require(p, character.only = TRUE)
  }
}

# --- 2. DEFINE PATHS AND FIND FILES ---
## Naming your directories. These will be consistent throughout the script and will prevent each user from having to change paths for every file
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
data_root <- file.path(script.dir, "Tracks")
wrkdir  <- file.path(data_root, "/") # working directory (folder where data files are located)

# Check if directory exists before proceeding (Safety check)
if (!dir.exists(wrkdir)) {
  stop(paste("The directory", wrkdir, "does not exist in your working directory."))
}

# Get list of all .csv files recursively (in all subfolders)
file_list <- list.files(path = wrkdir, 
                        pattern = "\\.csv$", 
                        recursive = TRUE, 
                        full.names = TRUE)

message(paste("Found", length(file_list), "CSV files to process."))

# --- 3. DATA INGESTION & PROCESSING ---

# Function to read a single file and extract Column J (10th Column)
read_and_extract_col_j <- function(file_path) {
  
  # Read CSV without headers to ensure we strictly get the 10th column position
  # We read as character ('c') to avoid type conflicts between files
  tryCatch({
    data <- read_csv(file_path, skip=2, col_names = FALSE, col_types = cols(.default = "c"), show_col_types = FALSE)
    
    # Check if the file actually has a 10th column (Column J)
    if (ncol(data) >= 10) {
      col_j_data <- data[[10]]
      # create tibble
      clean_tbl <- tibble(
        source_file = basename(file_path),
        classification = col_j_data
      ) %>%
        # --- MODIFICATION: FILTER BLANKS ---
        # 1. Remove NAs
        # 2. Remove empty strings ""
        # 3. Remove strings that are just spaces " " using trimws()
        filter(!is.na(classification), 
               trimws(classification) != "")
      
      return(clean_tbl)
    } else {
      warning(paste("File skipped (fewer than 10 columns):", basename(file_path)))
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("Error reading file:", basename(file_path)))
    return(NULL)
  })
}

# Use purrr::map_dfr to loop through files and combine them into one dataframe
all_data <- map_dfr(file_list, read_and_extract_col_j) %>%
  dplyr::mutate(station_name = gsub("_|-", "", substr(source_file, 0, 12)))

# --- 4. ANALYSIS & OUTPUT ---

if (nrow(all_data) > 0) {
  
  # A. Count per Classification per File
  counts_per_file <- all_data %>%
    group_by(source_file, classification, station_name) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(source_file, desc(count))
  
  # B. Total Count across ALL files
  counts_total <- all_data %>%
    group_by(classification) %>%
    summarise(total_count = n(), .groups = "drop") %>%
    arrange(desc(total_count))
  
  # --- 5. DISPLAY RESULTS ---
  
  print("--- SUMMARY: Counts Per File ---")
  print(kable(head(counts_per_file, 20), caption = "Top 20 Counts (Per File)"))
  
  print("--- SUMMARY: Total Counts Across All Files ---")
  print(kable(counts_total, caption = "Total Classification Counts"))
  
  # Optional: Visualizing the top 10 classifications totals
  plot_summary <- counts_total %>%
    slice_max(total_count, n = 10) %>%
    ggplot(aes(x = reorder(classification, total_count), y = total_count)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 10 Classifications (Column J)", x = "Classification", y = "Count")
  
  print(plot_summary)
  
} else {
  message("No data found in Column J across the provided files.")
}
# --- 6. EXPORT TO EXCEL ---

# distinct package for Excel writing (not in original list)
excel_pkg <- "openxlsx"
if(!require(excel_pkg, character.only = TRUE)) {
  install.packages(excel_pkg)
  require(excel_pkg, character.only = TRUE)
}

# Define your filename
file_name <- paste0("Annotation Library Tracking ", Sys.Date(), ".xlsx")

# Create the full path
full_path <- file.path(script.dir, file_name)

# Create a new Workbook object
wb <- createWorkbook()

# Add Sheet 1: Total Counts (The high-level summary)
addWorksheet(wb, "Total Counts")
writeData(wb, "Total Counts", counts_total)

# Add Sheet 2: Counts by File (The granular detail)
addWorksheet(wb, "Counts by File")
writeData(wb, "Counts by File", counts_per_file)

# Save the workbook
saveWorkbook(wb, full_path, overwrite = TRUE)

message(paste("Successfully exported results to Library Analysis"))
