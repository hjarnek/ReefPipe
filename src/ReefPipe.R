#!/usr/bin/env Rscript

#############
## MESSAGE ##
#############

cat(' ____  _____ _____ _____ ____ ___ ____  _____
|  _ \\| ____| ____|  ___|  _ |_ _|  _ \\| ____|
| |_) |  _| |  _| | |_  | |_) | || |_) |  _|   
|  _ <| |___| |___|  _| |  __/| ||  __/| |___  
|_| \\_|_____|_____|_|   |_|  |___|_|   |_____|\n\n')


###############
## FUNCTIONS ##
###############

# Print header
print_header <- function(number){
  headers <- c("\n _ _ _     \n| (_) |__  \n| | | \'_ \\ \n| | | |_) |\n|_|_|_.__/ \n\nInstalling and loading all packages\n\n",                                                      # lib
               " _ __  _ __ ___ _ __\n| \'_ \\| \'__/ _ \\ \'_ \\\n| |_) | | |  __/ |_) |\n| .__/|_|  \\___| .__/\n| |            | |\n|_|            |_|\n\nSetup & checks\n\n",     # prep
               "\n  __ _ _____   __\n / _` / __\\ \\ / /\n| (_| \\__ \\\\ V / \n \\__,_|___/ \\_/  \n",                                                                               # asv
               "\n      _ _           \n     | (_)          \n  ___| |_ _ __ ___  \n / _ \\ | | '_ ` _ \\ \n|  __/ | | | | | | |\n \\___|_|_|_| |_| |_|\n\n",                         # elim
               "\n       _       _   \n      | |     | |  \n _ __ | | ___ | |_ \n| '_ \\| |/ _ \\| __|\n| |_) | | (_) | |_ \n| .__/|_|\\___/ \\__|\n| |                \n|_|\n\n",    # plot
               " _             \n| |_ __ ___  __\n| __/ _` \\ \\/ /\n| || (_| |>  < \n \\__\\__,_/_/\\_\\\n\n")                                                                       # tax
  
  cat(headers[number])
}

# Download and load R packages
install_and_load_packages <- function(packages) {
  installed_packages <- installed.packages()[, 'Package']                   # Summarize all installed packages
  
  for (package in packages) {
    if (!(package %in% installed_packages)) {                               # If package not installed
      if (package %in% c('dada2', 'Biostrings', 'ShortRead')) {     
        options("BioC_mirror" = "https://bioconductor.org")                 # Set the Bioconductor image
        BiocManager::install(package)                                       # Install package with BiocManager
      } else {                                                      
        options(repos = "https://cran.rstudio.com/")                        # Set the cran mirror
        install.packages(package)                                           # Install package
      }
    }
    
    suppressPackageStartupMessages(library(package, character.only = TRUE)) # Silently load package
  }
}


####################################
## Parsing command line arguments ##
####################################

install_and_load_packages('argparse')                                       # Install and load command line argument parser

parser <- ArgumentParser(description = 'Reefpipe command line arguments')   # Initialize command line argument parser

# General command line arguments
parser$add_argument('-b', '--base_dir', metavar = 'BASEDIR', type = 'character', required = TRUE, help = 'Specify the base directory path for the analysis.')
parser$add_argument('-d', '--download', metavar = 'FILENAME', required = FALSE, help = 'Download data from ENA for the given accessions listed in FILENAME.')
parser$add_argument('-r', '--run_mode', choices = c('single', 'multi'), required = TRUE, help = 'Specify the run mode: \'single\' for analyzing a single sequencing run or \'multi\' for analyzing multiple runs.')
parser$add_argument('-g', '--gene', choices = c('COI', 'ITS', '18S'), required = TRUE, help = 'Specify the gene to analyze.')

# Primer removal arguments
parser$add_argument('-t', '--trim_primers', action = 'store_true', help = 'Enable primer trimming using cutadapt. By default, primers are not trimmed.')
parser$add_argument('-p', '--primers', nargs=2, type='character', metavar=c('Fwd_Primer', 'Rev_Primer'), help='Specify the forward and reverse primer sequences for trimming with cutadapt.')
parser$add_argument('-e', '--error_rate', type='numeric', default = 0.1, help = 'Specify the maximum allowed error rate for cutadapt (if 0 <= E < 1) or the absolute number of errors for full-length adapter match (if E is an integer >= 1)')

# Command line arguments for filtering and trimming
parser$add_argument('-m', '--minlen', type='integer', metavar='minlen', default=50, help='Specify the minimum length threshold for trimmed reads. Default is 50.')
parser$add_argument('-l', '--trunclen', nargs=2, type='integer', metavar=c('Fwd', 'Rev'), default=c(200,140), help='Specify the maximum length for trimmed reads. Default is 200 bases for forward reads and 140 bases for reverse reads.')
parser$add_argument('-n', '--max_ambiguous', type='integer', default=0, help='Specify the maximum number of ambiguous reads allowed. Default is 0.')
parser$add_argument('-E', '--max_error_rates', nargs=2, type='numeric', metavar=c('Fwd', 'Rev'), default=c(2, 4), help='Specify the maximum error rates for forward and reverse reads. Default is 2 for forward reads and 4 for reverse reads.')
parser$add_argument('-q', '--min_quality_score', type='integer', default=2, help='Specify the minimum quality score that each base should have. Default is 2.')
parser$add_argument('-x', '--contaminants', action='store_true', help='Disable the removal of contaminant reads from the PhiX genome during filtering and trimming.')
parser$add_argument('-u', '--uncompressed', action='store_true', help='Disable the compression of the output files (yields non-gzipped fastq files).')

# Command line arguments for sample inference
parser$add_argument('-N', '--nbases', type='integer', default=1e+08, help='The minimum number of total bases to use for error rate learning during sample inference.')
parser$add_argument('-o', '--min_overlap', type='integer', default=10, help='Specify the minimum overlap length required for merging pairs. Default is 10.')
parser$add_argument('-i', '--max_mismatch', type='integer', default=1, help='Specify the maximum number of mismatches allowed during merging. Default is 1.')

# ASV fine-tuning
parser$add_argument('-s', '--singletons', action = 'store_true', help = 'Keep singletons in the sequence table. By default, singletons are removed (except when only one sample is analyzed).')

# Command line arguments for taxonomic classification
parser$add_argument('-B', '--BOLD', nargs=2, type='character', metavar=c('USERNAME', 'PASSWORD'), help = 'Perform taxonomic classification using BOLD. Specify BOLDsystems.org username and password.')
parser$add_argument('-S', '--batch_size', type = 'numeric', required = FALSE, default = 50, help = 'Specify the BOLD batch size.')
parser$add_argument('-R', '--reference', action = 'store_true', help = 'Perform taxonomic classification with RDPclassifier using custom reference databases.')
parser$add_argument('-M', '--minBoot', type = 'numeric', default = 80, help = 'Specify the minimal bootstrap value for taxonomic classification with RDPclassifier. Default is 80.')

# Command line arguments for taxonomic table fusing
parser$add_argument('-F', '--fuse', action = 'store_true', help = 'Fuse the information from all taxonomic tables.')
parser$add_argument('-f', '--fuseLevels', type = 'character', default = 'Phylum,Class,Order,Family,Genus,Species', help = 'Specify the taxonomic levels to be used for fusing all taxonomic tables. Default levels are Phylum,Class,Order,Family,Genus,Species.')

# Parse the arguments
args <- parser$parse_args()

# Access the argument values
mainpath <- normalizePath(args$base_dir, winslash = '/')
download <- args$download
run_mode <- args$run_mode
GOI <- args$gene
trim_primers <- args$trim_primers
primers <- args$primers
cutadapt_error_rate <- args$error_rate
minlen <- args$minlen
trunclen <- args$trunclen
singletons <- args$singletons
bold <- args$BOLD
reference <- args$reference
minBoot <- args$minBoot
fuse <- args$fuse
fuseLevels <- args$fuseLevels
max_ambiguous <- args$max_ambiguous
max_error_rates <- args$max_error_rates
min_quality_score <- args$min_quality_score
contaminants <- args$contaminants
uncompressed <- args$uncompressed
nbases <- args$nbases
min_overlap <- args$min_overlap
max_mismatch <- args$max_mismatch
batch_size <- args$batch_size


########################
## MAIN PIPELINE PATH ##
########################

args <- commandArgs()
pipeline_path <- NULL
for (arg in args) {
  if (startsWith(arg, "--file=")) {
    pipeline_path <- sub("^--file=", "", arg)
    break
  }
}

if (!is.null(pipeline_path)) {
  pipeline_path <- normalizePath(pipeline_path, winslash = '/')
} else {
  stop("Pipeline path not found.")
}


#################
## R PACKAGES ##
################

print_header(1)

# Install and load
install_and_load_packages(c('BiocManager','openxlsx','dada2','ggplot2','stats', 'Biostrings','ShortRead','vegan','readxl','stringr','purrr','dplyr','gtools','gplots','tidyr'))


########################################
## COMMAND LINE CONDITIONS AND CHECKS ##
########################################

print_header(2)
cat('Checking command parameters... ')

  ##########################
  ## BASE DIRECTORY CHECK ##
  ##########################
  
  file_info <- file.info(mainpath)
  if(is.na(file_info$isdir)){
    stop(paste(mainpath, 'does not exist.'))
  } else if(!file_info$isdir){
    stop(paste(mainpath, 'is not a directory.'))
  }
  
  ################################
  ## TRIMMING AND PRIMERS CHECK ##
  ################################
  
  if(trim_primers) {
    
    # Define forward and reverse primer sequences
    if(!is.null(primers)) {
      FWD <- primers[1] 
      REV <- primers[2]
    } else {
      default_primers <- read.csv2(file.path(dirname(dirname(pipeline_path)), 'default_primers.csv'))
      FWD <- default_primers[which(default_primers[,"Region_name"] == GOI), "Fwd_primer"]
      REV <- default_primers[which(default_primers[,"Region_name"] == GOI), "Rev_primer"]
    }
    
    # Check that primer sequences only contain valid characters
    primercheck <- c(grepl("^[ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]+$", FWD), grepl("^[ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn]+$", REV))
    if(primercheck[1] == F & primercheck[2] == F) {
      stop(paste0('Both forward and reverse primers contain invalid characters. If no primers were specified, check the primers given in the default_primers.csv file, and that the corresponding Region_name agrees with the specified gene of interest (currently \"', GOI, '\")'))
    } else if(primercheck[1] == F & primercheck[2] == T) {
      stop(paste0('Forward primer contains invalid characters. If no primers were specified, check the primers given in the default_primers.csv file, and that the corresponding Region_name agrees with the specified gene of interest (currently \"', GOI, '\")'))
    } else if(primercheck[1] == T & primercheck[2] == F) {
      stop(paste0('Reverse primer contains invalid characters. If no primers were specified, check the primers given in the default_primers.csv file, and that the corresponding Region_name agrees with the specified gene of interest (currently \"', GOI, '\")'))
    }
  } else {
    FWD <- character()
    REV <- character()
  }

  if(GOI %in% c('ITS', '18S')){
    trunclen = 0
    cat('trunclen has been disabled.\n')
  }
    
  ###############################
  ## REFERENCE DATABASE CHECKS ##
  ###############################
  
  if(reference){
    
    ref_files <- list.files(path = file.path(dirname(dirname(pipeline_path)), 'reference'), full.names = T)  # Get list of all references
    config_file <- ref_files[basename(ref_files) == 'config.txt']                                            # Get the configuration file
    ref_files <- ref_files[ref_files != config_file]                                                         # Remove configuration file from reference list
    
    # Check if reference databases and config file exist when --reference is selected
    if(!file.exists(config_file)){
      stop(paste("The file", config_file,
                 "does not exist, which is necessary for the process of taxonomic classification using locally stored reference databases. It appears that the file may have been (re)moved or is missing from its expected location."))
    }
    if(length(ref_files) == 0){
      stop('The reference database folder is empty. Please make sure there are reference databases available in the directory.')
    }
    
    config_df <- read.csv(config_file)                                                                            # Read in the configuration file
    config_refs <- config_df[,'Database']                                                                         # Get the database names
    
    # Loop over every database in the reference directory
    for(ref_file in ref_files){
      
      # Check if the names of the reference databases match with those in the config file
      if(!basename(ref_file) %in% config_refs){
        stop(paste('The reference database config file is missing information or contains a typo related to the file', basename(ref_file), '. Please review the config file and try again.'))
      } 
      
      # Check if the names of the reference databases occur only once in the config file
      else if(sum(grepl(pattern = basename(ref_file), x = config_refs)) > 1){
        stop(paste('The reference database config file has multiple entries for the', basename(ref_file), 'database. Please review the config file and try again.'))
      }
    }
    
    #########################
    ## TABLE FUSION CHECKS ##
    #########################
    
    # Check if --fuseLevels is valid
    if(fuse){
      
      # Extract the individual levels from --fuseLevels
      individual_levels <- gsub(' ', '', str_to_title(strsplit(fuseLevels, ",")[[1]]))
      
      for(i in 1:nrow(config_df)){
        
        # Extract the taxonomic levels per reference database from config_df and convert them to 1 string
        ref_levels <- gsub('[\t ]', '', str_to_title(strsplit(config_df[i, "TaxonomicLevel"], ";")[[1]]))
        
        for(individual_level in individual_levels){
          if(!individual_level %in% ref_levels){
            stop(paste('Invalid taxonomic levels for --f argument: The level of', 
                       individual_level, 
                       'cannot be used for merging taxonomic tables as it is not present in all reference databases.'))
          }
        }
      }
    }
  # Check if fusing taxonomic tables is possible
  } else if((bold == T) & (fuse == T)){
    warning('You cannot merge taxonomic tables if only the BOLDSYSTEMS table is generated.')
    fuse = F    # Change fuse to false
  }
  
  cat('OK\n')
  
  #################
  ## BOLD CHECKS ##
  #################
  
  if(!is.null(bold)) {
    username <- bold[1]
    password <- bold[2]
    bold <- TRUE
  } else {
    bold <- FALSE
  }
  
  if(bold){
    if(GOI == '18S'){
      bold = F
      cat('BOLD has been disabled.\n')
    } else {
      
      cat('Checking BOLD access... ')
      
      # Check if the BOLD batch_size is between 1 and 50
      if(batch_size < 1 | batch_size > 50){
        stop('The BOLD batch size must be between 0 and 50.') 
      }
      
      # Check if connection to BOLDSYSTEMS is possible 
      loginresult <- system2(command = 'python3', args = c(paste0("\"", file.path(dirname(pipeline_path), 'BOLD_logincheck.py'), "\""), 
                                                           '-u', username, 
                                                           '-p', password),
                             stdout = TRUE)
      if('Login successful.' %in% loginresult){
        cat('OK')
      } else {
        cat(loginresult)
        stop('Could not connect to BOLDSYSTEMS.')
      }
    }
  }


#############################
## DOWNLOAD ENA ACCESSIONS ##
#############################

if(!is.null(download)){
  
  cat('\n\nFetching fastq files from ENA\n')
  
  # Source the R-script
  source(file.path(dirname(pipeline_path), 'ENAdownload.R'))
  
  ENAdownload(download)
  
  run_mode = 'multi'  # Set run_mode to multi
}


#########################################
## METABARCODE DATA BASE FOLDER CHECKS ##
#########################################

cat('\nChecking metabarcoding data...')

  ###############################
  ## LOCATE METABARCODING DATA ##
  ###############################
  
  # Determine location metabarcoding data folders
  if (run_mode == 'single' & is.null(download)){
    metabarcode_dirs <- mainpath
  } else if(run_mode == 'multi'){
    # Get all subdirectories in the base directory
    metabarcode_dirs <- normalizePath(list.dirs(path = mainpath, full.names = TRUE, recursive = FALSE), winslash = '/')
    
    # Remove any *_trial[0-9]+$ directories
    metabarcode_dirs <- metabarcode_dirs[!grepl(pattern = '_trial[0-9]+$', x = metabarcode_dirs)]
  }
  
  
  ##############################
  ## CHECK METABARCODING DATA ##
  ##############################
  
    ########################
    ## MULTI RUN ANALYSIS ##
    ########################
    
    if(run_mode == 'multi'){
        
      # Iterate over all metabarcode folders
      for(metabarcode_dir in metabarcode_dirs){
        
        # List contents of metabarcode folder
        files <- list.files(path = metabarcode_dir, full.names = TRUE)
        is_match <- grepl(pattern = "_1.fastq$|_2.fastq$|_1.fastq.gz$|_2.fastq.gz$", x = files) 
        
        # Check if metabarcode folder contains non-fastq files
        if(any(!is_match)){
          cat('\n')
          stop(paste("The sequencing run folders include files with prefixes that do not match the patterns '_1.fastq', '_2.fastq', '_1.fastq.gz', or '_2.fastq.gz.'\n",
                     "      Please check the following directories and remove any files that do not match the aforementioned patterns:\n",
                     paste("     ", metabarcode_dir)))
        }
      }
    } else if(run_mode == 'single'){
    
      
    #########################
    ## SINGLE RUN ANALYSIS ##
    #########################
    
      # List contents of metabarcode folder (i.e. base folder)
      files <- list.files(path = metabarcode_dirs, full.names = TRUE)
      is_match <- grepl(pattern = "_1.fastq$|_2.fastq$|_1.fastq.gz$|_2.fastq.gz$", x = files)
      
      # Check presence of fastq files
      if(!any(is_match)){
        stop(paste('No metabarcoding data could be found in', mainpath,
                   "\n       If the metabarcoding data is in a subfolder of the base directory, use '-r multi' instead of '-r single' for single sequencing run analysis.\n"))
      }
    }

##################
## TRIAL FOLDER ##
##################

  #########################
  ## CREATE TRIAL FOLDER ##
  #########################
  
  foldernumber <- 1
  
  while(TRUE){
    
    # Generate trial folder name
    trialpath <- file.path(mainpath, paste0(basename(mainpath), '_trial', foldernumber))
    
    # Check if the trial folder exists
    if(dir.exists(trialpath)){
      foldernumber <- foldernumber + 1
      next
    }
    
    else{
      
      # Create the trial folder
      dir.create(trialpath)
      
      # Copy the metabarcode data in base folder to trial folder
      metabarcode_content <- list.files(mainpath, full.names = T)
      metabarcode_content <- metabarcode_content[!grepl(pattern = '_trial[0-9]+$', x = metabarcode_content)]
      
      file.copy(metabarcode_content, 
                trialpath, 
                recursive = TRUE)
      
      # Reconstruct the command
      command_arguments <- commandArgs(trailingOnly = T)
      
      if(bold){
        command_arguments[command_arguments == username] <- "username_is_confidential"
        command_arguments[command_arguments == password] <- "password_is_confidential"
      }
      
      command <- paste('Rscript ReefPipe.R', paste(command_arguments, collapse = ' '))
      
      # Construct a parameters file
      parameters <- c(
        'Command\n-------',
        command,
        
        '\n\nGene\n----',
        paste0("GOI:", GOI),
        
        '\n\nPrimer Removal\n--------------',
        paste0("trim_primers:", trim_primers),
        paste0("Fwd_Primer:", FWD),
        paste0("Rev_Primer:", REV),
        
        '\n\nFiltering & Trimming\n--------------------',
        paste0("minlen:", minlen),
        paste0("trunclen:", trunclen),
        paste0("max_ambiguous:", max_ambiguous),
        paste0("max_error_rates:", max_error_rates),
        paste0("min_quality_score:", min_quality_score),
        paste0("contaminants:", contaminants),
        paste0("uncompressed:", uncompressed),
        
        '\n\nSample inference\n------------------',
        paste0("nbases:", nbases),
        paste0("min_overlap:", min_overlap),
        paste0("max_mismatch:", max_mismatch),
        
        '\n\nASV Fine-Tuning\n---------------',
        paste0("singletons:", singletons),
        
        '\n\nTaxonomic classification\n------------------------',
        paste0("BOLD:", bold),
        paste0("reference:", reference),
        paste0("minBoot:", minBoot),
        
        '\n\nTaxonomic Table fusion\n----------------------',
        paste0("fuse:", fuse),
        paste0("fuseLevels:", fuseLevels)
      )
      
      writeLines(parameters, con = file.path(trialpath, 'parameters.txt'))
  
      break
    }
  }
  
  ###############################################
  ## LOCATE METABARCODING DATA IN TRIAL FOLDER ##
  ###############################################
  
  if (run_mode == 'single'){
    paths = trialpath
  } else if(run_mode == 'multi'){
    paths = normalizePath(list.dirs(path = trialpath, full.names = TRUE, recursive = FALSE), winslash = '/')
  }
  
  cat(' OK\n\n')

#############################################
## ASV GENERATION FOR EVERY SEQUENCING RUN ##
#############################################

print_header(3)

# Construct comprehensive sequence table
main.seqtab <- NULL

# Construct a comprehensive read track table
main.track <- NULL

for(iter in 1:length(paths)){
  
  
  ##################################
  ## ITERATION AND STEP VARIABLES ##
  ##################################
  
  # Iteration message
  cat(paste0('\nIteration ', iter, ' out of ', length(paths), ': ', basename(paths[iter])))
  
  # Iteration label
  label <- paste0('\n[', iter, '/', length(paths), ' - Step ')
  
  # Step counter
  step <- 0
  
  
  ##########################
  ## FETCHING FASTQ FILES ##
  ##########################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Fetching fastq file paths\n'))
  
  # Forward and reverse read paths
  FwdRead <- sort(list.files(paths[iter], pattern="_1.fastq", full.names = TRUE))
  RevRead <- sort(list.files(paths[iter], pattern="_2.fastq", full.names = TRUE))
  
  # Extract sample names
  sample.names <- sapply(strsplit(basename(FwdRead), split = "_(?=[^_]*$)", perl = TRUE), `[`, 1)
  sample.names.check <- sapply(strsplit(basename(RevRead), split = "_(?=[^_]*$)", perl = TRUE), `[`, 1)
  
  # Check if forward and reverse sample names match
  if (length(sample.names)!= length(sample.names.check)){
    stop('Amount of forward and reverse files do not match!')
  } else {
    for (i in 1:length(sample.names)){
      if(sample.names[i] != sample.names.check[i]){
        stop(paste0('The forward read file of ', 
                    sample.names[i], 
                    ' matched with the reverse read file of ', 
                    sample.names.check[i], '. Something in the main sample folder is wrong and needs to be fixed first.'))
      }
    }
  }
  
  
  ##################################
  ## REMOVE PRIMERS WITH CUTADAPT ##
  ##################################
  
  if (trim_primers == T){
    
    # Message
    step <- step + 1
    cat(paste0(label, step, '] Removing primers with cutadapt and prefiltering the reads\n'))
    
    # Make a new directory to store prefiltered sequences
    path.cut <- file.path(paths[iter], '01.Prefiltered')
    
    if(!dir.exists(path.cut)){
      cat(paste('Creating output directory:', path.cut, '\n'))
      dir.create(path.cut)
    }
    
    # Make files to store the prefiltered sequences
    FwdRead.cut <- file.path(path.cut, basename(FwdRead))
    RevRead.cut <- file.path(path.cut, basename(RevRead))
    
    # If the gene of interest is COI
    if(GOI == 'COI'){
      FWD.argument <- paste0('-g', ' ^', FWD)
      REV.argument <- paste0('-G', ' ^', REV)
      
      # Run cutadapt
      for(i in seq_along(FwdRead)){
        system2('cutadapt', args = c(FWD.argument,                                                                      # Define the forward primer
                                     REV.argument,                                                                      # Define the reverse primer
                                     '-m 1',                                                                            # Only keep reads with a minimal length of 1,
                                     '-e', cutadapt_error_rate,
                                     '--discard-untrimmed',                                                             # Discard reads that were not trimmed
                                     '-o', paste0("\"", FwdRead.cut[i], "\""), '-p', paste0("\"",RevRead.cut[i], "\""), # Output files
                                     paste0("\"", FwdRead[i], "\""), paste0("\"", RevRead[i], "\""),                    # Input files
                                     '-j 0'))                                                                           # Number of cores to run on. 0 for auto-detect
      }
    }
    
    # If the gene of interest is ITS or 18S rRNA
    else if (GOI %in% c('ITS', '18S')){
      FWD.RC <- dada2::rc(FWD)
      REV.RC <- dada2::rc(REV)
      
      FWD.argument <- paste0('-g', ' ^', FWD, ' -a ', REV.RC)
      REV.argument <- paste0('-G', ' ^', REV, ' -A ', FWD.RC)
      
      # Run cutadapt
      for(i in seq_along(FwdRead)){
        system2('cutadapt', args = c(FWD.argument,                                                                      # Define the forward arguments
                                     REV.argument,                                                                      # Define the reverse arguments
                                     '-m 1',                                                                            # Only keep reads with a minimal length of 1,
                                     '-n 2',                                                                            # Required to remove FWD and REV from reads
                                     '-e', cutadapt_error_rate,
                                     '--discard-untrimmed',                                                             # Discard reads that were not trimmed
                                     '-o', paste0("\"", FwdRead.cut[i], "\""), '-p', paste0("\"",RevRead.cut[i], "\""), # Output files
                                     paste0("\"", FwdRead[i], "\""), paste0("\"", RevRead[i], "\""),                    # Input files
                                     '-j 0'))                                                                           # Number of cores to run on. 0 for auto-detect
      }
    }
  }
  
  
  ##################################
  ## FILTER READS BASED ON LENGTH ##
  ##################################
  
  if (trim_primers == F){
    
    # Message
    step <- step + 1
    cat(paste0(label, step, '] Prefiltering the reads\n'))
    
    # Make directory path to store prefiltered sequences
    path.cut <- file.path(paths[iter], '01.Prefiltered')
    
    # Make a new directory to store prefiltered sequences
    if(!dir.exists(path.cut)){
      cat(paste('Creating output directory:', path.cut, '\n'))
      dir.create(path.cut)
    }
    
    # Make files to store the prefiltered sequences
    FwdRead.cut <- file.path(path.cut, basename(FwdRead))
    RevRead.cut <- file.path(path.cut, basename(RevRead))
    
    length_filtered <- filterAndTrim(FwdRead, FwdRead.cut, RevRead, RevRead.cut, minLen = 1, multithread = TRUE)
    
    # Get the path of the forward and reverse trimmed fastq files
    FwdRead.cut.check <- sort(list.files(path.cut, pattern="_1.fastq", full.names = TRUE))
    RevRead.cut.check <- sort(list.files(path.cut, pattern="_2.fastq", full.names = TRUE))
    
    # Check if forward and reverse files match
    if(length(FwdRead.cut.check) != length(RevRead.cut.check)){
      stop('Something went wrong! Forward and reverse files do not match anymore!\n')
    }
  }
  
  
  ###################################
  ## INSPECT READ QUALITY PROFILES ##
  ###################################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Plotting quality of prefiltered reads\n'))
  
  # Make PDF with read quality profiles
  suppressWarnings({
    
    pdf(file = file.path(path.cut, 'QualityProfile.pdf'))
      
      # Forward reads
      if(length(FwdRead.cut) <= 10){
        
        show(plotQualityProfile(FwdRead.cut) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      } else {
        
        show(plotQualityProfile(FwdRead.cut[1:10]) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      }
      
      # Reverse reads
      if(length(RevRead.cut) <= 10){
        
        show(plotQualityProfile(RevRead.cut) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      } else {
        
        show(plotQualityProfile(RevRead.cut[1:10]) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      }
    
    dev.off()
  })
  
  
  #####################
  ## FILTER AND TRIM ##
  #####################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Trimming the prefiltered reads\n'))
  
  # Make directory path to store prefiltered sequences
  path.filt <- file.path(paths[iter], '02.Filtered_Trimmed')
  
  # Make files to store the trimmed sequences
  FwdRead.filt <- file.path(path.filt, paste0(sample.names, "_1_filt.fastq.gz"))
  RevRead.filt <- file.path(path.filt, paste0(sample.names, "_2_filt.fastq.gz"))
  
  names(FwdRead.filt) <- sample.names
  names(RevRead.filt) <- sample.names
  
  # Trim the prefiltered reads
  out <- filterAndTrim(FwdRead.cut,                            # Input files forward reads
                       FwdRead.filt,                           # Output files forward reads
                       RevRead.cut,                            # Input files reverse reads
                       RevRead.filt,                           # Output files reverse reads
                       truncLen= trunclen,                     # Max length of the reads (F will be truncated to 200, reverse to 140)
                       maxN= max_ambiguous,                    # Amount of ambiguous bases per read that are allowed
                       maxEE= max_error_rates,                 # Maximum error rates for F and R read
                       truncQ= min_quality_score,              # Minimum quality score that each base should have
                       minLen = minlen,                        # Minimum length of the reads after trimming
                       rm.phix= !contaminants,                 # Remove contaminant reads of the PhiX genome
                       compress= !uncompressed,                # Output files are compressed
                       multithread = T)                        # On Windows set multithread = FALSE
  
  saveRDS(out, file.path(path.filt, 'Filtered_Trimmed_Logfile.rds'))
  
  ####################################################
  ## INSPECT READ QUALITY PROFILES OF TRIMMED READS ##
  ####################################################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Plotting quality of trimmed reads\n'))
  
  # Make PDF with read quality profiles
  suppressWarnings({
    pdf(file = file.path(path.filt, 'QualityProfilesFilteredTrimmed.pdf'))
    
      # Forward reads
      if(length(FwdRead.filt) <= 10){
        
        show(plotQualityProfile(FwdRead.filt) +
          geom_hline(yintercept = 30) +
            scale_color_manual(guide = "none"))
        
      } else {
        
        show(plotQualityProfile(FwdRead.filt[1:10]) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      }
      
      # Reverse reads
      if(length(RevRead.filt) <= 10){
        
        show(plotQualityProfile(RevRead.filt) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      } else {
        
        show(plotQualityProfile(RevRead.filt[1:10]) +
          geom_hline(yintercept = 30) +
          scale_color_manual(guide = "none"))
        
      }
    
    dev.off()
    
  })
  
  
  ###########################
  ## LEARN THE ERROR RATES ##
  ###########################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Learning error rates\n'))
  
  # Make directory path to store the error plots
  path.error <- file.path(paths[iter], '03.Error_Rates')
  
  if(!dir.exists(path.error)){
    cat(paste('Creating output directory:', path.error,'\n'))
    dir.create(path.error)
  }
  
  # Learn the error rates
  set.seed(100)
  
  errF <- learnErrors(FwdRead.filt, nbases = nbases, multithread=TRUE)
  errR <- learnErrors(RevRead.filt, nbases = nbases, multithread=TRUE)
  
  # Construct and store the error plots
  suppressWarnings({
    pdf(file = file.path(path.error, paste0('errorRates', '.pdf')))
      
      # Error plots for forward reads
      show(plotErrors(errF, nominalQ=TRUE) +
             ggtitle('Forward') + 
             theme(plot.title = element_text(hjust = 0.5)))
      
      # Error plots for reverse reads
      show(plotErrors(errR, nominalQ=TRUE)+
             ggtitle('Reverse') + 
             theme(plot.title = element_text(hjust = 0.5)))
    
    dev.off()
  })
  
    
  ######################
  ## SAMPLE INFERENCE ##
  ######################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Inferring sample\n'))
  
  # Make directory path to store the dada objects
  path.infer <- file.path(paths[iter], '04.Sample_Inference')
  
  if(!dir.exists(path.infer)){
    cat(paste('Creating output directory:', path.infer, '\n'))
    dir.create(path.infer)
  }
  
  # Infer the samples
  dadaFwd <- dada(FwdRead.filt, err=errF, multithread=TRUE, pool = "pseudo")
  cat('\n')
  dadaRev <- dada(RevRead.filt, err=errR, multithread=TRUE, pool = "pseudo")
  cat('\n')
  
  # Store the dada objects
  saveRDS(dadaFwd, file.path(path.infer, 'dadaFwd.rds'))
  saveRDS(dadaRev, file.path(path.infer, 'dadaRev.rds'))
  

  ########################
  ## MERGE PAIRED READS ##
  ########################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Merging paired reads\n'))
  
  # Make directory path to store the merger object
  path.merge <- file.path(paths[iter], '05.Merged_Reads')
  
  if(!dir.exists(path.merge)){
    cat(paste('Creating output directory:', path.merge, '\n'))
    dir.create(path.merge)
  }
  
  # Merge the paired reads
  mergers <- mergePairs(dadaFwd, 
                        FwdRead.filt, 
                        dadaRev, 
                        RevRead.filt, 
                        minOverlap = min_overlap, 
                        maxMismatch = max_mismatch, 
                        verbose=T)

  # Write mergers object to an rds file
  saveRDS(mergers, file.path(path.merge, 'mergers.rds'))
  
  
  ##############################
  ## CONSTRUCT SEQUENCE TABLE ##
  ##############################
  
  # Message
  step <- step + 1
  cat(paste0(label, step, '] Constructing sequence tables\n'))
  
  # Construct the sequence table
  
  seqtab <- makeSequenceTable(mergers)
  
  
  ##########################
  ## CHECK SEQUENCE TABLE ##
  ##########################
  
  # If there are no sequences in the sequence table, continue
  if(dim(seqtab)[2] == 0){
    
    cat(paste('\nASVS could not be generated for', basename(paths[iter]), 
              '\nExcluding', paths[iter], 'from paths to take into consideration.\n'))
    
    # Find index of path to remove
    index_to_remove <- grep(pattern = paths[iter], paths)
    
    # Remove path from paths variable
    paths <- paths[-index_to_remove]
    
    # Continue
    next
  }
  
  # If there is only one row in the sequence table, add the sample name manually
  if(dim(seqtab)[1] == 1){
    rownames(seqtab) <- sample.names
  }
  
  # Add the seqtab to the main seqtab
  if(is.null(main.seqtab)){
    main.seqtab <- seqtab
  } else {
    main.seqtab <- mergeSequenceTables(main.seqtab, seqtab)
  }
  
  # Make a read output file
  getN <- function(x){sum(getUniques(x))}
  
  if(dim(seqtab)[1] == 1){
    track <- cbind(out, getN(dadaFwd), getN(dadaRev), getN(mergers), rowSums(seqtab))
  } else{
    track <- cbind(out, sapply(dadaFwd, getN), sapply(dadaRev, getN), sapply(mergers, getN), rowSums(seqtab))
  }
  
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab")
  rownames(track) <- sample.names
  
  
  # Add the read track to the main read track 
  if(is.null(main.track)){
    main.track <- track
  } else{
    main.track <- rbind(main.track, track)
  }
  
  # Remove the unfiltered .fastq.gz files
  unlink(c(FwdRead, RevRead))
}


#######################
## STORE MAIN SEQTAB ##
#######################

# Make directory path to store all results
path.result <- file.path(trialpath, 'Results')

if(!dir.exists(path.result)){
  cat(paste('Creating output directory:', path.result, '\n'))
  dir.create(path.result)
}

# Make directory to store all ASV sequence information
path.asv_sequence <- file.path(path.result, '01.ASV')

if(!dir.exists(path.asv_sequence)){
  cat(paste('Creating output directory:', path.asv_sequence, '\n'))
  dir.create(path.asv_sequence)
}

# Store main seqtab table
saveRDS(main.seqtab, file.path(path.asv_sequence, 'seqtab.rds'))


#####################
## PROCESSING ASVs ##
#####################

print_header(4)

  #####################
  ## REMOVE CHIMERAS ##
  #####################

  # Message
  step <- 1
  cat(paste0('[', 'step ', step, '] Removing chimeras\n'))
  
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(main.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  
  # Store seqtab.nochim
  saveRDS(seqtab.nochim, file.path(path.asv_sequence, 'seqtab_nochim.rds'))
  
  # Update the main.track table
  main.track <- cbind(main.track, rowSums(seqtab.nochim))
  colnames(main.track)[ncol(main.track)] <- "seqtab_nochim"
  
  
  #######################
  ## Remove singletons ##
  #######################
  
  if (singletons == F & dim(seqtab.nochim)[1] > 1){
    
    # Message
    step <- step + 1
    cat(paste0('\n[', 'step ', step, '] Removing singletons\n'))
    
    mode(seqtab.nochim) = 'numeric'
    seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1, drop = F]
    
    # Store seqtab.nochim without singletons
    saveRDS(seqtab.nochim, file.path(path.asv_sequence, 'seqtab_nochim_nosingle.rds'))
    
    # Update the main.track table
    main.track <- cbind(main.track, rowSums(seqtab.nochim))
    colnames(main.track)[ncol(main.track)] <- "seqtab_nochim_nosingle"
  }
  
  write.table(x = main.track, file = file.path(path.asv_sequence, 'summary.txt'))
  
  # If there are no sequences in the sequence table, stop
  if(dim(seqtab.nochim)[2] == 0){
    stop('\nASV sequences could not be generated.')
  }
  
  # Get ASV sequences
  asv_seqs <- colnames(seqtab.nochim)
  
  # Produce ASV headers
  asv_headers <- vector(dim(seqtab.nochim)[2], mode = 'character')
  
  for (i in 1:dim(seqtab.nochim)[2]){
    asv_headers[i] <- paste0('>ASV', i)
  }
  
  # Combine ASV sequences and headers
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  
  # Write ASV sequences and headers to a text (i.e. multifasta) and rds file
  write(asv_fasta, file.path(path.asv_sequence, paste0(GOI, '_ASV.fasta')))
  saveRDS(seqtab.nochim, file.path(path.asv_sequence, 'seqtab.rds'))
  
  if(bold){
    # Write second ASV multifasta file (-> For compatibility with Windows)
    write(asv_fasta, file.path(path.asv_sequence, paste0(GOI, '_ASV2.fasta')))
  }

  
#################
## ASV FIGURES ##
#################

print_header(5)
  
  #####################
  ## ASV RAREFACTION ##
  #####################
  
  # Message
  step <- 1
  cat(paste0('[', 'step ', step, '] Plotting ASV rarefaction\n'))
  
  # Make directory to store all ASV plots
  path.asv_plot <- file.path(path.result, '02.Plots')
  
  if(!dir.exists(path.asv_plot)){
    cat(paste('Creating output directory:', path.asv_plot, '\n\n'))
    dir.create(path.asv_plot)
  }
  
  if(nrow(seqtab.nochim) > 1){
    
    # Make the rarefaction plot
    suppressWarnings({
      # Remove rows where the sum of the observed count data across entire row is smaller than 2 
      seqtab.nochim <- seqtab.nochim[!rowSums(seqtab.nochim) <= 1,]
      
      S <- specnumber(seqtab.nochim) # observed number of ASV sequences, same as apply(seqtab.nochim, 1, function(x){sum(x != 0)})
      
      (raremax <- min(rowSums(seqtab.nochim))) # rarefaction threshold
      
      # The rarefy() function is used to estimate the expected number of species 
      # (or ASV sequences) in random subsamples of a specified size from a community 
      # or dataset. It is a common method for comparing species richness or diversity 
      # across different samples or communities. In the line below, rarefy(seqtab.nochim, raremax) 
      # calculates the rarefied species richness by simulating the expected number of 
      # species in random subsamples of size raremax from the seqtab.nochim dataset. 
      # The result is assigned to the variable Srare. The rarefaction process involves 
      # randomly selecting a subset of individuals (or ASV sequences) from the dataset 
      # without replacement, and then calculating the number of unique species (or ASV 
      # sequences) present in each subsample. This process is repeated multiple times 
      # to estimate the expected species richness for each subsample size.
      
      Srare <- rarefy(seqtab.nochim, raremax) # rarefy the data
      
      pdf(file = file.path(path.asv_plot, 'ASVRarefaction.pdf')) ## Write to pdf
      
      plot(S, Srare, xlab = "Observed No. of ASVs", ylab = "Rarefied No. of ASVs")
      abline(0, 1) # reference line
      
      rarecurve(seqtab.nochim, step = 20, sample = raremax, col = "blue", cex = 0.6, ylab = 'ASVs') # rarefaction curve
      
      dev.off()
    })
    
  } else{
      cat('Rarefaction not possible due to the sequence table having only 1 row.\n')
  }
  
  
  ##################################
  ## SEQUENCE LENGTH DISTRIBUTION ##
  ##################################
  
  # Message
  step <- step + 1
  cat(paste0('[', 'step ', step, '] Plotting sequence length distribution\n'))
  
  # Calculate the distribution
  distribution <- nchar(getSequences(seqtab.nochim))
  
  # Construct plots
  if(length(distribution) > 1){
    
    # Plot density line
    suppressWarnings({
      pdf(file = file.path(path.asv_plot, 'SequenceLengthDistribution.pdf')) ## Write to pdf
      
      # Plot histogram
      hist(distribution, main = 'Sequence length distribution',
           xlab = 'Sequence length',
           ylab = 'Frequency')
      
      plot(density(distribution), main = 'Sequence length distribution', 
           xlab = 'Sequence length', 
           ylab = 'Density')
      
      invisible(dev.off())
      
    })
  } else{
    cat('Unable to plot sequence length distribution due to inadequate number of ASV sequences.\n')
  }
  
  
  #################################
  ## ASV DISTRIBUTION PER SAMPLE ##
  #################################
  
  # Transpose the matrix to have samples as columns and sequences as rows
  abundance_matrix <- t(seqtab.nochim)
  
  # Convert the matrix to a data frame
  abundance_df <- as.data.frame(abundance_matrix)
  
  # Reshape the data into long format
  abundance_long <- pivot_longer(abundance_df, everything(), names_to = "Sample", values_to = "Abundance")
  
  # Divide the samples into groups of 5
  num_samples <- length(unique(abundance_long$Sample))
  group_size <- 5
  num_groups <- ceiling(num_samples / group_size)
  group_names <- rep(seq_len(num_groups), each = group_size)
  group_names <- rep(paste('Group', group_names), length.out = num_samples)
  groups <- setNames(group_names, colnames(abundance_df))
  
  abundance_long$Group <- groups[abundance_long$Sample]
  
  # Remove 0 values
  abundance_long <- abundance_long[abundance_long$Abundance != 0,]
  
  suppressWarnings({
    pdf(file = file.path(path.asv_plot, 'ASVAbundanceDistribution.pdf')) ## Write to pdf
      
      # Create the violin swarm plot
      show(ggplot(abundance_long, aes(x = Sample, y = Abundance)) +
        geom_violin(trim = TRUE) +
        geom_point(aes(color = Group), size = 1.5, position = position_jitter(width = 0.1, height = 0), show.legend = F) +
        labs(x = "Sample", y = "Abundance") +
        facet_wrap(~ Group, ncol = 2, scales = "free_x") +
        ggtitle(label = 'Distribution of ASV Abundance per sample') +
          theme(plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))
      
    invisible(dev.off())
  })
  
  cat('\n')


##############################
## TAXONOMIC CLASSIFICATION ##
##############################

# Check if taxonomic classification needs to be performed
if((reference == T | bold == T) & length(paths) > 0){
  
  # Taxonomy message
  print_header(6)
  
  # Get paths to ASV multifasta files
  path.ASV <- file.path(path.asv_sequence, paste0(GOI, '_ASV.fasta'))        # To avoid potential errors in Windows, 
  path.ASV2 <- file.path(path.asv_sequence, paste0(GOI,'_ASV2.fasta'))       # 2 identical multifasta files are used.
  
  # Get paths to transformed ASV multifasta files (BOLDigger renames the used multifasta files to GOI_ASVS2_done.fasta)
  bold.path.ASV <- file.path(path.asv_sequence, paste0(GOI, '_ASV2_done.fasta'))
  
  # Construct paths to directories where taxonomy files will be stored
  path.taxon <- file.path(path.result, '03.Taxonomy')
  
  # Create directory where taxonomy is stored
  if(!dir.exists(path.taxon)){
    cat(paste('Creating output directories:', path.taxon, '\n'))
    dir.create(path.taxon)
  }
  
  
  #############################################################
  ## TAXONOMIC CLASSIFICATION WITH LOCAL REFERENCE DATABASES ##
  #############################################################
  
  if(reference){
    
    cat('\n[Reference libraries]\n')
    
    # Read in the input multifasta file
    seqs <- readDNAStringSet(path.ASV)
    IDs <- names(seqs)
    
    # Loop through each reference database file
    for (ref_file in ref_files) {
      
      # Get the reference database base and name
      ref_base <- basename(ref_file)
      ref_name <- (basename(fs::path_ext_remove(ref_file)))
      
      # Extract the taxonomic levels per reference database from config_df
      ref_levels <- gsub(' ', '', str_to_title(strsplit(config_df[config_df$Database == ref_base, "TaxonomicLevel"], ";")[[1]]))
      
      cat(paste('Reference database:', ref_name, '\n'))
      
      if(!dir.exists(file.path(path.taxon, ref_name))){
        cat(paste('Creating output directories:', file.path(path.taxon, ref_name), '\n'))
        dir.create(file.path(path.taxon, ref_name))
      }
      
      # Create the taxonomy table using the DADA2 assignTaxonomy function
      tax_table <- assignTaxonomy(seqs, ref_file, outputBootstraps = TRUE, taxLevels = ref_levels, minBoot = minBoot, multithread = T)
      tax_table <- cbind(ID = IDs, tax_table)
      
      # Print the taxonomy table to an Excel file
      library(openxlsx)
      write.xlsx(x = as.data.frame(tax_table), file = file.path(path.taxon, ref_name, paste0(ref_name,'_tax_classification', '.xlsx')), asTable = T, sheetName = 'Taxonomy')
    }
  }

  
  ########################################
  ## TAXONOMIC CLASSIFICATION WITH BOLD ##
  ########################################
  
  if(bold){
    
    cat('\n[BOLD] ')
    
    # Specify the GOI BOLDigger command line argument
    if(GOI == 'COI'){
      bold_argument <- 'ie_coi'
    } else if(GOI == 'ITS'){
      bold_argument <- 'ie_its'}
    
    # Execute the BOLDigger command line tool: find top 20 hits
    system2(command = 'boldigger-cline', args = c(bold_argument, 
                                                  paste0("\"", username, "\""), 
                                                  paste0("\"", password, "\""), 
                                                  paste0("\"", path.ASV2, "\""), 
                                                  paste0("\"", path.taxon, "\""),
                                                  batch_size))
    
    # Get the file names in the Taxonomy directory (not recursive!)
    files <- list.files(path = path.taxon, full.names = T)
    
    non_dirs <- c()
    
    for(file in files){
      if(dir.exists(paths = file)){next} 
      else{non_dirs <- c(non_dirs, file)}
    }
    
    # Search for the BOLDigger output Excel file using regular expressions
    excel_file <- files[grepl(".xlsx$", files)]
    
    # Execute the BOLDigger command line tool: find first hit (top-scoring match)
    system2(command = 'boldigger-cline', args = c('first_hit', paste0("\"", excel_file, "\"")))
    
    # Read in the second sheet of the BOLDigger output excel file
    bold_taxonomy <- read.xlsx(xlsxFile = excel_file, sheet = 'First hit')
    
    # Change the > symbol from the ASV ID (Windows compatibility)
    bold_taxonomy$ID <- gsub(pattern = '^&gt;', replacement = '>', bold_taxonomy$ID)
    
    # Create directory to store only the first hits
    if(!dir.exists(file.path(path.taxon, 'BOLDSYSTEMS'))){
      cat(paste('Saving first hit outputs to:', file.path(path.taxon, 'BOLDSYSTEMS')), '\n')
      dir.create(file.path(path.taxon, 'BOLDSYSTEMS'))
    }
    
    # Write the contents of the second sheet to a separate excel file
    write.xlsx(x = as.data.frame(bold_taxonomy), file = file.path(path.taxon, 'BOLDSYSTEMS', 'BOLD_first_hit.xlsx'), asTable = T, sheetName = 'Sheet1')
    
    # Remove the original BOLDigger output files
    for(non_dir in non_dirs){unlink(x = non_dir)}
    
    # Remove the duplicate multifasta file with ASV sequences
    unlink(x = bold.path.ASV)
  }
  
  if(fuse){
    cat('\n[Merging]\n\n')
      
      # Source the taxonomic table merging script
      source(file.path(dirname(pipeline_path), 'TaxTableMerger.R'))
  }
}
