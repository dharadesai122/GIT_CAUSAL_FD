rm(list=ls())

#install.packages("tools")
#install.packages("pracma")


# Install the following libraries if you are missing any of them.
library(XML)
library(pcalg)
library(pracma)
library(tools)
library(graph)
library(pls)
library(permute)
library(C50)

# Set the path where functions.R file is located
home <- 'G:\\NC\\summer\\Mandar\\fwddeliverables\\climateData\\sahel'
setwd(home)
source("functions.R")
source("ida.pcr.pvalue.R")

# Climate Indices.csv file should contain the monthly values of indices
# For e.g. in case of Sahle it should contain the first six months of all climate indices
# Rainfall.csv should contain twelve months of rainfall.
data <- read.csv('Climate Raw Data Sahel.csv', row.names=1)
raw_rainfall_data <- read.csv('Sahel_Rainfall_Data.csv', header = TRUE)

# Enter the months for the rainfall season. 
season <- 7:9
col_names <- paste('Rainfall', season, sep='_')
seasonal_rainfall_data <- raw_rainfall_data[, col_names]
colnames(seasonal_rainfall_data) <- col_names
raw_data <- cbind(data, seasonal_rainfall_data)

# Generate training and testing data sets.
# The following function call will generate 10 folders, each containing sub-folders with training and testing data.
# For eg.
# CV_Run_1 -> CV_1, CV_2,...., CV_8; where
# CV_1 -> TrainData_CV_1, TestData_CV_1
#kfold_cv(raw_data, dim(raw_data)[1]) #leave one out

# Build causal networks using the Tetrad IV software.
# You will have to build a causal network for each TrainData_CV_1 and output them in CV_1->Causal_Graphs folder.

# Perform CDFD.
# The frequent set of features discovered will be output in a csv file titled "Feature Frequency Count.csv"
num_years <- nrow(raw_data)

years <- matrix(1:57, nrow=num_years, ncol=1, dimnames=list(1951:2007))

build_c50 <- predict_c50 <- result_c50 <- list()

disc_org_rainfall_data <-matrix(0, nrow=num_years, ncol=1)
disc_org_rainfall_data[, 1] <- getDiscretizedRainfall(raw_data, col_names)

cv_run <- 1:10
all_cdfd_features <- list()
unique_cdfd_features <- c()
k<-1
avg_num_features <- 0
avg_acc_per_run <- 0
for(i in 1:length(cv_run)){
  cv_run_dir <- paste('CV_Run', cv_run[i], sep='_')
  setwd(paste(home,  cv_run_dir, sep='/'))
  cv_dirs <- grep('CV_[1-8]', list.dirs(path='.', recursive=FALSE), value=TRUE)
  num_folds <- length(cv_dirs)
  print(cv_run_dir)
  avg_acc_per_cv <- 0
  for(j in 1:num_folds){
    train_data_file <- list.files(path=cv_dirs[j], pattern='^TrainData*.*.txt')
    test_data_file <- list.files(path=cv_dirs[j], pattern='TestData*.*.txt')
    
    print(train_data_file)
    
    train_data <- read.table(paste(cv_dirs[j], train_data_file, sep='/'), header=TRUE, sep='\t')
    test_data <- read.table(paste(cv_dirs[j], test_data_file, sep='/'), header=TRUE, row.names=1, sep='\t')
    cg_file_path <- file.path(paste(cv_dirs[j], 'Causal_Graphs', sep='/'))
    
    # Read Causal Graph from training data
    pc_graph_file <- list.files(path=cg_file_path, pattern='^ClimData_*.*.xml')
    
    cg_filename <- paste(file_path_sans_ext(pc_graph_file), 'xml', sep='.')
    
    print(cat('Reading files from ...', cg_filename))
    
    test_years <- array(0, dim=nrow(test_data))
    test_years <- years[row.names(test_data), 1]
    true_class <- array(0, dim=length(test_years))
    
    cv_no <- na.omit(as.numeric(unlist(strsplit(pc_graph_file, "[^0-9]+"))))[1]
    pc_graph <- tetradToR(paste(cg_file_path, cg_filename, sep="/"))    # Get the pc_graph object from the xml file
    
    pc_wm <- wgtMatrix(pc_graph)
    pc_wmD <- t(pc_wm - t(pc_wm))
    
    dir_pc <- which(pc_wmD==1, arr.ind=TRUE)  # Get list of directed edges from the graph
    dir_edge_rank <- matrix(0, nrow=nrow(dir_pc), ncol=3) # Store the edges and their causal effect on rainfall
    dir_edge_rank[, 1:2] <- dir_pc
    causal_effects <- matrix(0, nrow=length(pc_graph@nodes), ncol=1)   # Store the causal effect of each index on rainfall in anomalous phases
    rainfall_idx <- which(pc_graph@nodes=='Rainfall')
    
    disc_rainfall_data <- discVariable(as.matrix(train_data[, 'Rainfall']))   # Discretize the normalized rainfall 
    colnames(disc_rainfall_data) <- 'Rainfall'
    true_class <- disc_org_rainfall_data[test_years, 1]
    
    constructed_data <- cdfd(train_data, test_data, pc_graph, dir_edge_rank, rainfall_idx, option=1)
    if(is.null(constructed_data)){
      print(cat('Graph ', pc_graph_file[i], ' does not find any directed edges as features'))
      next
    }
    train_edge_features <- constructed_data[[1]]
    test_edge_features <- constructed_data[[2]]
    
    train_edge_features$Rainfall <- as.factor(disc_rainfall_data)
    
    model <- buildModel(train_edge_features, test_edge_features, response=disc_rainfall_data)
    build_c50[[cv_no]]     <- model[[1]];      predict_c50[[cv_no]]     <- model[[2]]
    avg_num_features <- avg_num_features + nrow(as.matrix(C5imp(build_c50[[cv_no]])))
    avg_acc_per_cv <- avg_acc_per_cv + calcAccuracy(true_class, model[[2]])
    cdfd_features <- as.matrix(C5imp(build_c50[[cv_no]]))
    all_cdfd_features[[k]] <- row.names(cdfd_features)[cdfd_features[, 1]>0]
    unique_cdfd_features <- union(unique_cdfd_features, all_cdfd_features[[k]])
  }
  print((avg_acc_per_cv / num_folds))
  avg_acc_per_run <- avg_acc_per_run + (avg_acc_per_cv / num_folds) #add avg accuracy for this run
  avg_num_features <- avg_num_features / num_folds
}

overall_acc <- avg_acc_per_run / length(cv_run)
freq_feature_mat <- matrix(0, nrow=num_folds*length(cv_run), ncol=length(unique_cdfd_features), dimnames=list(c(), unique_cdfd_features))
for(i in 1:length(all_cdfd_features)){
  freq_feature_mat[, all_cdfd_features[[i]]] <- 1
}

feature_count <- as.matrix(apply(freq_feature_mat, 2, sum))
sorted_feature_count<-as.matrix(feature_count[order(feature_count[,1], decreasing=TRUE),])
colnames(sorted_feature_count) <- 'Frequency'
setwd(home)
write.csv(sorted_feature_count, 'Feature Frequency Count.csv', row.names=TRUE)
