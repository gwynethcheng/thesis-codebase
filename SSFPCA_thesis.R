# testing with example data from author
load("/Users/gwyneth/Desktop/Masters_Thesis/Codebase/SSFPCA/Taro_Chl_trimmed.RData")
# overview of all objects
ls()
# "index.lk"  "lat.trim"  "lon.trim"  "mask.mat"  "Taro.data" "Taro.trim" "Time"

# check structure of each object
str(Taro.data)  # full lake Chl data
str(Taro.trim)  # trimmed lake data, used in the SS-FPCA
str(mask.mat)   # binary mask matrix — defines the lake boundary
str(index.lk)   # pixel indices identifying valid lake pixels
str(lat.trim)   # latitude coordinates of trimmed grid
str(lon.trim)   # longitude coordinates of trimmed grid
str(Time)       # time index of 100 monthly images (June 2002-April 2012)

# check dimensions of each object
dim(Taro.data)  # 2D array, 4635 x 100 (spatial extent x monthly)
dim(Taro.trim)  # 3D array, 103 x 45 x 101 (spatial grid x monthly)
dim(mask.mat)   # 2D array, 103 x 45 (spatial grid)

# first six values
head(Time)
head(lat.trim)
head(lon.trim)

# To load the data matrix in R:
  # data <- read.csv('/Users/gwyneth/Desktop/Masters_Thesis/LEVEL_2_SSFPCA_PREPPED/SEN3TEST/ssfpca_data_chl_oc4me_data_matrix.csv')
# data is a data.frame with n rows and N columns
# Convert to matrix for SS-FPCA:
  # data.mat <- as.matrix(data)   # n x N, NAs for missing

# data.mat corresponds directly to Taro.data in the author's R script.
# Note: data are already log-transformed and need centering
# (subtract per-time-point mean) before running SS-FPCA,
# matching the 'Taro.mean' centering step in the R script.

# some considerations about run time of the ss-fpca:
# might need to reduce resolution (multiples of 300m)
# and/or make it monthly composite instead of daily composite
# or see if Claude has a more efficient way, but don't edit the original ss-fpca code too much