# NOTE: ANYTHING BETWEEN THE DOUBLE LINED NOTE SYMBOLS WILL NEED CHANGED PRIOR TO RUN
#       #################
#
#       #################

# TO RUN CODE - QUICK KEY
# CTRL + SHIFT + ENTER

# Establishing Packages and Working Directory
library (dplyr)
library (xlsx2dfs)
library (readxl)
library (utils)
library (forecast)

#___________________________________________________________________________________________________
# FOLLOW INSTRUCTIONS BELOW
#     You must set your working directory prior to running your file through the following script
#     Your Working directory should be in the same file location as your raw Ponemah excel temp data
#     This only needs done once - just to establish the location

#####     Session -> Set Working Directory -> To Source File Location
#___________________________________________________________________________________________________

# NOTE: You cannot have the excel file open while you're extracting data from it
# Data Trimming Notes

#     Set Baseline Length
#     Your baseline length is the number of rows of temp data you have prior to the day of infection 
#     Follow this key as a reference:
#     3 days = (288) / 4 days = (384) / 5 days = (480)
####################
Base_Length <- (384)
####################

# Variables
# Enter the excel file name where your data is located
###########################
filename <- "TestData.xlsx"
###########################

# Enter the sheet in the excel file your data is located
# Example: 'Sheet1' / 'VEEV' / 'TempData'
################################################
filename <- read_excel(filename, sheet='Sheet3')
################################################


# Change the column_number variable to the column your temp data is located
# If it's in the first column it will be (1)
# The 5th column would be (5)
####################
column_number <- (3)
####################

# Enter the name you would like to give to your new excel workbook
# Example: "M100-10.xlsx" / "F20-23.xlsx" / etc.
###############################
animal_number <- "MBaselineTest-23.xlsx"
###############################

# varname.xlsx
RealTemp <- (filename[column_number])
wb <- createWorkbook()
variable <- animal_number

# Importing Temp Baseline and Data Set ------------------------------------
colnames(RealTemp) <- ("RealTemp")
Baseline <- RealTemp[[1]][1:Base_Length]
BaseLength <- length(Baseline)
BaseModel <- Baseline[1:Base_Length]
TotalLength <- lengths(RealTemp)


# Arima Model and Forecaster Conditions ------------------------------------
EstTempMdl <- arima(BaseModel, order=c(0,0,10), seasonal=list(order=c(0,1,1), period=96), method="ML")


# Generate Points of Forecasted Data
Predicted <- forecast(EstTempMdl, h=TotalLength, level=c(95))


# Baseline Averages -------------------------------------------------------

avtemp <- mean(Baseline, na.rm = TRUE)
sttemp <- sd(Baseline, na.rm = TRUE)
maxbase <- max(Baseline, na.rm = TRUE)
minbase <- min(Baseline, na.rm = TRUE)


# Residual Calculations ---------------------------------------------------

# Calculate the Residuals (Difference between actual temp and predicted)
Residual <- data.frame()  #Empty Vector
for (x in 1:BaseLength){
  T1 <- Predicted[[4]][x]
  T2 <- RealTemp[[1]][x]
  T3 <- T2-T1
  if(is.na(T3)){
    Residual <- rbind(Residual,0)
  }else{
    Residual <- rbind(Residual,T3)
  }
}
colnames(Residual) <- ("Residual")


# Calculate the Max and Min Temperature and Residuals
Rmax <- max(Residual, na.rm=TRUE)
Tmax <- max(RealTemp, na.rm=TRUE)
Rmin <- min(Residual, na.rm=TRUE)
Tmin <- min(RealTemp, na.rm=TRUE)

# Calculate Upper and Lower Limits for significant deviations in temperature
ResidualSq <- data.frame()
for (x in 2:length(Baseline)){
  S1 <- Baseline[x]
  S2 <- Baseline[x-1]
  S3 <- ((S1-S2)^2)
  if(is.na(S3)){
    ResidualSq <- rbind(ResidualSq,0)
  }else{
    ResidualSq <- rbind(ResidualSq,S3)
  }
}
colnames(ResidualSq) <- ("ResidualSq")

# Calculate the Residual Sum of Squares
ResidualSumSq <- sum(ResidualSq[1], na.rm = TRUE)
ResidualSumSqVal <- (ResidualSumSq/(BaseLength-1))

# Upper and Lower Residuals
ResidualUpper <- data.frame()
for (x in 1:BaseLength){
  V1 <- (3*(sqrt(ResidualSumSqVal)))
  ResidualUpper <- rbind(ResidualUpper,V1)
}
colnames(ResidualUpper) <- ("Residual Upper")

# Residual Lower
ResidualLower <- data.frame()
for (x in 1:BaseLength){
  V1 <- (-1*ResidualUpper[[1]][x])
  ResidualLower <- rbind(ResidualLower,V1)
}
colnames(ResidualLower) <- ("Residual Lower")


# Time Calculations -------------------------------------------------------

TimePre <- data.frame("Time"=numeric())
StartTime <- (-1*(BaseLength)/96)
for (x in 1:BaseLength){
  ST1 <- StartTime + (x/96)
  TimePre <- rbind(TimePre,ST1)
}
colnames(TimePre) <- ("Time")


#Determining Hour Fever / Hypothermia Determinations
BaseHrLength <- 24
BaseHrMedianTemp <- data.frame("BaseHrMedianTemp"=numeric())
FeverThreshold <- data.frame("FeverThreshold"=numeric())
HypoThreshold <- data.frame("HypoThreshold"=numeric())
BaseHrTime <- data.frame("BaseHrTime"=numeric())



# Hour Increments ---------------------------------------------------------
A <- 1
for (x in 1:BaseHrLength){
  B <- A+3
  BaseHrMedianTemp[x,1] <- median((Predicted[[4]])[A:B])
  FeverThreshold[x,1] <- BaseHrMedianTemp[x,1] + ResidualUpper[x,1]
  HypoThreshold[x,1] <- BaseHrMedianTemp[x,1] - ResidualUpper [x,1]
  BaseHrTime[x,1] <- x
  A <- A+4
}


# Excel Workbook Creation and Data Filing ---------------------------------
# Export Data - Create a Workbook
#wb <- createWorkbook()

Label1 <- list("Degrees of Freedom","Average Baseline Temp","Standard Deviation Baseline","Max Temp Baseline","Min Temp Baseline")
Label2 <- list("Max Temp","Max Residual","Duration","Fever-Hours","Hypothermia Duration","Hypothermia Hours")



#Create Worksheet
addWorksheet(wb,"Sheet1")
writeData(wb, "Sheet1", TimePre, startCol = 1, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet1", Predicted, startCol = 2, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet1", RealTemp, startCol = 3, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet1", Residual, startCol = 4, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet1", ResidualUpper, startCol = 5, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet1", ResidualLower, startCol = 6, startRow = 1,colNames = TRUE)

writeData(wb, "Sheet1", "Residual Sum of Squares", startCol = 8, startRow = 2)
writeData(wb, "Sheet1", ResidualSumSq, startCol = 8, startRow = 3)

writeData(wb, "Sheet1", Label1, startCol = 8, startRow = 4)
writeData(wb, "Sheet1", BaseLength, startCol = 8, startRow = 5)
writeData(wb, "Sheet1", avtemp, startCol = 9, startRow = 5)
writeData(wb, "Sheet1", sttemp, startCol = 10, startRow = 5)
writeData(wb, "Sheet1", maxbase, startCol = 11, startRow = 5)
writeData(wb, "Sheet1", minbase, startCol = 12, startRow = 5)

#Baseline Data
addWorksheet(wb,"Sheet2")
writeData(wb, "Sheet2", BaseHrTime, startCol = 1, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", HypoThreshold, startCol = 2, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", BaseHrMedianTemp, startCol = 3, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", FeverThreshold, startCol = 4, startRow = 1,colNames = TRUE)


saveWorkbook(wb, file=variable, overwrite = TRUE)



