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

#     After you have trimmed your data to midnight and added in any missing time chunks
#     you will need to check how many baseline data points you have and change the Base_Length variable

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
#####################
column_number <- (3)
#####################

# Enter the name you would like to give to your new excel workbook
# Example: "M100-10.xlsx" / "F20-23.xlsx" / etc.
################################
animal_number <- "MBiphasicTest2.xlsx"
################################

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
for (x in 1:TotalLength){
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
ResidualSumSq <- sum(ResidualSq[1])
ResidualSumSqVal <- (ResidualSumSq/(length(Baseline)-1))

# Upper and Lower Residuals
ResidualUpper <- data.frame()
for (x in 1:TotalLength){
  V1 <- (3*(sqrt(ResidualSumSqVal)))
  ResidualUpper <- rbind(ResidualUpper,V1)
}
colnames(ResidualUpper) <- ("Residual Upper")

ResidualLower <- data.frame()
for (x in 1:TotalLength){
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

TimePost <- data.frame()
PostLength <- TotalLength - BaseLength
for (x in 1:PostLength){
  ST2 <- (x/96)
  TimePost <- rbind(TimePost,ST2)
}
colnames(TimePost) <- ("Time")

# Concatenate the two lists together
TimeTotal <- bind_rows(TimePre,TimePost)
colnames(TimeTotal) <- ("TimeTotal")


# 6 Hour Temp Derivations -------------------------------------------------
MedianLength <- floor(TotalLength/24)
SixhrFeverHours <- data.frame("SixhrFeverHours"=numeric())
SixhrFeverDuration <- data.frame("SixhrFeverDuration"=numeric())
SixhrTempMedian <- data.frame("SixhrTempMedian"=numeric())
SixhrMax <- data.frame("SixhrMax"=numeric())
TimeMedian <- data.frame("TimeMedian"=numeric())

A <- 1
for (x in 1:MedianLength){
  B <- A+23
  SixhrFeverSum <- 0
  SixhrFeverPoints <- 0
  for (i in A:B){
    if((Residual[i,1]>ResidualUpper[i,1])&(!is.na(Residual[i,1]))) {
      SixhrFeverSum <- SixhrFeverSum+Residual[i,1]
      SixhrFeverPoints <- SixhrFeverPoints+1
    }
    
    SixhrFeverHours[x,1] <- SixhrFeverSum/4
    SixhrFeverDuration[x,1] <- SixhrFeverPoints/4
    
    SixhrTempMedian[x,1] <- median(Residual[[1]][A:B])
    SixhrMax[x,1] <- max(Residual[[1]][A:B])
  }
  TimeMedian[x,1] <- TimeTotal[A,1]
  A <- A+24;
}
colnames(SixhrFeverHours) <- ("SixhrFeverHours")
colnames(SixhrFeverDuration) <- ("SixhrFeverDuration")
colnames(SixhrTempMedian) <- ("SixhrTempMedian")
colnames(SixhrMax) <- ("SixhrMax")
colnames(TimeMedian) <- ("TimeMedian")


# Daily Temp Derivations --------------------------------------------------
DailyFeverLength <- floor(TotalLength/96)
DailyFeverHours <- data.frame("DailyFeverHours"=numeric())
DailyFeverDuration <- data.frame("DailyFeverDuration"=numeric())
DailyRmax <- data.frame("DailyRmax"=numeric())
DailyMedian <- data.frame("DailyMedian"=numeric())
TimeDays <- data.frame("TimeDays"=numeric())



A <- 1
for (x in 1:DailyFeverLength){
  B <- A+23
  DailyFeverSum <- 0
  DailyFeverPoints <- 0
  B <- A+95
  for (z in A:B){
    if((Residual[z,1]>ResidualUpper[z,1])){
      
      DailyFeverSum <- DailyFeverSum+Residual[z,1]
      
      DailyFeverPoints <- DailyFeverPoints+1
    }
    DailyFeverHours[x,1] <- DailyFeverSum/4
    DailyFeverDuration[x,1] <- DailyFeverPoints/4
    DailyRmax[x,1] <- max(Residual[[1]][A:B])
    DailyMedian[x,1] <- median(Residual[[1]][A:B])
  }
  TimeDays[x,1] <- TimeTotal[A,1]
  A <- A+96;
}
colnames(DailyFeverHours) <- ("DailyFeverHours")
colnames(DailyFeverDuration) <- ("DailyFeverDuration")
colnames(DailyRmax) <- ("DailyRmax")
colnames(DailyMedian) <- ("DailyMedian")
colnames(TimeDays) <- ("TimeDays")


# Fever Duration and Fever-Hours Calculations -----------------------------
FeverSum <- 0
FeverPoints <- 0
for (x in BaseLength:TotalLength){
  if(Residual[x,1]>ResidualUpper[x,1]){
    FeverSum <- FeverSum + Residual[x,1]
    FeverPoints <- FeverPoints+1
  }
}
FeverHours <- FeverSum/4
FeverDuration <- FeverPoints/4


# Hypothermia Calculations ------------------------------------------------
HypoThermSum <- 0
HypoThermPoints <- 0
for(i in BaseLength:TotalLength){
  if(Residual[i,1]<ResidualLower[i,1]){
    HypoThermSum <- HypoThermSum + Residual[x,1]
    HypoThermPoints <- HypoThermPoints+1
  }                                     
}
HypoThermHours <- abs(HypoThermSum)/4
HypothermDuration <- (HypoThermPoints/4)


# Febrile Periods ---------------------------------------------------------

# Determine Febrile Periods
# Incubation = 0-0.5d dpi
# First Febrile = 0.5-2.5 dpi
# Second Febrile = 2.5-8 dpi
# Recovery period = 8+ dpi
PeakSum <- 
  PeakPoints <- 0
IncubationStart <- BaseLength
IncubationEnd <- BaseLength + 48
for(x in IncubationStart:IncubationEnd){
  if(is.na(Residual[x,1])){
  }
  else if (Residual[x,1]>ResidualUpper[x,1]){
    PeakSum <- PeakSum + Residual[x,1]
    PeakPoints <- PeakPoints + 1
  }
}
IncubationDuration <- PeakPoints/4
IncubationSeverity <- PeakSum/4
IncubationMax <- max(na.exclude(Residual[[1]][IncubationStart:IncubationEnd]))


# First Period
PeakSum <- 0
PeakPoints <- 0
FirstStart <- BaseLength+49
FirstEnd <- BaseLength+240
for(x in FirstStart:FirstEnd){
  if(is.na(Residual[x,1])){
  }
  else if (Residual[x,1]>ResidualUpper[[1]][5]){
    PeakSum <- PeakSum + Residual[x,1]
    PeakPoints <- PeakPoints + 1
  }
}
FirstDuration <- PeakPoints/4
FirstSeverity <- PeakSum/4
FirstMax <- max(na.exclude(Residual[[1]][FirstStart:FirstEnd]))

# Second Period
PeakSum <- 0
PeakPoints <- 0
SecondStart <- BaseLength + 241
SecondEnd <- BaseLength + 768
for(x in SecondStart:SecondEnd){
  if(is.na(Residual[x,1])){
  }
  else if(Residual[x,1]>ResidualUpper[[1]][5]){
    PeakSum <- PeakSum + Residual[x,1]
    PeakPoints <- PeakPoints + 1
  }
}
PeakPoints
PeakSum
SecondDuration <- PeakPoints/4
SecondSeverity <- PeakSum/4
SecondMax <- max(na.exclude(Residual[[1]][SecondStart:SecondEnd]))

# Recovery Period
PeakSum <- 0
PeakPoints <- 0
RecoveryStart <- BaseLength + 768

#There has to be data for the recovery period to check if there isn't then we will return zero for the values we need
if (RecoveryStart<TotalLength){
  for(x in RecoveryStart:TotalLength){
    if(is.na(Residual[x,1])){
    }
    else if (Residual[x,1]>ResidualUpper[[1]][5]){
      PeakSum <- PeakSum + Residual[x,1]
      PeakPoints <- PeakPoints + 1
    }
  }
}else{
  RecoveryDuration <- 0
  RecoverySeverity <- 0
}

RecoveryDuration <- PeakPoints/4
RecoverySeverity <- PeakSum/4

if (is.na(RecoveryStart)| RecoveryStart>TotalLength){
  RecoveryMax <- 0
}else{
  RecoveryMax <- max(na.exclude(Residual[[1]][RecoveryStart:TotalLength]))
}


# Excel Workbook Creation and Data Filing ---------------------------------
#wb <- createWorkbook()


Label1 <- list("Degrees of Freedom","Baseline Mean","St Dev","Maximum","Minimum")
Label2 <- list("Max Temp","Max Residual","Duration","Fever-Hours","Hypothermia Duration","Hypothermia Hours")
Label3 <- list("Min Temp","Min Residual")
Label4 <- list("Period", "Incubation", "First Febrile Period", "Second Febrile Period", "Recovery Period")

#Create Worksheet
addWorksheet(wb,"Sheet1")
writeData(wb, "Sheet1", TimeTotal, startCol = 1, startRow = 1,colNames = TRUE)
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

writeData(wb, "Sheet1", Label2, startCol = 8, startRow = 7)
writeData(wb, "Sheet1", Tmax, startCol = 8, startRow = 8)
writeData(wb, "Sheet1", Rmax, startCol = 9, startRow = 8)
writeData(wb, "Sheet1", FeverDuration, startCol = 10, startRow = 8)
writeData(wb, "Sheet1", FeverHours, startCol = 11, startRow = 8)
writeData(wb, "Sheet1", HypothermDuration, startCol = 12, startRow = 8)
writeData(wb, "Sheet1", HypoThermHours, startCol = 13, startRow = 8)

writeData(wb, "Sheet1", Label3, startCol = 8, startRow = 10)
writeData(wb, "Sheet1", Tmin, startCol = 8, startRow = 11)
writeData(wb, "Sheet1", Rmin, startCol = 9, startRow = 11)

writeData(wb, "Sheet1", Label4, startCol = 8, startRow = 13)
writeData(wb, "Sheet1", "RMax", startCol = 8, startRow = 14)
writeData(wb, "Sheet1", IncubationMax, startCol = 9, startRow = 14)
writeData(wb, "Sheet1", FirstMax, startCol = 10, startRow = 14)
writeData(wb, "Sheet1", SecondMax, startCol = 11, startRow = 14)
writeData(wb, "Sheet1", RecoveryMax, startCol = 12, startRow = 14)

writeData(wb, "Sheet1", "Duration", startCol = 8, startRow = 15)
writeData(wb, "Sheet1", IncubationDuration, startCol = 9, startRow = 15)
writeData(wb, "Sheet1", FirstDuration, startCol = 10, startRow = 15)
writeData(wb, "Sheet1", SecondDuration, startCol = 11, startRow = 15)
writeData(wb, "Sheet1", RecoveryDuration, startCol = 12, startRow = 15)

writeData(wb, "Sheet1", "Fever-Hours", startCol = 8, startRow = 16)
writeData(wb, "Sheet1", IncubationSeverity, startCol = 9, startRow = 16)
writeData(wb, "Sheet1", FirstSeverity, startCol = 10, startRow = 16)
writeData(wb, "Sheet1", SecondSeverity, startCol = 11, startRow = 16)
writeData(wb, "Sheet1", RecoverySeverity, startCol = 12, startRow = 16)

addWorksheet(wb,"Sheet2")
writeData(wb, "Sheet2", TimeMedian, startCol = 1, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", SixhrTempMedian, startCol = 2, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", SixhrMax, startCol = 3, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", SixhrFeverHours, startCol = 4, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet2", SixhrFeverDuration, startCol = 5, startRow = 1,colNames = TRUE)

addWorksheet(wb,"Sheet3")
writeData(wb, "Sheet3", TimeDays, startCol = 1, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet3", DailyMedian, startCol = 2, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet3", DailyFeverHours, startCol = 3, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet3", DailyFeverDuration, startCol = 4, startRow = 1,colNames = TRUE)
writeData(wb, "Sheet3", DailyRmax, startCol = 5, startRow = 1,colNames = TRUE)

saveWorkbook(wb, file=variable, overwrite = TRUE)




