%Temperature Forecast (ARIMA)
%assumes 15 minute averages
%data should start at midnight; be sure baseline is full days
%load data
clear;
load('CyCMV_312.mat'); %This is the data.
RealTemp = CyCMVOHSU2220220127M31221S1.M31221;
BaseLine = RealTemp(1:384);
BaseLength = length(BaseLine);% length of vector of baseline values (number of points)
BaseModel = BaseLine(1:288);
TotalLength = length(RealTemp);% length of vector for all temperature values (number of points)TotalLength = length(RealTemp);% length of vector for all temperature values (number of points)
filename = 'M312-21_fever.xlsx';

%Estimate ARIMA Model
TempMdl = arima('Constant',0,'D',0,'Seasonality',96,'MALags',10,'SMALags',96);
        
EstTempMdl = estimate(TempMdl,BaseModel); % Generation of predictive model from baseline data with specified parameters.
[Predicted,yMSE] = forecast(EstTempMdl,TotalLength,'Y0',BaseModel); % Generates points of forecasted data.
avtemp = nanmean(BaseLine); % Average temperature, for reference.
sttemp = nanstd(BaseLine); % Standard deviation of baseline temp, for reference.
maxbase = nanmax(BaseLine); % maximum temperature during baseline
minbase = nanmin(BaseLine); % minimum temperature during baseline

[E,V,logL] = infer(EstTempMdl, BaseLine);
% Calculates goodness-of-fit statistics.   
 
%Calculate Residuals (difference between actual temp and predicted)
Residual = zeros(TotalLength,1);close all  
% This is an empty vector that will eventually hold the values of the residuals.

for i=1:TotalLength                       
    Residual(i) = RealTemp(i)-Predicted(i);
    % This for-loop populates the Residual vector with residual values.
end

% calculate max and min temperature and residuals
Rmax = max(Residual(BaseLength:TotalLength));
Tmax = nanmax(RealTemp(BaseLength:TotalLength));
Rmin = min(Residual(BaseLength:TotalLength));
Tmin = nanmin(RealTemp(BaseLength:TotalLength));
ResidualSq = zeros(BaseLength,1);close all

% calculate upper and lower limits for significant deviations in temperature
for i=2:BaseLength
    if isnan(BaseLine(i)) || isnan(BaseLine(i-1))%make sure NaN isn't included in RSS
        ResidualSq(i) = 0;    
    else
        ResidualSq(i) = (BaseLine(i)-BaseLine(i-1))^2;
        % calculates squared value for residuals in the baseline period
    end
end

ResidualSumSq = sum(ResidualSq(1:BaseLength));
% calculates residual sum of squares
ResidualSumSqVal = ResidualSumSq/(BaseLength-1);

ResidualUpper = zeros(TotalLength,1);
for i = 1:length(ResidualUpper)                  
    ResidualUpper(i,1) = 3*sqrt(ResidualSumSqVal); 
    % This for-loop generates an upper bound line for residuals.
end

ResidualLower = zeros(TotalLength,1);
for i = 1:length(ResidualLower)
    ResidualLower(i,1) = -1* ResidualUpper(i,1);
    % this for loop creates the lower bound line for residuals
end

% calculate time
TimePre = zeros(BaseLength,1);
StartTime = -1*(BaseLength/96);
for i = 1:BaseLength
    TimePre(i,1) = StartTime + (i/96);
end

PostLength = TotalLength - BaseLength;
TimePost = zeros(PostLength,1);
for i = 1:PostLength
    TimePost(i,1) = i/96;
end
TimeTotal = vertcat(TimePre,TimePost);

% Calculate 6 hour deviations in temperature
MedianLength = floor(TotalLength/24);
SixhrFeverHours = zeros(MedianLength,1);
SixhrFeverDuration = zeros(MedianLength,1);
SixhrTempMedian = zeros(MedianLength,1);
SixhrMax = zeros(MedianLength,1);
TimeMedian = zeros(MedianLength,1);
A = 1; % use this to count
for i = 1:MedianLength
    B = A+23;
    SixhrFeverSum = 0; % initialize to 0
    SixhrFeverPoints = 0; % initialize to 0
    for z = A:B
        if Residual(z,1)> ResidualUpper(z,1)
            SixhrFeverSum = SixhrFeverSum + Residual(z,1);
            SixhrFeverPoints = SixhrFeverPoints + 1;
        end
    SixhrFeverHours(i,1) = SixhrFeverSum/4;
    SixhrFeverDuration(i,1) = SixhrFeverPoints/4;
    SixhrTempMedian(i,1) = nanmedian(Residual(A:B));
    SixhrMax(i,1) = nanmax(Residual(A:B));
    end
    TimeMedian(i,1) = TimeTotal(A,1);
    A = A+24;
end

%daily deviations in temperature
DailyFeverLength = floor(TotalLength/96);
DailyFeverHours = zeros(DailyFeverLength,1);
DailyFeverDuration = zeros(DailyFeverLength,1);
DailyRmax = zeros(DailyFeverLength,1);
DailyMedian = zeros(DailyFeverLength,1);
TimeDays = zeros(DailyFeverLength,1);
A = 1;
for i = 1:DailyFeverLength
    DailyFeverSum = 0; % initialize to 0
    DailyFeverPoints = 0; % initialize to 0
    B = A+95;
    for z = A:B
        if Residual(z,1)> ResidualUpper(z,1)
            DailyFeverSum = DailyFeverSum + Residual(z,1);
            DailyFeverPoints = DailyFeverPoints + 1;
        end
        DailyFeverHours(i,1) = DailyFeverSum/4;
        DailyFeverDuration(i,1) = DailyFeverPoints/4;
        DailyRmax(i,1) = nanmax(Residual(A:B));
        DailyMedian(i,1) = nanmedian(Residual(A:B));
    end
    TimeDays(i,1) = TimeTotal(A,1);
    A = A+96;
end

%calculate fever duration and fever-hours
FeverSum = 0;
FeverPoints = 0;
for i = BaseLength:TotalLength
    if Residual(i,1) > ResidualUpper(i,1)
        FeverSum = FeverSum + Residual(i,1);
        FeverPoints = FeverPoints + 1;
    end
end

FeverHours = FeverSum/4;
FeverDuration = FeverPoints/4;

HypoThermSum = 0;
HypoThermPoints = 0;
for i = BaseLength:TotalLength
    if Residual(i,1) < ResidualLower(i,1)
        HypoThermSum = HypoThermSum + Residual(i,1);
        HypoThermPoints = HypoThermPoints + 1;
    end
end

HypoThermHours = abs(HypoThermSum)/4;
HypoThermDuration = (HypoThermPoints/4);

T = table(TimeTotal, Predicted, RealTemp, yMSE, Residual, ResidualUpper, ResidualLower);
writetable(T,filename)
RSS = {'Residual Sum of Squares'; ResidualSumSq};
xlswrite(filename, RSS, 1, 'H2')
DegF = {'Degrees of Freedom'; BaseLength};
xlswrite(filename, DegF, 1, 'H4')
BaseAve = {'Baseline Mean'; avtemp};
xlswrite(filename,BaseAve,1,'I4')
BaseStd = {'St Dev'; sttemp};
xlswrite(filename,BaseStd,1,'J4')
BaseMax = {'Maximum'; maxbase};
xlswrite(filename,BaseMax,1,'K4')
BaseMin = {'Minimum'; minbase};
xlswrite(filename,BaseMin,1,'L4')
TempMax = {'Max Temp'; Tmax};
xlswrite(filename,TempMax,1,'H7')
ResMax = {'Max Residual'; Rmax};
xlswrite(filename,ResMax,1,'I7')
FevDur = {'Duration'; FeverDuration};
xlswrite(filename,FevDur,1,'J7')
FevH = {'Fever-Hours'; FeverHours};
xlswrite(filename,FevH,1,'K7')
TempMin = {'Min Temp'; Tmin};
xlswrite(filename,TempMin,1,'H10')
ResMin = {'Min Residual'; Rmin};
xlswrite(filename,ResMin,1,'I10')
HypoTD = {'Duration'; HypoThermDuration};
xlswrite(filename,HypoTD,1,'J10')
HypoSev = {'Severity'; HypoThermHours};
xlswrite(filename,HypoSev,1,'K10')
ST = table(TimeMedian, SixhrTempMedian, SixhrMax, SixhrFeverHours, SixhrFeverDuration);
writetable(ST, filename, 'Sheet', 2)
RT = table(TimeDays, DailyMedian, DailyFeverHours, DailyFeverDuration, DailyRmax);
writetable(RT, filename, 'Sheet', 3)