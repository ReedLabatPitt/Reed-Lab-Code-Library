clear all

cd '~\Dropbox\amy\'

clear args AEG
args.rawDir = '~\Dropbox\amy\EEG\';
args.exp    =  '16026773';
args.animal =   '123';
args.chan   =   'ch1';

args.extn   = '16';
args.rawDir = '~\Dropbox\amy\Export\';
args.exp    =  '16026773';

targetSR = 4; % 4 samples per hour

clear d dirList
%d = dir([args.rawDir 'M' args.animal  '-16_EEG\*.edf']);
d = dir( [args.rawDir 'M' args.animal  '-' args.extn '_EEG\*.edf'] );
[dirList{1:length(d),1}] = deal(d.name);

args.setDir = ['~\Dropbox\amy\plots\' args.exp '_' args.animal '\'];
clear d setList
d = dir([args.setDir  '*fft15*.set']);
[setList{1:length(d),1}] = deal(d.name);

clear d setListMA
d = dir([args.setDir  '*ma15*.set']);
[setListMA{1:length(d),1}] = deal(d.name);


startTimeList = zeros(20,6) ;
durList       = zeros(20,1);

if strcmp(args.animal,'112')
    infection = [5 2 9 15 0]; %% 112
end

if strcmp(args.animal,'113')
    infection = [5 2 9 45 0]; %% 113
end

if strcmp(args.animal,'120')
    infection = [5 12 9 00 0]; %% 120
end

if strcmp(args.animal,'121')
    infection = [5 4 9 10 0]; %% 121
end

if strcmp(args.animal,'122')
    infection = [5 4 9 40 0]; %% 122
end

if strcmp(args.animal,'123')
    infection = [1 13 10 50 0]; %% 123
end


if strcmp(args.animal,'053')
    infection = [6 26 09 10 0]; %% 053
    args.extn   = '17';
end

if strcmp(args.animal,'054')
    infection = [6 26 09 40 0]; %% 054
    args.extn   = '17';
end

if strcmp(args.animal,'057')
    infection = [7 28 08 40 0]; %% 057
    args.extn   = '17';
end

if strcmp(args.animal,'058')
    infection = [8 04 09 10 0]; %% 057
    args.extn   = '17';
end


%% read all data sets that belong to this animal
for (i = 1:length(setList))
    
    % read the set
    tdat     = pop_loadset( setList{i}, args.setDir );
    mdat     = pop_loadset( setListMA{i}, args.setDir );
    %plot(tdat.data)
    
    if size(mdat.data,2)~= size(tdat.data,2)
        keyboard
    end
    
    %
    args.file  = dirList{i};
    tmp1       = strsplit(dirList{i},'_');
    tmp2       = strsplit(tmp1{2},'.');
    args.date  = tmp2{1};
    
    header = edfread( [args.rawDir 'M' args.animal  '-' args.extn '_EEG\' dirList{i}], 'AssignToVariables',true );
    
    tmpTime              = strsplit(header.starttime, '.');
    tmpDate              = strsplit(header.startdate, '.');
    
    if i == 1
        absoluteZero         = [ str2num( tmpDate{2} ) str2num( tmpDate{1} ) 0 0 0 ];
        % prepare the full data set
        AEG                  = tdat;
        AEG.data             = NaN( size(tdat.data,1), length(setList)*24*targetSR  );
        AEG.times            = (1:size(AEG.data,2))/(targetSR);
        keepTrack            =  NaN( length(setList), length(setList)*24*targetSR  );
        
        MA15                  = mdat;
        MA15.data             = NaN( size(mdat.data,1), length(setList)*24*targetSR  );
        MA15.times            = (1:size(MA15.data,2))/(targetSR);
    end
    
    thisTime             = [ str2num( tmpDate{2} ) str2num( tmpDate{1} ) str2num( tmpTime{1} ) 60*mod(tdat.times(1),1) 0];
    startTimeList(i,1:5) = thisTime;
    startTimeList(i,6)   = size(tdat.data,2)/(targetSR);
    
    daysInMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    dayMultiplier = daysInMonth(absoluteZero(1));
    
    startIndex = ...
        targetSR * dayMultiplier*24*(thisTime(1)-absoluteZero(1)) + ...
        targetSR *               24*(thisTime(2)-absoluteZero(2)) + ...
        targetSR *                  (thisTime(3)-absoluteZero(3)) + ...
        targetSR *             1/60*(thisTime(4)-absoluteZero(4));
    
    if(size(tdat.data,2)==1)
        %keyboard
    end
    
    if( size(tdat.data,2)>1 )
        toInd             = startIndex:(startIndex+size(tdat.data,2)-1);
        AEG.data(:,toInd) = NaN;
        if min( sum(tdat.data,1) )<1
            plot( sum(tdat.data,1) )
            %keyboard
        end
        
        valInd = find( sum(tdat.data,1)>1 );
        AEG.data( :,toInd(valInd)) = tdat.data(:,valInd);
        MA15.data(:,toInd(valInd)) = mdat.data(:,valInd);
        
        
        keepTrack(i,toInd(valInd)) = 1;
    end
    
    if size(AEG.data,2)~= size(MA15.data,2)
        keyboard
    end
    
end

%%

AEG.times            = (1:size(AEG.data,2))/(targetSR);

InfectionIndex = ...
    targetSR * dayMultiplier*24*(infection(1)-absoluteZero(1)) + ...
    targetSR *               24*(infection(2)-absoluteZero(2)) + ...
    targetSR *                  (infection(3)-absoluteZero(3)) + ...
    targetSR *             1/60*(infection(4)-absoluteZero(4));

InfectionIndexDay = ...
    targetSR * dayMultiplier*24*(infection(1)-absoluteZero(1)) + ...
    targetSR *               24*(infection(2)-absoluteZero(2)) ;

InfectionIndex    = floor(InfectionIndex);


AEG.times   = (AEG.times - AEG.times(InfectionIndexDay))/24;
AEG.setname = 'AEG';
AEG.pnts    = size(AEG.data,2);
AEG.filename = ['ALL_' args.animal '_AEG'];
pop_saveset(AEG, 'filename', AEG.filename, 'filepath',args.setDir );


MA15.times  = AEG.times;
MA15.setname = 'MA15';
MA15.pnts    = size(MA15.data,2);
MA15.filename = ['ALL_' args.animal '_MA15'];
pop_saveset(MA15, 'filename', MA15.filename, 'filepath',args.setDir );

plot(mean(AEG.data,1))
plot(mean(MA15.data,1))
%% create data set with four frequency bands
NFrq        = 500;
binSec      = 4;
freqs       = (0:1:(NFrq-1))/binSec ;

deltaPower = mean(AEG.data(freqs>0.4 & freqs<3 ,:),1);
thetaPower = mean(AEG.data(freqs>4   & freqs<8 ,:),1);
alphaPower = mean(AEG.data(freqs>8   & freqs<12,:),1);
betaPower  = mean(AEG.data(freqs>12  & freqs<30,:),1);

BND         = AEG;
BND.setname = 'BND';
BND.pnts    = size(BND.data,2);
BND.filename = ['ALL_' args.animal '_BND'];
BND.data   = [deltaPower; thetaPower; alphaPower; betaPower; MA15.data(1,:); MA15.data(2,:)];
BND.nbchan = size(BND.data,1);

pop_saveset(BND, 'filename', BND.filename, 'filepath',args.setDir );

%% create normalized data set
NRM = BND;
NRM.setname = 'NRM';
NRM.filename = ['ALL_',args.animal '_NRM'];

for j = 1:size(NRM.data,1)
    useInd                           = find( ~isnan(BND.data(j,:)) );
    NRM.data(j,:)                    = (BND.data(j,:)-mean(BND.data(j, useInd )) )/std(BND.data(j, useInd ));
end

% ============ low-pass filter the data
hicutoff   =  .5;
filtorder  =  8;  %

[NRMflt, com, b] = pop_firws(NRM, 'fcutoff', hicutoff, 'forder',filtorder,'ftype','lowpass');


plot(NRM.times,NRM.data)
plot(NRM.times,NRMflt.data)
pop_saveset(NRM, 'filename', NRM.filename, 'filepath',args.setDir );

% ============================================================
%% SOME EXAMPLES OF HOW TO DISPLAY THE DATA
% ============================================================

%% -------------------------- time-series
cd( args.setDir )
    
figure;
hold on
for d = min(NRM.times):max(NRM.times)
    fill( (d) + [18 18 30 30]/24, [-2 5 5 -2 ], [0.7 0.7 0.7] )
end

showChan = [1 2 3 4 5];
for j = 1:length(showChan)
    plot(NRM.times,(j-1)+0* (NRM.data(1,:)), 'k')
    plot(NRM.times,(j-1)+.2*NRMflt.data(showChan(j),:),'Linewidth',1)
end
plot(NRM.times,(-1)+0* (NRM.data(1,:)), 'k')
plot(NRM.times,(-1)+.1* (NRMflt.data(1,:) - NRMflt.data(4,:)),'k','Linewidth',1 )
plot( [ NRM.times(InfectionIndex) NRM.times(InfectionIndex)], [-2 5], 'k', 'Linewidth', 2  )

axis([-4 6 -1.5 5])
saveas(gcf,'GA_timeSeries', 'pdf');

axis([0 10 -1.5 5])
saveas(gcf,'GA_timeSeries2', 'pdf');

axis([-4 15 -1.5 5])
saveas(gcf,'GA_timeSeries3', 'pdf');


close all

%% ==============================
%    daytime vs night time during baseline
figure
AEG.hours = mod(AEG.times,1);
dayind = find( ~isnan(AEG.data(1,:)) & AEG.times<0 &  AEG.hours>0.25 & AEG.hours<0.75);
nitind = find( ~isnan(AEG.data(1,:)) & AEG.times<0 & (AEG.hours<0.25 | AEG.hours>0.75));

daySpc = mean(AEG.data(:,dayind),2);
nitSpc = mean(AEG.data(:,nitind),2);

plot(log10(freqs(1,1:160)'),log10(daySpc), 'g','Linewidth',2)
hold on
plot(log10(freqs(1,1:160)'),log10(nitSpc), 'k','Linewidth',2)

saveas(gcf,'GA_DayNight_baseline', 'pdf');
close all

%% ========================
% pre vs post
preind = find( ~isnan(AEG.data(1,:)) & AEG.times<0 );
pstind = find( ~isnan(AEG.data(1,:)) & AEG.times>2 );

preSpc = mean(AEG.data(:,preind),2);
pstSpc = mean(AEG.data(:,pstind),2);

figure
plot(log10(freqs(1,1:160)'),log10(preSpc), 'k','Linewidth',2)
hold on
plot(log10(freqs(1,1:160)'),log10(pstSpc), 'r','Linewidth',2)

saveas(gcf,'GA_PrePost', 'pdf');
close all


%% ========================
% day-time pre vs post
preind = find( ~isnan(AEG.data(1,:)) & AEG.times<0 & AEG.hours>0.25 & AEG.hours<0.75);
pstind = find( ~isnan(AEG.data(1,:)) & AEG.times>2 & AEG.hours>0.25 & AEG.hours<0.75);

preSpc = mean(AEG.data(:,preind),2);
pstSpc = mean(AEG.data(:,pstind),2);

figure
plot(log10(freqs(1,1:160)'),log10(preSpc), 'k','Linewidth',2)
hold on
plot(log10(freqs(1,1:160)'),log10(pstSpc), 'r','Linewidth',2)

saveas(gcf,'GA_PrePost_day', 'pdf');
close all


%% ========================
% night-time pre vs post
preind = find( ~isnan(AEG.data(1,:)) & AEG.times<0 & (AEG.hours<0.25 | AEG.hours>0.75));
pstind = find( ~isnan(AEG.data(1,:)) & AEG.times>2 & (AEG.hours<0.25 | AEG.hours>0.75));

preSpc = mean(AEG.data(:,preind),2);
pstSpc = mean(AEG.data(:,pstind),2);

figure
plot(log10(freqs(1,1:160)'),log10(preSpc), 'k','Linewidth',2)
hold on
plot(log10(freqs(1,1:160)'),log10(pstSpc), 'r','Linewidth',2)

saveas(gcf,'GA_PrePost_night', 'pdf');
close all

%% =======================
% spectra averaged by day

figure

for dx = -5:0 % pre-infection days
    dayind = find( ~isnan(AEG.data(1,:)) & floor(AEG.times)==dx);
    if(length(dayind)>24)
        daySpc = mean(AEG.data(:,dayind),2);

        plot(log10(freqs(1,1:160)'),log10(daySpc ), 'k','Linewidth',3, 'Color', [.8+dx/8 .8+dx/8 .8+dx/8])
        hold on
    end
end

for dx = 1:7 % post-infection days
    dayind = find( ~isnan(AEG.data(1,:)) & floor(AEG.times)==dx);
    if(length(dayind)>24)
        daySpc = mean(AEG.data(:,dayind),2);

        plot(log10(freqs(1,1:160)'),log10(daySpc ), 'k','Linewidth',2, 'Color', [dx/7 .05 1-dx/7])
        hold on
    end
end

saveas(gcf,'GA_SpactraByDay', 'pdf');
close all