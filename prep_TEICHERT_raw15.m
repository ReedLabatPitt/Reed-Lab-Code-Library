function [ output_args ] = prep_TEICHERT_raw15( args )
%%
% ultimate goal: create a set sampled at 1/15 minutes that includes:
%     delta, theta, alpha, beta and gamma power
%     activity
%     EEG envelope


FFTWinSec  = 4;              % 4 second long FFT windows

binDurMin  = 15;             % 15 minutes bin width
targetSR   = 1/(binDurMin*60); % average in bins of 15 minutes

% EEGLab structures created during the process
% EEG : channel 1: EEG trace at 500 Hz 50 Hz low-pass filtered
%       channel 2: chewing artefact detection channel. high-pass, rect, low-pass
%       channel 3: up-sampled activity monitor
%
% FFT15 : Spectra averaged in bins of 15 minutes
% MA15  : High-frequency envelope and activity in bins of 15 minute

%% load the data
[header, recorddata] = edfread( [args.rawDir 'M' args.animal '-' args.extn '_EEG\' args.file] );

%% identify all EEG channels
eegChanInd = strcmp(header.label, 'EEG') | strcmp(header.label, 'ECG');
NEEG       = sum(eegChanInd);
EEGchan    = find(eegChanInd);

minEEG     = 1000*header.physicalMin(EEGchan)*.99;          % droped values default to -0.0142

actChanInd = strcmp(header.label, 'Activity');
NACT       = sum(actChanInd);
ACTchan    = find(actChanInd);

tmpChanInd = strcmp(header.label, 'Temp');
NTMP       = sum(tmpChanInd);
TMPchan    = find(tmpChanInd);

tmp             = strsplit(header.starttime, '.');
args.onsetTime  = 3600 * str2num(tmp{1} ) + 60 * str2num(tmp{2} ) + str2num(tmp{3});


for chnd = 1: NEEG
    
    args.picDir = ['~\Dropbox\amy\plots\' args.exp '_' args.animal '\' ]; %'_ch' num2str(chnd) '\'];
    if ~exist(args.picDir,'dir')
        mkdir(args.picDir);
    end
    cd( args.picDir )
    
    args.fileName1 =  [args.date, '_sleepindex_' num2str(chnd)];
    args.fileName2 =  [ args.date, '_spectrumByActivity_' num2str(chnd)];
    
    %% put EEG data into a EEGlab structure, and set dropped data to 0; sampled at 500 Hz
    EEG                   = struct();
    taxis                 = (0:(size(recorddata,2)-1)) / header.frequency(EEGchan(chnd));
    
    tst                   = max(find(abs(recorddata(EEGchan(chnd),:)) > 1e-6));
    rawEEG                = recorddata(EEGchan(chnd),1:tst);
    taxis                 = taxis(1:max(tst));
    
    EEG.txS               = taxis + args.onsetTime;
    EEG.txH               = EEG.txS/3600;
    EEG.txD               = EEG.txH/24;
    
    EEG.data              = 1000*recorddata(EEGchan(chnd),1:tst);% + 10*sin(2*pi*taxis*10);
    
    dropEEGBin            = EEG.data < minEEG;
    dropEEGInd            = find(dropEEGBin);
    EEG.data(dropEEGInd)  = 0;
    
    EEG.srate             = 500;
    EEG.trials            = 1;
    EEG.event             = {};
    EEG.pnts              = size(EEG.data,2);
    EEG.nchan             = NEEG;
    
    % ==================================
    %    highpass filter to isolate epochs with strong muscle artifacts
    locutoff   =  50;
    filtorder  = 128;  %
    
    [tmp, com, b] = pop_firws(EEG, 'fcutoff', locutoff, 'forder',filtorder,'ftype','highpass');
    tmp.data = abs(tmp.data);
    [MUA, com, b] = pop_firws(tmp, 'fcutoff', 20, 'forder',filtorder,'ftype','lowpass');
    
    
    % ============ low-pass filter the data
    hicutoff   =  50;
    locutoff   = 1/50;
    filtorder  = 128;  %
    
    %[EEG, com, b] = pop_firws(EEG, 'fcutoff', [locutoff hicutoff], 'forder',filtorder,'ftype','bandpass');
    [EEG, com, b] = pop_firws(EEG, 'fcutoff', hicutoff, 'forder',filtorder,'ftype','lowpass');
    
    %
    EEG.data(2,:) = MUA.data;
    EEG.nchan     = NEEG + 1;
    
    %% put activty data into a EEGlab structure sampled at 10 HZ
    ACT                   = struct();
    taxis                 = (0:(size(recorddata,2)-1)) / header.frequency(ACTchan(chnd));
    
    tst                   = find(recorddata(ACTchan(chnd),:) > 0);
    %rawACT                = recorddata(ACTchan(chnd),1:max(tst));
    ACT.data              = recorddata(ACTchan(chnd),1:max(tst));% + 10*sin(2*pi*taxis*10);
    taxis                 = taxis(1:max(tst));
    
    ACT.txS               = taxis + args.onsetTime;
    ACT.times             = ACT.txS;
    ACT.txH               = ACT.txS/3600;
    ACT.txD               = ACT.txH/24;
    
    mnACT                 = mean(ACT.data);
    ACT.data              = ACT.data;% - mnACT;
    ACT.srate             = header.frequency(ACTchan(chnd));
    ACT.trials            = 1;
    ACT.event             = {};
    ACT.pnts              = size(ACT.data,2);
    ACT.nchan             = NACT;
    
    upACT = resample(ACT.data,header.frequency(EEGchan(chnd)),header.frequency(ACTchan(chnd)));
    
    tmx = max( size(EEG.data,2),size(upACT,2) );
    EEG.data(3,1:tmx) = upACT(1:tmx);
    EEG.nchan = 3;
    
    %% =============================================
    % calculate FFT in bins of 4 seconds length in non-overlapping steps of 4 second
    
    tmp         = strsplit(header.starttime, '.');
    startSec    = str2num(tmp{3});
    
    % find the first full 4 second period
    skipSec     = mod( FFTWinSec - mod(startSec,FFTWinSec), FFTWinSec);
    startInd    = skipSec * EEG.srate;
    binSec      = FFTWinSec;                           % length of each fft bin in seconds
    Nbin        = floor( (size(EEG.txS,2) - startInd)/(binSec*EEG.srate)); % number of bins
    
    %Nseconds    = floor(EEG.pnts/EEG.srate)-binSec-1;
    NFrq        = EEG.srate * binSec;                     % number of FFT bins
    freqs       = (0:1:(NFrq-1))/binSec ;
    valFrqInd   = find(freqs<40);
    NvalFrq     = length(valFrqInd);
    valfreqs    = freqs(valFrqInd);
    
    fftdat      = zeros(NvalFrq,Nbin);
    MAmax       = zeros(3,Nbin);
    fftaxS      = (args.onsetTime + startSec + FFTWinSec*(0:(Nbin-1)));
    fftaxH      = fftaxS/3600;
    fftaxD      = fftaxH/24;
    fftme       = mod(fftaxD,24);
    taper       = ([1:200 200*ones(1,NFrq-400) 200:-1:1]/200)-1;
    taper       = (1+cos(pi*taper))/2;
    
    for (i = 1:Nbin)
        thisStartInd  = startInd + binSec*(i-1)*EEG.srate;
        thisInd       = thisStartInd + (1:NFrq);
        
        thisfft       = fft(taper .* EEG.data(1, thisInd) );
        fftdat(:,i)   = thisfft(valFrqInd);
        
        MAmax(1,i)    = max(EEG.data(2,thisInd));
        MAmax(2,i)    = max(upACT(1,thisInd));
        MAmax(3,i)    = max(dropEEGBin(1,thisInd));
    end
    
    powfft    = abs(fftdat);%.^2;
    mnpow     = mean(powfft,1);
    
    %% ================================
    % split day-time spectra in four groups based on activity monitor
    % is there a relation between level of activity and specific frequency
    % bands?
    
    
    nightTimeIndx = find( (fftme<0.25 | fftme>0.75));
    dayTimeIndx   = find(  fftme>0.25 & fftme<0.75);
    brks          = quantile(MAmax(2,dayTimeIndx), [0 1/4 2/4 3/4 1] );
    
    %if length(dayTimeIndx)>40 && length(nightTimeIndx)>40
    %    figure
    %    plot(sort( MAmax(2,dayTimeIndx) ), (1:length(dayTimeIndx))/length(dayTimeIndx), 'r', 'Linewidth',3)
    %    hold on
    %    plot(sort( MAmax(2,nightTimeIndx) ), (1:length(nightTimeIndx))/length(nightTimeIndx),'k', 'Linewidth',3)
    %end
    
    if length(dayTimeIndx)>40
        figure
        plot(log10(valfreqs), log10(mean(abs(fftdat(:,dayTimeIndx)),2 )), 'r', 'Linewidth',3 );
        hold on
        
        for (qx = 1:4)
            qxind =  fftme>0.25 & fftme<0.75 & MAmax(2,:)> brks(qx) & MAmax(2,:)< brks(qx+1);
            plot(log10(valfreqs), log10(mean(abs(fftdat(:,qxind)),2 )),'Linewidth',1, 'Color',[qx/4,.5,1-qx/4]  );
            hold on
        end
    end
    
    if length(nightTimeIndx)>40
        plot(log10(valfreqs), log10(mean(abs(fftdat(:,nightTimeIndx)),2 )), 'k', 'Linewidth',3 );
        hold on
        
         for (qx = 1:4) 
            qxind =  (fftme<0.25 | fftme>0.75) & MAmax(2,:)> brks(qx) & MAmax(2,:)< brks(qx+1);
            plot(log10(valfreqs), log10(mean(abs(fftdat(:,qxind)),2 )),'Linewidth',1, 'Color',[qx/4,.5,1-qx/4], 'Linestyle','--'  );
            hold on
         end
    end
    
    saveas(gcf,args.fileName2, 'pdf');
    close all
    %% ================================================ downsample to 1 / 15 minutes
    % find the first 15 minute bin
    tmp         = strsplit(header.starttime, '.');
    startHr     = str2num(tmp{1});
    startMin    = str2num(tmp{2});
    startSec    = str2num(tmp{3});
    
    startMinBin = startMin - mod(startMin,binDurMin);
    startBinS   = 3600 * startHr + 60 * startMinBin + 0;
    startBinM   = startBinS/60;
    startBinH   = startBinM/60;
    
    % number of bins to be created
    N15bin     = ceil( FFTWinSec * (max(EEG.txS)-startBinS)/3600);
    
    fftdat15      = zeros(NvalFrq,N15bin);
    ma15dat       = zeros(5,N15bin);
    Nval15        = zeros(1,N15bin);
    valFFTaxSM    = MAmax(1,:) < .1;
    valFFTaxSA    = MAmax(2,:) < 1000000; %126000;
    valFFTaxSD    = MAmax(3,:) < .5;
    
    for bx = 1:N15bin
        % range of times in seconds that contribute to this 15 minute bin
        trngS = startBinS + [bx-1, bx] * binDurMin*60;
        
        % indeces of all spectra in the 15 minute bin
        thisFFTaxSbin    = find(fftaxS>trngS(1) & fftaxS<=trngS(2));
        ma15dat(1,bx)    = mean(MAmax(1,thisFFTaxSbin),2);
        ma15dat(2,bx)    = mean(MAmax(2,thisFFTaxSbin),2);
        
        %thisFFTaxSbin = find(TMP.txS>trngS(1) & TMP.txS<=trngS(2));
        %ma15dat(3,bx)    = mean(TMP.data(1,thisFFTaxSbin),2);
        
        % indeces of all valid spectra in the 15 minute bin
        thisFFTaxSbin = find(fftaxS>trngS(1) & fftaxS<=trngS(2) & valFFTaxSM & valFFTaxSA & valFFTaxSD);
        Nval15(bx)    = length(thisFFTaxSbin);
        if (length(thisFFTaxSbin)>20)
            fftdat15(:,bx)       = mean( powfft(:,thisFFTaxSbin),2);
            ma15dat(3,bx)        = mean(  MAmax(1,thisFFTaxSbin),2);
            ma15dat(4,bx)        = mean(  MAmax(2,thisFFTaxSbin),2);
        end
    end
    
    %% =============================
    FFT15            = struct();
    FFT15.srate      = 60/binDurMin;
    FFT15.data       = fftdat15;
    
    FFT15.setname    = 'FFT15';
    FFT15.icawinv    = [];
    FFT15.icaweights = [];
    FFT15.icasphere  = [];
    FFT15.icaact     = [];
    FFT15.chanlocs   = [];
    FFT15.trials     = [];
    FFT15.nbchan     = size(FFT15.data,1);
    FFT15.pnts       = size(FFT15.data,2);
    FFT15.filename   = [args.date, '_fft15_ch_' num2str(chnd)];
    FFT15.xmin       = 0;
    FFT15.xmax       = FFT15.pnts/FFT15.srate;
    FFT15.times      = (0:(FFT15.pnts-1))/FFT15.srate  + startBinH;
    pop_saveset(FFT15, 'filename', FFT15.filename, 'filepath',args.picDir );
     
    MA15             = FFT15;
    MA15.filename    = [args.date, '_ma15_ch_' num2str(chnd)];
    MA15.data        = ma15dat;
    MA15.nbchan      = size(MA15.data,1);
    MA15.pnts        = size(MA15.data,2);
    pop_saveset(MA15, 'filename', MA15.filename, 'filepath',args.picDir );
    
end

output_args = -1;
end

