%% ======================================
%  batch processing of all edf files
%

args.rawDir = '~\Dropbox\amy\Export\';
args.exp    =  '16026773';

expNumber = '16';
%expNumber = '17';

animalList16 = {'112','113','120','121','122','123'}; % two each:   WEEV(112 113) VEEV (121 122) EEEV (120 123)
animalList17 = {'053','054','057','058'};   % all eastern

if strcmp(expNumber,'16')
    animalList = animalList16;
end

if strcmp(expNumber,'17')
    animalList = animalList17;
end


for ax = 1:length(animalList)
    args.animal = animalList{ax};
    args.extn   = expNumber;
    
    clear d dir dirList
    d = dir( [args.rawDir 'M' args.animal  '-' args.extn '_EEG\*.edf'] );
    [dirList{1:length(d),1}] = deal(d.name);
    
    dirList
    
    %clear d xlsList
    %d = dir(           [args.rawDir '\Ponemah\' args.exp  '*.xls']);
    %[xlsList{1:length(d),1}] = deal(d.name);
    
    %%
    for i = 1:length(dirList)
        i
        args.file  = dirList{i};
        tmp1       = strsplit(dirList{i},'_');
        tmp2       = strsplit(tmp1{2},'.');
        args.date  = tmp2{1};
        args
        
        %disc = prep_amy( args );
        disc = prep_amy_raw15( args );
    end
end
