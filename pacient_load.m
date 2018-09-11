function [ E ] = pacient_load( nick,testname,filename ,frequencies,classname,channels,loadall)
%PACIENT_LOAD nacte soubor pacienta, bud ze zdrojvych data (pokud neni udano filename) nebo uz ulozeny soubor
%   zdrojova data rovnou zpracuje, rozepochuje atd 
%   E= pacient_load( nick,testname,filename ,frequencies,classname,channels)
if ~exist('frequencies','var') || isempty(frequencies), frequencies = []; end
if ~exist('classname','var') || isempty(classname), classname = []; end
if ~exist('channels','var') || isempty(channels), channels = []; end
if ~exist('loadall','var') || isempty(loadall), loadall = 1; end
if strcmp(testname,'aedist')
    pacienti = pacienti_aedist(); %nactu celou strukturu pacientu
    setup = setup_aedist(0); %nactu nastaveni aedist  
elseif strcmp(testname,'menrot')
    pacienti = pacienti_menrot(); %nactu celou strukturu pacientu
    setup = setup_menrot(0); %nactu nastaveni aedist  
elseif strcmp(testname,'ppa')
    pacienti = pacienti_ppa(); %nactu celou strukturu pacientu
    setup = setup_ppa(0); %nactu nastaveni ppa  
else
    error('zatim pracuji jen s aedist, menrot a ppa');
end
nalezen = false;
for p = 1:numel(pacienti)
    if strfind(pacienti(p).folder,nick)
        nalezen = true;
        break; %nasel jsem pacienta
    end        
end
if ~nalezen
    disp(['pacient nenalezen: ' nick]);    
    E = [];
    return;
end
disp(['loading pacient ' pacienti(p).folder ]);

if exist('filename','var') && ~isempty(filename) %pokud filename ~= []
    %nactu si existujici CHilbert nebo CiEEGdata
    fullfilename = [setup.basedir pacienti(p).folder '\' setup.subfolder '\' filename];
    if exist(fullfilename,'file')==2
        if strfind(filename,'CHilbert') 
            E = CHilbert(fullfilename,loadall);
        elseif strfind(filename,'CMorlet') 
            E = CMorlet(fullfilename,loadall);
        else
            E = CiEEGData(fullfilename);
        end      
    else
        E = [];
        disp(['soubor neexistuje: ' fullfilename]);
    end
else    
    %vytvorim novy objekt CiEEGData z existujicich EEG dat
    load([ setup.basedir pacienti(p).folder '\' pacienti(p).data]);
    if ~exist('mults','var'), mults = [];  end
    if ~isempty(frequencies) 
        if exist('classname','var') && strcmp(classname,'cmorlet')
            E = CMorlet(d,tabs,fs,mults);
        else
            E = CHilbert(d,tabs,fs,mults);
        end
    else
        E = CiEEGData(d,tabs,fs,mults);
    end

    %header
    load([ setup.basedir pacienti(p).folder '\' pacienti(p).header]);
    E.GetHHeader(H);

    %epievents
    filename = [ setup.basedir pacienti(p).folder '\' setup.subfolder '\' pacienti(p).epievents];   
    E.GetEpiEvents(getepievents(filename));       

    %vyradim kanaly
    if ~isfield(pacienti(p),'rjch'), return; end
    E.RejectChannels(pacienti(p).rjch);
    
    %frekvencni analyza
    if  ~isempty(frequencies) 
        if ~exist('channels','var'), channels = 1:E.channels; end
        E.ChangeReference('b');
        E.PasmoFrekvence(frequencies,channels);
    end
    %vytvorim epochy
    filename = [ setup.basedir pacienti(p).folder '\' setup.subfolder '\' pacienti(p).psychopy];    
    E.ExtractEpochs(getpsychopydata(filename,testname),setup.epochtime,setup.baseline);

    %vyradim epochy
    if ~isempty(pacienti(p).rjepoch)
        load([setup.basedir pacienti(p).folder '\' setup.subfolder '\' pacienti(p).rjepoch]);
    end
    if ~exist('RjEpochCh','var'), RjEpochCh = []; end %pokud neexistuje RjEpochCh - starsi data
    if ~exist('RjEpoch','var'), RjEpoch = []; end %pokud neexistuje RjEpoch - starsi data
    E.RejectEpochs(RjEpoch, RjEpochCh); %uz drive ulozene vyrazene epochy
    %E.RjEpochsEpi([],0);
    %E.RjEpochsEpi(30);
    if ~isempty(frequencies)  %defaultne delam statistiku jen pro frekvencni data
        E.ResponseSearch(0.1,setup.stat_kats); %vypocitam statistiku
    end
end
end
function DE = getepievents(filename)
    if exist(filename,'file')~=2
            disp(['no epievents: ' filename]);    
            DE = [];
    else
            load(filename);        
    end
end
function psychopy = getpsychopydata(filename,testname)
    if exist(filename,'file')~=2
        disp(['no psychopy data:' pacienti(p).psychopy]);
        psychopy = []; 
    else
        load(filename);
        if strcmp(testname,'aedist')
            psychopy = aedist; %matrix s psychopy daty
        elseif strcmp(testname,'menrot')
            psychopy = menrot; %matrix s psychopy daty
        else
            psychopy = ppa;
        end
    end
end


