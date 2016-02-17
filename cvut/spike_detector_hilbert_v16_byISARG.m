function [out,MARKER,envelope,background,discharges,envelope_pdf]=spike_detector_hilbert_v16_byISARG(d,fs,settings)

% Mutichannel spike detector using statistical description of bacground
% activity to finding of spike suspected section. Algorithm resample signal
% to 200 Hz, detect fast discharges (spike, poly-spike) in defined
% frequenci band. Algorithm contains main hum filtering (50 Hz DEFAULT).
% Parallel computing is supported (matlabpool). Detailed description is in
% paper: DOI: 10.1007/s10548-014-0379-1
% "Detection of Interictal Epileptiform Discharges Using Signal Envelope Distribution Modelling: Application to Epileptic and Non-Epileptic Intracranial Recordings"
%
% _________________________________________________________________________
% recomended calling:
% [...]=spike_detector_hilbert_v16_byISARG(d,fs);
% 
% or
% [...]=spike_detector_hilbert_v16_byISARG(d,fs,settings);
% _________________________________________________________________________
%
%
% inputs: -----------------------------------------------------------------
% d ... signal (time x channel)
% fs ... sampling frequency (Hz)
% settings ... string options (exmple: settings='-bl 10 -bh 60 -k1 3.0')
%   -bl low frequency of filtering ('-bl 10' DEFAULT)
%   -bh high frequency of filtering ('-bl 10' DEFAULT)
%   -ft filter type: 1-Chebyshev II, 2-Butterworth, 3-FIR ('-ft 1' DEFAULT)
%   -k1 threshold value for obvious spike decision ('-k1 3' DEFAULT)
%   -k2 defines ambiguous spike treshold. Ambiguous 
%       spike is accepted, when simultaneous obvious detection is in other 
%       channel k1 > k2 (not used in DEFAULT) 
%   -w winsize - size of segment in sample around spike for background
%                   definition (5*fs DEFAULT)
%   -n noverlap - overlap of segment in sample or absolute (novelap=0.9 ~ 90%)
%                   (4*fs DEFAULT)
%   -buf buffer ... value in second of length subsection of signal for batch
%               analyse. High value with combination of parallel computing, 
%               many analyzed channels is memory-comsuming. Minimaly 10x 
%               winsize time is recomended. ('-buf 300' DEFAULT)
%   -h main_hum_freq ... main hum frequency in Hz ('-h 50' DEFAULT)
%   -b low hrequency boundary of beta activity detection beta-25 Hz ('-b 15' recomended)
%       Inf - beta-activity detector off ('-b Inf' DEFAULT)
%   -bw length of tested segment on beta-activity in second
%   -br autoregresive model order of beta-activity detector ('-br 12' DEAFAULT)
%   -dt time between spikes in events
%   -pt polyspike union time: spike in same channel nearest then time will
%        be united
%
%         
%
% outputs:-----------------------------------------------------------------
% MARKER.d ... resampled signal "d" to 200 Hz
% MARKER.fs ... new sampling frequency
% MARKER.M ... zero matrix with same size as signal "d", where 1 indicate 
%              obviuos spike, 0.5 ambiguous 
% envelope ... instant hilbert envelope of filtered signal "d"
% bacground ... threshold curves with same size as signal "d"
%               background(:,:,1) ... for high threshold (k~k1)
%               background(:,:,2) ... for low threshol (k~k2)
% disharges ... structure of multichannel event describes occurence of
%               spikes suddenly 
%           discharges.MV ... matrix of type (n x ch) 
%                             (1-obvious spike, 0.5-ambiguous)
%           discharges.MA ... matrix of max. amplitude of envelope above 
%                             backround (n x ch)
%           discharges.MP ... start position of multichannel event 
%           discharges.MD ... duration of event
%           discharges.MW ... statistical significance "CDF"
%           discharges.MPDF ... probability of occurence
% out ... structure describing each detected spike. Usefull for export to
%         GDF files using biosig toolbox (biosig.sourceforge.net)
%           out.pos ... position (second)
%           out.dur ... duration (second) - fixed 1/fs 
%           out.chan ... channel
%           out.con ... type (1-obvious 0.5-ambiguous)
%           out.weigth ... statistical significance "CDF"
%           out.pdf ... statistical significance "PDF"
%
% required toolboxes:
%   Signal Processing Toolbox, Statistic Toolbox, 
%
% recomended toolboxes:
%   Parallel Computing Toolbox
%
% MAKING BY ISARG - Radek Janca, 19.9.2013



% settings defaults:------------------------------------------------------
band_low=10; % (-fl)
band_high=60; % (-fh)
k1=3.65; %(-k1)
k2=[]; % (-k2)
winsize=5*fs;   % (-w)
noverlap=4*fs;  % (-n)
buffering=300; % (-buf)
main_hum_freq=50; % (-h)
beta=Inf;    % (-b)
beta_win=20; % (-bw)
beta_AR=12; % (-br)
f_type=1; % 1-cheby2, 2-but, 3-fir (-ft)
discharge_tol=0.005; % (-dt)
polyspike_uniom_time=0.12; % (-pt)

if nargin>2; 
    poz=strfind(settings,'  '); % multi spaces removing
    settings(poz)=[];
    SET=textscan(settings,'%s %f','delimiter',' '); % paremeters reading
else    
    SET{1,1}=[]; 
end




for i=1:size(SET{1,1},1)
    prefix=SET{1,1}{i,:};
    prefix(1)=[];
    
    switch prefix
        case 'fl'
            band_low=SET{1,2}(i);
        case 'fh'
            band_high=SET{1,2}(i);
        case 'k1'
            k1=SET{1,2}(i);
        case 'k2'
            k2=SET{1,2}(i);
        case 'w'
            winsize=SET{1,2}(i);
        case 'n'
            noverlap=SET{1,2}(i);
        case 'buf'
            buffering=SET{1,2}(i);
        case 'h'
            main_hum_freq=SET{1,2}(i);
        case 'b'
            beta=SET{1,2}(i);
        case 'bw'
            beta_win=SET{1,2}(i);
        case 'br'
            beta_AR=SET{1,2}(i);
        case 'ft'
            f_type=SET{1,2}(i);
        case 'dt'
            discharge_tol=SET{1,2}(i);
        case 'pt'
            polyspike_uniom_time=SET{1,2}(i);
    end
end


% signal resampling to 200 Hz -------------------------------------------
if fs>200
    for ch=1:size(d,2)
        d_res(:,ch)=resample(d(:,ch),200,fs,100);
    end
    winsize=winsize/fs;
    noverlap=noverlap/fs;
    fs=200; d=d_res; clear d_res;
    winsize=round(winsize*fs);
    noverlap=round(noverlap*fs);
end

MARKER.fs=fs;

% signal buffering ---------------------------------------------------
N_seg=floor(size(d,1)/(buffering*fs));
if N_seg<1; N_seg=1; end
T_seg=round(size(d,1)/N_seg/fs);

% indexs of segments with two-side overlap
index_start=1:T_seg*fs:size(d,1);
if length(index_start)>1
    index_start(2:end)=index_start(2:end)-(3*winsize);
    
    index_stop=index_start+T_seg*fs+2*(3*winsize)-1;
    index_stop(1)=index_stop(1)-(3*winsize);
    index_stop(end)=size(d,1);
    
    if index_stop(end)-index_start(end)<T_seg*fs
        index_start(end)=[];
        index_stop(end-1)=[];
    end
else
    index_stop=size(d,1);
end


% detection calling ---------------------------------------------------
out.pos=[]; % spike position
out.dur=[]; % spike duration - fix value 5 ms
out.chan=[]; % spike in channel
out.con=[]; % spike condition
out.weight=[]; % spike weight "CDF"
out.pdf=[]; % spike probability "PDF"

discharges.MV=[]; % spike type 1-obvious, 0.5- ambiguous
discharges.MA=[]; % max. amplitude of envelope above backround
discharges.MP=[]; % event start position
discharges.MD=[]; % event duration 
discharges.MW=[]; % statistical weight
discharges.MPDF=[]; % statistical weight

fprintf(1,'progress:   0 %%')
for i=1:length(index_stop)
    % subsection signal spike detection ----------------
    [sub_MARKER,sub_envelope,sub_background,sub_discharges,sub_out,sub_envelope_pdf]=spike_detector(d(index_start(i):index_stop(i),:),fs,[band_low band_high],k1,k2,winsize,noverlap,f_type,main_hum_freq,beta,beta_win,beta_AR,polyspike_uniom_time,discharge_tol);
    
    % progress stat
    percent=num2str(100*i/length(index_stop),'%3.0f');
    percent=[repmat(' ',1,3-length(percent)) , percent,' %%'];
    fprintf(1,['\b\b\b\b\b' percent])
    
    
    % conecting of subsection ----------------
    MARKER.M(index_start(i)+(i>1)*(3*winsize):index_stop(i)-(i<length(index_stop))*(3*winsize),:)=sub_MARKER.M(1+(i>1)*(3*winsize):end-(i<length(index_stop))*(3*winsize),:);
    MARKER.d(index_start(i)+(i>1)*(3*winsize):index_stop(i)-(i<length(index_stop))*(3*winsize),:)=sub_MARKER.d(1+(i>1)*(3*winsize):end-(i<length(index_stop))*(3*winsize),:);
    envelope(index_start(i)+(i>1)*(3*winsize):index_stop(i)-(i<length(index_stop))*(3*winsize),:)=sub_envelope(1+(i>1)*(3*winsize):end-(i<length(index_stop))*(3*winsize),:);
    background(index_start(i)+(i>1)*(3*winsize):index_stop(i)-(i<length(index_stop))*(3*winsize),:,:)=sub_background(1+(i>1)*(3*winsize):end-(i<length(index_stop))*(3*winsize),:,:);
    envelope_pdf(index_start(i)+(i>1)*(3*winsize):index_stop(i)-(i<length(index_stop))*(3*winsize),:)=sub_envelope_pdf(1+(i>1)*(3*winsize):end-(i<length(index_stop))*(3*winsize),:);
    
    % removing of two side overlap detections
    
    if ~isempty(sub_out.pos)
        if length(index_stop)>1
 
%             if i==1
%                 idx_evt=sub_out.pos>T_seg;
%                 idx_disch=sub_discharges.MP(:,1)>T_seg;
%             elseif i>1 && i<length(index_stop)
%                 idx_evt=sub_out.pos<=(3*winsize)/fs | sub_out.pos>T_seg+(3*winsize)/fs;
%                 idx_disch=sub_discharges.MP(:,1)<=(3*winsize)/fs | sub_discharges.MP(:,1)>T_seg+(3*winsize)/fs;
%             elseif i==length(index_stop)
%                 idx_evt=sub_out.pos<=(3*winsize)/fs;
%                 idx_disch=sub_discharges.MP(:,1)<=(3*winsize)/fs;
%             end
            
            idx_evt=sub_out.pos<((i>1)*(3*winsize)/fs) | sub_out.pos>((index_stop(i)-index_start(i))-(i<length(index_stop))*(3*winsize))/fs;
            idx_disch=min(sub_discharges.MP,[],2)<((i>1)*(3*winsize)/fs) | min(sub_discharges.MP,[],2)>((index_stop(i)-index_start(i))-(i<length(index_stop))*(3*winsize))/fs;
            
            sub_out.pos(idx_evt)=[];
            sub_out.dur(idx_evt)=[];
            sub_out.chan(idx_evt)=[];
            sub_out.con(idx_evt)=[];
            sub_out.weight(idx_evt)=[];
            sub_out.pdf(idx_evt)=[];
            
            sub_discharges.MV(idx_disch,:)=[];
            sub_discharges.MA(idx_disch,:)=[];
            sub_discharges.MP(idx_disch,:)=[];
            sub_discharges.MD(idx_disch,:)=[];
            sub_discharges.MW(idx_disch,:)=[];
            sub_discharges.MPDF(idx_disch,:)=[];
            
        end
    end
    
    out.pos=[out.pos; sub_out.pos+index_start(i)/fs-1/fs];
    out.dur=[out.dur; sub_out.dur];
    out.chan=[out.chan; sub_out.chan];
    out.con=[out.con; sub_out.con];
    out.weight=[out.weight; sub_out.weight];
    out.pdf=[out.pdf; sub_out.pdf];
    
    discharges.MV=[discharges.MV; sub_discharges.MV];
    discharges.MA=[discharges.MA; sub_discharges.MA];
    discharges.MP=[discharges.MP; sub_discharges.MP+index_start(i)/fs-1/fs];
    discharges.MD=[discharges.MD; sub_discharges.MD];
    discharges.MW=[discharges.MW; sub_discharges.MW];
    discharges.MPDF=[discharges.MPDF; sub_discharges.MPDF];
    
end
fprintf(1,'\n')

end







function [MARKER,envelope,background,discharges,out,envelope_pdf,markers_high,markers_low]=spike_detector(d,fs,bandwidth,k1,k2,winsize,noverlap,f_type,main_hum_freq,beta,beta_win,beta_AR,polyspike_uniom_time,discharge_tol)


%--------------------------------------------------------------------------
% segmentation index
%--------------------------------------------------------------------------
if noverlap<1
    index=1:round(winsize*(1-noverlap)):size(d,1)-winsize+1; % indexy zaèátkù oken
else
    index=1:winsize-noverlap:size(d,1)-winsize+1;
end


df=filtering(bandwidth,d,fs,f_type);
df=filt50Hz(df,fs,main_hum_freq);

if exist('parfor')==5 % If you don't have "Parallel Computing Toolbox", only standard for-cycle will be performed
        
    parfor ch=1:size(df,2)
        [envelope(:,ch),markers_high(:,ch),markers_low(:,ch),background(:,ch,:),envelope_cdf(:,ch),envelope_pdf(:,ch)]=one_channel_detect(df(:,ch),fs,index,winsize,k1,k2,polyspike_uniom_time); % ~ = prah_int(:,ch,:)
    end
    
else % If you do not have "Parallel Computing Toolbox", only standard for-cycle will be performed
        
    for ch=1:size(df,2)
        [envelope(:,ch),markers_high(:,ch),markers_low(:,ch),background(:,ch,:),envelope_cdf(:,ch),envelope_pdf(:,ch)]=one_channel_detect(df(:,ch),fs,index,winsize,k1,k2,polyspike_uniom_time); % ~ = prah_int(:,ch,:)
    end  
end

% first and last second is not analyzed (filter time response etc.)
markers_high([1:fs, end-fs+1:end],:)=false;
markers_low([1:fs, end-fs+1:end],:)=false;

% beta activity detection
if beta<fs/2 && beta_win>0
    M_beta=beta_detect(d,fs,beta,beta_win,beta_AR);
    markers_high(M_beta)=false;
    markers_low(M_beta)=false;
end

ovious_M=sum(markers_high,2)>0;

% obvious spike events output
out.pos=[];
out.dur=[];
out.chan=[];
out.con=[];
out.weight=[];
out.pdf=[];

t_dur=0.005;
for ch=1:size(markers_high,2)
    if sum(markers_high(:,ch))>0
        idx=find(markers_high(:,ch));
        for i=1:length(idx)
            out.pos=[out.pos; idx(i)/fs];
            out.dur=[out.dur; t_dur]; % fix 50 ms
            out.chan=[out.chan; ch];
            out.con=[out.con; 1];
            out.weight=[out.weight; envelope_cdf(idx(i),ch)];
            out.pdf=[out.pdf; envelope_pdf(idx(i),ch)];
        end
    end
end

if ~(k2==k1)
    % ambiguous spike events output
    for ch=1:size(markers_low,2)
        if sum(markers_low(:,ch))>0
            idx=find(markers_low(:,ch));
            for i=1:length(idx)
                if markers_high(idx(i),ch) % is true
                    continue
                elseif sum(ovious_M(round(idx(i)-0.01*fs:idx(i)-0.01*fs)))>0
                    out.pos=[out.pos; idx(i)/fs];
%                     out.dur=[out.dur; 1/fs];
                    out.dur=[out.dur; t_dur]; %
                    out.chan=[out.chan; ch];
                    out.con=[out.con; 0.5];
                    out.weight=[out.weight; envelope_cdf(idx(i),ch)];
                    out.pdf=[out.pdf; envelope_pdf(idx(i),ch)];
                end
            end
        end
    end
end



%--------------------------------------------------------------------------
% making M stack pointer of events
%--------------------------------------------------------------------------
M=zeros(size(df));
for k=1:size(out.pos,1)
    M(round(out.pos(k)*fs:out.pos(k)*fs+discharge_tol*fs),out.chan(k))=out.con(k);
end

%--------------------------------------------------------------------------
% definition of multichannel events vectors
%--------------------------------------------------------------------------

point(:,1)=find(diff([0;(sum(M,2))>0])>0);
point(:,2)=find(diff([(sum(M,2))>0;0])<0);

discharges.MV=[]; % spike type 1-obvious, 0.5- ambiguous
discharges.MA=[]; % max. amplitude of envelope above backround
discharges.MP=[]; % event start position
discharges.MD=[]; % event duration 
discharges.MW=[]; % statistical weight
discharges.MPDF=[]; % statistical weight

for k=1:size(point,1)
    seg=M(point(k,1):point(k,2),:);
    mv=max(seg,[],1);
    
%     seg=df(point(k,1):point(k,2),:);
    seg=envelope(point(k,1):point(k,2),:)-(background(point(k,1):point(k,2),:,1)/k1);
    ma=max(abs(seg),[],1);
    
%     seg=1-envelope_pdf(point(k,1):point(k,2),:);
    seg=envelope_cdf(point(k,1):point(k,2),:);
%     mw=max(seg.*(M(point(k,1):point(k,2),:)>0),[],1);
    mw=max(seg,[],1);
    
    seg=envelope_pdf(point(k,1):point(k,2),:);
    mpdf=max(seg.*(M(point(k,1):point(k,2),:)>0),[],1);
    
    discharges.MV=[discharges.MV;mv];
    discharges.MA=[discharges.MA;ma];
    discharges.MW=[discharges.MW;mw];
    discharges.MPDF=[discharges.MPDF;mpdf];
    % discharges.MP=[discharges.MP; repmat(point(k,1),1,size(d,2))/fs];
    discharges.MD=[discharges.MD; repmat(point(k,2)-point(k,1),1,size(d,2))/fs];
    % ------ precision position in discharges.MP -------
    [row,col]=find(M(point(k,1):point(k,2),:)>0);
    mp=nan(1,size(d,2));
    mp(col)=row+point(k,1)-1;
    discharges.MP=[discharges.MP;mp/fs];
    
end

MARKER.M=M;
MARKER.d=d;
end






%% -------------------------------------------------------------------------
% subfunctions
%--------------------------------------------------------------------------


function d=filt50Hz(d,fs,hum_fs)

if nargin<3
    hum_fs=50;
end

if min(size(d))==1
   d=d(:); 
end

R = 1; r = 0.985;


f0 = hum_fs:hum_fs:fs/2; % Hz


for i=1:length(f0)
    b = [1 -2*R*cos(2*pi*f0(i)/fs) R*R];
    a = [1 -2*r*cos(2*pi*f0(i)/fs) r*r];
    for ch=1:size(d,2)
        d(:,ch)=filtfilt(b,a,d(:,ch));
    end
end
end


%--------------------------------------------------------------------
function df=filtering(bandwidth,d,fs,type)

% bandpass filtering
switch type
    case 1
        % IIR-cheby2
        % low pass
        Wp = 2*bandwidth(2)/fs; Ws = 2*bandwidth(2)/fs+ 0.1;
        Rp = 6; Rs = 60;
        [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
        [bl,al] = cheby2(n,Rs,Ws);
        
        % high pass
        Wp = 2*bandwidth(1)/fs; Ws = 2*bandwidth(1)/fs- 0.05;
        Rp = 6; Rs = 60;
        [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
        [bh,ah] = cheby2(n,Rs,Ws,'high');
    case 2
        % IIR butterworth
        [bl,al]=butter(5,2*bandwidth(2)/fs);
        [bh,ah]=butter(5,2*bandwidth(1)/fs,'high');
        
    case 3
        % FIR
        % low pass
        bl = fir1(fs/2,2*bandwidth(2)/fs); al=1;
        
        % high pass
        bh = fir1(fs/2,2*bandwidth(1)/fs,'high'); ah=1;
end


if exist('parfor')==5 % If you don't have "Parallel Computing Toolbox", only standard for-cycle will be performed
    parfor ch=1:size(d,2); df(:,ch)=filtfilt(bh,ah,d(:,ch)); end
    if bandwidth(2)==fs/2; return; end
    parfor ch=1:size(d,2); df(:,ch)=filtfilt(bl,al,df(:,ch)); end
    
else % If you do not have "Parallel Computing Toolbox", only standard for-cycle will be performed
    for ch=1:size(d,2); df(:,ch)=filtfilt(bh,ah,d(:,ch)); end
    if bandwidth(2)==fs/2; return; end
    for ch=1:size(d,2); df(:,ch)=filtfilt(bl,al,df(:,ch)); end
    
end
end

%--------------------------------------------------------------------
function [envelope,markers_high,markers_low,prah_int,envelope_cdf,envelope_pdf]=one_channel_detect(d,fs,index,winsize,k1,k2,polyspike_uniom_time)


envelope=abs(hilbert(d)); % Hilbert's envelope (intense envelope)

for k=1:length(index) % for each segment
    
    segment=envelope(index(k):index(k)+winsize-1);
    segment(segment<=0)=[];
    
    % estimation of segment's distribution using MLE 
%     [phat(k,:),~] = mle(segment(:)','distribution','logn'); % paremeters estimation 
    phat(k,1)=mean(log(segment)); % median
    phat(k,2)=std(log(segment));
%     
end

r=size(envelope,1)/length(index);

n_average=winsize/fs;
phat(:,1)=filtfilt(ones(round(n_average*fs/r),1)/(round(n_average*fs/r)),1,phat(:,1));
phat(:,2)=filtfilt(ones(round(n_average*fs/r),1)/(round(n_average*fs/r)),1,phat(:,2));



% interpolation of thresholds value to threshold curve (like backround)
phat_int=[];
if size(phat,1)>1
    phat_int(:,1) = interp1(index+round(winsize/2),phat(:,1),(index(1):index(end))+round(winsize/2),'spline');
    phat_int(:,2) = interp1(index+round(winsize/2),phat(:,2),(index(1):index(end))+round(winsize/2),'spline');
    
    phat_int=[ones(floor(winsize/2),size(phat,2)).*repmat(phat_int(1,:),floor(winsize/2),1); phat_int ; ones(size(envelope,1)-(length(phat_int)+floor(winsize/2)),size(phat,2)).*repmat(phat_int(end,:),size(envelope,1)-(length(phat_int)+floor(winsize/2)),1)];   
else
    phat_int=phat.*ones(size(d,1),2);
end


lognormal_mode= exp(phat_int(:,1)-phat_int(:,2).^2);
lognormal_median=exp(phat_int(:,1));
prah_int(:,1)=k1*(lognormal_mode+lognormal_median);
if ~(k2==k1)
    prah_int(:,2)=k2*(lognormal_mode+lognormal_median);
end



envelope_cdf=0.5+0.5*erf((log(envelope)-phat_int(:,1))./sqrt(2*phat_int(:,2).^2)); % CDF of lognormal distribution
envelope_pdf=exp(-0.5 .* ((log(envelope) - phat_int(:,1))./phat_int(:,2)).^2) ./ (envelope  .* phat_int(:,2) .* sqrt(2*pi)); % PDF of lognormal distribution

% -------------------------------------------------------------------------
% detection of obvious and ambiguous spike
% -------------------------------------------------------------------------
markers_high=local_maxima_detection(envelope,prah_int(:,1),fs,polyspike_uniom_time);
if ~(k2==k1)
    markers_low=local_maxima_detection(envelope,prah_int(:,2),fs,polyspike_uniom_time);
else
    markers_low=markers_high;
end

end

%--------------------------------------------------------------------
function marker1=local_maxima_detection(envelope,prah_int,fs,polyspike_uniom_time)


marker1=zeros(size(envelope));
marker1(envelope(:)>prah_int(:))=1; % crossing of high threshold

point=[];
point(:,1)=find(diff([0;marker1])>0); % strat crossing
point(:,2)=find(diff([marker1;0])<0); % end crossing

marker1=false(size(envelope));
for k=1:size(point,1)
    
    % detection of local maxima in section which crossed threshold curve
    if point(k,2)-point(k,1)>2
        seg=envelope(point(k,1):point(k,2));
        seg_s=diff(seg);
        seg_s=sign(seg_s); 
        seg_s=find(diff([0;seg_s])<0); % positions of local maxima in the section
        
        marker1(point(k,1)+seg_s-1)=true;
    elseif point(k,2)-point(k,1)<=2
        seg=envelope(point(k,1):point(k,2));
        [~,s_max]=max(seg); % positions of local maxima in the section
        marker1(point(k,1)+s_max-1)=true;
    end
end


% union of section, where local maxima are close together <(1/f_low + 0.02 sec.)~ 120 ms
pointer=find(marker1==true); % index of local maxima
state_previous=false;
for k=1:length(pointer)
    if ceil(pointer(k)+polyspike_uniom_time*fs)>size(marker1,1)
        seg=marker1(pointer(k)+1:end);
    else
        seg=marker1(pointer(k)+1:ceil(pointer(k)+polyspike_uniom_time*fs));
    end
    
    if state_previous
        if sum(seg)>0
            state_previous=true;
        else
            state_previous=false;
            marker1(start:pointer(k))=true;
        end
        
    else
        
        if sum(seg)>0
            state_previous=true;
            start=pointer(k);
        end
    end
end

% finding of the highes maxima of the section with local maxima
point=[];
point(:,1)=find(diff([0;marker1])>0); % start 
point(:,2)=find(diff([marker1;0])<0); % end 

% local maxima with gradient in souroundings
for k=1:size(point,1)
    if point(k,2)-point(k,1)>1
        lokal_max=pointer(pointer>=point(k,1) & pointer<=point(k,2)); % index of local maxima
        
        marker1(point(k,1):point(k,2))=false;
        
        lokal_max_val=envelope(lokal_max); % envelope magnitude in local maxima
        lokal_max_poz=(diff(sign(diff([0;lokal_max_val;0]))<0)>0);
        
        marker1(lokal_max(lokal_max_poz))=true;
    end
end


% % local maxima with gradient in souroundings
% for k=1:size(point,1)
%     if point(k,2)-point(k,1)>1
%         lokal_max=pointer(pointer>=point(k,1) & pointer<=point(k,2)); % index of local maxima
%         
%         marker1(point(k,1):point(k,2))=false;
%         
%         lokal_max_val=envelope(lokal_max); % envelope magnitude in local maxima
%         lokal_max(lokal_max_val<0.50*max(lokal_max_val))=[];
%         lokal_max_val=envelope(lokal_max);
%         
%         
%         lokal_max_poz=(diff(sign(diff([0;lokal_max_val;0]))<0)>0);
%         
%         
%         
%         
%         marker1(lokal_max(lokal_max_poz))=true;
%     end
% end

% % local maxima 
% for k=1:size(point,1)
%     if point(k,2)-point(k,1)>1
%         [~,max_poz]=max(envelope(point(k,1):point(k,2)));
%         marker1(point(k,1):point(k,2))=false;
%         
%         marker1(point(k,1)+max_poz-1)=true;
%     end
% end

end


function M=beta_detect(d,fs,beta,winsize,beta_AR)


M=[];
winsize=winsize*fs;
noverlap=round(0.5*winsize);

index=1:winsize-noverlap:size(d,1)-winsize+1;
if isempty(index)
   index=1;
   winsize=size(d,2);
end

[bb,aa]=butter(4,2*30/fs);
parfor ch=1:size(d,2)
    MM=[];
%     winsize=5*fs;
%     noverlap=2.5*fs;
    
    for i=1:length(index)
        seg=d(index(i):index(i)+winsize-1,ch);
        seg=filtfilt(bb,aa,seg);
        
        a=lpc(seg-mean(seg),beta_AR);
        [h,f]=freqz(1,a,[],fs);
        h=abs(h);
        poz=f(diff([0; sign(diff([0; h]))])<0);
        
        MM(i)=sum(poz<25 & poz>beta)>0;
    end
    MM=[MM MM(end)];
    M(:,ch)=interp1([index size(d,1)],MM,1:size(d,1),'nearest');
end
M=logical(M);
end

