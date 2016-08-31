%% Tutorial of spike detector using
clear all; close all; clc;

%% loading a patient mat-file
load('d:\eeg\motol\pacienti\p083 Kol VT10\VT10_2015-05-19_10-00_001_X_ppa.mat');
% ----signal----
% d... iEEG records: channels are in columns
% fs... sampling frequency
%
% ----labels----
% DE... output structure of detected events (recalculated below)
%       DE.pos ... vector of events positions in seconds
%       DE.chan ... vector of events channel-positions corresponding to DE.pos
%       DE.con ... vector of conditions (obvious 1, ambiguous 0.5) corresponding to DE.pos
% REV_01... output structure of first reviewer
% REV_02... output structure of second reviewer
% REV_03... output structure of third reviewer
% REV_GS... output structure of "gold standard" - at least two reviewers
%           intersection, where
%           obvious+ambiguous spike = obvious gold stadard (con 1)
%           ambiguous+ambiguous spike = ambiguous gold standard (con 0.5)


%% spike detection in default mode, for more setting and outputs see help:
% signal decimating to fs=200 Hz
% signal filtering 10-60 Hz (ChebyII)
% segmentation by 5 second window with 80% noverlap
% sensitivity koef k=3;
% 50 Hz main hum removing 

%--------------------------------------------------------------------------
[DE]=spike_detector_hilbert_v16_byISARG(d,fs);
%--------------------------------------------------------------------------
% DE.pos ... vector of events positions in second
% DE.chan ... vector of events channel-positions corresponding to DE.pos
% DE.* ... other variables is not yet used
 

%%
% Example of visualization - channel 5

M=zeros(size(d)); % making of signal markers
for i=1:length(DE.pos)
   M(round(DE.pos(i)*fs):round((DE.pos(i)+0.02)*fs),DE.chan(i))=DE.weight(i);
end

figure(1)
subplot(211)
plot(d(:,5),'k')
hold on
stairs(M(:,5)*max(d(:,5)),'c')
title('all detected events of channel 5')

%% comparison to "gold standard"
% kamil 31.8.2016 - dal to nefunguje, nemame REV promenne s rucni analyzou
%--------------------------------------------------------------------------
[total_stat,chan_stat,TP,FP,FN]=label_stat(DE,REV_GS,0.1);
%--------------------------------------------------------------------------
% see command window
% outputs can be used for visualization of detection




% Example of visualization TP - channel 5
M_TP=zeros(size(d));
for i=1:length(TP.pos)
   M_TP(round(TP.pos(i)*fs):round((TP.pos(i)+0.02)*fs),TP.chan(i))=1; 
end

subplot(212)
plot(d(:,5),'k')
hold on
stairs(M_TP(:,5)*max(d(:,5)),'g')

%%
% Example of visualization FP - channel 5
M_FP=zeros(size(d));
for i=1:length(FP.pos)
   M_FP(round(FP.pos(i)*fs):round((FP.pos(i)+0.02)*fs),FP.chan(i))=1; 
end

stairs(M_FP(:,5)*max(d(:,5)),'r')

%%
% Example of visualization FN - channel 5
M_FN=zeros(size(d));
for i=1:length(FN.pos)
   M_FN(round(FN.pos(i)*fs):round((FN.pos(i)+0.02)*fs),FN.chan(i))=1; 
end

stairs(M_FN(:,5)*max(d(:,5)),'b')
title('statistic of channel 5')

legend('iEEG','TP','FP','FN')