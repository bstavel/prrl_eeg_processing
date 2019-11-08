clear all
%addpath(genpath('~/Dropbox/Sarah/MATLAB/EEG Data/eeglab13_6_5b/'));
addpath(genpath('/eeglab11_0_5_4b/')); % older version

%%
subno=15;
load(['../eeg_data/processed_data/Oddball_',num2str(subno),'/Oddball_',num2str(subno),'_interpolated_rereferenced_ica_filtered_FINAL'])

%%

data=EEG.data;
times=EEG.times;
epoch=EEG.epoch;

%% trim the epoch

T = find(times>-200 & times<800);
times=times(T);
data=data(:,T,:);

%% subsample the data
times=times(1:4:end);
data=data(:,1:4:end,:);

%% baseline correct the data

baseline = find(times>-200 & times<0);

data = data-repmat(mean(data(:,baseline,:),2),[1,length(times),1]);

%% extract event types

for e=1:length(epoch)
    thisepoch=epoch(e);
    
    latency=cell2mat(thisepoch.eventlatency);
    type= cell2mat(thisepoch.eventtype(find(latency==0)));
    eventtype(e)=type;
end


%%
types=unique(eventtype);
%% sanity check
num=[];
k=0;
for t=types
    k=k+1;
    num(k,:)=[t,sum(eventtype==t)];
    
    indices{k}=find(eventtype==t);
end


%% make erps


for k =1:length(indices)
    
    ERP(:,:,k)=mean(data(:,:,indices{k}),3);
    
%     figure;
%     plot(times,squeeze(data(47,:,indices{k})))
%     
end


%% plot results


elec=47;

figure;
hold on
for k=1:length(indices)
    plot(times, ERP(elec,:,k),'linewidth',2)
end
title(EEG.chanlocs(elec).labels)
xlabel('stim-locked time (ms)')
ylabel('\mu V')
legend('high (target)','low (baseline)','distractor (oddball)')

chanlocs = EEG.chanlocs;

%% oddball effect
% plot all channels of a data epoch on the same axis and map its 
% scalp map(s) at selected latencies.

diff = ERP(:,:,2)-mean(ERP(:,:,3),3);

figure
timtopo(diff,chanlocs,'limits',[min(times),max(times)])

%% topoplot
% plot a topographic map of a scalp data field in a 2-D circular view 
%(looking down at the top of the head) using interpolation on a fine 
%cartesian grid. 

T = find(times>300 &times<400);
P300=squeeze(mean(ERP(:,T,:),2));
figure
topoplot_vect = reshape(P300(:,2),1,[]);
topoplot(topoplot_vect,chanlocs); % first arg must be vector

