function plot_bad_epochs(subject, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  this function plots all the epochs for one channel, then helps you
%%  identify the epochs with max min values by channel
%%
%%  INPUTS:
%%    subject: the subject structure, with EEG and PRRL specific metadata
%%    k: number of min/max pairs to display
%%
%% OUTPUTS:
%%    datatable: table with k max and mins pairs for worrisome channels
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Pull EEG out of subject structure
    EEG = subject.EEG;

  %% Identify channels with abberant ERPs for epoch ejection
    time = EEG.times/1000; %% Q: Why divide by 1000?

    t1 = find(time == -.5);
    t2 = find(time <1.5);t2 = t2(end);

    FILT_STIMULI = EEG.data;
    subject.badchans = {};
    elec = 1;
    for electrode = 1:64
        figure(3)
        clf
        subplot(2,1,1)
        plot(time(t1:t2),squeeze(FILT_STIMULI(electrode,t1:t2,:)),'k')
        subplot(2,1,2)
        plot(time(t1:t2),squeeze(mean(FILT_STIMULI(electrode,t1:t2,:),3)),'r','linewidth',2)

        title(EEG.chanlocs(electrode).labels)
        spit = input('Investigate? y or n ', 's');
        if strcmp(spit,'y')
          subject.badchans{elec} = {EEG.chanlocs(electrode).labels}
          elec = elec + 1;
        end
    end
    pause(2);

  %takes bad channel names and finds their index in the EEG struct based on
  %chan_names
    chan_names = [];
    for i = 1:64 %cycle through channels
        var = EEG.chanlocs(1,i).labels; % save channel names in table
        temp = compose(var);
        chan_names = [chan_names; temp];
    end

    if ~isempty(subject.badchans)
      for i = 1:length(subject.badchans )
          ENTER_CHANNEL(i) = find(strcmp(subject.badchans{i},chan_names));
      end
    end

    if isempty(subject.badchans)
        ENTER_CHANNEL = [];
    end

  %% plot worrisome channels
    for ii = ENTER_CHANNEL
        chan = squeeze(FILT_STIMULI(ii,t1:t2,:));
        Chanmax = max(chan);
        Chanmin = min(chan);
        Chanstd = std(chan);
        Chanmean = mean(chan);

        [sorted,index]=sort(Chanmax);
        maxes = [index(end-k+1:end)];
        maxes = fliplr(maxes);
        [sorted,index]=sort(Chanmin);
        mins = [index(1:k)];

        figure
        plot(time(t1:t2),chan,'k')
        hold on
        plot(time(t1:t2),chan(:,[maxes,mins]),'linewidth',2)
        title(chan_names(ii))

        dist_from_mean_maxes = chan(ii,maxes,:) - Chanmean(maxes);
        dist_from_mean_mins = chan(ii,mins,:) - Chanmean(mins);

        data_table = table;
        data_table.chan = repmat(chan_names(ii),k,1);
        data_table.Max = maxes';
        data_table.Min = mins';
        data_table.STDfromMean_Max = (dist_from_mean_maxes./Chanstd(maxes))';
        data_table.STDfromMean_Min = (dist_from_mean_mins./Chanstd(mins))';
        data_table
        
      end

return
