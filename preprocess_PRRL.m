%% Quesstions
%
% For analysis: Should data be re-sampled at 250 Hz? At 500 currently. Yes
% - after cleaning.
% add save points, rename files if doing a second pass, modularize the interp and remove epochs

clearvars
close all

try

%% Create a subject structure that will be saved as a .mat file %%
 subject = {};
 script_home = fileparts(mfilename('fullpath'));
 subject.datapath =  fullfile(script_home, '../eeg_data');
 subject.subject_num = input('Enter Subject #:');
 subject_num = subject.subject_num; % rename just to keep it short throughout the script

%% Pick the type of data cleaning. Epochs should be a cell of strings. See read.me for details. %%
  % triggers is a structure with the name and event triggers for each type of processing
   triggers{1}.name = 'goal_point_triggers'
   triggers{1}.cell_string = {'100', '110', '115'};
   triggers{2}.name = 'bandit_triggers';
   triggers{2}.cell_string = {'120', '121', '122', '123'};
   triggers{3}.name = 'choice_triggers';
   triggers{3}.cell_string = {'124', '125'};
   triggers{4}.name = 'feed_back_triggers';
   triggers{4}.cell_string = {'130', '131', '132', '133'};
   triggers{5}.name = 'all_prrl_eeg_triggers';
   triggers{5}.cell_string = {'100', '110', '115', '119', '120', '121', '122', '123', '124', '125', '129', '130', '131', '132', '133', '139', '199', '200', '210', '211', '212', '213', '220', '221', '230', '239', '244', '245', '249', '253', '254', '255'};
 % get input for which triggers to use
  input_trigger = input('Enter one of the following: \n goal_point_triggers OR\n bandit_triggers OR\n choice_triggers OR\n  feed_back_triggers\n :', 's');
  % get index based on inputted response
  num_names = 1:size(triggers, 2);
  index = arrayfun(@(num) strcmp(triggers{num}.name, input_trigger), num_names);
  % put it in the subject structure
  subject.triggers_name = triggers{num_names(index)}.name;
  subject.triggers = triggers{num_names(index)}.cell_string;

%% add the eeg lab functions
  locpath=('../../Sarah/MATLAB/EEG Data/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
  % addpath('../../Sarah/MATLAB/EEG Data/eeglab13_6_5b/');
  addpath(genpath('../../Sarah/MATLAB/EEG Data/eeglab13_6_5b/'))

%% initialize eeglab
  [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% are we running this from the right place?
  try
    cd(subject.datapath); % will crash if we are not where we should be
  catch
    'Please open matlab from the directory ~/Dropbox (Personal)/Brooke/eeg_processing_scripts!'
  end

%% make new directory
  subject_string = sprintf('PRRL_%s_%d', subject.triggers_name, subject_num);
  status = mkdir(sprintf('PRRL/%s', subject_string));
  folder = sprintf('PRRL/%s', subject_string);
  subject_string = sprintf('PRRL_%d', subject_num); % remove trigger name for all output files

%% File structure. This helps us to determine what files in the dir correspond to the next step %%
  output_files{1}.name = sprintf('%s_interpolated_rereferenced.mat', subject_string);
  output_files{2}.name = sprintf('%s_interpolated_rereferenced_ica.mat', subject_string);
  output_files{3}.name = sprintf('%s_interpolated_rereferenced_ica_filtered.mat', subject_string);
  output_files{1}.name = sprintf('%s_interpolated_rereferenced.mat', subject_string);

%%% Important Note! %%%%
% If you're trying to pick up where someone left off but MATLAB got closed:
% Run up to here and then double-click on file (e.g. ICA file) to re-load data

% Load data - must be done by hand for now.
  eeglab redraw % redraw updates the eeglab gui

step = 0;
steps = [step];
num_files = 1:size(output_files, 2);
already_created_files = dir(folder);
for i = 1:length(already_created_files)
  index = arrayfun(@(num) strcmp(output_files{num}.name, already_created_files(i).name), num_files);
  steps = [steps num_files(index)] ;
end
step = max(steps);

if step == 0;
%% Load the data
  if exist(sprintf('%s/%s.bdf', subject.datapath , subject_string))
    EEG = pop_biosig(sprintf('%s/%s.bdf', subject.datapath , subject_string));
    EEG.setname= subject_string;
    pop_saveset(EEG, subject_string ,folder);
    EEG = eeg_checkset(EEG);
    subject.EEG = EEG % save the EEG structure in the frame
  else
    'EEG DATA DOES NOT EXIST'
    sprintf('Looking here: %s', sprintf('%s/%s.bdf', subject.datapath , subject_string) )
    return
  end

%% Only run the following data cleaning steps if the second pass information is nonexistant. Else skip to interpolation and epoch rejection. Note, this will get confusing once you clean with different trigger sets

  %% add electrode cap location
    EEG = pop_chanedit(EEG,  'lookup', locpath); % if this fails, use loc path as above but from local EEGLab
  %% Extract epochs from the time series.
    % 3rd arg is x seconds before and y seconds after epoch
    EEG = pop_epoch( EEG, subject.triggers, [-1    2]);
    % makes sure things are still coherent
    EEG = eeg_checkset( EEG );

  %% Remove VEOG - gets rid of externals
    VEOG=squeeze(EEG.data(64,:,:));
    HEOG=squeeze(EEG.data(65,:,:));
    EEG.data=EEG.data(1:63,:,:);
    EEG.chanlocs = EEG.chanlocs(1,1:63);
    EEG.nbchan=63;

  %% Re add CPz  - but don't re-ref yet (CPz is for re-referencing stuff, an average of other channels)
    EEG = pop_chanedit(EEG, 'append',63,'changefield',{64 'CPz' 180 0.1266 -32.9279 [-4.0325e-15] 78.3630 -180 67.2080 85});
    EEG = pop_chanedit(EEG,  'lookup', locpath);
    % Remove mean
    EEG = pop_rmbase(EEG,[],[]); % removes baseline activity
  %% Update gui and plot
    eeglab redraw
    pop_eegplot( EEG, 1, 1, 0, [], 'color','on');

  %% Grab channel names
    chan_names = [];
    subject.badchans  = [];

    for i = 1:EEG.nbchan %cycle through channels
        var = EEG.chanlocs(1,i).labels; % save channel names in table
        temp = compose(var);
        chan_names = [chan_names; temp];
    end


  %% --- IN EEG LAB (BLUE POP-UP) ---
  % First pass: (If second pass, skip)
  % Go to "Tools > Reject data epochs > Reject by inspection"
  % Click on messy epochs to highlight them. Scroll through all epochs.
  % When finished, press "update marks."
  % Then press control enter.
    subject.first_rejected_epochs = input('Please enter any rejected epochs as a vector, ie [39 42 54]:');
    subject.rejected_epochs = subject.first_rejected_epochs;

    subject.first_interp = input('Please enter the channels that require interpolation as a cell of strings:');
    subject.interp = subject.first_interp ;


  %%% Reject epochs, Interpolate, and Rereference %%%

  epochNumbers = [1:EEG.trials];
  epochNumbers(subject.rejected_epochs)=[];
  df = 63 - length(subject.interp);
  for i = 1:length(subject.interp)
      ENTER_CHANNEL(i) = find(strcmp(subject.interp(i),chan_names)); %% fix structural issues
  end
  if isempty(subject.interp)
      ENTER_CHANNEL = [];
  end

  %% Interpolate
    if ~isempty(subject.interp)
        for xi=1:size(subject.interp,2)
            for ei=1:63
                if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                    subject.badchans (xi)=ei;
                end
            end
        end
        EEG.data=double(EEG.data);
        EEG = pop_interp(EEG,subject.badchans ,'spherical');
    end


  %% Reject bad epochs
    if ~isempty(subject.rejected_epochs)
        binarized=zeros(1,EEG.trials);
        binarized(subject.rejected_epochs)=1;
        EEG = pop_rejepoch(EEG,binarized,0);
    end


  %% Rereference
    EEG = pop_reref(EEG, [],'refloc',struct('labels',{'CPz'},'type',{[]},'theta',{180},'radius',{0.12662},'X',{-32.9279},'Y',{-4.0325e-15},'Z',{78.363},'sph_theta',{-180},'sph_phi',{67.208},'sph_radius',{85},'urchan',{64},'ref',{''}));

  subject.EEG = EEG % save the EEG structure in the frame
  save(sprintf('%s/%s_interpolated_rereferenced.mat', folder, subject_string),'-struct', 'subject');
  step = 1;
end


if step == 1
%% ICA
  %% Load saved data
  input('Ready to start ICA?');
  subject = load(sprintf('%s/%s_interpolated_rereferenced.mat', folder, subject_string));
  EEG = subject.EEG;

  if isempty(subject.interp)
      ENTER_CHANNEL = [];
  end

  channels = [1:63];
  idx = ismember(channels,ENTER_CHANNEL); %%% can channels have the same name?
  channels(idx) = [];

  %Remove channels from list
  tic
  EEG = pop_runica(EEG, 'chanind',channels);
  toc

  % EEG = pop_runica(EEG,'icatype','runica');
  disp('ICA finished. Saving now.')
  subject.EEG = EEG % save the EEG structure in the frame
  save(sprintf('%s/%s_interpolated_rereferenced_ica.mat', folder, subject_string),'-struct', 'subject');
  disp('Saving finished.')
  step = 2  ;
end
%% Data reset, display the ICA components selected

clearvars -except subject_string folder step

if step == 2
  subject = load(sprintf('%s/%s_interpolated_rereferenced_ica.mat', folder, subject_string));
  EEG = subject.EEG;

  %% Display ICA components %%
  pop_selectcomps(EEG, [1:10] ); %Displays topoplots
  pop_eegplot( EEG, 0, 1, 1); %Displays component scroll

  %% Reject marked ICA components according to subno

  subject.ICAreject = input('Please enter any bad component as a vector, ie [1 2]:');

  EEG = pop_subcomp(EEG, subject.ICAreject, 0); %Removes marked component

  % Re-interpolate channels after removing bad components
  if ~isempty(subject.interp)
      for xi=1:size(subject.interp,2)
          for ei=1:63
              if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                  subject.badchans (xi)=ei;
              end
          end
      end
      EEG.data=double(EEG.data);
      EEG = pop_interp(EEG,subject.badchans ,'spherical');
  end

  EEG = eeg_checkset( EEG );
  eeglab redraw
  ALLEEG=[];
  pop_eegplot( EEG, 1, 1, 1);

  %% Filter
    input('Ready to start filtering?');
    dims = size(EEG.data);
    FILT=eegfilt(EEG.data,EEG.srate,[],15); % low pass filter below 15.
    FILT=eegfilt(FILT,EEG.srate,.5,[]); % high pass filter above .5.   eegfilt
    EEG.data = reshape(FILT,dims(1),dims(2),dims(3));
    pop_eegplot( EEG, 1, 1, 1);

    subject.EEG = EEG;
    save(sprintf('%s/%s_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
    disp('Filtering and saving complete.')
    step = 3;
end

%%% Second Pass cleaning %%%

clearvars -except subject_string folder step;

if step == 3
  subject = load(sprintf('%s/%s_interpolated_rereferenced_ica_filtered.mat', folder, subject_string));
  EEG = subject.EEG;

  eeglab redraw
  pop_eegplot( EEG, 1, 1, 0,[],'color','on');

  subject.second_rejected_epochs = input('Please enter any NEW rejected epochs as a vector, ie [39 42 54]:');
  subject.second_interp = input('Please enter any NEW channels that require interpolation as a cell of strings:');

  input('Press Enter to Continue:');

  time = EEG.times/1000; %% Q: Why divide by 1000?

  t1 = find(time == -.5);
  t2 = find(time <1.5);t2 = t2(end);

  FILT_STIMULI = EEG.data;

  for electrode = 1:64
      figure(3)
      clf
      subplot(2,1,1)
      plot(time(t1:t2),squeeze(FILT_STIMULI(electrode,t1:t2,:)),'k')
      subplot(2,1,2)
      plot(time(t1:t2),squeeze(mean(FILT_STIMULI(electrode,t1:t2,:),3)),'r','linewidth',2)

      title(EEG.chanlocs(electrode).labels)
      pause(6)
  end

  %% Identify channels with abberant ERPs for epoch ejection
  clear ENTER_CHANNEL subject.badchans  var data_table
  chan_names = [];

  for i = 1:64 %cycle through channels
      var = EEG.chanlocs(1,i).labels; % save channel names in table
      temp = compose(var);
      chan_names = [chan_names; temp];
  end

  % Enter channels (can enter multiple, will evaluate one at a time)
  % Make sure to scroll all the way up in command window to see ALL bad eps

  subject.badchans = input('Please enter the bad channels as a cell of strings:');


  %takes bad channel names and finds their index in the EEG struct based on
  %chan_names
  if ~isempty(subject.badchans)
    for i = 1:length(subject.badchans )
        ENTER_CHANNEL(i) = find(strcmp(subject.badchans (i),chan_names));
    end
  end

  if isempty(subject.badchans)
      ENTER_CHANNEL = [];
  end

  k = 1; %k is the number of max and min pairs
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
      %[maxFp1 badep1] = max(Chanmax)
      %[minFp1 badep2] = min(Chanmin)

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

      disp(['Epochs to take out - enter these below.  Might only want max or min. Max is first #, min is 2nd #.']);
      %epochNumbers([badep1,badep2])
      data_table
  end

  input('Press Enter to Continue:');

  eeglab redraw
  pop_eegplot( EEG, 1, 1, 0,[],'color','on');


  subject.second_rejected_epochs = input('Please enter any rejected epochs as a vector, ie [39 42 54]:');
  subject.second_interp = input('Please enter the channels that require interpolation as a cell of strings:');

  disp(['do not double select: ' num2str(sort(subject.second_rejected_epochs)) ]) % Q: should this be the first epochs
end
%% Second pass?

% to_reject_second = [to_reject_second, find(EEG.reject.rejmanual)]; %add to the bad epochs from vis inspection should I use this?

spit = input('Are you going to do a second pass? y or n ', 's');
if strcmp(spit,'y')
    subject.second_pass = 'TRUE';
    save(sprintf('%s/%s_first_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
else
  step = 5;
end

if step == 4

  %% reload dataset
  EEG = pop_biosig(sprintf('%s/%s.bdf', subject.datapath , subject_string));
  EEG.setname= subject_string;
  pop_saveset(EEG, subject_string ,folder);
  EEG = eeg_checkset(EEG);
  subject.EEG = EEG % save the EEG structure in the frame

  %% add electrode cap location
    EEG = pop_chanedit(EEG,  'lookup', locpath); % if this fails, use loc path as above but from local EEGLab
  %% Extract epochs from the time series.
    % 3rd arg is x seconds before and y seconds after epoch
    EEG = pop_epoch( EEG, subject.triggers, [-1    2]);
    % makes sure things are still coherent
    EEG = eeg_checkset( EEG );

  %% Remove VEOG - gets rid of externals
    VEOG=squeeze(EEG.data(64,:,:));
    HEOG=squeeze(EEG.data(65,:,:));
    EEG.data=EEG.data(1:63,:,:);
    EEG.chanlocs = EEG.chanlocs(1,1:63);
    EEG.nbchan=63;

  %% Re add CPz  - but don't re-ref yet (CPz is for re-referencing stuff, an average of other channels)
    EEG = pop_chanedit(EEG, 'append',63,'changefield',{64 'CPz' 180 0.1266 -32.9279 [-4.0325e-15] 78.3630 -180 67.2080 85});
    EEG = pop_chanedit(EEG,  'lookup', locpath);
    % Remove mean
    EEG = pop_rmbase(EEG,[],[]); % removes baseline activity
  %% Update gui and plot
    eeglab redraw
    pop_eegplot( EEG, 1, 1, 0, [], 'color','on');

  %% Grab channel names
    chan_names = [];
    subject.badchans  = [];

    for i = 1:EEG.nbchan %cycle through channels
        var = EEG.chanlocs(1,i).labels; % save channel names in table
        temp = compose(var);
        chan_names = [chan_names; temp];
    end

    epochNumbers = [1:EEG.trials];
    epochNumbers(subject.rejected_epochs)=[];
    df = 63 - length(subject.interp);
    for i = 1:length(subject.interp)
        ENTER_CHANNEL(i) = find(strcmp(subject.interp(i),chan_names)); %% fix structural issues
    end
    if isempty(subject.interp)
        ENTER_CHANNEL = [];
    end

    %% Interpolate
      if ~isempty(subject.interp)
          for xi=1:size(subject.interp,2)
              for ei=1:63
                  if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                      subject.badchans (xi)=ei;
                  end
              end
          end
          EEG.data=double(EEG.data);
          EEG = pop_interp(EEG,subject.badchans ,'spherical');
      end


    %% Reject bad epochs
      if ~isempty(subject.rejected_epochs)
          binarized=zeros(1,EEG.trials);
          binarized(subject.rejected_epochs)=1;
          EEG = pop_rejepoch(EEG,binarized,0);
      end


    %% Rereference
      EEG = pop_reref(EEG, [],'refloc',struct('labels',{'CPz'},'type',{[]},'theta',{180},'radius',{0.12662},'X',{-32.9279},'Y',{-4.0325e-15},'Z',{78.363},'sph_theta',{-180},'sph_phi',{67.208},'sph_radius',{85},'urchan',{64},'ref',{''}));

    subject.EEG = EEG % save the EEG structure in the frame

    %% ICA
      %% Load saved data
      input('Ready to start ICA?');

      if isempty(subject.interp)
          ENTER_CHANNEL = [];
      end

      channels = [1:63];
      idx = ismember(channels,ENTER_CHANNEL); %%% can channels have the same name?
      channels(idx) = [];

      %Remove channels from list
      tic
      EEG = pop_runica(EEG, 'chanind',channels);
      toc

      % EEG = pop_runica(EEG,'icatype','runica');
      disp('ICA finished. Saving now.')
      subject.EEG = EEG % save the EEG structure in the frame
      save(sprintf('%s/%s_second_interpolated_rereferenced_ica.mat', folder, subject_string),'-struct', 'subject');
      disp('Saving finished.')

    %% Data reset, display the ICA components selected

    clearvars -except subject_string folder step

    subject = load(sprintf('%s/%s_second_interpolated_rereferenced_ica.mat', folder, subject_string));
    EEG = subject.EEG;

    %% Display ICA components %%
    pop_selectcomps(EEG, [1:10] ); %Displays topoplots
    pop_eegplot( EEG, 0, 1, 1); %Displays component scroll

    %% Reject marked ICA components according to subno

    subject.ICAreject2 = input('Please enter any bad component as a vector, ie [1 2]:');

    EEG = pop_subcomp(EEG, subject.ICAreject2, 0); %Removes marked component

    % Re-interpolate channels after removing bad components
    if ~isempty(subject.interp)
        for xi=1:size(subject.interp,2)
            for ei=1:63
                if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                    subject.badchans (xi)=ei;
                end
            end
        end
        EEG.data=double(EEG.data);
        EEG = pop_interp(EEG,subject.badchans ,'spherical');
    end

    EEG = eeg_checkset( EEG );
    eeglab redraw
    ALLEEG=[];
    pop_eegplot( EEG, 1, 1, 1);

    %% Filter
      input('Ready to start filtering?');
      dims = size(EEG.data);
      FILT=eegfilt(EEG.data,EEG.srate,[],15); % low pass filter below 15.
      FILT=eegfilt(FILT,EEG.srate,.5,[]); % high pass filter above .5.   eegfilt
      EEG.data = reshape(FILT,dims(1),dims(2),dims(3));
      pop_eegplot( EEG, 1, 1, 1);

      subject.EEG = EEG;
      save(sprintf('%s/%s_second_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
      disp('Filtering and saving complete.')
      step = 5;
end

if step == 6
  subject = load(sprintf('%s/%s_second_interpolated_rereferenced_ica_filtered.mat', folder, subject_string));
  EEG = subject.EEG;

  time = EEG.times/1000; %% Q: Why divide by 1000?

  t1 = find(time == -.5);
  t2 = find(time <1.5);t2 = t2(end);

  FILT_STIMULI = EEG.data;

  for electrode = 1:64
      figure(3)
      clf
      subplot(2,1,1)
      plot(time(t1:t2),squeeze(FILT_STIMULI(electrode,t1:t2,:)),'k')
      subplot(2,1,2)
      plot(time(t1:t2),squeeze(mean(FILT_STIMULI(electrode,t1:t2,:),3)),'r','linewidth',2)

      title(EEG.chanlocs(electrode).labels)
      pause(6)
  end

  %% Identify channels with abberant ERPs for epoch ejection
  clear ENTER_CHANNEL subject.badchans  var data_table
  chan_names = [];

  for i = 1:64 %cycle through channels
      var = EEG.chanlocs(1,i).labels; % save channel names in table
      temp = compose(var);
      chan_names = [chan_names; temp];
  end

  % Enter channels (can enter multiple, will evaluate one at a time)
  % Make sure to scroll all the way up in command window to see ALL bad eps

  subject.badchans = input('Please enter the bad channels as a cell of strings:');


  %takes bad channel names and finds their index in the EEG struct based on
  %chan_names
  if ~isempty(subject.badchans)
    for i = 1:length(subject.badchans )
        ENTER_CHANNEL(i) = find(strcmp(subject.badchans (i),chan_names));
    end
  end

  if isempty(subject.badchans)
      ENTER_CHANNEL = [];
  end

  k = 1; %k is the number of max and min pairs
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
      %[maxFp1 badep1] = max(Chanmax)
      %[minFp1 badep2] = min(Chanmin)

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

      disp(['Epochs to take out - enter these below.  Might only want max or min. Max is first #, min is 2nd #.']);
      %epochNumbers([badep1,badep2])
      data_table
  end

end

%% save - DON'T FORGET THIS PART
subject.EEG = EEG;
save(sprintf('%s/%s_interpolated_rereferenced_ica_filtered_final.mat', folder, subject_string),'-struct', 'subject');

disp('Cleaning stats:')
%% I took out the if statement on second pass, because in both cases it is length(reject_epochs) + lenght(second_rejected_epochs), just if second_pass == TRUE, the second_rejected is technically a third pass
subject.total_num_rej_epochs = length(subject.rejected_epochs)+length(subject.second_rejected_epochs);
subject.total_num_interps = (length(subject.interp)/epochNumbers(end))*100;
disp(['Number of epochs rejected, note in lab notebook: ' num2str(subject.total_num_rej_epochs )]) % in this case
disp(['Percentage of epochs rejected, note in lab notebook: ' num2str(subject.total_num_interps ) '%'])


catch err
  subject.EEG = EEG;
  save(sprintf('%s/%s_interrupted_during_step%d.mat', folder, subject_string, step),'-struct', 'subject');
  err
end
