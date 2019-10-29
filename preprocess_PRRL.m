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
 cd(script_home);
 subject.datapath = input('What is the path to the data?\n ../eeg_data\n ../../Brooke/eeg_data\n other\n:', 's');
 subject.subject_num = input('Enter Subject #:');
 subject_num = subject.subject_num; % rename just to keep it short throughout the script
 subject.first_rejected_epochs = [];
 subject.second_rejected_epochs = [];
 subject.third_rejected_epochs = [];
 subject.second_pass = 'FALSE';
 subject.third_pass = 'FALSE';

%% Pick the type of data cleaning. Epochs should be a cell of strings. See read.me for details. %%
  % triggers is a structure with the name and event triggers for each type of processing
   triggers{1}.name = 'goal_point_triggers';
   triggers{1}.cell_string = {'100', '110', '115'};
   triggers{2}.name = 'bandit_triggers';
   triggers{2}.cell_string = {'120', '121', '122', '123'};
   triggers{3}.name = 'choice_triggers';
   triggers{3}.cell_string = {'124', '125'};
   triggers{4}.name = 'feed_back_triggers';
   triggers{4}.cell_string = {'130', '131', '132', '133'};
   triggers{5}.name = 'nback_ons_triggers';
   triggers{5}.cell_string = {'20', '25', '30', '35', '40', '45'};
   triggers{6}.name = 'all_prrl_eeg_triggers';
   triggers{6}.cell_string = {'100', '110', '115', '119', '120', '121', '122', '123', '124', '125', '129', '130', '131', '132', '133', '139', '199', '200', '210', '211', '212', '213', '220', '221', '230', '239', '244', '245', '249', '253', '254', '255'};
   triggers{7}.name = 'nback_off_triggers';
   triggers{7}.cell_string = {'59'};
 % get input for which triggers to use
 input_trigger = input('Enter one of the following: \n goal_point_triggers OR\n bandit_triggers OR\n choice_triggers OR\n feed_back_triggers OR\n nback_ons_triggers\n:', 's');
  % get index based on inputted response
  num_names = 1:size(triggers, 2);
  index = arrayfun(@(num) strcmp(triggers{num}.name, input_trigger), num_names);
  % put it in the subject structure
  subject.triggers_name = triggers{num_names(index)}.name;
  subject.triggers = triggers{num_names(index)}.cell_string;

%% add the eeg lab functions
  locpath=sprintf('%s/EEGLAB/eeglab13_6_5b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp', script_home);
  addpath('./EEGLAB/eeglab13_6_5b/');	  addpath(sprintf('%s/EEGLAB/eeglab13_6_5b/', script_home));

%% initialize eeglab
  [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% are we running this from the right place?
  try
    cd(subject.datapath); % will crash if we are not where we should be
  catch
    'Please open matlab from the directory ~/Dropbox (Personal)/Brooke/eeg_processing_scripts!'
  end

%% make new directory
  subject_string = sprintf('Oddball_%s_%d', subject.triggers_name, subject_num);
  status = mkdir(sprintf('preprocessed_data/%s', subject_string));
  folder = sprintf('preprocessed_data/%s', subject_string);
  subject_string = sprintf('Oddball_%d', subject_num); % remove trigger name for all output files

%% File structure. This helps us to determine what files in the dir correspond to the next step %%
  output_files{1}.name = sprintf('%s_interpolated_rereferenced.mat', subject_string);
  output_files{2}.name = sprintf('%s_interpolated_rereferenced_ica.mat', subject_string);
  output_files{3}.name = sprintf('%s_interpolated_rereferenced_ica_filtered.mat', subject_string);
  output_files{4}.name = sprintf('%s_first_full_interpolated_rereferenced_ica_filtered.mat', subject_string);
  output_files{5}.name = sprintf('%s_second_interpolated_rereferenced_ica.mat', subject_string);
  output_files{6}.name = sprintf('%s_second_interpolated_rereferenced_ica_filtered.mat', subject_string);
  output_files{7}.name = sprintf('%s_second_full_interpolated_rereferenced_ica_filtered.mat', subject_string);
  output_files{8}.name = sprintf('%s_third_interpolated_rereferenced_ica.mat', subject_string);
  output_files{9}.name = sprintf('%s_third_interpolated_rereferenced_ica_filtered.mat', subject_string);

%% See if data files exist for this subject/trigger set, and hop to correct step
  step = 0;
  steps = [step];
  num_files = 1:size(output_files, 2);
  already_created_files = dir(folder);
  for i = 1:length(already_created_files)
    index = arrayfun(@(num) strcmp(output_files{num}.name, already_created_files(i).name), num_files);
    steps = [steps num_files(index)] ;
  end
  step = max(steps);

%%% First pass, load data, identify bad channels/epochs %%%
if step == 0;

%% Load the data
  subject = load_eeg_data(subject.datapath, folder, subject_string, locpath, subject);
  subject.num_epochs = subject.EEG.trials;
  EEG = subject.EEG;
  %% Update gui and plot
    % eeglab redraw

    pop_eegplot(EEG, 1, 1, 0, [],'color','on');


  %% Identify bad epochs and channels
    input('Are you done rejecting epochs? Remember to press reject on the plot');
    idx = size(ALLEEG, 2);
    subject.first_rejected_epochs = ALLEEG(idx).reject.rejmanual ;
    subject.first_interp = input('Please enter the channels that require interpolation as a cell of strings:');
    subject.interp = subject.first_interp ;

  %% Reject epochs, Interpolate, and Rereference
    subject = reject_interp_reref(subject);
    EEG = subject.EEG;
  %% Save
    subject.EEG = EEG;
    save(sprintf('%s/%s_interpolated_rereferenced.mat', folder, subject_string),'-struct', 'subject');
    step = 1;

end

%%% First ICA %%%
if step == 1

  %% Load saved data
  load_string = input('Do you need to reload the structure? If so, press y', 's');
  if strcmp(load_string,'y')
    subject = load(sprintf('%s/%s_interpolated_rereferenced.mat', folder, subject_string));
    EEG = subject.EEG;
  end

  input('Ready to start ICA?');

  %% Grab channel names and prep for ICA
    chan_names = [];
    for i = 1:EEG.nbchan %cycle through channels
        var = EEG.chanlocs(1,i).labels; % save channel names in table
        temp = compose(var);
        chan_names = [chan_names; temp];
    end

    for i = 1:length(subject.interp)
        ENTER_CHANNEL(i) = find(strcmp(subject.interp(i),chan_names)); %% fix structural issues
    end
    if isempty(subject.interp)
        ENTER_CHANNEL = [];
    end

    channels = [1:64];
    idx = ismember(channels,ENTER_CHANNEL);
    channels(idx) = [];

  %% run the ica
  tic
  EEG = pop_runica(EEG, 'chanind',channels);
  toc
  disp('ICA finished. Saving now.')

  %% save the data
  subject.EEG = EEG;
  save(sprintf('%s/%s_interpolated_rereferenced_ica.mat', folder, subject_string),'-struct', 'subject');
  disp('Saving finished.')
  step = 2  ;

end

%%% First Pass, Identify/Remove bad ICA components, Filter %%%

if step == 2
  load_string = input('Do you need to reload the structure? If so, press y', 's');
  if strcmp(load_string,'y')
    subject = load(sprintf('%s/%s_interpolated_rereferenced_ica.mat', folder, subject_string));
    EEG = subject.EEG;
  end

  %% Display ICA components %%
  pop_selectcomps(EEG, [1:10] ); %Displays topoplots
  pop_eegplot( EEG, 0, 1, 1); %Displays component scroll

  %% Reject marked ICA components according to subno

  subject.ICAreject = input('Have you finished rejecting components?:');
  subject.ICAreject = find(EEG.reject.gcompreject == 1);

  EEG = pop_subcomp(EEG, subject.ICAreject, 0); %Removes marked component

  % Re-interpolate channels after removing bad components
  if ~isempty(subject.interp)
      for xi=1:size(subject.interp,2)
          for ei=1:64
              if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                  subject.badchans (xi)=ei;
              end
          end
      end
      EEG.data=double(EEG.data);
      EEG = pop_interp(EEG,subject.badchans ,'spherical');
  end

  EEG = eeg_checkset( EEG )
  % eeglab redraw
  pop_eegplot( EEG, 1, 1, 1);

  %% Filter
    input('Ready to start filtering?');
    dims = size(EEG.data);
    FILT=eegfilt(EEG.data,EEG.srate,[],20); % low pass filter below 20.
    FILT=eegfilt(FILT,EEG.srate,.5,[]); % high pass filter above .5.
    EEG.data = reshape(FILT,dims(1),dims(2),dims(3));
    pop_eegplot( EEG, 1, 1, 1);

 %% Save data
    subject.EEG = EEG;
    save(sprintf('%s/%s_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
    disp('Filtering and saving complete.')
    step = 3;
end

%%% Second pass, Identify Bad Channels/Epochs %%%

if step == 3
  %% load data
    load_string = input('Do you need to reload the structure? If so, press y', 's');
    if strcmp(load_string,'y')
      subject = load(sprintf('%s/%s_interpolated_rereferenced_ica_filtered.mat', folder, subject_string));
      EEG = subject.EEG;
    end

  %% Plot channels to find abberant ERPs for epoch ejection
    plot_bad_epochs(subject, 3)

  %% Mark epochs for rejection
    pop_eegplot(EEG, 1, 1, 0, [],'color','on');


    input('Look at the datatable above and determine if there are any more epochs ot delete. If there are, go to the epoch number on the plot and select the epochs for deletion. Press enter to continue')
    idx = size(ALLEEG, 2);
    subject.second_rejected_epochs = ALLEEG(idx).reject.rejmanual ;
    subject.second_interp = input('Please enter any NEW channels that require interpolation as a cell of strings:');
    subject.interp = [subject.first_interp subject.second_interp];

  %% Will you do a second ICA and filter? Probably always yes %%
    spit = input('Are you going to do a second ica & filter ? y or n ', 's');
    if strcmp(spit,'y')
        subject.second_pass = 'TRUE';
        subject.EEG = EEG;
        save(sprintf('%s/%s_first_full_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
        step = 4;
    else
      step = 9;
    end

end

%%% Second ICA %%%
if step == 4
  %% clear vars
    clear ENTER_CHANNEL

  %% reload the subject structure with accurate metadata
    load_string = input('Do you need to reload the structure? If so, press y', 's');
    if strcmp(load_string,'y')
      subject = load(sprintf('%s/%s_first_full_interpolated_rereferenced_ica_filtered.mat', folder, subject_string));
      EEG = subject.EEG;
    end

  %% load raw eeg data
    new_eeg_data = load_eeg_data(subject.datapath, folder, subject_string, locpath, subject)
    subject.EEG = new_eeg_data.EEG;
    EEG = subject.EEG;

  %% Reject epochs, Interpolate, and Rereference
    subject = reject_interp_reref(subject);
    EEG = subject.EEG;

  %% ICA
    input('Ready to start ICA?');

  %% Grab channel names and prep for ICA
    chan_names = [];
    for i = 1:EEG.nbchan %cycle through channels
        var = EEG.chanlocs(1,i).labels; % save channel names in table
        temp = compose(var);
        chan_names = [chan_names; temp];
    end

    for i = 1:length(subject.interp)
        ENTER_CHANNEL(i) = find(strcmp(subject.interp(i),chan_names)); %% fix structural issues
    end
    if isempty(subject.interp)
        ENTER_CHANNEL = [];
    end

    channels = [1:64];
    idx = ismember(channels,ENTER_CHANNEL);
    channels(idx) = [];

    % Remove channels from list
    tic
    EEG = pop_runica(EEG, 'chanind',channels);
    toc
    disp('ICA finished. Saving now.')

  %% save the data
    subject.EEG = EEG;
    save(sprintf('%s/%s_second_interpolated_rereferenced_ica.mat', folder, subject_string),'-struct', 'subject');
    disp('Saving finished.')
    step = 5;

end

%%% Second Pass, Identify/Remove bad ICA components, Filter %%%
if step == 5
    %% Reload the structure
      load_string = input('Do you need to reload the structure? If so, press y', 's');
      if strcmp(load_string,'y')
        subject = load(sprintf('%s/%s_second_interpolated_rereferenced_ica.mat', folder, subject_string));
        EEG = subject.EEG;
      end

    %% Display ICA components %%
      pop_selectcomps(EEG, [1:10] ); %Displays topoplots
      pop_eegplot( EEG, 0, 1, 1); %Displays component scroll

    %% Reject Bad ICA components
      subject.ICAreject2 = input('Have you finished rejecting components?:');
      subject.ICAreject2 = find(EEG.reject.gcompreject == 1);
      EEG = pop_subcomp(EEG, subject.ICAreject2, 0); %Removes marked component

      % Re-interpolate channels after removing bad components
      if ~isempty(subject.interp)
          for xi=1:size(subject.interp,2)
              for ei=1:64
                  if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                      subject.badchans (xi)=ei;
                  end
              end
          end
          EEG.data=double(EEG.data);
          EEG = pop_interp(EEG,subject.badchans ,'spherical');
      end

      EEG = eeg_checkset( EEG )
      % eeglab redraw
      pop_eegplot( EEG, 1, 1, 1);

    %% Filter
      input('Ready to start filtering?');
      dims = size(EEG.data);
      FILT=eegfilt(EEG.data,EEG.srate,[],20); % low pass filter below 20.
      FILT=eegfilt(FILT,EEG.srate,.5,[]); % high pass filter above .5.
      EEG.data = reshape(FILT,dims(1),dims(2),dims(3));
      pop_eegplot( EEG, 1, 1, 1);

    %% Save dataset
      subject.EEG = EEG;
      save(sprintf('%s/%s_second_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
      disp('Filtering and saving complete.')
      step = 6;
end

%%% Third Pass, Identify Bad Channels/Epochs %%%

if step == 6

  %% Reload the structure
    load_string = input('Do you need to reload the structure? If so, press y', 's');
    if strcmp(load_string,'y')
      subject = load(sprintf('%s/%s_second_interpolated_rereferenced_ica_filtered.mat', folder, subject_string));
      EEG = subject.EEG;
    end

  %% Plot channels to find abberant ERPs for epoch ejection
    plot_bad_epochs(subject, 3)

  %% Mark epochs for rejection
    % eeglab redraw

    pop_eegplot(EEG, 1, 1, 0, [],'color','on');


    input('Look at the datatable above and determine if there are any more epochs ot delete. If there are, go to the epoch number on the plot and select the epochs for deletion. Press enter to continue')
    idx = size(ALLEEG, 2);
    subject.third_rejected_epochs = ALLEEG(idx).reject.rejmanual ;
    subject.third_interp = input('Please enter any NEW channels that require interpolation as a cell of strings:');
    subject.interp = [subject.first_interp subject.second_interp subject.third_interp];

    %% Third ICA & Filter? %%
      spit = input('Are you going to do a third ica & filter? y or n ', 's');
      if strcmp(spit,'y')
          subject.third_pass = 'TRUE';
          subject.EEG = EEG;
          save(sprintf('%s/%s_second_full_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
          step = 7;
      else
        step = 9;
      end

end

%%% Third ICA %%%
if step == 7
  %% clear vars
    clear ENTER_CHANNEL

  %% reload the subject structure with accurate metadata
    load_string = input('Do you need to reload the structure? If so, press y', 's');
    if strcmp(load_string,'y')
      subject = load(sprintf('%s/%s_second_full_interpolated_rereferenced_ica_filtered.mat', folder, subject_string));
      EEG = subject.EEG;
    end

  %% load raw eeg data
    new_eeg_data = load_eeg_data(subject.datapath, folder, subject_string, locpath, subject)
    subject.EEG = new_eeg_data.EEG;
    EEG = subject.EEG;

  %% Reject epochs, Interpolate, and Rereference
    subject = reject_interp_reref(subject);
    EEG = subject.EEG;

  %% ICA
    input('Ready to start ICA?');

  %% Grab channel names and prep for ICA
    chan_names = [];
    for i = 1:EEG.nbchan %cycle through channels
        var = EEG.chanlocs(1,i).labels; % save channel names in table
        temp = compose(var);
        chan_names = [chan_names; temp];
    end

    for i = 1:length(subject.interp)
        ENTER_CHANNEL(i) = find(strcmp(subject.interp(i),chan_names)); %% fix structural issues
    end
    if isempty(subject.interp)
        ENTER_CHANNEL = [];
    end

    channels = [1:64];
    idx = ismember(channels,ENTER_CHANNEL);
    channels(idx) = [];

    % Remove channels from list
    tic
    EEG = pop_runica(EEG, 'chanind',channels);
    toc
    disp('ICA finished. Saving now.')

  %% save the data
    subject.EEG = EEG;
    save(sprintf('%s/%s_third_interpolated_rereferenced_ica.mat', folder, subject_string),'-struct', 'subject');
    disp('Saving finished.')
    step = 8;

end

%%% Third Pass, Identify/Remove bad ICA components, Filter %%%

if step == 8
    %% Reload the structure
      load_string = input('Do you need to reload the structure? If so, press y', 's');
      if strcmp(load_string,'y')
        subject = load(sprintf('%s/%s_third_interpolated_rereferenced_ica.mat', folder, subject_string));
        EEG = subject.EEG;
      end

    %% Display ICA components %%
      pop_selectcomps(EEG, [1:10] ); %Displays topoplots
      pop_eegplot( EEG, 0, 1, 1); %Displays component scroll

    %% Reject Bad ICA components
      subject.ICAreject3 = input('Have you finished rejecting components?:');
      subject.ICAreject3 = find(EEG.reject.gcompreject == 1);
      EEG = pop_subcomp(EEG, subject.ICAreject3, 0); %Removes marked component

      % Re-interpolate channels after removing bad components
      if ~isempty(subject.interp)
          for xi=1:size(subject.interp,2)
              for ei=1:64
                  if strmatch(EEG.chanlocs(ei).labels,subject.interp{xi});
                      subject.badchans (xi)=ei;
                  end
              end
          end
          EEG.data=double(EEG.data);
          EEG = pop_interp(EEG,subject.badchans ,'spherical');
      end

      EEG = eeg_checkset( EEG )
      pop_eegplot(EEG, 1, 1, 0, [],'color','on');


    %% Filter
      input('Ready to start filtering?');
      dims = size(EEG.data);
      FILT=eegfilt(EEG.data,EEG.srate,[],20); % low pass filter below 20.
      FILT=eegfilt(FILT,EEG.srate,.5,[]); % high pass filter above .5.
      EEG.data = reshape(FILT,dims(1),dims(2),dims(3));
      pop_eegplot( EEG, 1, 1, 1, [],'color','on');

    %% Save data
      subject.EEG = EEG;
      save(sprintf('%s/%s_third_interpolated_rereferenced_ica_filtered.mat', folder, subject_string),'-struct', 'subject');
      disp('Filtering and saving complete.')
      step = 9;

end

if step == 9
  %% extract trial type info and save in human usable format %%
  event_list = NaN(length(subject.EEG.epoch), 1);
  for i = 1:length(subject.EEG.epoch)
    tmp_cell = subject.EEG.epoch(i).eventtype;
    trigNum = str2num(str2mat(subject.triggers));
    cellIdx = cellfun(@(cell) sum(cell == trigNum), tmp_cell);
    event_list(i) = tmp_cell{find(cellIdx == 1)} ;
  end
  subject.final_event_list = event_list ;

  %% any additional notes?
  subject.notes = input('Please enter any additional notes you may for this subject:', 's');

  %% save - DON'T FORGET THIS PART
  subject.EEG = EEG;
  save(sprintf('%s/%s_interpolated_rereferenced_ica_filtered_FINAL.mat', folder, subject_string),'-struct', 'subject');
  %% clean up
  spit = input('Are you ready to delete intermediary files? If so, press y', 's');
  if strcmp(spit,'y')
    cd(folder)
    interrupted_files = dir(['*interrupted*']);
    all_files = dir();
    not_match = cellfun(@isempty, regexp({all_files.name},'filtered'));
    extra_files = all_files(not_match);
    if ~isempty(interrupted_files)
        delete(interrupted_files.name);
    end
    delete(extra_files.name);
  end

  disp('Cleaning stats:')
  if strcmp(subject.third_pass, 'TRUE')
    subject.total_num_rej_epochs = subject.num_epochs - sum(subject.third_rejected_epochs == 0);
  elseif strcmp(subject.second_pass, 'TRUE')
      subject.total_num_rej_epochs = subject.num_epochs - sum(subject.second_rejected_epochs == 0);
  else
    subject.total_num_rej_epochs = subject.num_epochs - sum(subject.first_rejected_epochs == 0);
  end
  subject.percentage_rejected = (subject.total_num_rej_epochs/subject.num_epochs) * 100;
  subject.total_num_interps = length(subject.interp);
  disp(['Number of epochs rejected, note in lab notebook: ' num2str(subject.total_num_rej_epochs )])
  disp(['Percentage of epochs rejected, note in lab notebook: ' num2str(subject.percentage_rejected ) '%'])
end


catch err
  subject.EEG = EEG;
  save(sprintf('%s/%s_interrupted_during_step%d.mat', folder, subject_string, step),'-struct', 'subject');
  err
  'On lines'
  err.stack.line
end
