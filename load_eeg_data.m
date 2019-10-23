function subject = load_eeg_data(datapath, folder, subject_string, locpath, subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  this function loads in the .bdf data, adds channel location data,
%%  removes externals, and removes baseline activity
%%
%%  INPUTS:
%%    datapath: path to the .bdf file, stored in subject.datapath
%%    folder: output folder
%%    subject_string: string name of the bdf file, with the trigger set name
%%    locpath: file with the channel location data
%%    subject: the subject structure, with EEG and PRRL specific metadata
%%
%% OUTPUTS:
%%    subject: modified eeg data within subject structure
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Load data
    if exist(sprintf('%s/%s.bdf', datapath , subject_string))
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

  %% add electrode cap location
    EEG = pop_chanedit(EEG,  'lookup', locpath);
  %% Extract epochs from the time series.
    EEG = pop_epoch( EEG, subject.triggers, [-1 2]); % 3rd arg is x secs before and y secs after epoch
  %% makes sure things are still coherent
    EEG = eeg_checkset( EEG )

  %% Remove externals
    VEOG=squeeze(EEG.data(65,:,:));
    HEOG=squeeze(EEG.data(66,:,:));
    EEG.data=EEG.data(1:64,:,:);
    EEG.chanlocs = EEG.chanlocs(1,1:64);
    EEG.nbchan=64;

  %% Remove mean
    EEG = pop_rmbase(EEG,[],[]); % removes baseline activity

  %% save to structure
    subject.EEG = EEG;

return
