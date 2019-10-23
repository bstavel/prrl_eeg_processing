function subject = reject_interp_reref(subject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  this function loads rejects the given epochs, interpolates bad channels
%%  and the re - rereferences the data
%%
%%  INPUTS:
%%    subject: the subject structure, with EEG and PRRL specific metadata
%%
%% OUTPUTS:
%%    subject: cleaned eeg data within subject structure
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Pull EEG out of subject structure
    EEG = subject.EEG;

  %% Interpolate
    subject.badchans  = [];
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


  %% Reject bad epochs
    % round 1
    if ~isempty(subject.first_rejected_epochs)
        binarized=zeros(1,EEG.trials);
        binarized(subject.first_rejected_epochs)=1;
        EEG = pop_rejepoch(EEG,binarized,0);
    end
    % round 2
    if ~isempty(subject.second_rejected_epochs)
        binarized=zeros(1,EEG.trials);
        binarized(subject.second_rejected_epochs)=1;
        EEG = pop_rejepoch(EEG,binarized,0);
    end
    % round 3
    if ~isempty(subject.third_rejected_epochs)
        binarized=zeros(1,EEG.trials);
        binarized(subject.second_rejected_epochs)=1;
        EEG = pop_rejepoch(EEG,binarized,0);
    end

  %% Rereference
    EEG = pop_reref(EEG, []);

  %% save to structure
    subject.EEG = EEG;

return
