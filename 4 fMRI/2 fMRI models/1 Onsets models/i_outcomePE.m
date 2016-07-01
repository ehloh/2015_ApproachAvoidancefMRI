function [variables c] = i_outcomePE(data, prefix, WhichVar, variables, col, c )
% [variables c] = i_outcomePE(data, WhichVar, variables, c )
%   Create variables for (of-interest) outcome regressors
%
%   Note: ONLY valid (correct-response) trials are counted. If >1 outcome
%   variable/pmod is attached to outcome onset, they should be set up to
%   compete (i.e. Duplicate 'Outcome_onset' regressor needs to be removed
%   after model specification, after model estimation).
% ------------------------------------------------------------------------------------

%% Format data

% WHICH trials?
alldata=data;
data= data(data(:, col.OutcomePresented)==1,:);  % Outcome presented 
data= data(data(:, col.TrialValid)==1,:);  % Exclude errors

% Outcome timings
d.OutcomeOnset=data(:, col.Onset_Outcome);
d.OutcomeDuration=data(:, col.Duration_Outcome);
d.OutcomeDuration=0; % Reset to 0

% Calculate quantities
d.OutcomeMagnitude=data(:, col.OutcomeMagnitude);
d.OutcomeValence=d.OutcomeMagnitude;
d.OutcomeValence(d.OutcomeValence>0)=1;
d.OutcomeValence(d.OutcomeValence<0)=-1;

% Outcome Prediction errors (Putative)
d.OutcomePE=data(:, col.PEoutcome);

% Cue values at outcome
d.OutcomevGamble=data(:, col.vGamble);
d.OutcomevModalchoice=data(:, col.vModalchoice);
d.OutcomevChosen=data(:, col.vChosen);

%% Construct regressors
% They are set up to compete!

for p=1:length(WhichVar)
    if strcmp(WhichVar{p}, 'Presented')==1;
        variables.names{c}=[prefix 'OutcomePresentationOnset'];
        variables.onsets{c}=d.OutcomeOnset;
        variables.durations{c}=d.OutcomeDuration;
        
    else
        variables.names{c}=[prefix 'Outcome_onset'];
        variables.onsets{c}=d.OutcomeOnset;
        variables.durations{c}=d.OutcomeDuration;
        
        % Which outcome variable?
        variables.pmod(c).name{1}=[prefix 'Outcome' WhichVar{p}];
        eval(['variables.pmod(c).param{1}=d.Outcome' WhichVar{p} ';'])
        variables.pmod(c).poly{1}=1;
    end
    
    c=c+1; 
end

% Add error feedback here if you'd like, but it's probably not necessary (error trials are already coded)

end

