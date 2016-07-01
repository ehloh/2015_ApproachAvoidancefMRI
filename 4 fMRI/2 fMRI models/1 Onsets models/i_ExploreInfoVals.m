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

% Identifying qualifying trials:  Explored, Outcome shown, no errors
data= data(data(:, col.TrialValid)==1,:);  % Exclude errors
data=data(   data(:, col.Resp1)==3 & data(:, col.OutcomePresented)==1, :);

% Timings
d.InfoOnset=data(:, col.Onset_Info);
d.InfoDuration=0;

% Quantities
d.vExplore=data(:, col.vExplore);
d.vExploreInfoVal=data(:, col.vExploreInfo);
d.ExploreInfoPE=data(:, col.ExploreInfoPE);
d.ExploreInfoPEsign=data(:, col.ExploreInfoPEsign);


%% Construct regressors (set up to compete)

for p=1:length(WhichVar)    
    if strcmp(WhichVar{p}, 'Presented')
        variables.names{c}=[prefix 'ExplorePresentedInfo'];
        variables.onsets{c}=d.InfoOnset;
        variables.durations{c}=d.InfoDuration;         
    else 
        variables.names{c}=[prefix 'ExploreInfo_onset'];
        variables.onsets{c}=d.InfoOnset;
        variables.durations{c}=d.InfoDuration; 

        % Quantity
        variables.pmod(c).name{1}=[prefix WhichVar{p}];
        eval(['variables.pmod(c).param{1}=d.' WhichVar{p} ';'])
        variables.pmod(c).poly{1}=1;
    end
    c=c+1; 
end

% Add error feedback here if you'd like, but it's probably not necessary (error trials are already coded)

end

