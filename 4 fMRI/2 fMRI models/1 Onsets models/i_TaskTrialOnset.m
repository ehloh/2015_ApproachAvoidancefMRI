function [variables c] =  i_TaskTrialOnset(data, variables, prefix,  varnames, col, c )
% Create trial event (all trials in data treated as the same type of event)
%   alone. Other onsets (created for the pmods) will be removed, but this
%   onset will stay in.
%
% ---------------------------------------------------------------

%%

variables.names{c}=[prefix 'trialonset'];
variables.onsets{c}=data(:,col.Onset_Offer);
variables.durations{c}=data(:,col.Duration_Offer);
variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0

c=c+1;

end

