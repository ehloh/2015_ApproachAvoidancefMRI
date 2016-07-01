function [variables c] = i_Trialtype(data, prefix, variables, col, c )
% [variables c] = i_Trialtype(data, prefix, variables, col, c )
% 'Trial Type' first-level model - allows for flexible analysis at 2nd level
%
% ------------------------------------------------------------------------------------------

% Evaluate to debug: data=datastruc{1}, prefix='cF_'; 

%% Create regressor variables

for i=1:36
    datatrial=data(data(:,col.TrialType)==i,:);
    
    variables.names{c}=[prefix 't' num2str( ceil(i/6)  ) '-' num2str(i-  (ceil(i/6)-1)*6   )];
    variables.onsets{c}=datatrial(:, col.Onset_Offer);
    variables.durations{c}=datatrial(:, col.Duration_Offer);
    variables.durations{c}=zeros(size(variables.onsets{c})); % Reset durations to 0
    
    %
    c=c+1;
end

end
