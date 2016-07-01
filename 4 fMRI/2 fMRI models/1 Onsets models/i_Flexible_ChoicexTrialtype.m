function [variables c] = i_Flexible_ChoicexTrialtype(data, prefix, choicenames,variables, col, c )
% [variables c] = i_Flexible_ChoicexTrialtype(data, prefix, choicenames,variables, col, c )
% 'Trial Type' first-level model - allows for flexible analysis at 2nd level
%
% ------------------------------------------------------------------------------------------

% Evaluate to debug: data=datastruc{1}, prefix='cF_'; choicenames={'Accept';'Reject';'Explore'}; 

%% Create regressor variables

for j=1:3 % Choices
    for i=1:36
        datatrial=data(data(:,col.TrialType)==i & data(:,col.Resp1)==j,:);
        
        if isempty(datatrial)==0
            variables.names{c}=[prefix choicenames{j} '_t' num2str( ceil(i/6)  ) '-' num2str(i-  (ceil(i/6)-1)*6   )];
            variables.onsets{c}=datatrial(:, col.Onset_Offer);
            variables.durations{c}=datatrial(:, col.Duration_Offer);
            variables.durations{c}=zeros(size(variables.onsets{c})); % Reset durations to 0
           c=c+1;
        end
    end
end

end

