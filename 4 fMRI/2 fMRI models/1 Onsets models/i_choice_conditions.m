function [variables c] =  i_choice_conditions(data, variables, prefix,  choicenames, col, c )
% Create regressors for Choice, each modelled as a condition
%
%   choicenames: cell with names, should correspond to code for choice
%                           (1=Accept/NoBomb, 2=Reject/Bomb, 3=Explore)
% ---------------------------------------------------------------

% Evaluate to debug: datastruc=ws.d; data=datastruc{1}; prefix='cF_'; choicenames={'Accept';'Reject';'Explore'};

for p=1:length(choicenames)
    choiced=data(data(:,col.Resp1)==p,:);
    
    variables.names{c}=[prefix choicenames{p}];
    variables.onsets{c}=choiced(:,col.Onset_Motor1); % Onset= Locked to time of decision
    variables.durations{c}=choiced(:,col.Duration_Choice); 
    variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0
    
    c=c+1;
end


end

