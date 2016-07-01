function [variables c] = m_t1_Trialtype(datastruc, col, c)
% 'Trial Type' first-level model - allows for flexible analysis at 2nd level
% Choice regressors (as events) included right at start
%
% 	datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
% Evaluate to debug: datastruc=ws.d;

% Exclude errors?
%     datastruc{1}=datastruc{1}(datastruc{1}(:,col.TrialValid)==1,:);
%     datastruc{2}=datastruc{2}(datastruc{2}(:,col.TrialValid)==1,:); 

%% 

% CHOICE regressors
[variables c] =  i_choice_conditions(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
[variables c] =  i_choice_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control

% RL VARIABLE regressors
[variables c] = i_Trialtype(datastruc{1}, 'cF_', variables, col, c ); % Conflict
[variables c] = i_Trialtype(datastruc{2}, 'ct_',  variables, col, c ); % Control

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'; 'Info'}, variables, col, c );


end

