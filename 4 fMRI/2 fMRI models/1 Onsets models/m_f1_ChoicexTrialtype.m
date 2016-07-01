function [variables c] = m_f1_ChoicexTrialtype(datastruc, col, c)
% Choice x Trialtype first-level model - allows for fully flexible analysis at 2nd level
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
% [variables c] =  i_choice_conditions(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
% [variables c] =  i_choice_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control

% RL VARIABLE regressors
% [variables c] = i_Flexible_ChoicexTrialtype(datastruc{1}, 'cF_', {'Accept';'Reject';'Explore'}, variables, col, c ); % Conflict
[variables c] = i_Flexible_ChoicexTrialtype(datastruc{1}, 'cF_', {'Accept';'Reject';'Explore'}, [], col, c ); % Conflict
[variables c] = i_Flexible_ChoicexTrialtype(datastruc{2}, 'ct_', {'NoBomb';'Bomb';'Explore'}, variables, col, c ); % Control

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'; 'Info'}, variables, col, c );


end

