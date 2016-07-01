function [variables c] = m_t2_TrialtypeNc(datastruc, col, c)
% 'Trial Type No-Choice' first-level model - allows for flexible analysis at 2nd level
% NO choice regressors included at all.
%
% 	datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
% Evaluate to debug: datastruc=ws.d;

% Exclude errors?
%     datastruc{1}=datastruc{1}(datastruc{1}(:,col.TrialValid)==1,:);
%     datastruc{2}=datastruc{2}(datastruc{2}(:,col.TrialValid)==1,:); 

%% 

% RL VARIABLE regressors
[variables c] = i_Trialtype(datastruc{1}, 'cF_', [] , col, c ); % Conflict
[variables c] = i_Trialtype(datastruc{2}, 'ct_',  variables, col, c ); % Control

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'; 'Info'}, variables, col, c );


end

