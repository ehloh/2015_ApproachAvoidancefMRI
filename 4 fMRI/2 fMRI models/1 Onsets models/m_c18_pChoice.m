function [variables c] = m_c3_ChoiceFull_ULPEN(datastruc, col, c)
% Full model: Entropy (U), EV (L), pLoss (P), EnvThreat (E), NTokens (N)
% Competing: Each variable has its own onset + pmod
%                   Duplicate onsets are ignored in first-level contrasts
%
%       datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
% Evaluate to debug: datastruc=ws.d;

% RLvariables={'EnvThreat'; 'NTokens';'pLoss'; 'Entropy'; 'EV'; };

% Exclude errors?
%     datastruc{1}=datastruc{1}(datastruc{1}(:,col.TrialValid)==1,:);
%     datastruc{2}=datastruc{2}(datastruc{2}(:,col.TrialValid)==1,:); 

%% 

variables=[];

% BASE trial events (onsets only)
[variables c] =  i_TaskTrialOnset(datastruc{1}, variables, 'cF_', {}, col, c ); % Conflict
[variables c] =  i_TaskTrialOnset(datastruc{2}, variables, 'ct_', {}, col, c ); % Control

% % Observed Probability(Choice) regressors
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), variables, 'cF_', {'pReject';'pExplore'}, col, c ); % Conflict
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), variables, 'ct_', {'pReject';'pExplore'}, col, c ); % Control (column-choice names still= Accept/Reject/Explore)


% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );

end

