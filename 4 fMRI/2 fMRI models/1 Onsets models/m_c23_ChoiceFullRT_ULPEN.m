function [variables c] = m_c10_ChCluster6FullRT_OULPEN(datastruc, col, c)
% Full model, but with Choice split as a function of exploration cluster or not
%
% Competing: Each variable has its own onset + pmod
%                   Duplicate onsets are ignored in first-level contrasts
%
% Variables:  EntropyNTok (O), Entropy (U), EV (L), pLoss (P), EnvThreat (E), NTokens (N)
%
%       datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
%
% Evaluate to debug: datastruc=ws.d;

RLvariables={'Entropy';  'EV'; 'pLoss'; 'EnvThreat'; 'NTokens';}; 

% Exclude errors?
%     datastruc{1}=datastruc{1}(datastruc{1}(:,col.TrialValid)==1,:);
%     datastruc{2}=datastruc{2}(datastruc{2}(:,col.TrialValid)==1,:); 

%% 
 

% CHOICE regressors
[variables c] =  i_choiceRT_conditions(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c); % Conflict
[variables c] =  i_choiceRT_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c); % Control


% RL VARIABLE regressors
[variables c] = i_RLvariables_Competing(datastruc{1}, 'cF_', RLvariables, variables, col, c ); % Conflict
[variables c] = i_RLvariables_Competing(datastruc{2}, 'ct_', RLvariables, variables,col,  c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );

end

