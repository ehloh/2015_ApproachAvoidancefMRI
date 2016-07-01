function [variables c] = m_c7_ChCluster6Full_OULPEN(datastruc, col, c)
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

% Define inner cluster: 6 cells
InnerCluster=[3 4; 3 5; 3 6; 4 4; 4 5; 4 6];

% CHOICE regressors
[variables c] =  i_choice_conditions_cluster(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c,InnerCluster ); % Conflict
[variables c] =  i_choice_conditions_cluster(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c,InnerCluster ); % Control

% RL VARIABLE regressors
[variables c] = i_RLvariables_Competing(datastruc{1}, 'cF_', RLvariables, variables, col, c ); % Conflict
[variables c] = i_RLvariables_Competing(datastruc{2}, 'ct_', RLvariables, variables,col,  c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end
