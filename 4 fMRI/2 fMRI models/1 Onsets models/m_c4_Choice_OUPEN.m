function [variables c] = m_c4_Choice_OUPEN(datastruc, col, c)
% Full -1 model: EntropyNTok (O), Entropy (U), pLoss (P), EnvThreat (E), NTokens (N)
% Competing: Each variable has its own onset + pmod
%                   Duplicate onsets are ignored in first-level contrasts
%
%       datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
% Evaluate to debug: datastruc=ws.d;

RLvariables={'EntropyNTok'; 'Entropy'; 'pLoss'; 'EnvThreat'; 'NTokens'}; % XUPEN

% Exclude errors?
%     datastruc{1}=datastruc{1}(datastruc{1}(:,col.TrialValid)==1,:);
%     datastruc{2}=datastruc{2}(datastruc{2}(:,col.TrialValid)==1,:); 

%% 

% CHOICE regressors
[variables c] =  i_choice_conditions(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
[variables c] =  i_choice_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control

% RL VARIABLE regressors
[variables c] = i_RLvariables_Competing(datastruc{1}, 'cF_', RLvariables, variables, col, c ); % Conflict
[variables c] = i_RLvariables_Competing(datastruc{2}, 'ct_', RLvariables, variables,col,  c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end
