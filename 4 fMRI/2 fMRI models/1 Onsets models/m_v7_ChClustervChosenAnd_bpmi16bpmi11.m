function [variables c] = m_v7_ChClustervChosenAnd_bpm16bpm11(datastruc, col, c)
% Choice split as a function of exploration cluster or not  + V(Chosen) and V(Best Unchosen)
%       RL model: bpm16 for cF, bpm11 for ct 
%
% Value pmods compete with each other to explain variance
% i.e. Each variable has its own onset + pmod, onset-events are removed prior to model estimate
%
% datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------

%% 

% Define inner cluster: 6 cells
InnerCluster=[3 4; 3 5; 3 6; 4 4; 4 5; 4 6];

% CHOICE regressors
[variables c] =  i_choice_conditions_cluster(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c,InnerCluster ); % Conflict
[variables c] =  i_choice_conditions_cluster(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c,InnerCluster ); % Control

% VALUE regressors (cF/ct trials, with pmods describing values)
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), variables, 'cF_', {'vChosen'; 'vBestUnchosen'}, col, c ); % Conflict
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), variables, 'ct_', {'vChosen'; 'vBestUnchosen'}, col, c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end

