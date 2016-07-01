function [variables c] = m_v9c_vChosenAndposneg_bpm16bpm11(datastruc, col, c)
% Choice regressors + V(Chosen) and V(Best Unchosen), for each subject
%   Value regressors are split as a function of task (but not task x choice)
%       RL model: bpm16 for cF, bpm11 for ct
%       
% Value pmods compete with each other to explain variance
% i.e. Each variable has its own onset + pmod, onset-events are removed prior to model estimate
%
% datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------

data_vbup{1}=datastruc{1}(datastruc{1}(:, col.TrialValid)==1 & datastruc{1}(:, col.vBestUnchosen)>=0, :);
data_vbun{1}=datastruc{1}(datastruc{1}(:, col.TrialValid)==1 & datastruc{1}(:, col.vBestUnchosen)<0, :);
data_vbup{2}=datastruc{2}(datastruc{2}(:, col.TrialValid)==1 & datastruc{2}(:, col.vBestUnchosen)>=0, :);

%% 

variables=[];

% CHOICE regressors
% [variables c] =  i_choice_conditions(datastruc{1}, variables, 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
% [variables c] =  i_choice_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control

% BASE trial events (onsets only)
[variables c] =  i_TaskTrialOnset(datastruc{1}, variables, 'cF_', {}, col, c ); % Conflict
[variables c] =  i_TaskTrialOnset(datastruc{2}, variables, 'ct_', {}, col, c ); % Control

% VALUE regressors (cF/ct trials, with pmods describing values)
[variables c] =  i_Tasktrial_pmod_Competing(data_vbup{1}, variables, 'cF_', {'vBUpos_vMargChosen'}, col, c ); % Conflict 
[variables c] =  i_Tasktrial_pmod_Competing(data_vbun{1}, variables, 'cF_', {'vBUneg_vMargChosen'}, col, c ); 
[variables c] =  i_Tasktrial_pmod_Competing(data_vbup{2}, variables, 'ct_', {'vBUpos_vMargChosen'}, col, c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end

