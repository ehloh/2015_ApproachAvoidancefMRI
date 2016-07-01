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

%% 

data_r{1}=datastruc{1}(datastruc{1}(:, col.TrialValid)==1 & datastruc{1}(:, col.Resp1)==2, :);
data_nr{1}=datastruc{1}(datastruc{1}(:, col.TrialValid)==1 & datastruc{1}(:, col.Resp1)~=2, :);
data_r{2}=datastruc{2}(datastruc{2}(:, col.TrialValid)==1 & datastruc{2}(:, col.Resp1)==2, :);
data_nr{2}=datastruc{2}(datastruc{2}(:, col.TrialValid)==1 & datastruc{2}(:, col.Resp1)~=2, :);

%%
variables=[];

% CHOICE regressors
% [variables c] =  i_choice_conditions(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
% [variables c] =  i_choice_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control

% BASE trial events (onsets only)
[variables c] =  i_TaskTrialOnset(datastruc{1}, variables, 'cF_', {}, col, c ); % Conflict
[variables c] =  i_TaskTrialOnset(datastruc{2}, variables, 'ct_', {}, col, c ); % Control

% cF 
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:, col.TrialValid)==1,:), variables, 'cF_', {'vChosen'}, col, c );                          % vChosen (not split by choice)
[variables c] =  i_Tasktrial_pmod_Competing( data_r{1}, variables, 'cF_', {'Rej_vBestUnchosen';}, col, c );             % vBU split by RejOr
[variables c] =  i_Tasktrial_pmod_Competing( data_nr{1}, variables, 'cF_', {'NonRej_vBestUnchosen';}, col, c );   

% ct
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:, col.TrialValid)==1,:), variables, 'ct_', {'vChosen'}, col, c );                          % vChosen (not split by choice)
[variables c] =  i_Tasktrial_pmod_Competing( data_r{2}, variables, 'ct_', {'Rej_vBestUnchosen';}, col, c );             % vBU split by RejOr
[variables c] =  i_Tasktrial_pmod_Competing( data_nr{2}, variables, 'ct_', {'NonRej_vBestUnchosen';}, col, c );   

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end

