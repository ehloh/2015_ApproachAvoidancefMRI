function [variables c] = m_v2_vBestAnd_bpm16bpm11(datastruc, col, c)
% [variables c] = m_v2_BestAnd_bpm16bpm11(datastruc, col, c)
% V(Best choice) and V(2nd Best choice), for each subject
%       RL model: bpm16 for cF, bpm11 for ct 
%
% Value pmods compete with each other to explain variance
% i.e. Each variable has its own onset + pmod, onset-events are removed prior to model estimate
%
% datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------

%% 

% BASE trial events (onsets only)
[variables c] =  i_TaskTrialOnset(datastruc{1}, [], 'cF_', {}, col, c ); % Conflict
[variables c] =  i_TaskTrialOnset(datastruc{2}, variables, 'ct_', {}, col, c ); % Control

% VALUE regressors (cF/ct trials, with pmods describing values)
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}, variables, 'cF_', {'vBest';'vSecondBest'}, col, c ); % Conflict
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}, variables, 'ct_', {'vBest';'vSecondBest'}, col, c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );

end

