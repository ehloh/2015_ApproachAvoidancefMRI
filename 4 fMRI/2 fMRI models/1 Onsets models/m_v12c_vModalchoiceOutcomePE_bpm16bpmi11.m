function [variables c] = model(datastruc, col, c)
% Choice regressors + Immediate value of the gamble (i.e. without exploring)
%   EV for cF: pLoss*SubjectiveLoss + (1-pLoss)*NTok     [ V(Accept) ]
%   EV for ct: (1-Entropy)*NTok                                        [ V(Not Exploring) ]
%   Value regressors are split by task. All values are subjective according
%   to the model specified (i.e. optimal if model is b010b1)
%       RL model: bpm16 for cF, bpmi11 for ct
%       
% Value pmods compete with each other to explain variance
% i.e. Each variable has its own onset + pmod, onset-events are removed prior to model estimate
%
% datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------

%% 

% VALUE regressors (cF/ct trials, with pmods describing values)
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), [], 'cF_', {'vModalchoice'}, col, c ); % Conflict
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), variables, 'ct_', {'vModalchoice'}, col, c ); % Control

% % EXPLORE information regressors
% [variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% OUTCOME regressors (by task)
[variables c] = i_outcomePE(datastruc{1}, 'cF_', {'PE'}, variables, col, c );
[variables c] = i_outcomePE(datastruc{2}, 'ct_', {'PE'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'}, variables, col, c );

end

