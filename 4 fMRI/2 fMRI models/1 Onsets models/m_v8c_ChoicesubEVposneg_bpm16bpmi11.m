function [variables c] = model(datastruc, col, c)
% Choice regressors + subjective EV (Value of the gamble), split into positive and negative, for each subject
%   EV for cF: pLoss*SubjectiveLoss + (1-pLoss)*NTok     [ V(Accept) ]
%   EV for ct: (1-Entropy)*NTok                                        [ V(Not Exploring) ]
%   Value regressors are split by task
%       RL model: bpm16 for cF, bpmi11 for ct
%       
% Value pmods compete with each other to explain variance
% i.e. Each variable has its own onset + pmod, onset-events are removed prior to model estimate
%
% datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------

%% 

% % CHOICE regressors (these are kept in)
[variables c] =  i_choice_conditions(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
[variables c] =  i_choice_conditions(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control


% NOTE: existing v8 (models already run) have regressors named 'subEV'.
% This notation is WRONG because all the value models are subjective - use
% of objective/optimal values is just a variant on each model is reflected
% in the MODEL used (e.g. model v6e uses values from model b010b1, i.e.
% optimal).
% The already-run v8 models are left as is because those models were not
% productive anyway. BUT future notation (AND the onsets script) refers to
% this valuation [V(Accept) for cF, V(Accept/Reject) for ct] as vGamble
%
% [ To use this script as is (with bad notation): ] ########################
% - Comment out the readme error command below - but re-insert into code as soon as done running!
% - May need to alter script that exports model values to the fMRI analysis:
%       D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\a6_exportmodelpars4fmri.m
%       scol.vGamble --> scol.subEV & scol.vGamblepos/neg --> scol.subEVpos/neg
error('See model code in v8c!!'); % Comment this out to use this script anyway!
%
col.subEVpos=col.vGamblepos;
col.subEVneg=col.vGambleneg;

% Really, future models shoould rename 'subEVpos' --> vGamblepos etc. If
% you think this model thread is going to be important later on, just
% fecking RERUN v8c with the right names (vGamble). 

% VALUE regressors (cF/ct trials, with pmods describing values)
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), variables, 'cF_', {'subEVpos'; 'subEVneg'}, col, c ); % Conflict
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), variables, 'ct_', {'subEVpos'; 'subEVneg'}, col, c ); % Control



% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end

