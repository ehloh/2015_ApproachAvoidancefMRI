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

% Separate outcome-shown vs unshown trials
d_shown{1}=datastruc{1}(  datastruc{1}(:, col.OutcomePresented)==1 & datastruc{1}(:, col.TrialValid)==1, :);
d_shown{2}=datastruc{2}(  datastruc{2}(:, col.OutcomePresented)==1 & datastruc{1}(:, col.TrialValid)==1, :);
d_unshown{1}=datastruc{1}(  datastruc{1}(:, col.OutcomePresented)==0 & datastruc{1}(:, col.TrialValid)==1, :);
d_unshown{2}=datastruc{2}(  datastruc{2}(:, col.OutcomePresented)==0 & datastruc{1}(:, col.TrialValid)==1, :);

%% 

% BASE trial events (onsets only)
[variables c] =  i_TaskTrialOnset(datastruc{1}, [], 'cF_', {}, col, c ); % Conflict
[variables c] =  i_TaskTrialOnset(datastruc{2}, variables, 'ct_', {}, col, c ); % Control

% VALUE regressors (cF/ct trials, with pmods describing values)
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), variables, 'cF_', {'vExplore'}, col, c ); % Conflict
[variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), variables, 'ct_', {'vExplore'}, col, c ); % Control (column-choice names still= Accept/Reject/Explore)
% [variables c] =  i_Tasktrial_pmod_Competing(datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), variables, 'cF_', {'vAccept';'vReject';'vExplore'}, col, c ); % Conflict
% [variables c] =  i_Tasktrial_pmod_Competing(datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), variables, 'ct_', {'vAccept';'vReject';'vExplore'}, col, c ); % Control (column-choice names still= Accept/Reject/Explore)

% % EXPLORE information regressors
[variables c] = i_ExploreInfoVals( datastruc{1}(datastruc{1}(:,col.TrialValid)==1, :), 'cF_', {'Presented'; 'ExploreInfoPE'}, variables, col, c );
[variables c] = i_ExploreInfoVals( datastruc{2}(datastruc{2}(:,col.TrialValid)==1, :), 'ct_', {'Presented'; 'ExploreInfoPE'}, variables, col, c );

% Outcome-SHOWN regressors: outcome PEs
[variables c] = i_outcomePE(datastruc{1}, 'cF_', {'Presented'}, variables, col, c ); % Outcome onsets only
[variables c] = i_outcomePE(datastruc{2}, 'ct_', {'Presented'}, variables, col, c );
[variables c] =  i_Outcome_pmod_Competing( d_shown{1}( d_shown{1}(:, col.Resp1)==3,:), variables, 'cF_', {'ExploreInfotoOutcomePE'}, col, c ); 
[variables c] =  i_Outcome_pmod_Competing( d_shown{2}( d_shown{2}(:, col.Resp1)==3,:), variables, 'ct_', {'ExploreInfotoOutcomePE'}, col, c ); 

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'}, variables, col, c );

end

