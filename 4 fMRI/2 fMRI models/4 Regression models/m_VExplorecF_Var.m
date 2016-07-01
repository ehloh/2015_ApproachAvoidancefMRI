function [ matlabbatch] = m_VExplorecF_Var(log,where,behvar)
% [ matlabbatch] = m_VExplorecF_Var(log,where,behvar)
% In-going contrasts: VExplore contrast from the t1_TrialType model
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%%

TargetConName='VExplore_cF';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end
