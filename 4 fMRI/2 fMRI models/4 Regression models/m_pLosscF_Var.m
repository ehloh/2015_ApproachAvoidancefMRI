function [ matlabbatch] = m_pLosscF_Var(log,where,behvar)
% [ batch] = pLoss_Var(log,where)
% In-going contrasts: pLoss contrast from the t1_TrialType model
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%%

TargetConName='pLoss_cF';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end