function [ matlabbatch] = m_RejMExp_Var(log,where,behvar)
% [ matlabbatch] = m_RejMExp_Var(log,where,behvar)
% In-going contrasts: cF_Reject-Explore 
%                           (Simple effects contrast set up at 1st level in par models)
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%%

TargetConName='Reject-Explore';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end