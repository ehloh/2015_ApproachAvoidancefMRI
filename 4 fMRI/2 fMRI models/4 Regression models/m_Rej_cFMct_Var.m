function [ matlabbatch] = m_Rej_cFMct_Var(log,where,behvar)
% [ matlabbatch] = m_Rej_cFMct_Var(log,where,behvar)
% In-going contrasts: cF_Reject-Explore 
%                           (Simple effects contrast set up at 1st level in par models)
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%%

TargetConName='Reject_cF-ct';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end