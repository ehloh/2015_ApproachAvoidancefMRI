function [ matlabbatch] = m_ppi_2_Rej_cFMct_Var(log,where,behvar)
% [ matlabbatch] = m_ppi_2_Rej_cFMct_Var(log,where,behvar)
% In-going contrasts: cF_Reject-Explore 
%                           (Simple effects contrast set up at 1st level in par models)
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%% Checks

if strcmp(log.ppimodel(length(log.ppimodel)-8:end), 'Rej_cF-ct')==0
    error('Invalid first-level PPI model chosen for this type of regression model.')
end

%%

TargetConName='PPI Pos';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end