function [ matlabbatch] = m_ppi_4_cF_RejMExp_Var(log,where,behvar)
% [ matlabbatch] = m_ppi_4_cF_RejMExp_Var(log,where,behvar)
% In-going contrast: Positive PPI contrast from specified PPI model
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%% Checks

if strcmp(log.ppimodel(length(log.ppimodel)-9:end), 'cF_Rej-Exp')==0;
    error('Invalid first-level PPI model chosen for this type of regression model.')
end

%%

TargetConName='PPI Pos';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end

