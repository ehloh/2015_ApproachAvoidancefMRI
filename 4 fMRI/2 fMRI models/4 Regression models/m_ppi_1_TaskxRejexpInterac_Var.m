function [ matlabbatch] = m_ppi_1_TaskxRejexpInterac_Var(log,where,behvar)
% [ matlabbatch] = m_ppi_1_TaskxRejexpInterac_Var(log,where,behvar)
% In-going contrast: Positive PPI contrast from specified PPI model
%
% -----------------------------------------------------------------------------

% Execute:   behvar=log.behvar;

%% Checks

if strcmp(log.ppimodel(length(log.ppimodel)-10:end), 'TaskxChoice')==0
    error('Invalid first-level PPI model chosen for this type of regression model.')
end

%%

TargetConName='PPI Pos';

% Set up batch
[ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName);

end
