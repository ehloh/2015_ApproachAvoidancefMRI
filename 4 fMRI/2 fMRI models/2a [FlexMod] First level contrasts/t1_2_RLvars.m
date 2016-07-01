function [Con log] = t1_2_RLvars(where, log)
% Contrasts for all RL variables (cF, ct, both tasks, cF-ct, ct-cF)
%       See function tf_Trialtype_RLvars for details
%       Assumes regressors are labelled e.g. t1-1
% --------------------------------------------------------------------

% Set up RL var contrasts for a Trialtype model
[Con log] = tf_Trialtype_RLvars(where, log);

end

