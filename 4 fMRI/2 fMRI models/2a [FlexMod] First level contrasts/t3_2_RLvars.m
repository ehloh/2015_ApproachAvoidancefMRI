function [Con log] = t3_2_RLvars(where, log)
% Contrasts for all RL variables (cF, ct, both tasks, cF-ct, ct-cF)
%       See function tf_Trialtype_RLvars for details
%       Assumes regressors are labelled e.g. cF_e1-n1
% --------------------------------------------------------------------

% Set up RL var contrasts for a Trialtype model
[Con log] = tf_Chunk4Trialtype_RLvars(where, log);

end

