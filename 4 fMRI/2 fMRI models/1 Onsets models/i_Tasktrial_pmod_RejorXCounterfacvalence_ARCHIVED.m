function [variables c] =  i_Tasktrial_pmod_RejorXCounterfacvalence(data, variables, prefix,  varnames, col, c )
% Create trial event (all trials in data treated as the same type of event) + pmods describing requested parameters
%       Pmods/variables  are allowed to compete for variance 
%
%   varnames: cell with names of the variables to attach as pmods
%         should correspond to the parameter column names
%
% ---------------------------------------------------------------

%%

% vChosen regressors (not split by choice)
[variables c] =  i_Tasktrial_pmod_Competing( data, variables, prefix, {'vChosen'}, col, c );

% vBestUnchosen: Reject trials
data_r=data(data(:, col.Resp1)==2, :);
[variables c] =  i_Tasktrial_pmod_Competing( data_r, variables, prefix, {'Rej_vBestUnchosen_pos';'Rej_vBestUnchosen_neg'}, col, c );

% vBestUnchosen: Non-Reject trials
data_or=data(data(:, col.Resp1)~=2, :);
[variables c] =  i_Tasktrial_pmod_Competing( data_or, variables, prefix, {'NonRej_vBestUnchosen_pos';'NonRej_vBestUnchosen_neg'}, col, c );

end

