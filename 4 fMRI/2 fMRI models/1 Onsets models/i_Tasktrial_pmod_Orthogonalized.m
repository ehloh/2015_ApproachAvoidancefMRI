function [variables c] =  i_tasktrial_pmod_Orthogonalized(data, variables, prefix,  varnames, col, c )
% Create trial event (all trials in data treated as the same type of event) + pmods describing requested parameters
%       Pmods/variables  are orthogonalized by SPM (i.e. order of entry matters)
%
%   varnames: cell with names of the variables to attach as pmods
%         should correspond to the parameter column names
%
% ---------------------------------------------------------------

%%

% Trial event
variables.names{c}=[ prefix 'onset'];
variables.onsets{c}=data(:,col.Onset_Offer);
variables.durations{c}=data(:,col.Duration_Offer);
variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0

% Variables as pmods, attached to the same event
for p=1:length(varnames)
    eval(['colnum=col.' varnames{p} ';']) % Which variable column?
    
    % Parametric modulators - Choice RT
    variables.pmod(c).name{p}=[prefix varnames{p}];
    variables.pmod(c).param{p}=data(:,colnum);
    variables.pmod(c).poly{p}=1;
end

c=c+1;

end

