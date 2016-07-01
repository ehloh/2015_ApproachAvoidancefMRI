function [variables c] =  i_Tasktrial_pmod_Competing(data, variables, prefix,  varnames, col, c )
% Create trial event (all trials in data treated as the same type of event) + pmods describing requested parameters
%       Pmods/variables  are allowed to compete for variance 
%
%   varnames: cell with names of the variables to attach as pmods
%         should correspond to the parameter column names
%
% ---------------------------------------------------------------

%%

% keyboard

for p=1:length(varnames)
    eval(['colnum=col.' varnames{p} ';']) % Which variable column?
    
    variables.names{c}=[prefix 'onset'];
    variables.onsets{c}=data(:,col.Onset_Offer);
    variables.durations{c}=data(:,col.Duration_Offer); 
    variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0
    
    % Parametric modulators - Choice RT
    variables.pmod(c).name{1}=[prefix varnames{p}];
    variables.pmod(c).param{1}=data(:,colnum);
    variables.pmod(c).poly{1}=1;
    
    c=c+1;
end

end

