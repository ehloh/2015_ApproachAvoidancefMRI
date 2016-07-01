function [variables c] =  i_Choice_pmod_Competing(data, variables, prefix,  varnames, col, c )
% Create task x choice events + pmods describing requested parameters
%       Pmods/variables  are allowed to compete for variance 
%
%   varnames: cell with names of the variables to attach as pmods
%         should correspond to the parameter column names
%
% ---------------------------------------------------------------

%%

% Choice names
if strcmp(prefix(1:2), 'cF')==1; choicenames={'Acc';'Rej';'Exp'};
elseif strcmp(prefix(1:2), 'ct')==1; choicenames={'NoB';'Bom';'Exp'};
else error('Choice names unspecified!')
end

for ch=1:3
    cdata=data(data(:,col.Resp1)==ch,:);
    
    for p=1:length(varnames)
        eval(['colnum=col.' varnames{p} ';']) % Which variable column?
        
        variables.names{c}=[prefix choicenames{ch}];
        variables.onsets{c}=cdata(:,col.Onset_Offer);
        variables.durations{c}=cdata(:,col.Duration_Offer);
        variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0
        
        % Parametric modulators
        variables.pmod(c).name{1}=[prefix  choicenames{ch} '_' varnames{p}];
        variables.pmod(c).param{1}=cdata(:,colnum);
        variables.pmod(c).poly{1}=1;
        
        c=c+1;
    end
    
end

end

