function [variables c] = i_ExploredInfo(data, ExploreVars, variables, col, c )
% [variables c] = i_ExploredInfo(data, ExploreVars, variables, col, c )
%   Create variables for Explored-info regressors (ExploredBomb vs ExploredNobomb)
%       with parametric modulators for requested RL variables
%
%   All variables: {'Error'; 'Motor'; 'OutcomePresented'; 'Info'}
% ------------------------------------------------------------------------------------

% Execute to debug: data=datastruc{3}; ExploreVars={'Error'; 'Motor'; 'OutcomePresented'; 'Info'}; 

task_prefixes={'cF';'ct'};
info_prefixes={'ExploredBomb' 1 ;'ExploredNobomb' 0};  % Name, specification in col.ShowedExploredBomb

for t=1:2
    for b=1:2  % Explored information types
                
        % Trials of the correst task where subject explores + certain bomb/nobomb info is actively shown
        wt=data(data(:, col.Task)==t & data(:,col.Resp1)==3 & data(:,col.OutcomePresented)==1 & data(:, col.ShowedExploredBomb)==info_prefixes{b, 2}, :);
        variables.names{c}=[task_prefixes{t} '_' info_prefixes{b} 'Info'];
        variables.onsets{c}=wt(:,col.Onset_Info);
        variables.durations{c}=0;
        
        % Variables as parametric modulators
        for p=1:length(ExploreVars)
            eval(['cn=col.' ExploreVars{p} ';'])
            variables.pmod(c).name{p}=[task_prefixes{t} '_' info_prefixes{b} 'x' ExploreVars{p}];
            variables.pmod(c).param{p}=wt(:, cn);
            variables.pmod(c).poly{p}=1;
        end
        %
        c=c+1;
    end
end

end

