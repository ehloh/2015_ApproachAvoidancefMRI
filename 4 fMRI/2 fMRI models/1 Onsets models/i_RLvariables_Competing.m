function [variables c] = i_RLvariables_Competing(data, prefix, RLvariables, variables, col, c )
% [variables c] = i_RLvariables_Competing(data, prefix, RLvariables, variables, col, c )
% Competing model: Each variable has its own onset + pmod
%                              Duplicate onsets are ignored in first-level contrasts
%
%     All variables:
%         EnvThreat (E)
%         NTokens (N)
%         pLoss (P)
%         Uncertainty (U)
%         EV (V)
%         VExplore (X)
%         Conflict (C)
% 
% ------------------------------------------------------------------------------------------

% Evaluate to debug: data=datastruc{1}, prefix='cF_'; 

%%

for p=1:length(RLvariables)
    variables.names{c}=['o' prefix RLvariables{p}];
    variables.onsets{c}=data(:,col.Onset_Offer);
    variables.durations{c}=data(:,col.Duration_Offer);
    variables.durations{c}=zeros(size(variables.onsets{c})); % Reset duration to 0
    
    % Parametric modulators
    variables.pmod(c).name{1}=['p' prefix RLvariables{p}];
    eval(['variables.pmod(c).param{1}=data(:,col.' RLvariables{p} ');'])  % parametric modulators
    variables.pmod(c).poly{1}=1;
    
    %
    c=c+1;
end

end

