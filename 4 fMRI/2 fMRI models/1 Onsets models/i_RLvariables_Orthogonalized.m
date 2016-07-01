function [variables c] = i_RLvariables_Orthogonalized(data, prefix, RLvariables, variables, col, c )
% [variables c] = i_RLvariables_Orthogonalized(data, prefix, RLvariables, variables, col, c )
% Orthogonalized: Each variable is applied as a parametric modulator to the trial onset 
%                           (Note: Order of variables matters!)
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

% Onsets (Main condition regressor)
variables.names{c}=[prefix 'onset'];
variables.onsets{c}=data(:,col.Onset_Offer);
variables.durations{c}=data(:,col.Duration_Offer);
variables.durations{c}=zeros(size(variables.onsets{c})); % Reset duration to 0
    
for p=1:length(RLvariables)
  
    % Parametric modulators
    variables.pmod(c).name{p}=['p' prefix RLvariables{p}];
    eval(['variables.pmod(c).param{p}=data(:,col.' RLvariables{p} ');'])  % parametric modulators
    variables.pmod(c).poly{p}=1;
    
end

c=c+1;

end

