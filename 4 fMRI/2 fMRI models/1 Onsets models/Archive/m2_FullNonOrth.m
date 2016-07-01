function [ws c log] = m2_FullNonOrth(ws, c)
% Full model: EnvThreat, NTokens, pLoss, Entropy, Conflict, OutcomeMean, OutcomeVariance
% Non-orthogonalized: Each variable has its own onset + pmod, Duplicate
%                                 onsets are ignored in first-level contrasts, but included in the model specification
%                                  
%  -----------------------------------------------

log.pmods={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'OutcomeMean'; 'OutcomeVariance'};

for t=1:2 % For both tasks - each task's regressors are kept separate
    for p=1:length(log.pmods)
        switch t
            case 1
                ws.v.names{c}=['cF_' log.pmods{p}];
            case 2
                ws.v.names{c}=['ct_' log.pmods{p}];
        end
        ws.v.onsets{c}=ws.ts.onsets{t};
        ws.v.durations{c}=ws.ts.durations{t};
        
        eval(['ws.v.pmod(c).param{1}=ws.ts.' log.pmods{p} '{t};'])  % parametric modulators
        ws.v.pmod(c).name{1}=['p' ws.v.names{c}];
        ws.v.pmod(c).poly{1}=1;
        c=c+1;
    end
end

end

