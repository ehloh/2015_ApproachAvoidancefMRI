function [ws c log] = m3_FullOrthog(ws, c)
% Full model (in order): EnvThreat, NTokens, pLoss, Entropy, Conflict, OutcomeMean, OutcomeVariance
% Orthogonalized: Each variable is applied as a parametric modulator to the trial onset 
% Note: Order of variables matters!
%                                  
%  -----------------------------------------------

log.pmods={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'OutcomeMean'; 'OutcomeVariance'};

for t=1:2 % For both tasks - each task's regressors are kept separate
    ws.v.onsets{c}=ws.ts.onsets{t};  % Onsets
    switch t
        case 1
            ws.v.names{c}='cFonset';
        case 2
            ws.v.names{c}='ctonset';
    end
    ws.v.durations{c}=ws.ts.durations{t};
    
    for p=1:length(log.pmods)
        switch t
            case 1
                ws.v.pmod(c).name{p}=['pcF_' log.pmods{p}];
            case 2
                ws.v.pmod(c).name{p}=['pct_' log.pmods{p}];
        end
        eval(['ws.v.pmod(c).param{p}=ws.ts.' log.pmods{p} '{t};'])  % parametric modulators
        ws.v.pmod(c).poly{p}=1;
    end
    c=c+1;
end

end

