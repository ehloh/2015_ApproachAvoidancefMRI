function [ws c log] = m4_FullNonOrth_combtaskonsets(ws, c)
% Full model: EnvThreat, NTokens, pLoss, Entropy, Conflict, OutcomeMean, OutcomeVariance
% Non-orthogonalized: Each variable has its own onset + pmod, Duplicate
%                                 onsets are ignored in first-level contrasts, but included in the model specification
% Combined tasks: Onsets & Durations for both tasks are combined, but PMods are specified according to task
% - The reason for combining both task onsets is to decorrelate the trial
%   onsets with Accepting (since, in the Conflict task, majority of
%   resposnes = Accept).
%
%  -----------------------------------------------

log.pmods={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'OutcomeMean'; 'OutcomeVariance'};

for t=1:2 % For both tasks
    for p=1:length(log.pmods)
        
        % Duplicate onsets + durations (includes both tasks)
        ws.v.onsets{c}=ws.tt.onsets;
        ws.v.durations{c}=ws.tt.durations;
        
        
        switch t
            case 1
                ws.v.names{c}=['allonsets_cF_p' log.pmods{p}];
                ws.v.pmod(c).name{1}=['pcF_' log.pmods{p}];
                eval(['ws.v.pmod(c).param{1}=ws.tt.ct_' log.pmods{p} ';'])  % parameter vector
            case 2
                ws.v.names{c}=['allonsets_ct_p' log.pmods{p}];
                ws.v.pmod(c).name{1}=['pct_' log.pmods{p}];
                eval(['ws.v.pmod(c).param{1}=ws.tt.cF_' log.pmods{p} ';'])  % parameter vector
        end
        ws.v.pmod(c).poly{1}=1;
        
        %
        c=c+1;
    end
end

end

