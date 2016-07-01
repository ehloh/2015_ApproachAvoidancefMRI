function [ws c log] = m5_FullOrthog_combtaskonsets(ws, c)
% Full model (in order): EnvThreat, NTokens, pLoss, Entropy, Conflict, OutcomeMean, OutcomeVariance
% Orthogonalized: Each variable is applied as a parametric modulator to the trial onset 
%                           (Note: Order of variables matters!)
% Combined tasks: Onsets & PMods include trials from both tasks (but 0s for
%                           trials in the task-specific pmods that do not belong to that task)
% - The reason for combining both task onsets is to decorrelate the trial
%   onsets with Accepting (since, in the Conflict task, majority of
%   resposnes = Accept).
%
% - Duplicate onsets (one for each task - though each onset vector includes
%   onsets for all tasks), and attach the pMods for each task to separate set
%   of onsets (i.e. AllOnsets_ControltaskallPmods). 
% - Onsets are duplicated so as to avoid the complaint that pMods for both tasks are
%   orthoganlized by SPM - even though technically this isn't a problem
%   since the pMods for variables are orthogonal across tasks anyway (i.e.
%   cF_Envthreat is orthogonal to ct_EnvThreat, ct_NTokens...  etc).
%
%  -----------------------------------------------

log.pmods={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'OutcomeMean'; 'OutcomeVariance'};

%%

for t=1:2
    % One set of onsets + durations PER task (all pmods for that task are attached to it)
    ws.v.onsets{c}=ws.tt.onsets;
    ws.v.durations{c}=ws.tt.durations;
    switch t
        case 1
            ws.v.names{c}='allonsets_cF';
        case 2
            ws.v.names{c}='allonsets_ct';
    end
    
    % Parametric modulators
    for p=1:length(log.pmods)
        switch t
            case 1
                ws.v.pmod(c).name{p}=['pcF_' log.pmods{p}];
                eval(['ws.v.pmod(c).param{p}=ws.tt.ct_' log.pmods{p} ';'])  % parameter vector
            case 2
                ws.v.pmod(c).name{p}=['pct_' log.pmods{p}];
                eval(['ws.v.pmod(c).param{p}=ws.tt.cF_' log.pmods{p} ';'])  % parameter vector
        end
        ws.v.pmod(c).poly{p}=1;
        
        %
    end
    
    %
    c=c+1;
end

end

