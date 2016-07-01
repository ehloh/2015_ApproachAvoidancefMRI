function [ws c] = i7_choice_pmodstogether(ws, c)
% Create regressors for Choice, each modelled as a pmod
% to the same onset.
%
% ---------------------------------------------------------------

    for t=1:2 
        switch t
            case 1
                ws.v.names{c}='cF_Choice';
            case 2
                ws.v.names{c}='ct_Choice';
        end
        ws.v.onsets{c}=ws.ts.Choice_Onsets{t};
        ws.v.durations{c}=ws.ts.Choice_Durations{t};
        
        % Choice - parametric modulator
        switch t
            case 1
                ws.v.pmod(c).name='pARE';
            case 2
                ws.v.pmod(c).name='pNBE';
        end
        ws.v.pmod(c).param{1}=ws.ts.Choice{t};
        ws.v.pmod(c).poly{1}=1;
       
        %
        c=c+1;
    end
    
end

