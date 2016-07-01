clear all; close all hidden, clc

% Load linked
cd('D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\6 Archived approaches\2 Link cFct params 14April15')
load('res_fitmodels_cFct (14-Apr-2015) bpjm_fe_owyq L22726.mat'), lr=r_res{1,2}; ld=details.models; 

% Load separate cF
cd('D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs')
load('res_fitmodels_cF (09-Apr-2015) all.mat') ; fcf=r_res{strcmp( r_res(:,1), 'bpjm16_feow'),2}; cfd=details.models(strcmp( details.models(:,1), 'bpjm16_feow'), :);

% Load separate ct
load('res_fitmodels_ct (09-Apr-2015) all.mat') ; fct=r_res{strcmp( r_res(:,1), 'bpjm19_eyw'),2}; ctd=details.models(strcmp( details.models(:,1), 'bpjm19_eyw'), :);

%% Comparing linked and separate runs: some sensible sanity checks

% For params that apply to cF and ct, does the linked fit correspond to the average of the separate?
for parnum=1:4
    
    subplot(1,4,parnum)
    %
    ws.l=lr(:,3+parnum);
    ws.cf=fcf(:,3+parnum) ;
    ws.ct= fct(:,3+parnum);
    [ws.r ws.p]=corr(ws.l,     mean([ws.cf ws.ct ],2));
    scatter(ws.l,     mean([ws.cf ws.ct ],2)), lsline, axis square
    title(['[' ld{3}{parnum}   ' - ' cfd{3}{parnum} ' - '  ctd{3}{parnum} ']   r= '  num2str(ws.r) ' , p=' num2str(ws.p)])
    
    % Is the linked an average of the two?
    ws.d=sortrows([ws.l ws.cf ws.ct],2);
    for s=1:20
        if ws.d(s,1)>ws.d(s,3) && ws.d(s,1)<ws.d(s,2);
            ws.d(s,4)=1;
        elseif ws.d(s,1)<ws.d(s,3) && ws.d(s,1)>ws.d(s,2);
            ws.d(s,4)=1;
        else  ws.d(s,4)=0;
        end
    end
    xlabel(['Linked is between cF and ct: ' num2str(sum(ws.d(:,4)))])
    
    xlim([min( [ws.l;  mean([ws.cf ws.ct ],2)] ) max( [ws.l;  mean([ws.cf ws.ct ],2)] )])
    ylim([min( [ws.l;  mean([ws.cf ws.ct ],2)] ) max( [ws.l;  mean([ws.cf ws.ct ],2)] )])
    
    ws=[];
    
end

% Likelihood penalty?        % BIC: 2*nll+ K*ln(n_trials)


ws.nl=[lr(:,2) sum([fcf(:,2) fct(:,2)],2)];
ws.bic=[lr(:,1) sum([fcf(:,1) fct(:,1)],2)];
ws.bicpen= ws.bic- ws.nl*2;
ws.plot = ws.bicpen; 
ws.plot = ws.nl; 
%
barwitherr(std(ws.plot)/sqrt(20), mean(ws.plot))
subplot(1,3,1), barwitherr(std(ws.bic)/sqrt(20), mean(ws.bic), 'y'); title('BIC', 'FontSize', 20), set(gca, 'xticklabel', {'Linked';'Separate'}), axis square
subplot(1,3,2), barwitherr(std(ws.nl)/sqrt(20), mean(ws.nl), 'y'); title('nLL', 'FontSize', 20), set(gca, 'xticklabel', {'Linked';'Separate'}), axis square
subplot(1,3,3), barwitherr(std(ws.bicpen)/sqrt(20), mean(ws.bicpen), 'y'); title('BIC penalty', 'FontSize', 20), set(gca, 'xticklabel', {'Linked';'Separate'}), axis square
ws=[]; 

 

