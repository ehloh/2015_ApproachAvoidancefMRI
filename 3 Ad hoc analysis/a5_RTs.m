% RT plots (both and ctrl)
clear all; close all hidden; clc
% clear all;clc




% c17 subjects included
% log.subjects={'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};
log.subjects={}; 
request.modversion='bpji08bpji11';   % cF, ct
% request.modversion='b01b01';   % cF, ct
request.rescale_IVs=1; 

for o1=1:1  % Settings
    
    w.w=pwd;
    if strcmp(w.w(1), '/')==1; 
        where.where='/Users/EleanorL/Dropbox/SCRIPPS/4 Explore experiment/3 Analysis'; where.scripts=[where.where '/4 Fit computational models']; where.data_beh='/Users/EleanorL/Dropbox/SCRIPPS/5 Explore fMRI/1 Behavioural data'; 
        
        where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/1 Brain data';
        where.moderes='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/2 Second level results/2 Behavioural model details';
    else where.moderes='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\2 Behavioural model details';
        where.where='D:\Dropbox\SCRIPPS\4 Explore experiment\3 Analysis';  where.data_beh='D:\Dropbox\SCRIPPS\5 Explore fMRI\1 Behavioural data'; where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data';
    end; where.scripts=[where.where '/4 Fit computational models'];

    % Load data
    %     Specific subject selection not implemented yet. If you want to, alter it from the start here. 
    f=spm_select('List', [where.moderes filesep request.modversion], 'Model values');  if size(f,1)~=1; error('No. modelres files ~=1!'); end
    d= load([where.moderes filesep request.modversion filesep f]);
    log.subjects= d.details.subjects;  log.n_subjs=length(log.subjects); 
    col=d.scol; log.moddetails={d.details.cf_moddetails d.details.ct_moddetails}; d_modpar=d.d_par.partable; 
    d_trialstats=d.d_trialstats;
    
    
    
    % Load more details here if you want 
    
    
end

%% GLM settings 

% Which IVs?
ivs={
% 'Choice';       % Choice 
% 'Choice_Accept';
'Choice_Reject';
'Choice_Explore';

% 
% 'EnvThreat';   % Psych
% 'NTokens';
% 'pLoss';
% 'Entropy';  
% 'EV';


'vChosen';          % Values 
'vBestUnchosen';
% 'vGamble';


% [ Exploration-related variables ] ---------------------------------------
% 'EntropyNTok';
% 'BinomVar';  % y
% 'StanDev';


};


col.Choice_Accept=51; 
col.Choice_Reject=52; 
col.Choice_Explore=53; 
request.ivs=nan(length(ivs),1);  for i=1:length(ivs), eval(['request.ivs(i)=col.' ivs{i} ';']), end; log.n_ivs=length(request.ivs);

%% GLM analysis

d_nerrors= nan(log.n_subjs,2); 
d_betas=[{nan(log.n_subjs, log.n_ivs)} {nan(log.n_subjs, log.n_ivs)}];
d_choice=nan(log.n_subjs,6); 
% First level analysis
for s=1:log.n_subjs
    for t=1:2
        ws.d=  d_trialstats{s,t+1}; 
        ws.d(:, [col.Choice_Accept col.Choice_Reject col.Choice_Explore])=0;
        ws.d(:, col.Choice_Accept) = ws.d(:, col.Choice)==1;
        ws.d(:, col.Choice_Reject) = ws.d(:, col.Choice)==2;
        ws.d(:, col.Choice_Explore) = ws.d(:, col.Choice)==3;
        ws.ivs=ws.d(:, request.ivs);

        for c=1:3
            d_choice(s,(t-1)*3 +c)= mean(ws.d(:, col.Choice)==c); 
        end
        
        if request.rescale_IVs
            ws.ivs = ws.ivs -  repmat(min(ws.ivs), size(ws.ivs ,1),1); % Adjust floor
            ws.ivs=(ws.ivs)./repmat(max(ws.ivs), size(ws.ivs ,1),1);
        end
        
        [ws.b ws.dd ws.stats]=glmfit(ws.ivs, ws.d(:,col.RT1));
        d_betas{t}(s,:) = ws.b(2:end)';  
        d_nerrors(s,t)= sum(ws.d(:,col.TrialValid)==0);
        
        
        ws=[]; 
    end
end
disp(d_nerrors)
error

% Second level plot
figure('Name', 'Betas on RTs', 'NumberTitle', 'off', 'Position', [100 250 600 800], 'Color', 'w');    f.FontSize=15; f.FontName='PT Sans Caption';   
f.yrangeslack =0.8;
subplot(3,1,1); % exp
barwitherr(std(d_betas{1})./sqrt(log.n_subjs), mean(d_betas{1}), 'y'), title('Exp','FontSize', f.FontSize, 'FontName', f.FontName);   
set(gca,'xticklabel', ivs,'FontSize', f.FontSize, 'FontName', f.FontName);
ylim([min([d_betas{1}(:); d_betas{2}(:)])*f.yrangeslack     max([d_betas{1}(:); d_betas{2}(:)])*f.yrangeslack])
subplot(3,1,2); % ctrl
barwitherr(std(d_betas{2})./sqrt(log.n_subjs), mean(d_betas{2}), 'y'), title('Ctrl','FontSize', f.FontSize, 'FontName', f.FontName);   
set(gca,'xticklabel', ivs,'FontSize', f.FontSize, 'FontName', f.FontName);
ylim([min([d_betas{1}(:); d_betas{2}(:)])*f.yrangeslack     max([d_betas{1}(:); d_betas{2}(:)])*f.yrangeslack])
subplot(3,1,3); % ctrl
barwitherr(std(d_betas{1}-d_betas{2})./sqrt(log.n_subjs), mean(d_betas{1}-d_betas{2}), 'y'), title('Exp - Ctrl','FontSize', f.FontSize, 'FontName', f.FontName);   
set(gca,'xticklabel', ivs,'FontSize', f.FontSize, 'FontName', f.FontName);
ylim([min([d_betas{1}(:)- d_betas{2}(:)])*f.yrangeslack     max([d_betas{1}(:)- d_betas{2}(:)])*f.yrangeslack])

% Second level stats 
r_res= [{' ';'cF';'ct';'cF-ct'} [ivs'; cell(3,log.n_ivs)]];  k=2; openvar r_res
[wr.h wr.p wr.ci wr.stats]=ttest(d_betas{1});   % cf
r_res(k, 2:end) = num2cell(wr.h);k=k+1;
[wr.h wr.p wr.ci wr.stats]=ttest(d_betas{2});    % ct 
r_res(k, 2:end) = num2cell(wr.h);k=k+1;
[wr.h wr.p wr.ci wr.stats]=ttest(d_betas{1} - d_betas{2});    % cf - ct
r_res(k, 2:end) = num2cell(wr.h);k=k+1;


%% OTHER analysis

% How many explore-info trials?
d_ei=log.subjects; 
for s=1:log.n_subjs
    for t=1:2
        wt.d=d_trialstats{s,t+1};
        wt.d= wt.d(wt.d(:, col.Choice)==3,:); 
        wt.d= wt.d(wt.d(:, col.OutcomePres)==1,:); 
        
        d_ei{s,1+(t-1)*2+1}= sum(wt.d(:, col.ExploredBomb));
        d_ei{s,1+(t-1)*2+2}= sum(1-wt.d(:, col.ExploredBomb));
        
        
%         d_ei{s,t+1}= sum(wt.d(:, col.ExploredBomb));
    end
end
d_ei


% Compare % choice across cF vs ct

figure, barwitherr(std(d_choice(:,1:6))./sqrt(log.n_subjs), mean(d_choice(:,1:6)))
[choice_anova]=teg_repeated_measures_ANOVA(d_choice(:,1:6),[2 3], {'Task' 'Choice'});
[{' ' 'F' 'df1' 'df2' 'p' ' ' ' '};  [choice_anova.labels' num2cell(choice_anova.R)]]

% Combined cF vs ct 
for c=1:3
    d_choice(:, 6+c) = mean([d_choice(:,c) d_choice(:,3+c)], 2);
end
% figure, barwitherr(std(d_choice(:,7:9))./sqrt(log.n_subjs), mean(d_choice(:,7:9)))

[h p ci stats]=ttest( d_choice(:,7) - d_choice(:,8));
disp(['A vs R: t(' num2str(stats.df) ')= '  num2str(stats.tstat) '    p=' num2str(p)])
[h p ci stats]=ttest( d_choice(:,7) - d_choice(:,9));
disp(['A vs E: t(' num2str(stats.df) ')= '  num2str(stats.tstat) '    p=' num2str(p)])
[h p ci stats]=ttest( d_choice(:,8) - d_choice(:,9));
disp(['R vs E: t(' num2str(stats.df) ')= '  num2str(stats.tstat) '    p=' num2str(p)])

