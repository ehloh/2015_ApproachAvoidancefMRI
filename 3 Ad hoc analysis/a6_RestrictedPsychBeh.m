% Examine behaviour within the in-zone (where Reject and Explore trade off)  
clear all; close all hidden; clc
% clear all;clc




% c17 subjects included
log.subs_6cell={'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};
log.subs_4cell={'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};
log.specificsubjects = log.subs_6cell; 
% log.specificsubjects = log.subs_4cell; 
% log.specificsubjects={}; 

request.modversion='bpji08bpji11';   % cF, ct
% request.modversion='b01b01';   % cF, ct
request.rescale_IVs=0; 

for o1=1:1  % Settings
    
    w.w=pwd;
    if strcmp(w.w(1), '/')==1; where.where='/Users/EleanorL/Dropbox/SCRIPPS/4 Explore experiment/3 Analysis'; where.scripts=[where.where '/4 Fit computational models']; where.data_beh='/Users/EleanorL/Dropbox/SCRIPPS/5 Explore fMRI/1 Behavioural data'; where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/1 Brain data';
        
        where.moderes='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/2 Second level results/2 Behavioural model details';
        
    else where.moderes='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\2 Behavioural model details'; 
%         where.where='D:\Dropbox\SCRIPPS\4 Explore experiment\3 Analysis'; where.scripts=[where.where '/4 Fit computational models']; where.data_beh='D:\Dropbox\SCRIPPS\5 Explore fMRI\1 Behavioural data'; where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data';
    end

    % Load data
    %     Specific subject selection not implemented yet. If you want to, alter it from the start here. 
    f=spm_select('List', [where.moderes filesep request.modversion], 'Model values');  if size(f,1)~=1; error('No. modelres files ~=1!'); end
    d= load([where.moderes filesep request.modversion filesep f]);
    log.subjects= d.details.subjects;  log.n_subjs=length(log.subjects); 
    col=d.scol; log.moddetails={d.details.cf_moddetails d.details.ct_moddetails}; d_modpar=d.d_par.partable; 
    d_trialstats=d.d_trialstats;
    
    % Subject selections 
    [w.subjects w.n_subjs d_modpar] = f_selectsubjects(d_modpar, log.specificsubjects, [{'Sub' 'All'};  [log.subjects num2cell(ones(log.n_subjs,1))]], 'All');
    [log.subjects log.n_subjs d_trialstats] = f_selectsubjects([cell(1, size(d_trialstats,2)); d_trialstats] , log.specificsubjects, [{'Sub' 'All'};  [log.subjects num2cell(ones(log.n_subjs,1))]], 'All'); d_trialstats=d_trialstats(2:end, :); 
    
    
end

conds={'cF';'ct';'cFct'; 'cF>ct'}; 
col.pAcc=60;
col.pRej=col.pAcc+1;
col.pExp=col.pAcc+2;

%%  

% Compile inzone data
incell_6={3 4; 3 5; 3 6;  4 4; 4 5; 4 6}; 
incell_4={3 5; 3 6;  4 5; 4 6}; 
incell=  incell_6; disp('Selecting inzone for 6 cells!'); 
% incell=  incell_4; disp('Selecting inzone for 6 cells!'); 



% In zones for Acc-Rej
% incell={
% % 6 1; 6 2; 
% 5 2; 5 3; 5 4; 
% 4 3; 4 4; 4 5;
% };   disp('Selecting for Acc-Rej tradeoff')
n_inexplore=nan(log.n_subjs,3); 
d_ind =cell(log.n_subjs,3); 
for s=1:log.n_subjs
    for t=1:2
        wt.d=  d_trialstats{s,t+1}; 
        
        wt.ind=[]; 
        for ic= 1:size(incell,1)
%             disp('Deleting exploration!') 
%             wt.d(wt.d(:, col.EnvThreatOriginal).*6 == incell{ic,1}  & wt.d(:, col.NTokens)./2  & wt.d(:, col.Choice)==3) =nan;
            
            
            wi.d= wt.d(wt.d(:, col.EnvThreatOriginal).*6 == incell{ic,1}  & wt.d(:, col.NTokens)./2 == incell{ic,2}, :);
            
            
            wi.d(:, col.pAcc) = mean(wi.d(:, col.Choice)==1 );
            wi.d(:, col.pRej) = mean(wi.d(:, col.Choice)==2 );
            wi.d(:, col.pExp) = mean(wi.d(:, col.Choice)==3 );
            
            
            
            
            wt.ind= [wt.ind;  wi.d];
        end
        n_inexplore(s,3)=sum(wt.ind(:, col.Choice)==3); 
        
        d_ind{s,t}=  wt.ind; 
        wt=[]; 
    end    
    d_ind{s,3}= sortrows([d_ind{s,1};d_ind{s,2}], col.Trialnum); 
end


disp(n_inexplore)
error

%%

d_nchoice= repmat({nan(log.n_subjs,2)}, 1,3);  
dc_pLosspRej=nan(log.n_subjs,3); 
for t=1:3
    for s=1:log.n_subjs
        
        % n Reject, n Explore 
%         d_nchoice{t}(s,1)=  mean(d_ind{s,t}(:, col.Choice)==2); 
%         d_nchoice{t}(s,2)=  mean(d_ind{s,t}(:, col.Choice)==3); 
        
        
        % Acc Rej
        d_nchoice{t}(s,1)=  mean(d_ind{s,t}(:, col.Choice)==1);  if t==1 && s==1; disp('Acc Reject tradeofF!!! '); end
        d_nchoice{t}(s,2)=  mean(d_ind{s,t}(:, col.Choice)==2); 
        
        
        
        % Reject correlated w pLoss? 
        [r p]= corr( d_ind{s,t}(:, col.pRej) , d_ind{s,t}(:, col.pLoss)); 
%         dc_pLosspRej(s,t) =p ;   if t==1 && s==1; input('corr pLos w pRej: recording p'); end
        dc_pLosspRej(s,t) =r ;  if t==1 && s==1; input('corr pLos w pRej: recording r'); end
    end
end

% pReject correlated w pLoss? 
figure('color','w', 'name','Choice corr w vars'), k=1; 
[h p ci st] = ttest( dc_pLosspRej );
% subplot(1,3, k); k=k+1;
barwitherr((std(d_nchoice{t})./sqrt(log.n_subjs))', mean(d_nchoice{t})')
set(gca,'xticklabel', conds), ylabel('Mean r'), title('Reject correlated w pLoss?')
disp('Reject correlated w pLoss? --------------------- '), p

k=1; figure('color','w', 'name','N trials/choice')
for t=1:3 % Compare n choice 
    [h p ci st] = ttest( d_nchoice{t}(:,1) - d_nchoice{t}(:,2) ); 
    subplot(1,4, k); k=k+1; 
    barwitherr((std(d_nchoice{t})./sqrt(log.n_subjs))', mean(d_nchoice{t})')
    set(gca,'xticklabel', {'Reject','Explore'}) 
    title(['[' conds{t} ']  t(' num2str(st.df) ')=' num2str(st.tstat) ', p='  num2str(p)])
end
d_nchoice{4} = d_nchoice{1} - d_nchoice{2}; 
[h p ci st] = ttest(d_nchoice{4}); 
subplot(1,4, k); k=k+1; 
barwitherr((std(d_nchoice{4})./sqrt(log.n_subjs))', mean(d_nchoice{4})')
set(gca,'xticklabel', {'Reject','Explore'}), title('cF > ct')
disp('In-zone, cF>ct difference in % Rej & Explore?  --------------------- '), p

