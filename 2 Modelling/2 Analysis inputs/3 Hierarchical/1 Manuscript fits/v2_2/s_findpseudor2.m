clear all; clc

%% file from MRI data folder
for o1=1:1
% load('/Users/EleanorL/Dropbox/WorkPC/Conflict beta/Model values (15-Sep-2014).mat')
% res=d_fits.cf.r_res;
% 
% % ------
% b1=res{1,3};
% b2=res{2,3};
% b1=a{1}; b2=a{2};
% B=(b1-b2)*-0.5; B
end
% 
% % Formula: Page 22 from here:  http://www.cns.nyu.edu/~daw/d10.pdf
% %     If R is the log data likelihood under chance (e.g., for 100 trials of a two-choice task, 100 · log(.5)) and L is the log
% %     likelihood under the fit model, then pseudo-r2 is 1 ? L/R.
% 
% 
% clear all; clc; cd('/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/3 Hierarchical/Used in manuscirpt/v2_2/')
% load('res_hierarfitmodels_cF (21-Jul-2014) top10.mat')
% % load('res_hierarfitmodels_ct (21-Jul-2014) top10.mat')
% 
% modres=r_res{1,2};
% for s=1:20
%     
%     
%     R=log(1/3)* details.subj_ntrials(s);
%     L=-modres(s,2);
%     psur(s)=1-L/R;
% end
% mean(psur)


%% Mark BICs = 2*nll+ K*ln(n_trials)
% Difference in BICs, comparing tasks?

cf=load('res_hierarfitmodels_cF (21-Jul-2014) top10.mat');
ct=load('res_hierarfitmodels_ct (21-Jul-2014) top10.mat');

% Which model
mcf=1; mct=1;



k_cf=size(cf.r_res{mcf, 2},2)-3;  k_ct=size(ct.r_res{mct, 2},2)-3;
for s=1:20
    % cF
    bcf(s)=2*cf.r_res{mcf, 2}(s,2) + k_cf*log(cf.details.subj_ntrials(s));
    
    % ct
    bct(s)=2*ct.r_res{mct, 2}(s,2) + k_ct*log(ct.details.subj_ntrials(s));
end

[h p]=ttest(bcf, bct)

