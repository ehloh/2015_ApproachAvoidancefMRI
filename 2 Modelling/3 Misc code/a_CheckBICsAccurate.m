clear all; clc

where='/Users/EleanorL/Dropbox/sandisk/4 Explore experiment/3 Analysis/4 Fit computational models';
r1=load([where '/2 Analysis inputs/Det Jpower/res_fitmodels_cF (02-Dec-2014) all fams.mat']);
r=load([where '/2 Analysis inputs/Det Jpower/res_fitmodels_cF (18-Jul-2014) all20iter.mat']);
d=load([where '/2 Analysis inputs/All data (09-May-2014).mat']);
path(pathdef);  addpath(genpath([where filesep '1 Value functions Det Jpower'])); addpath(where);

m=58;  task=1;
% modpars={'b';'p';'j';'m';'f';'w'};
% modpars={'b';'p';'m';'i'; 'f'; 'e';'w'};
modpars={'b'; 'f';'e'};


%%

subpars=r.r_res{m,2}(:,4:end); mod=r.r_res{m,1};
allnll=0; ntrials=0; cumbic=0;
whichsubs=1:20; % whichsubs=[1:11 13:20];


for s=whichsubs
    [ transparval ] = f_transpar(modpars, subpars(s,:), 'from');
%     [nll pch]=f_nllsoftmax_lapse(transparval , {mod d.subjdata{s,task+1} r.details.fixedpar d.details.col});
    [nll pch]=f_nllsoftmax(transparval , {mod d.subjdata{s,task+1} r.details.fixedpar d.details.col});
    

%     ws.d=d.subjdata{s,task+1};
%     ws.VA=
    
    
    
    % Record details
    allnll=allnll+nll;
    ntrials=ntrials+size(d.subjdata{s,task+1},1);
    cumbic=cumbic+2*nll+ length(transparval ) * log(  size(d.subjdata{s,task+1},1)  );
end
cumbic

% Calculate BIC: 2*nll+ K*ln(n_trials)

%% Any consistency on which is the winining model?

r.r_res=[r.r_res; r1.r_res]; modfamsize=16;
r.r_res=sortrows(r.r_res,1);

r_bic=nan(size(r.r_res,1), 20);  % model, sub
for m=1:size(r.r_res,1)
    for s=1:20
       r_bic(m,s) =r.r_res{m,2}(s,1);
    end
end

% r_bic=r_bic./repmat(min(r_bic), size(r_bic,1),1); %  Subjectwise scaling. 1=best, > = worse
% for s=1:20; r_bic(s,:)=zscore(r_bic(s,:)); end  % Subject wise scaling (z score)

r_bic=r_bic==repmat(min(r_bic), size(r_bic,1),1); %  Subjectwise find best
imagesc(-r_bic); colormap(hot) %colorbar
allmods=r.r_res(:,1);


for f=1:length(allmods)/modfamsize-1 % Demarcate model families
hold on; plot(1:20, f*modfamsize);
end

% ytick
% (1:modfamsize:length(allmods)) + modfamsize/2
allmods(1:modfamsize:length(allmods))
