% Setup data for modelling
clear all; close all hidden; clc

% where.root= 'C:\Users\e.loh\Dropbox\SCRIPPS\';
where.root='/Users/EleanorL/Dropbox/SCRIPPS/';
where.data=[where.root '1 Explore fMRI' filesep '1 Behavioural data']; where.model_folder=[where.root '2 Explore experiment' filesep '3 Analysis' filesep '4 Fit computational models'];

req.subsample_thirds=3;    % 1, 2, 3; otherwise ignore
req.subsample_1st2nd=0;   % 0 = No subsampling, 1=First half, 2=Second half 

%% Load subject data + stick into column format for modelling 
%   Variable 'subjdata': Col 1=Subject, Col 2=Conflict task, Col  3= Control task, Col 4= Data together (cF & ct)

for o1=1:1 % Columns 
    
    % Columns from fMRI data
    col.f.NTokPairs=2;
    col.f.EnvThreat=3;
    col.f.Task=5;
    col.f.Trialnum=7;
    col.f.Choice=8;
    col.f.Choice2=10;
    col.f.OutcomeMag=15;
    col.f.TrialValid=13;
    
    % Columns for modelling data
    
    col.Choice=3;
    col.EnvThreat=6;
    col.NTokens=2;
    col.pLoss=1;
    col.Entropy=4;
    col.VExplore=5;
    col.EV=10;
    col.OutcomeMagnitude=11;
    col.OutcomeMean=12;
    col.OutcomeVariance=13;
    col.Task=9;
    col.Trialnum=8;
    col.EntropyNTok=7;
    col.EntropyEV=14;
%     col.Conflict=15;
    col.StanDev=16;
    col.vMeanVar=17;
    col.BinomVar=18; 
    % 
    details.col=col;
end

% Load data
cd(where.data); w.dir=dir('p*'); details.subjects=cellstr(char(w.dir(:).name)); details.n_subjs=length(details.subjects); cd(where.model_folder)
subjdata=cell(details.n_subjs, 4);
for s=1:details.n_subjs
    
    ws=load([where.data filesep details.subjects{s} filesep details.subjects{s} '_file_taskfMRI.mat']);
    ws.ts=ws.alldata(ws.alldata(:,col.f.TrialValid)==1, :);
    ws.d=nan*zeros(size(ws.ts,1), 7);
    ws.d(:,col.NTokens)=ws.ts(:, col.f.NTokPairs)*2;
    ws.d(:,col.EnvThreat)=ws.ts(:, col.f.EnvThreat)/6;
    ws.d(:,col.Choice)=ws.ts(:, col.f.Choice);
    ws.d(:,col.Trialnum)=  1:size(ws.ts,1);
    ws.d(:,col.Task)=ws.ts(:, col.f.Task);
    ws.d(:,col.OutcomeMagnitude)=ws.ts(:,col.f.OutcomeMag);
    
    % Write parameters
    ws.cf=ws.d(ws.d(:,col.Task)==1, :);
    ws.ct=ws.d(ws.d(:,col.Task)==2, :);
    [ ws.cf] = fpar_conflict( ws.cf, col);
    [ ws.ct] = fpar_control( ws.ct, col);
    
    % Subsample by split 
    switch req.subsample_1st2nd
        case 1
            ws.cf=sortrows(ws.cf, col.Trialnum); 
            ws.cf=  ws.cf(1:round(size(ws.cf,1)/2),:); 
            ws.ct=sortrows(ws.ct, col.Trialnum); 
            ws.ct=  ws.ct(1:round(size(ws.ct,1)/2),:); 
        case 2
            ws.cf=sortrows(ws.cf, -col.Trialnum); 
            ws.cf=  ws.cf(1:round(size(ws.cf,1)/2),:); 
            ws.ct=sortrows(ws.ct, -col.Trialnum); 
            ws.ct=  ws.ct(1:round(size(ws.ct,1)/2),:); 
            %
            ws.cf=sortrows(ws.cf, col.Trialnum); 
            ws.ct=sortrows(ws.ct, col.Trialnum); 
    end  
    
    
    
    
    switch req.subsample_thirds
        case 1
            ws.cf=sortrows(ws.cf, col.Trialnum); 
            ws.tn = 1: round(size(ws.cf,1)/3);
            ws.cf= ws.cf(ws.tn, :);
            ws.ct=sortrows(ws.ct, col.Trialnum); 
            ws.tn = 1: round(size(ws.ct,1)/3);
            ws.ct= ws.ct(ws.tn, :);
        case 2
            ws.cf=sortrows(ws.cf, col.Trialnum); 
            ws.tn = round(size(ws.cf,1)/3)+1:round(size(ws.cf,1)*2/3); 
            ws.cf= ws.cf(ws.tn, :);
            ws.ct=sortrows(ws.ct, col.Trialnum); 
            ws.tn = round(size(ws.ct,1)/3)+1:round(size(ws.ct,1)*2/3); 
            ws.ct= ws.ct(ws.tn, :);
        case 3
            ws.cf=sortrows(ws.cf, col.Trialnum); 
            ws.tn = round(size(ws.cf,1)*2/3)+1; 
            ws.cf= ws.cf(ws.tn:end, :);
            ws.ct=sortrows(ws.ct, col.Trialnum); 
            ws.tn = round(size(ws.ct,1)*2/3)+1; 
            ws.ct= ws.ct(ws.tn:end, :);
            
    end   
    
    subjdata{s, 1}=details.subjects{s};
    subjdata{s, 2}=ws.cf;
    subjdata{s, 3}=ws.ct;
    subjdata{s, 4}=sortrows([ws.cf; ws.ct], col.Trialnum);
    
    ws=[];
end


% SAVE 
save([where.model_folder filesep '2 Analysis inputs' filesep 'All data (' date ').mat'], 'subjdata', 'details')

%%

col.ChoiceReject=20; 
col.pReject=21; 

subs=1:10; 
f.plotcols=2; k=1;
for s=1:length(subs)

ws=subjdata{s,2};

ws(:, col.ChoiceReject)=ws(:, col.Choice)==2;
for e=1:6 % Mark p(Reject)
    for n=1:6
        wc=ws(ws(:, col.EnvThreat)*6==e & ws(:, col.NTokens)/2==n, :);
        ws(ws(:, col.EnvThreat)*6==e & ws(:, col.NTokens)/2==n, col.pReject)=mean(wc(:, col.ChoiceReject));
    end
end


ws=sortrows(ws, col.ChoiceReject);
% ws=sortrows(ws, col.Choice);

% x= ws(:, [col.EnvThreat col.NTokens col.Entropy col.EV]);
x= ws(:, [col.EnvThreat col.NTokens col.pLoss col.Entropy col.EV]);
% y=ws(:, col.Choice);
y=ws(:, col.ChoiceReject);
[b bint residuals1 rint stats1]=regress( y,  [ones(size(ws,1),1) x]);   


d_resid(s,1)=mean(residuals1( ws(:, col.ChoiceReject) ==1));
d_resid(s,2) = mean(residuals1( ws(:, col.ChoiceReject) ==0));
[h p ci stat]=ttest(residuals1( ws(:, col.ChoiceReject) ==1));
d_resid(s,3)=stat.tstat;
[h p ci stat]=ttest(residuals1( ws(:, col.ChoiceReject) ==0));
d_resid(s,4)=stat.tstat;

% [b bint residuals rint stats1]=regress( residuals,  [ones(size(ws,1),1) ws(:, col.pLoss)  ] ); stats
subplot(ceil(length(subs)/f.plotcols),  f.plotcols, k); k=k+1; 
plot(1:size(ws,1), y), hold on
plot(1:size(ws,1), residuals1,'r')
[r p]=corr(y, residuals1);

end


[h p]=ttest( d_resid(:,1) - d_resid(:,2));

figure('color','w')
barwitherr(std(d_resid(:, 3:4))/sqrt(details.n_subjs), mean(d_resid(:, 3:4)))














