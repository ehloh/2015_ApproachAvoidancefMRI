% Extract behavioural stats and parameters from behavioural model fitting
clear all; close all hidden; clc

where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data']; where.parameter_scripts='/Volumes/PENNYDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
% where.where='I:\5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data'];  where.parameter_scripts='I:\4 Explore experiment\3 Analysis\4 Fit computational models';

%
log.specificsubjects={};
log.modelpar_cF='res_fitmodels_conflict (23-Aug-2013)';  % Results from model fitting
log.modelpar_ct='res_fitmodels_control (23-Aug-2013)';

request.n_behmods=3;

for o1=1:1 % General setup  
    
    addpath(where.where); addpath(where.parameter_scripts)
    
    % Load subjects
    log.w=load([where.data_beh filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Behaviour folders
    where.subfol=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        where.subfol{s}=[where.data_beh filesep log.subjects{s} filesep];
    end
    
end

%% (1) Load all data (variable dd.)

disp('Loading data ###########################')
for o1=1:1 % Data specifications (columns) 
    
    col.t.Task=5; % Event classifiers
    col.t.Resp1=8; 
    col.t.Resp2=10;
    col.t.TrialValid=13; 
    col.t.OutcomePresented=14;
    col.t.OutcomeMagnitude=15;
    col.t.RT1=9;
    col.t.Trialnum=30;
    
    col.t.EnvThreat=3; % Design columns
    col.t.NTokens=2;
    col.t.pLoss=31;
    col.t.Entropy=32;
    col.t.VExplore=33;
    col.t.EV=34;
    col.t.OutcomeMean=35;
    col.t.OutcomeVariance=36;
    col.t.TrialType=1;
    
    col.t.Onset_Offer=22; % Event onsets
    col.t.Onset_Info=25;
    col.t.Onset_Outcome=28;
    col.t.Onset_TrialEnd=29;
    col.t.Onset_Motor1=23;
    col.t.Onset_Motor2=26;

    % Model fitting results
    col.m.modfullname=1;
    col.m.modname=13;
    col.m.mod_singlepars=5;
    col.m.subpar=2;
    col.m.subpar_startpar=4; % within subpar array
end

% Output measures
d.subjects=log.subjects;
d.o.counterhands=zeros(log.n_subjs,1);
d.o.countertask=zeros(log.n_subjs,1);

% (1) Write all task data to structure d
dd.rawMRI=cell(log.n_subjs,3); % Col 1= trialstats, Col 2=stats, Col 3= settings 
dd.task=cell(log.n_subjs,1); dd.task_cF=cell(log.n_subjs,1); dd.task_ct=cell(log.n_subjs,1);
dd.tasklearn=cell(log.n_subjs,3);
for s=1:log.n_subjs
    
    % Load & mark parameters for fMRI-task trials
    for b=1:6
        ws.b{b}=load([where.subfol{s} log.subjects{s} '_file_6integratedfMRI_b' num2str(b) '.mat']);
    end
    dd.rawMRI{s,2}=load([where.subfol{s} log.subjects{s} '_file_taskfMRI.mat']);
    dd.rawMRI{s,3}=ws.b{1}.integratedfMRI.settings;
    dd.rawMRI{s,1}=vertcat(ws.b{1}.integratedfMRI.rdata,ws.b{2}.integratedfMRI.rdata,ws.b{3}.integratedfMRI.rdata,ws.b{4}.integratedfMRI.rdata,ws.b{5}.integratedfMRI.rdata,ws.b{6}.integratedfMRI.rdata);
    dd.rawMRI{s,1}(:,col.t.Trialnum)=1:size(dd.rawMRI{s,1},1);
    [dd.task_cF{s}] = fpar_conflict(dd.rawMRI{s,1}(dd.rawMRI{s,1}(:,col.t.Task)==1,:), col.t);
    [dd.task_ct{s}] = fpar_control(dd.rawMRI{s,1}(dd.rawMRI{s,1}(:,col.t.Task)==2,:), col.t);
    dd.task{s}=sortrows(vertcat(dd.task_cF{s},dd.task_ct{s}), col.t.Trialnum);
    
    % Settings
    d.o.counterhands(s)=dd.rawMRI{s,2}.settings.counterbal_hands;
    d.o.countertask(s)=dd.rawMRI{s,2}.settings.whichtaskfirst;
%     d.o.countertask(s)=dd.rawMRI{s,2}.settings.counterbal_startask;
    
    % Learning stage
    dd.tasklearn{s}=load([where.subfol{s} log.subjects{s} '_file_1learnenv.mat']);
    dd.tasklearn{s}=dd.tasklearn{s}.learnenv;
    
    ws=[];
end

% (2) Load results of behavioural model fitting 
dd.modelpar_cF=load([where.parameter_scripts filesep '2 Analysis inputs' filesep log.modelpar_cF '.mat']);
dd.modelpar_ct=load([where.parameter_scripts filesep '2 Analysis inputs' filesep log.modelpar_ct '.mat']);
for m=1:request.n_behmods % Short names 
    dd.modelpar_cF.r_res{m,col.m.modname}=['mcF_' dd.modelpar_cF.r_res{m,col.m.modfullname}(4) dd.modelpar_cF.r_res{m,col.m.modfullname}(find(dd.modelpar_cF.r_res{m,col.m.modfullname}=='_',1, 'last')+1:end)];
    dd.modelpar_ct.r_res{m,col.m.modname}=['mct_' dd.modelpar_ct.r_res{m,col.m.modfullname}(4) dd.modelpar_ct.r_res{m,col.m.modfullname}(find(dd.modelpar_ct.r_res{m,col.m.modfullname}=='_',1, 'last')+1:end)];
end

% (3) Load personality scores
% [dd.w dd.ww dd.personality]=xlsread([where.data_beh filesep 'Group behaviour files' filesep 'Personality scores.xlsx']);
% for p=2:size(dd.personality,2)
%     eval(['d.p.' dd.personality{1,p} '(1:log.n_subjs,1)=cell2mat(dd.personality(2:end,p));'])
% end

%% (2) Extract behaviour stats (everything in variable d)

disp('Extracting behaviour ###########################')

% Define inner cluster
log.InnerCluster4=[3 5; 3 6; 4 5; 4 6];
log.InnerCluster6=[3 4; 3 5; 3 6; 4 4; 4 5; 4 6];
if isempty(dd.task{1}(:,col.t.NTokens)==12)==0;  log.InnerCluster4(:,2)=log.InnerCluster4(:,2)*2; log.InnerCluster6(:,2)=log.InnerCluster6(:,2)*2; end

% (1) Percentage choice (Overall) ########################
d.per.Accept=zeros(log.n_subjs,1); 
d.per.Reject=zeros(log.n_subjs,1); 
d.per.Explore=zeros(log.n_subjs,1); 
d.per.cF_Accept=zeros(log.n_subjs,1); d.per.cF_Reject=zeros(log.n_subjs,1);  d.per.cF_Explore=zeros(log.n_subjs,1);
d.per.ct_NoBomb=zeros(log.n_subjs,1); d.per.ct_Bomb=zeros(log.n_subjs,1);  d.per.ct_Explore=zeros(log.n_subjs,1);
for s=1:log.n_subjs
    
    
    d.per.Accept(s) = (sum(dd.task_cF{s}(:,col.t.Resp1)==1) + sum(dd.task_ct{s}(:,col.t.Resp1)==1)) /(size(dd.task_cF{s},1) + size(dd.task_ct{s},1));
    d.per.Reject(s) = (sum(dd.task_cF{s}(:,col.t.Resp1)==2) + sum(dd.task_ct{s}(:,col.t.Resp1)==2) )/(size(dd.task_cF{s},1) + size(dd.task_ct{s},1));
    d.per.Explore(s) =( sum(dd.task_cF{s}(:,col.t.Resp1)==3) + sum(dd.task_ct{s}(:,col.t.Resp1)==3)) /(size(dd.task_cF{s},1) + size(dd.task_ct{s},1));
    d.per.cF_Accept(s)=sum(dd.task_cF{s}(:,col.t.Resp1)==1)/size(dd.task_cF{s},1);
    d.per.cF_Reject(s)=sum(dd.task_cF{s}(:,col.t.Resp1)==2)/size(dd.task_cF{s},1);
    d.per.cF_Explore(s)=sum(dd.task_cF{s}(:,col.t.Resp1)==3)/size(dd.task_cF{s},1);
    d.per.ct_NoBomb(s)=sum(dd.task_ct{s}(:,col.t.Resp1)==1)/size(dd.task_ct{s},1);
    d.per.ct_Bomb(s)=sum(dd.task_ct{s}(:,col.t.Resp1)==2)/size(dd.task_ct{s},1);
    d.per.ct_Explore(s)=sum(dd.task_ct{s}(:,col.t.Resp1)==3)/size(dd.task_ct{s},1);
end

% (2) General performance indicators ########################
d.o.win_all=zeros(log.n_subjs,1); % Winnings
d.o.win_cF=zeros(log.n_subjs,1); d.o.win_ct=zeros(log.n_subjs,1);
for s=1:log.n_subjs
    d.o.win_all(s)=sum(dd.task{s}(:, col.t.OutcomeMagnitude));
    d.o.win_cF(s)=sum(dd.task_cF{s}(:, col.t.OutcomeMagnitude));
    d.o.win_ct(s)=sum(dd.task_ct{s}(:, col.t.OutcomeMagnitude));
end
d.o.learnenv=zeros(log.n_subjs,1); % Initial learning accuracy
for s=1:log.n_subjs
    d.o.learnenv(s)=dd.tasklearn{s}.outcomes.learningaccuracy;
end

% (3) Cluster ########################
d.per_i4.Reject=zeros(log.n_subjs,1); d.per_i4.Explore=zeros(log.n_subjs,1); 
d.per_i4.cF_Reject=zeros(log.n_subjs,1);  d.per_i4.cF_Explore=zeros(log.n_subjs,1);
d.per_i4.ct_Bomb=zeros(log.n_subjs,1);  d.per_i4.ct_Explore=zeros(log.n_subjs,1);
d.per_i6.Reject=zeros(log.n_subjs,1);  d.per_i6.Explore=zeros(log.n_subjs,1); 
d.per_i6.cF_Reject=zeros(log.n_subjs,1);  d.per_i6.cF_Explore=zeros(log.n_subjs,1);
d.per_i6.ct_Bomb=zeros(log.n_subjs,1);  d.per_i6.ct_Explore=zeros(log.n_subjs,1);
%
d.rt.all=zeros(log.n_subjs,1); d.rt.cF=zeros(log.n_subjs,1); d.rt.ct=zeros(log.n_subjs,1); d.rt.cF_Accept=zeros(log.n_subjs,1); d.rt.cF_Reject=zeros(log.n_subjs,1);d.rt.cF_Explore=zeros(log.n_subjs,1); d.rt.ct_NoBomb=zeros(log.n_subjs,1); d.rt.ct_Bomb=zeros(log.n_subjs,1); d.rt.ct_Explore=zeros(log.n_subjs,1);
d.rt_i4.all=zeros(log.n_subjs,1); d.rt_i4.cF=zeros(log.n_subjs,1); d.rt_i4.ct=zeros(log.n_subjs,1); d.rt_i4.cF_Accept=zeros(log.n_subjs,1);  d.rt_i4.cF_Reject=zeros(log.n_subjs,1);d.rt_i4.cF_Explore=zeros(log.n_subjs,1); d.rt_i4.ct_NoBomb=zeros(log.n_subjs,1); d.rt_i4.ct_Bomb=zeros(log.n_subjs,1); d.rt_i4.ct_Explore=zeros(log.n_subjs,1);
d.rt_i6.all=zeros(log.n_subjs,1); d.rt_i6.cF=zeros(log.n_subjs,1); d.rt_i6.ct=zeros(log.n_subjs,1); d.rt_i6.cF_Accept=zeros(log.n_subjs,1);  d.rt_i6.cF_Reject=zeros(log.n_subjs,1);d.rt_i6.cF_Explore=zeros(log.n_subjs,1); d.rt_i6.ct_NoBomb=zeros(log.n_subjs,1); d.rt_i6.ct_Bomb=zeros(log.n_subjs,1); d.rt_i6.ct_Explore=zeros(log.n_subjs,1);
for s=1:log.n_subjs
    
    % All data
    ws.d_cF=dd.task_cF{s}; ws.d_ct=dd.task_ct{s};
    d.rt.all(s)=mean(dd.task{s}(:,col.t.RT1));
    d.rt.cF(s)=mean(ws.d_cF(:,col.t.RT1));
    d.rt.ct(s)=mean(ws.d_ct(:,col.t.RT1));
    d.rt.cF_Accept(s)=mean(ws.d_cF(ws.d_cF(:,col.t.Resp1)==1, col.t.RT1));
    d.rt.cF_Reject(s)=mean(ws.d_cF(ws.d_cF(:,col.t.Resp1)==2, col.t.RT1));
    d.rt.cF_Explore(s)=mean(ws.d_cF(ws.d_cF(:,col.t.Resp1)==3, col.t.RT1));
    d.rt.ct_NoBomb(s)=mean(ws.d_ct(ws.d_ct(:,col.t.Resp1)==1, col.t.RT1));
    d.rt.ct_Bomb(s)=mean(ws.d_ct(ws.d_ct(:,col.t.Resp1)==2, col.t.RT1));
    d.rt.ct_Explore(s)=mean(ws.d_ct(ws.d_ct(:,col.t.Resp1)==3, col.t.RT1));
    
    % Only exploration cluster (4 cells)
    ws.i4_cF=[]; ws.i4_ct=[];
    for i=1:size(log.InnerCluster4,1)
        ws.i4_cF=vertcat(ws.i4_cF, ws.d_cF(ws.d_cF(:,col.t.EnvThreat)==log.InnerCluster4(i,1) &  ws.d_cF(:,col.t.NTokens)==log.InnerCluster4(i,2), :));
        ws.i4_ct=vertcat(ws.i4_ct, ws.d_ct(ws.d_ct(:,col.t.EnvThreat)==log.InnerCluster4(i,1) &  ws.d_ct(:,col.t.NTokens)==log.InnerCluster4(i,2), :));
    end
    d.per_i4.cF_Reject(s)=sum(ws.i4_cF(:,col.t.Resp1)==2)/size(ws.i4_cF,1);     % Choice
    d.per_i4.cF_Explore(s)=sum(ws.i4_cF(:,col.t.Resp1)==3)/size(ws.i4_cF,1);
    d.per_i4.ct_Bomb(s)=sum(ws.i4_ct(:,col.t.Resp1)==2)/size(ws.i4_ct,1);
    d.per_i4.ct_Explore(s)=sum(ws.i4_ct(:,col.t.Resp1)==3)/size(ws.i4_ct,1); 
    d.per_i4.Reject=(sum(ws.i4_cF(:,col.t.Resp1)==2)+sum(ws.i4_ct(:,col.t.Resp1)==2))     /    (size(ws.i4_cF,1) + size(ws.i4_ct,1));
    d.per_i4.Explore=(sum(ws.i4_cF(:,col.t.Resp1)==3)+sum(ws.i4_ct(:,col.t.Resp1)==3))     /    (size(ws.i4_cF,1) + size(ws.i4_ct,1));
    d.rt_i4.all=mean(vertcat(ws.i4_cF(:,col.t.RT1),ws.i4_ct(:,col.t.RT1)));             % RT
    d.rt_i4.cF=mean(ws.i4_cF(:,col.t.RT1));
    d.rt_i4.ct=mean(ws.i4_ct(:,col.t.RT1));
    d.rt_i4.cF_Accept(s)=mean(ws.i4_cF(ws.i4_cF(:,col.t.Resp1)==1, col.t.RT1));
    d.rt_i4.cF_Reject(s)=mean(ws.i4_cF(ws.i4_cF(:,col.t.Resp1)==2, col.t.RT1));
    d.rt_i4.cF_Explore(s)=mean(ws.i4_cF(ws.i4_cF(:,col.t.Resp1)==3, col.t.RT1));
    d.rt_i4.ct_NoBomb(s)=mean(ws.i4_ct(ws.i4_ct(:,col.t.Resp1)==1, col.t.RT1));
    d.rt_i4.ct_Bomb(s)=mean(ws.i4_ct(ws.i4_ct(:,col.t.Resp1)==2, col.t.RT1));
    d.rt_i4.ct_Explore(s)=mean(ws.i4_ct(ws.i4_ct(:,col.t.Resp1)==3, col.t.RT1));
    
        
    % Only exploration cluster (6 cells)
    ws.i6_cF=[]; ws.i6_ct=[];
    for i=1:size(log.InnerCluster6,1)
        ws.i6_cF=vertcat(ws.i6_cF, ws.d_cF(ws.d_cF(:,col.t.EnvThreat)==log.InnerCluster6(i,1) &  ws.d_cF(:,col.t.NTokens)==log.InnerCluster6(i,2), :));
        ws.i6_ct=vertcat(ws.i6_ct, ws.d_ct(ws.d_ct(:,col.t.EnvThreat)==log.InnerCluster6(i,1) &  ws.d_ct(:,col.t.NTokens)==log.InnerCluster6(i,2), :));
    end    
    d.per_i6.cF_Reject(s)=sum(ws.i6_cF(:,col.t.Resp1)==2)/size(ws.i6_cF,1);     % Choice
    d.per_i6.cF_Explore(s)=sum(ws.i6_cF(:,col.t.Resp1)==3)/size(ws.i6_cF,1);
    d.per_i6.ct_Bomb(s)=sum(ws.i6_ct(:,col.t.Resp1)==2)/size(ws.i6_ct,1);
    d.per_i6.ct_Explore(s)=sum(ws.i6_ct(:,col.t.Resp1)==3)/size(ws.i6_ct,1); 
    d.per_i6.Reject=(sum(ws.i6_cF(:,col.t.Resp1)==2)+sum(ws.i6_ct(:,col.t.Resp1)==2))     /    (size(ws.i6_cF,1) + size(ws.i6_ct,1));
    d.per_i6.Explore=(sum(ws.i6_cF(:,col.t.Resp1)==3)+sum(ws.i6_ct(:,col.t.Resp1)==3))     /    (size(ws.i6_cF,1) + size(ws.i6_ct,1));
    d.rt_i6.all=mean(vertcat(ws.i6_cF(:,col.t.RT1),ws.i6_ct(:,col.t.RT1)));              % RT
    d.rt_i6.cF=mean(ws.i6_cF(:,col.t.RT1));
    d.rt_i6.ct=mean(ws.i6_ct(:,col.t.RT1));
    d.rt_i6.cF_Accept(s)=mean(ws.i6_cF(ws.i6_cF(:,col.t.Resp1)==1, col.t.RT1));
    d.rt_i6.cF_Reject(s)=mean(ws.i6_cF(ws.i6_cF(:,col.t.Resp1)==2, col.t.RT1));
    d.rt_i6.cF_Explore(s)=mean(ws.i6_cF(ws.i6_cF(:,col.t.Resp1)==3, col.t.RT1));
    d.rt_i6.ct_NoBomb(s)=mean(ws.i6_ct(ws.i6_ct(:,col.t.Resp1)==1, col.t.RT1));
    d.rt_i6.ct_Bomb(s)=mean(ws.i6_ct(ws.i6_ct(:,col.t.Resp1)==2, col.t.RT1));
    d.rt_i6.ct_Explore(s)=mean(ws.i6_ct(ws.i6_ct(:,col.t.Resp1)==3, col.t.RT1));
    
    
    ws=[];
end

% (4) Parameters from model fitting ########################
for m=1:request.n_behmods % Conflict task parameters 
    wm.subpar_cF=sortrows([dd.modelpar_cF.details.subjects num2cell(dd.modelpar_cF.r_res{m,2}(:,col.m.subpar_startpar:end))],1);
    for p=1:size(wm.subpar_cF,2)-1
        eval(['d.' dd.modelpar_cF.r_res{m, col.m.modname} '.' dd.modelpar_cF.r_res{m, col.m.mod_singlepars}{p} '=cell2mat(wm.subpar_cF(:,1+p));']);
    end
    wm=[];
end
for m=1:request.n_behmods % Control task parameters 
    wm.subpar_ct=sortrows([dd.modelpar_ct.details.subjects num2cell(dd.modelpar_ct.r_res{m,2}(:,col.m.subpar_startpar:end))],1);
    for p=1:size(wm.subpar_ct,2)-1
        eval(['d.' dd.modelpar_ct.r_res{m, col.m.modname} '.' dd.modelpar_ct.r_res{m, col.m.mod_singlepars}{p} '=cell2mat(wm.subpar_ct(:,1+p));']);
    end
    wm=[];
end
    
%% (3) Write data to table

% Read titles
disp('Reading titles - recursions in structure ''d'' ----------------')
titles=cell(50,1); k=1; read.lev1=fieldnames(d);
for i1=1: length(read.lev1)   % 1st recursion
    wi1=[];
    eval(['wi1.c=d.' read.lev1{i1} ';']);
    wi1.name=read.lev1{i1};
    if isstruct(wi1.c)
        disp('1st recursion')
        read.lev2=fieldnames(wi1.c);
        for i2=1:length(read.lev2)   % 2nd recursion
            wi2=[];
            eval(['wi2.c=wi1.c.' read.lev2{i2} ';']);
            wi2.name=[wi1.name '.' read.lev2{i2}];
            if isstruct(wi2.c)
                disp('2nd recursion')
                read.lev3=fieldnames(wi2.c);
                for i3=1:length(read.lev3)   % 3rd recursion
                    wi3=[];
                    eval(['wi3.c=wi2.c.' read.lev3{i3} ';']);
                    wi3.name=[wi2.name '.' read.lev3{i3}];
                    if isstruct(wi3.c)
                        disp('3rd recursion')
                        read.lev4=fieldnames(wi3.c);
                        for i4=1:length(read.lev4)   % 4th recursion
                            wi4=[];
                            eval(['wi4.c=wi3.c.' read.lev4{i4} ';']);
                            wi4.name=[wi3.name '.' read.lev4{i4}];
                            if isstruct(wi4.c)
                                disp('4rd recursion')
                                read.lev5=fieldnames(wi4.c);
                                for i5=1:length(read.lev5)   % 5th recursion
                                    wi5=[];
                                    eval(['wi5.c=wi4.c.' read.lev5{i5} ';']);
                                    wi5.name=[wi4.name '.' read.lev5{i5}];
                                    if isstruct(wi5.c)
                                        error('Too many recursions! Either set up less, or adapt script for 1 more')
                                    else % end at 5th recursion
                                        titles{k}=wi5.name; k=k+1;
                                    end
                                end
                            else % end at 4th recursion
                                titles{k}=wi4.name; k=k+1;
                            end
                        end
                    else % end at 3rd recursion
                        titles{k}=wi3.name; k=k+1; 
                    end
                end
            else % end of 2nd recursion
                titles{k}=wi2.name; k=k+1; 
            end
        end
    else % end of 1st recursion
        titles{k}=wi1.name; k=k+1;
    end
end
disp(titles)

% Read data, on the basis of titles
dat=cell(log.n_subjs, size(titles,1));
for t=1:size(titles,1)
    eval(['wt=d.' titles{t} ';'])
    if iscell(wt)
        eval(['dat(:,t)=d.' titles{t} ';']); % Char
    else
        eval(['dat(:,t)=num2cell(d.' titles{t} ');']) % numerical
    end
end

% Save data
dat=vertcat(titles',dat);
% save([where.data_beh filesep 'Group behaviour files' filesep 'Behavioural profile (' date ').mat'], 'dat')
for i=1:size(dat,1)  % Correct nans for printing to text documents 
    for j=1:size(dat,2)
        if isnan(dat{i,j})==1 &  ischar({i,j})==0
            dat{i,j}=[];
        end
    end
end
print2txt([where.data_beh filesep 'Group behaviour files'], ['Behavioural stats modelpars (' date ')'], dat);
