% s7_Firstlevel_checkcorrelatedconditions
% clear all; close all hidden; clc; 
clear all; clc

where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\Brain data';
log.w=load([where.data_brain filesep 'datalog_test2.mat']); log.datalog=log.w.datalog;

% Requested analysis
log.specificsubjects={'p02_YY'}; % BLANK to process all subjects

% % Model details
log.model='m1_AllPmods';
log.pmods={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'OutcomeMean'; 'OutcomeVariance'};
% log.model='m2_Basic';
% log.pmods={'EnvThreat'; 'NTokens'; 'Entropy'};
% log.model='m3_BasicOutcome';
% log.pmods={'EnvThreat'; 'NTokens'; 'Entropy'; 'OutcomeMean'; 'OutcomeVariance'};

for o1=1:1 % General settings and specifications
    
    % Subjects
    for o2=1:1 % Specific subjects requested?
        if isempty(log.specificsubjects)==0
            w.log=cell(length(log.specificsubjects)+1,size(log.datalog,2));
            for i=1:size(log.datalog,2) % Headers
                w.log{1,i}=log.datalog{1,i};
            end
            for i=1:length(log.specificsubjects)
                w.sok=0;
                while w.sok==0
                    for j=2:size(log.datalog,1)
                        if strcmp(log.datalog{j,1}, log.specificsubjects{i})==1
                            for ii=1:size(log.datalog,2)            
                                w.log{1+i,ii}=log.datalog{j,ii};
                            end
                            w.sok=1;
                        end
                        if w.sok==0 && j==size(log.datalog,1)
                            input('Error: Could not find requested specific subject');
                        end
                    end
                end
            end
            log.allsubs=log.datalog;
            log.datalog=w.log;
%             log.subjects=[];
        end
    end
    log.n_subjs=size(log.datalog,1)-1;
    for  i=1:log.n_subjs 
        log.subjects{i,1}=log.datalog{i+1,1};
    end
    errorlog=cell(1,1); e=1;
   
    % Interface
    disp('=======================================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location (brain): ' where.data_brain])
    disp(' ')
    disp(['First level model: ' log.model])
    disp('pMods included in this model:' )
    disp(log.pmods)
    disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
%     where.spm='D:\My Documents\MATLAB\spm8'; 
    
end

%%

% Load SPM.mat files (subjdata column 2)
subjdata=cell(log.n_subjs,2);
for s=1:log.n_subjs
    subjdata{s,1}=log.subjects{s};
    ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep];
    f=spm_select('List', ws.where, ['^SPM_' log.model '.*.mat$']);
    if isempty(f)==1
        f=spm_select('List', [ws.where 'Estimated ' log.model], ['^SPM_' log.model '.*.mat$']);
        if isempty(f)==1
            input(['ERROR: Could not find the SPM.mat file for the requested model  ('  log.model '   -    ' log.subjects{s}  ')']);
        else
            subjdata{s,2}=load([ws.where 'Estimated ' log.model filesep f]);
        end
    else
        subjdata{s,2}=load([ws.where f]);
    end
end

% Collect stick functions for regressors (also in subjdata column 2)
for s=1:log.n_subjs
    
    % Identify regressors of interest
    ws.names=subjdata{s,2}.SPM.xX.name'; r=1;
    for t=1:2
        for i=1:length(log.pmods)
            switch t
                case 1
                    ws.regs{r,1}=['cF_' log.pmods{i}];
                case 2
                    ws.regs{r,1}=['ct_' log.pmods{i}];
            end
            subjdata{s,2}.regdetails{r,1}=ws.regs{r,1};
            wr.regfind=strfind(ws.names, ['xp' ws.regs{r,1} '^1*bf(1)']);
            for ii=1:length(wr.regfind); if isempty(wr.regfind{ii})==0; wr.regnum(ii,1)=1; else wr.regnum(ii,1)=0;end; end
             subjdata{s,2}.regdetails{r,2}=find(wr.regnum==1);
            wr=[]; r=r+1;
        end
    end
    
    % Collect stick functions
    ws.regs=nan*zeros(size(subjdata{s,2}.SPM.xX.X,1),length(log.pmods)*2);
    ws.regs_cF=nan*zeros(size(subjdata{s,2}.SPM.xX.X,1),length(log.pmods));
    ws.regs_ct=nan*zeros(size(subjdata{s,2}.SPM.xX.X,1),length(log.pmods));
    for r=1:length(log.pmods)*2
        ws.regs(:,r)=subjdata{s,2}.SPM.xX.X(:,subjdata{s,2}.regdetails{r,2});
        if r<length(log.pmods)+1
            ws.regs_cF(:,r)=subjdata{s,2}.SPM.xX.X(:,subjdata{s,2}.regdetails{r,2});
        else
            ws.regs_ct(:,r-length(log.pmods))=subjdata{s,2}.SPM.xX.X(:,subjdata{s,2}.regdetails{r,2});
        end
    end
    subjdata{s,2}.regs=ws.regs;
    subjdata{s,2}.regs_cF=ws.regs_cF;
    subjdata{s,2}.regs_ct=ws.regs_ct;
    %
    ws=[];
end

% Perform correlations
for s=1:log.n_subjs
    ws.cF_corrs=nan*zeros(length(log.pmods)*2); ws.cF_corr_p=nan*zeros(length(log.pmods)*2);
    ws.ct_corrs=nan*zeros(length(log.pmods)*2); ws.ct_corr_p=nan*zeros(length(log.pmods)*2);
    ws.regs_cF=subjdata{s,2}.regs_cF;
    ws.regs_ct=subjdata{s,2}.regs_ct;
    for i=1:length(log.pmods)
        for j=1:length(log.pmods)
            % Conflict task
            [r p]=corr(ws.regs_cF(:,i),ws.regs_cF(:,j));
            subjdata{s,2}.cF_corrs(i,j)=r;
            subjdata{s,2}.cF_corr_p(i,j)=p;
            % Control task
            [r p]=corr(ws.regs_ct(:,i),ws.regs_ct(:,j));
            subjdata{s,2}.ct_corrs(i,j)=r;
            subjdata{s,2}.ct_corr_p(i,j)=p;
        end
    end
end

% Plot correlations
figure('Name', 'Absolute Correlation coefficients for Conflict & Control regressors (scaled to maximum)', 'Position', [100 400 800 500])
labels=cell(length(log.pmods),1);
for i=1:length(log.pmods)
    labels{i}=[log.pmods{i} '   (' num2str(i) ')'];
end
for s=1:log.n_subjs
    % Conflict 
    subplot(log.n_subjs,2, (s-1)*2+1)
    imagesc(abs(subjdata{s,2}.cF_corrs), [0 0.5]); colorbar;  axis square
    set(gca,'YTick', 1:length(labels), 'YTickLabel', labels)
    set(gca,'XTick', 1:length(labels));
    % Control
    subplot(log.n_subjs,2, (s-1)*2+2)
    imagesc(abs(subjdata{s,2}.ct_corrs), [0 0.5]); colorbar; axis square
    set(gca,'YTick', 1:length(labels), 'YTickLabel', labels)
    set(gca,'XTick', 1:length(labels));
end







