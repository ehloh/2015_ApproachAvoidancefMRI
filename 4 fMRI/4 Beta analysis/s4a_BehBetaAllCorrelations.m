% Perform ALL correlations between betas and behaviour.
clear all; close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data']; where.parameter_scripts='/Volumes/PENNYDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
where.where='I:\5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data'];  where.parameter_scripts='I:\4 Explore experiment\3 Analysis\4 Fit computational models'; where.secondlevelresults='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results';

%
log.specificsubjects={};
log.behprofile='Behavioural profile (23-Aug-2013)';
% rlvar='VExplore';

% Beta values  *** Edit here ***
% log.where_betafile=[where.secondlevelresults filesep 'm_c3_CompeteFull_XUVPEN\choice_2x3\ROI'];
% log.where_betafile='/Volumes/PENNYDISK/5 Explore fMRI/7 Behavioural analysis';
% log.name_betatextfile='(20-Aug-2013) Extracted betas 4 SPSS TaskxChoice';
% log.name_betatextfile='(29-Aug-2013) Extracted betas 4 SPSS - c3 MEChoice';
% log.name_betatextfile='(29-Aug-2013) Extracted betas 4 SPSS c3 TxC';
% log.name_betatextfile='(29-Aug-2013) Extracted betas 4 SPSS c7';
% log.where_betafile='I:\5 Explore fMRI\7 Behavioural analysis';
log.where_betafile='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c3_CompeteFull_XUVPEN\choice_2x3\ROI\5 ME Choice striatum';
% 'C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c3_CompeteFull_XUVPEN\choice_2x3\ROI\4 Task x Choice Excluding MEChoice\masked by contrast';


% 'C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c7_Cluster6CompeteFull_XUVPEN\choice_cluster2x2\ROI\HPC only';
% 'C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c3_CompeteFull_XUVPEN\choice_2x3\ROI\4 Task x Choice Excluding MEChoice\masked by image';
% 'C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c3_CompeteFull_XUVPEN\choice_2x3\ROI';
log.name_betatextfile='(13-Sep-2013) Extracted betas';

for o1=1:1 % General setup
    
    % Load subjects
    addpath(where.where);
    log.w=load([where.data_beh filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects_all log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
end

%% Load beta file & behaviour file (data in variable d, other details in variable log)

% Load beta file
w.bet=importdata([log.where_betafile filesep log.name_betatextfile '.txt']);
d.beta_table=vertcat([{'Subject'} w.bet.textdata(1, 2:end)], sortrows([w.bet.textdata(2:end,1) num2cell(w.bet.data)],1));
log.beta_list=d.beta_table(1,2:end)'; log.n_betas=length(log.beta_list);

% Load behavioural profile
w.beh=importdata([where.data_beh filesep 'Group behaviour files' filesep log.behprofile '.txt']);
d.beh_table_all=vertcat([{'Subject'} w.beh.textdata(1, 2:end)], sortrows([w.beh.textdata(2:end,1) num2cell(w.beh.data)],1));
log.beh_list=d.beh_table_all(1,2:end)'; log.n_beh=length(log.beh_list);

% Apply subject selection according to beta table
k=size(d.beh_table_all,2)+1;
d.beh_table_all{1,k}='BetaInclude';
for s=1:size(d.beh_table_all,1)-1
    if sum(strcmp(d.beta_table(2:end,1), d.beh_table_all{s+1,1}))==1
        d.beh_table_all{s+1,k}=1;
    else
        d.beh_table_all{s+1,k}=0;
    end
end
[log.subjects log.n_subjs d.beh_table]=f_selectsubjects(d.beh_table_all, {}, d.beh_table_all, 'BetaInclude');

% Remove names and tables from data
d.beh_table(:,1)=[]; d.beh_table(1,:)=[];
d.beta_table(:,1)=[]; d.beta_table(1,:)=[];

%% (2) Execute correlations (grouped by ROI)

r_corr=cell(log.n_betas,2); % Col 2=correlation results (Col 1=r, Col 2=p, Col 3=Sig?)
for r=1:log.n_betas
    r_corr{r,1}=log.beta_list{r};
    r_corr{r,2}=nan*zeros(log.n_beh,3);
    for b=1:log.n_beh
        [r_corr{r,2}(b,1) r_corr{r,2}(b,2)]=corr(cell2mat(d.beh_table(:,b)),cell2mat(d.beta_table(:,r)));
        if r_corr{r,2}(b,2)<0.001
            r_corr{r,2}(b,3)=3;
        elseif r_corr{r,2}(b,2)<0.01
            r_corr{r,2}(b,3)=2;
        elseif r_corr{r,2}(b,2)<0.05
            r_corr{r,2}(b,3)=1;
%         elseif r_corr{r,2}(b,2)<0.1 % No trends
%             r_corr{r,2}(b,3)=0.5;
        else
            r_corr{r,2}(b,3)=0;
        end
    end
    r_corr{r,3}=sum(r_corr{r,2}(:,3)>0);
end

%% (3) Report correlation results

%  Group by ROIs
r_report=cell(sum(cell2mat(r_corr(:,3)))+size(r_corr,1)*1,1); k=1;
for r=1:log.n_betas
    if r_corr{r,3}>0
        k=k+1;r_report{k,1}=['---- ' log.beta_list{r} ' -----']; 
        k=k+2;
        for b=1:log.n_beh
            if r_corr{r,2}(b,3)>0.5 % Exclude trends
                switch r_corr{r,2}(b,3)
                    case 0.5; s='    (t)';
                    case 1; s='    *';
                    case 2; s='    **';
                    case 3; s='    ***';
                end
                r_report{k,1}=log.beh_list{b};
                r_report{k,2}=['r= ' num2str(r_corr{r,2}(b,1),2) ',  p= ' num2str(r_corr{r,2}(b,2),3)  s];
                k=k+1;
            end
        end
%         r_report{k,1}='      #####             ';
%         k=k+1;
if b==log.n_beh
    k=k+1;
end
    else
        r_report{k,1}=['---- ' log.beta_list{r} ' -----']; 
        k=k+1;
    end
end


% Print out 
printok=print2txt(log.where_betafile, ['(' date ') Beta correlations'], r_report);
disp(printok)

%%

a=find(strcmp(log.beh_list, 'per.cF_Reject '))
b=find(strcmp(log.beh_list, 'p.st.Trait '));
[r p]=corr(cell2mat(d.beh_table(:,a)),cell2mat(d.beh_table(:,b)))

