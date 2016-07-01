% Perform TARGETTED correlations between betas and behaviour.
clear all; close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data']; where.parameter_scripts='/Volumes/PENNYDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
where.where='I:\5 Explore fMRI';  where.expt_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.expt_folder filesep '1 Brain data']; where.data_beh=[where.where filesep '1 Behavioural data'];  where.parameter_scripts='I:\4 Explore experiment\3 Analysis\4 Fit computational models'; where.secondlevelresults='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results';

% Model details
log.factorial=[2 2]; % 2 x N Choices, cF then ct; Empty if not factorial
%
log.n_contrasts=4; % No. of contrasts per ROI
log.n_char4con1=13; % No characters for contrast name (cF_Accept=10, PPI Pos=8, in_cF_Reject=13)
log.clustermodel=1;


% Requested analysis
log.specificsubjects={};

% Beta file *** Edit here ***   
log.betatextfile_where='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c7_Cluster6CompeteFull_XUVPEN\choice_cluster2x2\ROI\2 Selected c7 rois';
% 'C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c7_Cluster6CompeteFull_XUVPEN\choice_cluster2x2\ROI\c3 ROI battery';
% log.betatextfile_where= 'C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c3_CompeteFull_XUVPEN\choice_2x3\ROI\6 c3 Battery betas\c3 ROI battery';




log.betatextfile_name='(03-Nov-2013) Extracted betas';

for o1=1:1 % General setup
    
    % General
    addpath(where.where); addpath([where.where filesep '4 Set up models'])
        
    % Load beta file
    w.bet=importdata([log.betatextfile_where filesep log.betatextfile_name '.txt']);
    d_beta_table=vertcat([{'Subject'} strtrim(w.bet.textdata(1, 2:end))], sortrows([strtrim(w.bet.textdata(2:end,1)) num2cell(w.bet.data)],1));
    [log.subjects_all log.n_subjs d_beta_table] = f_selectsubjects(d_beta_table, log.specificsubjects,  [d_beta_table [{'ok'}; num2cell(ones(size(d_beta_table,1)-1,1))]], 'ok');
    
   
    % Apply further subject selection for some models
    w.slashes=strfind(log.betatextfile_where, filesep);
    log.firstlevelmodel= log.betatextfile_where(strfind(log.betatextfile_where,'m_'):   w.slashes(find(w.slashes>strfind(log.betatextfile_where,'m_'),1,'first'))-1);
    disp(['First level model:    ' log.firstlevelmodel  '       (if incorrect, specify manually - for subject selection)']); input('Hit enter to continue       ');
    w.modelsneedingsubselect={'m_c6_Cluster4CompeteFull_XUVPEN';'m_c7_Cluster6CompeteFull_XUVPEN';'m_c8_Cluster4MovCompeteFull_XUVPEN';'m_c9_Cluster6MovCompeteFull_XUVPEN';'m_c10_Cluster6CompeteRT_XUVPEN'};
    if sum(strcmp(log.firstlevelmodel, w.modelsneedingsubselect))==1
        log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    end
     
    % Check data for nans
    openvar d_beta_table
    if sum(sum(isnan(cell2mat(d_beta_table(2:end, 2:end)))))~=0;  input('Check for NANS!!!!!'); else disp('No NaNs found in data table'); end; disp(' '); disp(' ')
    
    
    % Read details of betas (ROIs, beta list)
    log.beta_list=d_beta_table(1,2:end)'; log.n_betas=length(log.beta_list); % Get list of betas
    disp('(1) Full list of betas:'); disp(char(log.beta_list));
    disp(['(2) No. of contrasts per ROI: ' num2str(log.n_contrasts)]);
    log.n_rois=log.n_betas/log.n_contrasts;  % Get list of ROIs
    log.rois=log.beta_list([1:log.n_rois]*log.n_contrasts-log.n_contrasts+1);
    log.rois=cellfun(@(x)[x(1:(length(x)-log.n_char4con1))], log.rois, 'UniformOutput',0);
    disp(['(3) ROI list     (Assuming ' num2str(log.n_char4con1) ' char for contrast in 1st beta)']); disp(log.rois);
    log.contrasts=log.beta_list(1:log.n_contrasts);  % Get contrast names
    for i=1:log.n_contrasts
        log.contrasts{i}=log.contrasts{i}(length(log.rois{1})+2:end);
    end
    disp('(4) Contrasts for each ROI:'); disp(log.contrasts)
    input('Hit enter if OK   ');
    
end


%% Format beta txt file into table (row=ROI, col=Contrast)

% Specific subset of ROIs?
log.requestedROIs={};

% Format into table
d_subjbetas=cell(log.n_rois+1, log.n_contrasts+1); k=2; % k=first col of d_beta_table to read data out of
for r=1:log.n_rois
    d_subjbetas{r+1,1}=log.rois{r};
    for c=1:log.n_contrasts
        if r==1; d_subjbetas{1,c+1}=log.contrasts{c}; end
        d_subjbetas{r+1,c+1}=cell2mat(d_beta_table(2:end,k)); k=k+1;
    end
end
[log.rois log.n_rois d_subjbetas]=  f_selectsubjects(d_subjbetas, log.requestedROIs,  [d_subjbetas [{'ok'}; num2cell(ones(size(d_subjbetas,1)-1,1))]], 'ok');



openvar d_subjbetas
% % error
% % % Alter order!!
% d_subjbetas=[
%     
% d_subjbetas(1,:);    
% d_subjbetas(4,:);
% d_subjbetas(2,:);
% d_subjbetas(3,:);
% d_subjbetas(5,:);]
% 
% 
% openvar a
% 
% % d_subjbetas=[d_subjbetas(1,:);d_subjbetas(2,:);d_subjbetas(3,:);d_subjbetas(4,:);d_subjbetas(5,:);d_subjbetas(10,:);]; % c3 Battery, Frontal & Striatal
% % % d_subjbetas=[d_subjbetas(1,:);d_subjbetas(6,:);d_subjbetas(7,:);d_subjbetas(8,:);d_subjbetas(9,:)]; % c3 Battery, HPC
% log.rois=d_subjbetas(2:end,1); log.n_rois=length(log.rois);



% Mean centre (by subject x ROI)
request.meancentre=1;
if request.meancentre
    for s=1:log.n_subjs
        for r=1:log.n_rois
            wr.sum=0; % Calc mean
            for c=1:log.n_contrasts
                wr.sum=wr.sum+d_subjbetas{r+1,c+1}(s);
            end
            
            for c=1:log.n_contrasts % Mean centre
                d_subjbetas{r+1,c+1}(s)=d_subjbetas{r+1,c+1}(s)-wr.sum/log.n_contrasts;
            end
        end
    end
end


%% Calculate values for graphs + Stats

% (1) Cell means + error bars
d_graphs.means=cell(log.n_rois+1, log.n_contrasts+1); d_graphs.means(1,2:end)=cellstr(log.contrasts); d_graphs.means(2:end,1)=log.rois;
d_graphs.errorbar_se=cell(log.n_rois+1, log.n_contrasts+1);  d_graphs.errorbar_se(1,2:end)=cellstr(log.contrasts); d_graphs.errorbar_se(2:end,1)=log.rois;
for r=1:log.n_rois
    for c=1:log.n_contrasts
        d_graphs.means{r+1,c+1}=mean(d_subjbetas{r+1,c+1});
        d_graphs.errorbar_se{r+1,c+1}=std(d_subjbetas{r+1,c+1})/ sqrt(log.n_subjs);
    end
end


% (2) One-sample ttests
d_onesam.stats=cell(log.n_rois+1, log.n_contrasts+1); d_onesam.stats(1,2:end)=cellstr(log.contrasts); d_onesam.stats(2:end,1)=log.rois;
for r=1:log.n_rois
    for c=1:log.n_contrasts
        
        [wc.h wc.p wc.ci wc.stats]=ttest(d_subjbetas{r+1,c+1});
        if wc.p<0.001;   wc.sig=[num2str(wc.p,4) '  ***'];
        elseif wc.p<0.01;   wc.sig=[num2str(wc.p,3) '   **'];
        elseif wc.p<0.051;   wc.sig=[num2str(wc.p,3) '    *'];
        else wc.sig=num2str(wc.p,2);
        end
        d_onesam.stats{r+1,c+1}={wc.stats.tstat wc.stats.df  wc.p wc.sig};
        
        wc=[];
    end
end

% (3) Factorial ANOVA (only if factorial, 2 xN)
if isempty(log.factorial)==0
    log.fac_factors={'Task';'Choice'}; disp('Factorial factors:'); disp(log.fac_factors); input('OK?  ');
    r_anova.stats=cell(log.n_rois+1, 4);  % Stats: F, df1, df2, pval
    r_anova.stats{1,2}=log.fac_factors{1};  % Col 1= ME1, Col 2=ME 2, Col3=Interaction
    r_anova.stats{1,3}=log.fac_factors{2};
    r_anova.stats{1,4}=[log.fac_factors{1} 'x' log.fac_factors{1}];
    for r=1:log.n_rois
        r_anova.stats{r+1,1}=log.rois{r};
        %
        if log.factorial(2)==2; wr.d=[d_subjbetas{r+1,2} d_subjbetas{r+1,3} d_subjbetas{r+1,4} d_subjbetas{r+1,5}];
        elseif log.factorial(2)==3; wr.d=[d_subjbetas{r+1,2} d_subjbetas{r+1,3} d_subjbetas{r+1,4} d_subjbetas{r+1,5}  d_subjbetas{r+1,6} d_subjbetas{r+1,7} ];
        else error('Repeated measures ANOVA not yet specified for this number of levels for 2nd factor')
        end
        
        % Run repeated measures ANOVA
        [wr.anova.res]=teg_repeated_measures_ANOVA(wr.d,log.factorial, log.fac_factors); % R: Fstats, df1,df2, p val, ?, ?
        wr.sig=cell(3,1);
        for i=1:3
            if wr.anova.res.R(i,4)<0.001; wr.sig{i}=[num2str(wr.anova.res.R(i,4),4) '  ***'];
            elseif wr.anova.res.R(i,4)<0.01; wr.sig{i}=[num2str(wr.anova.res.R(i,4),3) '   **'];
            elseif wr.anova.res.R(i,4)<0.051; wr.sig{i}=[num2str(wr.anova.res.R(i,4),3) '    *'];
            else wr.sig{i}=num2str(wr.anova.res.R(i,4),2);
            end
        end
        r_anova.stats{r+1,2}= {wr.anova.res.R(1,1) wr.anova.res.R(1,2) wr.anova.res.R(1,3) wr.sig{1}};  % ME Fac 1
        r_anova.stats{r+1,3}= {wr.anova.res.R(2,1) wr.anova.res.R(2,2) wr.anova.res.R(2,3) wr.sig{2}};  % ME Fac 2
        r_anova.stats{r+1,4}= {wr.anova.res.R(3,1) wr.anova.res.R(3,2) wr.anova.res.R(3,3) wr.sig{3}};  % Interaction
        
        wr=[];
    end
end

% (4) Simple effects (Means, SE, Paired T)  - only if factorial, 2xN
if isempty(log.factorial)==0
    d_simplefx.means=cell(log.n_rois+1, log.factorial(2)+factorial(log.factorial(2))+1);
    d_simplefx.errorbar_se=cell(log.n_rois+1, log.factorial(2)+factorial(log.factorial(2))+1);
    r_simplefx.pairedt_stats=cell(log.n_rois+1, log.factorial(2)+factorial(log.factorial(2))+1);  % each cell: t, df, p, sigtype
    k=1;
    % Instructions for compiling simple effects (Name, Var 1, Transformation, Var 2, Short name for graph)
    instruc.simpfx={
        'Accept_cF-ct'      'cF_Accept'     '-'     'ct_NoBomb'         'A';    % Across task comparisons
        'Reject_cF-ct'      'cF_Reject'     '-'     'ct_Bomb'              'R';
        'Explore_cF-ct'    'cF_Explore'   '-'     'ct_Explore'            'E' ;
        'next'                  ' '                   ' '     ' '                           ' ';
        'cF_Acc-Rej'        'cF_Accept'     '-'     'cF_Reject'             'A - R';    % Within Conflict task
        'cF_Rej-Exp'        'cF_Reject'     '-'     'cF_Explore'           'R - E';
        'cF_Acc-Exp'        'cF_Accept'     '-'     'cF_Explore'           'A - E';
        'next'                  ' '                   ' '     ' '                           ' ';
        'ct_NoB-Bmb'      'ct_NoBomb'   '-'     'ct_Bomb'             'A - R';    % Within Control task
        'ct_Bmb-Exp'       'ct_Bomb'       '-'     'ct_Explore'          'R - E';
        'ct_NoB-Exp'       'ct_NoBomb'   '-'     'ct_Explore'           'A - E';
        };
    
    if log.clustermodel
        instruc.simpfx={
            'Reject_cF-ct'      'in_cF_Reject'     '-'     'in_ct_Bomb'              'R';
            'Explore_cF-ct'    'in_cF_Explore'   '-'     'in_ct_Explore'            'E' ;
            'next'                  ' '                   ' '     ' '                           ' ';
            'cF_Rej-Exp'        'in_cF_Reject'     '-'     'in_cF_Explore'           'R - E';    % Within Conflict task
            'ct_Bmb-Exp'       'in_ct_Bomb'       '-'     'in_ct_Explore'          'R - E';
            };
    end
    
    for i=1:size(instruc.simpfx,1)
        if strcmp(instruc.simpfx{k,1},'next')==1
            d_simplefx.means{1,k+1}=' ';
            d_simplefx.errorbar_se{1,k+1}=' ';
            r_simplefx.pairedt_stats{1,k+1}=' ';
            for r=1:log.n_rois
                d_simplefx.means{r+1,k+1}=nan;
                d_simplefx.errorbar_se{r+1,k+1}=nan;
                r_simplefx.pairedt_stats{r+1,k+1}=nan;
            end
            k=k+1;
        else
            try
                d_simplefx.means{1,k+1}=instruc.simpfx{k,1};
                d_simplefx.errorbar_se{1,k+1}=instruc.simpfx{k,1};
                r_simplefx.pairedt_stats{1,k+1}=instruc.simpfx{k,1};
                for r=1:log.n_rois
                    if k==1;
                        d_simplefx.means{r+1,1}=log.rois{r};
                        d_simplefx.errorbar_se{r+1,1}=log.rois{r};
                        r_simplefx.pairedt_stats{r+1,1}=log.rois{r};
                    end
                    
                    % Data
                    wr.d1=d_subjbetas{r+1, find(strcmp(log.contrasts, instruc.simpfx{k,2}))+1};
                    wr.d2=d_subjbetas{r+1, find(strcmp(log.contrasts, instruc.simpfx{k,4}))+1};
                    eval(['wr.d=wr.d1' instruc.simpfx{k,3} 'wr.d2;'])
                    
                    % Stats
                    d_simplefx.means{r+1,k+1}=mean(wr.d);
                    d_simplefx.errorbar_se{r+1,k+1}=std(wr.d)/sqrt(log.n_subjs);
                    [wr.h wr.p wr.ci wr.stats]=ttest(wr.d);
                    if wr.p<0.001;   wr.sig=[num2str(wr.p,4) '  ***'];
                    elseif wr.p<0.01;   wr.sig=[num2str(wr.p,3) '   **'];
                    elseif wr.p<0.051;   wr.sig=[num2str(wr.p,3) '    *'];
                    else wr.sig=num2str(wr.p,2);
                    end
                    r_simplefx.pairedt_stats{r+1,k+1}={wr.stats.tstat wr.stats.df wr.p wr.sig}; wr=[];
                    
                    %
                    wr=[];
                end
                disp(['Simple effects comparison: ' instruc.simpfx{k,1}]); k=k+1;
            catch
                disp(['Simple effects FAILED: ' instruc.simpfx{k,1}])
            end
        end
    end
    
    log.n_simpfx=size(d_simplefx.means,2)-1;
    log.simpfx=d_simplefx.means(1,2:end)';
end

%% Graphs: Factorial, Simpfx;  Print Anovastat, mark simpfx sig

if isempty(log.factorial)==0
    f.figheight=log.n_rois*250+120;
    ff=figure('Name', ['ROI group statistics (n=' num2str(log.n_subjs) ')'],'NumberTitle','off','Position',[200,00,1400,f.figheight]); set(gcf,'Color', 'w');
    if log.clustermodel; f.bargraph_errorbarslack=1.35; else f.bargraph_errorbarslack=0.65;  end
    
    f.subplotcols=4; f.marksig_shiftstarsleft=0.08;  f.marksig_positions=(1:size(instruc.simpfx,1))/size(instruc.simpfx,1)-f.marksig_shiftstarsleft;
    
    f.subplot_VerHorz=[0.07 0.07]; f.fig_BotTop=[0.08 0.08]; f.submw=[0.05 0.1];
    f.printbetafile_xy=[0.15 4.2];
    for r=1:log.n_rois
        wr.mean=cell2mat(d_graphs.means(r+1,2:end)); % formatting
        wr.SE=cell2mat(d_graphs.errorbar_se(r+1,2:end));
        
        % Group means graph
        if r==1
            subtightplot(log.n_rois,f.subplotcols,1,f.subplot_VerHorz,f.fig_BotTop, f.submw); % Print details to figure
            text(f.printbetafile_xy(1),f.printbetafile_xy(2),  ['Beta file: ' log.betatextfile_where filesep log.betatextfile_name])
            hold on;
        else
            subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols +1,f.subplot_VerHorz,f.fig_BotTop, f.submw);
        end
        if log.clustermodel
            barwitherr([d_graphs.errorbar_se{r+1,2} d_graphs.errorbar_se{r+1,3}; d_graphs.errorbar_se{r+1,4} d_graphs.errorbar_se{r+1,5} ],[d_graphs.means{r+1,2} d_graphs.means{r+1,3}; d_graphs.means{r+1,4} d_graphs.means{r+1,5} ])
        else
            barwitherr([d_graphs.errorbar_se{r+1,2} d_graphs.errorbar_se{r+1,3} d_graphs.errorbar_se{r+1,4}; d_graphs.errorbar_se{r+1,5} d_graphs.errorbar_se{r+1,6} d_graphs.errorbar_se{r+1,7}],[d_graphs.means{r+1,2} d_graphs.means{r+1,3} d_graphs.means{r+1,4}; d_graphs.means{r+1,5} d_graphs.means{r+1,6} d_graphs.means{r+1,7}])
        end        
        axis([0.5, log.n_contrasts/log.factorial(2)+0.5,  min(wr.mean)-max(wr.SE)-f.bargraph_errorbarslack,  max(wr.mean)+max(wr.SE)+f.bargraph_errorbarslack])
        title(log.rois{r}); ylabel('beta value');set(gca,'XTick',[1 2]); set(gca,'XTickLabel', {'cF' 'ct'});
        
        % Simple effects
        subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols +2,f.subplot_VerHorz,f.fig_BotTop, f.submw);
        barwitherr(cell2mat(d_simplefx.errorbar_se(r+1,2:end)), cell2mat(d_simplefx.means(r+1,2:end)), 'y')
        axis([0.5, log.n_simpfx+0.5,  min(wr.mean)-max(wr.SE)-f.bargraph_errorbarslack,  max(wr.mean)+max(wr.SE)+f.bargraph_errorbarslack])
        if log.clustermodel
            set(gca,'XTick',[1:2 4:5]); set(gca,'XTickLabel', instruc.simpfx([1:2 4:5],5)); rotateXLabels(gca,90); ylabel('difference (betas)');
            title(' cF - ct                         cF          ct')
        else
            set(gca,'XTick',[1:3 5:7 9:11]); set(gca,'XTickLabel', instruc.simpfx([1:3 5:7 9:11],5)); rotateXLabels(gca,90); ylabel('difference (betas)');
            title('            cF - ct            cF              ct            ')
        end
        for s=1:log.n_simpfx % Mark significance
            if size(r_simplefx.pairedt_stats{r+1, s+1},2)>1  && cell2mat(r_simplefx.pairedt_stats{r+1, s+1}(3))<0.051
                text(f.marksig_positions(s), 0.85, '*', 'Units','Normalized', 'color',[ 0.8 0 0], 'Fontsize', 23)
            end
        end
        
        
        % Write ANOVA stats
        f.lineanova=0.75:-0.15:-0.5; k=-0.25;
        subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols+3,f.subplot_VerHorz,f.fig_BotTop, f.submw); axis off;
        text(k,f.lineanova(1),  [log.fac_factors{1} ':   p='  char(r_anova.stats{r+1,2}(4))] )
        text(k,f.lineanova(2),  [log.fac_factors{2} ':   p='  char(r_anova.stats{r+1,3}(4))] )
        text(k,f.lineanova(3),  ['Interac:   p='  char(r_anova.stats{r+1,4}(4))] )
        
        % Write simple effects stats
        f.linett=0.9:-0.15:-0.85; k=0.55; kk=-0.1;
        if log.clustermodel
            subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols+4,f.subplot_VerHorz,f.fig_BotTop, f.submw); axis off;
            s=1; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=2; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=4; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF :   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=5; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   ct :   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        else
            subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols+4,f.subplot_VerHorz,f.fig_BotTop, f.submw); axis off;
            s=1; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=2; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=3; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            %
            s=5; text(-k,f.linett(5),  ['[cF]  '  instruc.simpfx{s,5}  ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=6; text(-k,f.linett(6),  ['[cF]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=7; text(-k,f.linett(7),  ['[cF] '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            %
            s=9; text(k,f.linett(5),  ['[ct]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=10; text(k,f.linett(6),  ['[ct]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
            s=11; text(k,f.linett(7),  ['[ct]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        end
    end
end

%% Non-factorial: One-sample ttest against baseline 0


% d_onesam.stats

if isempty(log.factorial)==1
    f.figheight=log.n_rois*230+120; f.figwidth=600;
    ff=figure('Name', ['ROI group statistics (n=' num2str(log.n_subjs) ')'],'NumberTitle','off','Position',[200,00,f.figwidth,f.figheight]); set(gcf,'Color', 'w');
    f.bargraph_errorbarslack=0.95; f.subplotcols=2; f.marksig_shiftstarsleft=0.08;  
    
    stopped here. rethinking. maybe use the combined script after all :S
    
    
    
    f.marksig_positions=(1:size(log.simpfx,1))/size(instruc.simpfx,1)-f.marksig_shiftstarsleft;   
    f.subplot_VerHorz=[0.07 0.07]; f.fig_BotTop=[0.08 0.08]; f.submw=[0.05 0.1];
    f.printbetafile_xy=[0.15 2.7];
    
    for r=1:log.n_rois
        wr.mean=cell2mat(d_graphs.means(r+1,2:end)); % formatting
        wr.SE=cell2mat(d_graphs.errorbar_se(r+1,2:end));
        
        % Group means graph
        if r==1
            subtightplot(log.n_rois,f.subplotcols,1,f.subplot_VerHorz,f.fig_BotTop, f.submw); % Print details to figure
            text(f.printbetafile_xy(1),f.printbetafile_xy(2),  ['Beta file: ' log.betatextfile_where filesep log.betatextfile_name])
            hold on;
        else
            subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols +1,f.subplot_VerHorz,f.fig_BotTop, f.submw);
        end
        barwitherr([d_graphs.errorbar_se{r+1,2} d_graphs.errorbar_se{r+1,3} d_graphs.errorbar_se{r+1,4}; d_graphs.errorbar_se{r+1,5} d_graphs.errorbar_se{r+1,6} d_graphs.errorbar_se{r+1,7}],[d_graphs.means{r+1,2} d_graphs.means{r+1,3} d_graphs.means{r+1,4}; d_graphs.means{r+1,5} d_graphs.means{r+1,6} d_graphs.means{r+1,7}])
        axis([0.5, log.n_contrasts/log.factorial(2)+0.5,  min(wr.mean)-max(wr.SE)-f.bargraph_errorbarslack,  max(wr.mean)+max(wr.SE)+f.bargraph_errorbarslack])
        title(log.rois{r}); ylabel('beta value');set(gca,'XTick',[1 2]); set(gca,'XTickLabel', {'cF' 'ct'});
        
        % Simple effects
        subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols +2,f.subplot_VerHorz,f.fig_BotTop, f.submw);
        barwitherr(cell2mat(d_simplefx.errorbar_se(r+1,2:end)), cell2mat(d_simplefx.means(r+1,2:end)), 'y')
        axis([0.5, log.n_simpfx+0.5,  min(wr.mean)-max(wr.SE)-f.bargraph_errorbarslack,  max(wr.mean)+max(wr.SE)+f.bargraph_errorbarslack])
        set(gca,'XTick',[1:3 5:7 9:11]); set(gca,'XTickLabel', instruc.simpfx([1:3 5:7 9:11],5)); rotateXLabels(gca,90); ylabel('difference (betas)');
        title('            cF - ct            cF              ct            ')
        for s=1:log.n_simpfx % Mark significance
            if size(r_simplefx.pairedt_stats{r+1, s+1},2)>1  && cell2mat(r_simplefx.pairedt_stats{r+1, s+1}(3))<0.051
                text(f.marksig_positions(s), 0.85, '*', 'Units','Normalized', 'color',[ 0.8 0 0], 'Fontsize', 23)
            end
        end
        
        % Write ANOVA stats
        f.lineanova=0.75:-0.15:-0.5; k=-0.25;
        subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols+3,f.subplot_VerHorz,f.fig_BotTop, f.submw); axis off;
        text(k,f.lineanova(1),  [log.fac_factors{1} ':   p='  char(r_anova.stats{r+1,2}(4))] )
        text(k,f.lineanova(2),  [log.fac_factors{2} ':   p='  char(r_anova.stats{r+1,3}(4))] )
        text(k,f.lineanova(3),  ['Interac:   p='  char(r_anova.stats{r+1,4}(4))] )
        
        % Write simple effects stats
        f.linett=0.9:-0.15:-0.85; k=0.55; kk=-0.1;
        subtightplot(log.n_rois,f.subplotcols,(r-1)*f.subplotcols+4,f.subplot_VerHorz,f.fig_BotTop, f.submw); axis off;
        s=1; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        s=2; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        s=3; text(kk,f.linett(s),  ['[' instruc.simpfx{s,5} ']   cF - ct:   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        %
        s=5; text(-k,f.linett(5),  ['[cF]  '  instruc.simpfx{s,5}  ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        s=6; text(-k,f.linett(6),  ['[cF]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        s=7; text(-k,f.linett(7),  ['[cF] '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        %
        s=9; text(k,f.linett(5),  ['[ct]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        s=10; text(k,f.linett(6),  ['[ct]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
        s=11; text(k,f.linett(7),  ['[ct]  '  instruc.simpfx{s,5} ':   p='  char(r_simplefx.pairedt_stats{r+1,s+1}(4))] )
    end





end




