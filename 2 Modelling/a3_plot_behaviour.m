% plot_behaviour (performance)
%       Choice contingencies for both individuals & overall group mean
clear all; close all hidden; clc
% clear all;clc

% Requested
request.tasktype=1; % 1=Conflict, 2=Control
request.sess_fMRIpractice=1;   % 1=fMRI, 2=Practice session
requests.subs_IncExcl=1;   % 1=Included subs, 2=Excluded subs
%
request.plot_individuals=0;
request.plot_mean=0;
request.plot_correlations=0;
request.plotrt_mean=1;
request.plotrt_individuals =0; 

% Specification details
request.onlylatterhalf=0;
request.outcome_includesexplorationcost=0;

for o1=1:1  % Settings
    
    w.w=pwd;
    if strcmp(w.w(1), '/')==1; 
        where.where='/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis'; 
        where.scripts=[where.where '/4 Fit computational models']; 
        where.where2='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI'; 
        where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI subjdata/1 Brain data';
    else  where.where='C:\Users\e.loh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis';
        where.scripts=[where.where fs '4 Fit computational models'];
        where.where2='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI'; 
        where.data_brain='G:\2 [Explore]\1 Brain data';
    end
    
    
    
    
    switch requests.subs_IncExcl
        case 1, where.data_beh= [where.where2 filesep '1 Behavioural data'];  disp('Subjects: Included (i.e. with fMRI');
        case 2,  where.data_beh=[where.where2 filesep '1 Behavioural data' fs 'Excluded subjects'];  disp('Subjects: Excluded (i.e. screened out');
            request.exclsubs.allwsomedata= {'p03_HA';'p05_PY';'p07_AM';'p09_MN';'p11_LB';'p12_AB';'p14_CW';'p16_AM';'p19_DM';'p20_GN';'p22_JL';'p25_AB';'p26_IK';'p28_RT';'p29_RS';'p31_MG';'p32_SP';'p39_TW';'p40_RB';}; % Excluded Subs W some Data
            request.exclsubs.allwallpracdata=  {'p03_HA';'p07_AM';'p09_MN';'p11_LB';'p12_AB';'p19_DM';'p22_JL';'p25_AB';'p26_IK';'p29_RS';'p32_SP';'p39_TW';};
    end
    cd(where.data_beh); w.s=dir('p*');  logg.subjects=cellstr(char(w.s.name)); logg.n_subjs=length(logg.subjects); cd(where.scripts)

    % Settings that don't change much
    switch request.sess_fMRIpractice
        case 1, request.taskfile='taskfMRI';
        case 2, request.taskfile='taskprac';
    end
    request.tasknum=[];
    
    % Settings for plots
    w.axison=1; % Axis of plots
    w.titleson=1;
    
    % Save variable?
    plots=cell(3,1);
    
    disp(request)
    disp('---------------------------------------------------------------------------------')
end

%%  Load data for all subjects
%   subjdata: col 1=subject, col 2= trialstats, col 3= 6x6 matrices 

for o1=1:1  % Columns from fMRI subjdata
    col.NTokPairs=2;
    col.EnvThreat=3;
    col.Task=5;
    col.Trialnum=7;
    col.Choice=8;
    col.Choice2=10;
    col.OutcomeMag=15;
    col.RT1=9;
    col.TrialValid=13;
    col.OutcomePres=14;
    col.ExploredSee=18;
end

% Load data
subjdata=cell(logg.n_subjs+1,2); subjdata{logg.n_subjs+1,3}.e_a=zeros(6,6); subjdata{logg.n_subjs+1,3}.e_r=zeros(6,6); subjdata{logg.n_subjs+1,2}=[];
for  s=1: logg.n_subjs 
    ws.w=load([where.data_beh filesep logg.subjects{s} filesep logg.subjects{s} '_file_' num2str(request.tasknum) request.taskfile]);
    switch request.taskfile
        case 'taskfMRI'   % fMRI task (cF and ct)
            switch request.tasktype
                case 1; ws.d.subjdata=ws.w.conflict; % Conflict
                case 2; ws.d.subjdata=ws.w.control; % Control
            end
            ws.d.settings=ws.w.settings; 
        case 'taskprac'  % Practice sessions for fMRI screening  (cF and ct)
            switch request.tasktype
                case 1; ws.d.subjdata=ws.w.conflict; % Conflict
                case 2; ws.d.subjdata=ws.w.control; % Control
            end
            ws.d.settings=ws.w.settings; 
%         case  'taskconflictcontrol'  % Integrated conflict & control  
%             ws.d=ws.w.taskconflictcontrol;
%             switch request.tasktype
%                 case 1; ws.d.subjdata=ws.d.conflictsubjdata; % Conflict
%                 case 2; ws.d.subjdata=ws.d.controlsubjdata; % Control
%             end
%         case 'taskcontrol', ws.d=ws.w.taskcontrol; % Control only
%         case 'taskconflict', eval(['ws.d=ws.w.' request.taskfile ';'])% Conflict only
        otherwise, error('Check sessions!'); 
    end
    if request.onlylatterhalf==1 % Latter half of session only (each task separtely)
        ws.d.subjdata= sortrows(ws.d.subjdata, col.Trialnum); 
        ws.d.subjdata= ws.d.subjdata(round(size(ws.d.subjdata,1)/2)+1:end, :);
    end    
    
    % RT
    ws.d.subjdata(ws.d.subjdata(:, col.RT1)<0)=nan; 
    ws.meanrt =  nanmean(ws.d.subjdata(:,col.RT1));  
    ws.d.subjdata(:, col.RT1) = log(ws.d.subjdata(:, col.RT1));  if s==1, disp('RTs logg'), end % Log RTs by subject
    % %     ws.d.subjdata(:, col.RT1) =wr.rts(wr.rts(:,1)==request.tasktype,2);  if s==1, disp('RTs z scored'), end
    %     ws.d.subjdata(:, col.RT1) = ws.d.subjdata(:, col.RT1)-ws.meanrt;  if s==1, disp('RTs mean centred'), end % Mean centre RTs by subject
    
    subjdata{s,1}=logg.subjects{s};
    subjdata{s,2}=ws.d.subjdata; 
    subjdata{s,2}=subjdata{s,2}(subjdata{s,2}(:,col.TrialValid)==1, :);
    subjdata{logg.n_subjs+1,2}=[subjdata{logg.n_subjs+1,2}; subjdata{s,2}];
    ws=[];
end

% Compile/format data
for s=1:logg.n_subjs+1
    ws.d=subjdata{s,2};
    ws.pa=zeros(6,6); ws.pr=zeros(6,6); ws.pe=zeros(6,6); % Set up data variables
    ws.rt_a=nan(6,6); ws.rt_r=nan(6,6); ws.rt_e=nan(6,6); ws.rt=nan(6,6);
    ws.e_a=nan(6,6); ws.e_r=ws.e_a; ws.ens_a=ws.e_a; ws.ens_r=ws.e_a;
    
    for e=1:6
        for n=1:6
            wc.d=ws.d(ws.d(:, col.EnvThreat)==e &  ws.d(:, col.NTokPairs)==n,:);
            wc.da=wc.d(wc.d(:,col.Choice)==1,:); wc.dr=wc.d(wc.d(:,col.Choice)==2,:); wc.de=wc.d(wc.d(:,col.Choice)==3,:);
            wc.dens=wc.de( wc.de(:,  col.ExploredSee)==0, :);
            
            % All requested variables
            ws.rt_a(7-e,n)=mean(wc.da(:,col.RT1)); 
            ws.rt_r(7-e,n)=mean(wc.dr(:,col.RT1));
            ws.rt_e(7-e,n)=mean(wc.de(:,col.RT1));
            ws.rt(7-e,n)=mean(wc.d(:,col.RT1));
            
            ws.pa(7-e,n)=mean(wc.d(:,col.Choice)==1);
            ws.pr(7-e,n)=mean(wc.d(:,col.Choice)==2);
            ws.pe(7-e,n)=mean(wc.d(:,col.Choice)==3);
%             if isnan(ws.pa(7-e,n))==1; ws.pa(7-e,n)=0; end            
%             if isnan(ws.pr(7-e,n))==1; ws.pr(7-e,n)=0; end
%             if isnan(ws.pe(7-e,n))==1; ws.p(7-e,n)=0; end
            if isempty(wc.de)==0 && sum(wc.de(:, col.OutcomePres))~=0
                ws.e_a(7-e,n)=mean(wc.de(:, col.Choice2)==1);
                ws.e_r(7-e,n)=mean(wc.de(:, col.Choice2)==2);
                if isempty(wc.dens)==0 && sum(wc.dens(:, col.OutcomePres))~=0
                    ws.ens_a(7-e,n)=mean(wc.dens(:, col.Choice2)==1);
                    ws.ens_r(7-e,n)=mean(wc.dens(:, col.Choice2)==2);
                end
            end
            wc=[];
        end  
    end
    
    % Record
    subjdata{s,3}.pa=ws.pa;
    subjdata{s,3}.pr=ws.pr;
    subjdata{s,3}.pe=ws.pe;
    subjdata{s,3}.rt=ws.rt;
    subjdata{s,3}.rt_a=ws.rt_a;
    subjdata{s,3}.rt_r=ws.rt_r;
    subjdata{s,3}.rt_e=ws.rt_e;
    subjdata{s,3}.e_a=ws.e_a;
    subjdata{s,3}.e_r=ws.e_r;
    subjdata{s,3}.ens_a=ws.ens_a;
    subjdata{s,3}.ens_r=ws.ens_r;
%     subjdata{logg.n_subjs+1,3}.e_a=subjdata{logg.n_subjs+1,3}.e_a+ws.e_a;
%     subjdata{logg.n_subjs+1,3}.e_r=subjdata{logg.n_subjs+1,3}.e_r+ws.e_r;
    
    ws=[];
end
% subjdata{logg.n_subjs+1,3}.e_a=subjdata{logg.n_subjs+1,3}.e_a/logg.n_subjs;
% subjdata{logg.n_subjs+1,3}.e_r=subjdata{logg.n_subjs+1,3}.e_r/logg.n_subjs;

sd=subjdata;

% Other calcs
d_pch=nan(logg.n_subjs,3);  openvar d_pch
for s=1:logg.n_subjs
    for c=1:3
        d_pch(s,c)=mean(subjdata{s,2}(:, col.Choice)==c);
        
    end
end

%% Plot individual choices 

% Figure details
if request.plot_individuals==1
    f.figwidth= 400; f.figheight=1600; f.subplotcols=4; f.subplot_VerHorz=[0.01 0.07]; f.fig_BotTop=[0.01 0.03]; f.fig_LeftRight=[0.2 0.3];
    figure('Color', 'w', 'Name', 'Individual choice performance (Accept, Reject, Explore)', 'NumberTitle', 'off', 'Position', [100 -250 f.figwidth f.figheight]); 
    
    for s=1:logg.n_subjs
        subtightplot(logg.n_subjs,f.subplotcols,(s-1)*f.subplotcols+1,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        text(0.4,0, logg.subjects{s}); axis off
        
        % Change here to plot some other variable 
        titles={'% Accept'; '% Reject'; '% Explore';};  ws.p{1}=subjdata{s,3}.pa;    ws.p{2}=subjdata{s,3}.pr; ws.p{3}=subjdata{s,3}.pe;
%         titles={'% Accept (post Explore)'; '% Reject  (post Explore)'; ' ';};   ws.p{1}=subjdata{s,3}.e_a; ws.p{2}=subjdata{s,3}.e_r; ws.p{3}=nan;
%         titles={'% Accept (post Explore ~See)'; '% Reject  (post Explore ~See)'; ' ';};    ws.p{1}=subjdata{s,3}.ens_a; ws.p{2}=subjdata{s,3}.ens_r; ws.p{3}=nan;
    
        
        for c=1:3 % Plot
            subtightplot(logg.n_subjs,f.subplotcols,(s-1)*f.subplotcols+c+1,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
            imagescnan(ws.p{c},'NanColor', [0 0 0]); axis square; axis off; % colorbar
            if s==1; title(titles{c}); end
            caxis([0 1])
        end
        
    end
    
end

%% Plot mean choices
% In variable subjdata, last row is pooled subjdata

if request.plot_mean==1;
    
%     Plot what?
    titles={'% Accept'; '% Reject'; '% Explore';}; plots{1}=subjdata{logg.n_subjs+1,3}.pa;   plots{2}=subjdata{logg.n_subjs+1,3}.pr;     plots{3}=subjdata{logg.n_subjs+1,3}.pe;
%     titles={'% Accept (post Explore)'; '% Reject  (post Explore)'; ' ';};   plots{1}=subjdata{logg.n_subjs+1,3}.e_a;   plots{2}=subjdata{logg.n_subjs+1,3}.e_r; plots{3}=nan(6,6); 
%     titles={'% Accept (post Explore ~See)'; '% Reject  (post Explore ~See)'; ' ';};   plots{1}=subjdata{logg.n_subjs+1,3}.ens_a;   plots{2}=subjdata{logg.n_subjs+1,3}.ens_r; plots{3}=nan(6,6); 
    
    % Execute
    f.FontSize=25; f.FontName='PT Sans Caption'; 
    
    figure('Color', 'w', 'Name',['Mean choice performance (n=' num2str(logg.n_subjs) ')' ],'NumberTitle','off','Position',[715,400,900,600]);
    for c=1:3
        subplot(1,3,c); imagesc(plots{c}); 
        axis square; title(titles{c},'FontSize',f.FontSize+10, 'FontName', f.FontName); 
        ylabel('EnvThreat','FontSize',f.FontSize, 'FontName', f.FontName);  xlabel('No. Tokens','FontSize',f.FontSize, 'FontName', f.FontName);  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
        caxis([0 1]),  colorbar
        set(gca,'FontSize',f.FontSize, 'FontName', f.FontName, 'LineWidth', 0.8);
        
    end
%     subplot(1,4,4); colorbar; axis off
    
end


%% Correlations
%     (Do design variables correlate with each other? correlation <0.4 is ok)

doplotvarcor=0;
if doplotvarcor
%     details.designvariables={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'OutcomeMean'; 'OutcomeVariance'};  % Full
    details.designvariables={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'}; %  'OutcomeMean'; 'OutcomeVariance'
    % details.designvariables={'OutcomeMean'; 'OutcomeVariance'};
    % details.designvariables={'Conflict'; 'OutcomeMean'};
    
    for o1=1:1 % subjdata specifications (columns)
        col.Task=5; % Event classifiers
        col.Resp1=8;
        col.Resp2=10;
        col.Error=13;
        col.OutcomePresented=14;
        col.OutcomeAmount=15;
        
        col.EnvThreat=3; % Design columns
        col.NTokens=2; % convert to n, not pairs
        col.pLoss=35;
        col.Entropy=36;
        col.Conflict=37;
        col.OutcomeMean=38;
        col.OutcomeVariance=39;
        
        for i=1:length(details.designvariables)
            details.labels{i}=[details.designvariables{i} '   (' num2str(i) ')'];
        end
        
    end
    
    % Calculate correlations
    for s=1:logg.n_subjs
        ws.d=subjdata{s,2};
        
        % Mark details - design
        ws.d(:,col.NTokens)=ws.d(:,col.NTokens).*2;  % nTokens, not nToken Pairs
        ws.d(:,col.pLoss)=(ws.d(:,3)/6).*(ws.d(:,col.NTokens)/12);
        for i=1:size(ws.d,1)
            if ws.d(i,col.pLoss)==1 % pLoss values  of 1 = corrected by -0.00001 for entropy calculation
                ws.d(i,col.Entropy)=-(ws.d(i,col.pLoss)-0.00001).*log(ws.d(i,col.pLoss)-0.00001)  -   (1-(ws.d(i,col.pLoss)-0.00001)).*log(1-(ws.d(i,col.pLoss)-0.00001));        % Uncertainty/entropy:  -p(loggp)-(1-p)logg(1-p) (Col 4)
            else
                ws.d(i,col.Entropy)=-(ws.d(i,col.pLoss)).*log(ws.d(i,col.pLoss))  -   (1-ws.d(i,col.pLoss)).*log(1-ws.d(i,col.pLoss));        % Uncertainty/entropy:  -p(loggp)-(1-p)logg(1-p) (Col 4)
            end
        end
        ws.d(:, col.Conflict)=ws.d(:,col.Entropy).*ws.d(:,col.NTokens); % Conflict: Entropy x N tokens
        ws.d(:, [col.OutcomeMean col.OutcomeVariance])=nan;
        if request.outcome_includesexplorationcost==0 % Include exploration cost in Outcome magnitude?
            ws.d(ws.d(:,col.Resp1)==3, col.OutcomeAmount)=ws.d(ws.d(:,col.Resp1)==3, col.OutcomeAmount)+2;
        end
        for i=1:36
            ws.d(ws.d(:,1)==i,col.OutcomeMean)=mean(ws.d(ws.d(:,1)==i,col.OutcomeAmount));
            ws.d(ws.d(:,1)==i,col.OutcomeVariance)=var(ws.d(ws.d(:,1)==i,col.OutcomeAmount));
        end
        
        % Perform correlations
        ws.corr=nan*zeros(length(details.designvariables)); ws.corr_pval=nan*zeros(length(details.designvariables));
        for i=1:length(details.designvariables)
            for j=1:length(details.designvariables)
                eval(['[r p]=corr(ws.d(:,col.' details.designvariables{i} '), ws.d(:,col.' details.designvariables{j} '));'])
                ws.corr(i,j)=r;
                ws.corr_pval(i,j)=p;
                if i==j; ws.corr(i,j)=nan; end
            end
        end
        
        %
        subjdata{s,3}.corr=ws.corr;
        subjdata{s,3}.corr_pval=ws.corr_pval;
        ws=[];
    end
    
    % Overall correlation between variables (averaging of subject-level correlations)
    for o1=1:1
        subjdata{logg.n_subjs +1,3}.corr=0; subjdata{logg.n_subjs +1,3}.corr_pval=0;
        for s=1:logg.n_subjs % Mean
            subjdata{logg.n_subjs +1,3}.corr=subjdata{logg.n_subjs +1,3}.corr+subjdata{s,3}.corr;
            subjdata{logg.n_subjs +1,3}.corr_pval=subjdata{logg.n_subjs +1,3}.corr_pval+subjdata{s,3}.corr_pval;
        end
        subjdata{logg.n_subjs +1,3}.corr=(subjdata{logg.n_subjs +1,3}.corr)./logg.n_subjs;
        subjdata{logg.n_subjs +1,3}.corr_pval=subjdata{logg.n_subjs +1,3}.corr_pval/logg.n_subjs;
        r_corrvals=[];
        for i=1:length(details.designvariables)
            r_corrvals=vertcat(r_corrvals,subjdata{logg.n_subjs +1,3}.corr(:,i));
        end
        r_corrvals=r_corrvals(isnan(r_corrvals(:))==0);
        disp(['Mean correlation (rho) between conditions: ' num2str(mean(r_corrvals))])
        disp(['Mean absolute correlation rho: ' num2str(mean(abs(r_corrvals)))])
        disp('(Mean correlations = averaging of subject-level correlations)')
    end
    
    % Plot correlations
    w.ol=700; w.ob=0; w.w=800; w.h=1300; % Size for figure (individual subjects)
    if request.plot_correlations==1
        disp('Design conditions:')
        disp(details.designvariables); disp(' ')
        
        figure('Name', 'Correlation coefficients (absolute)', 'NumberTitle', 'off', 'Position', [w.ol w.ob w.w w.h]);
        % Visual formatting
        side=0.2; x=0.3:side+0.05:1; y=0.005: 1/(logg.n_subjs) : 0.94;  y=sortrows(y',-1);
        if 1/(logg.n_subjs)<side; side=1/logg.n_subjs-0.0065; end;
        %
        k=2;
        for s=1:logg.n_subjs+1
            subplot(logg.n_subjs+1, k, (s-1)*k+1); axis 'off'  % Correlation coefficients (rho) -------
            %         imagesc(abs(subjdata{s,3}.corr)   ,[0 0.5]    ); axis square; colorbar
            imagesc(abs(subjdata{s,3}.corr)   ,[0 1]    ); axis square; colorbar
            set(gca,'YTick', 1:length(details.labels), 'YTickLabel', details.labels)
            set(gca,'XTick', 1:length(details.labels));
            subplot(logg.n_subjs+1, k, (s-1)*k+k); axis 'off' % P values -------------------------------
            imagesc(subjdata{s,3}.corr_pval  ,   [0 0.05]     ); axis square; colorbar
            set(gca,'YTick', 1:length(details.labels), 'YTickLabel', details.labels)
            set(gca,'XTick', 1:length(details.labels));
        end
    end
end

%% Plots RTs

if request.plotrt_mean
    f.subplotcols=5; f.subplot_VerHorz=[0.01 0.1]; f.fig_BotTop=[0.01 0.01]; f.fig_LeftRight=[0.01 0.01]; f.figwidth= 800; f.figheight=1600; 
    if request.plotrt_individuals , figure('Color', 'w', 'Position', [100 50 f.figwidth f.figheight]);  k=1;end
     f.FontSize=15; f.FontName='PT Sans Caption'; 
    
    
    subjdata{logg.n_subjs+1,3}.rt=zeros(6,6);  subjdata{logg.n_subjs+1,3}.rt_a=zeros(6,6);  subjdata{logg.n_subjs+1,3}.rt_r=zeros(6,6);  subjdata{s,3}.rt_e=zeros(6,6); 
    d_rts=nan(logg.n_subjs,4);   ch={'Acc';'Rej';'Expl';'All'}; 
    d_pcho=nan(logg.n_subjs,3);  
    for s=1:logg.n_subjs  % Plot and RECOMPILE
        d_rts(s,:)= [nanmean(subjdata{s,3}.rt_a(:)) nanmean(subjdata{s,3}.rt_r(:)) nanmean(subjdata{s,3}.rt_e(:)) nanmean(subjdata{s,3}.rt(:))];
        subjdata{logg.n_subjs+1,3}.rt=subjdata{logg.n_subjs+1,3}.rt+subjdata{s,3}.rt;
        subjdata{logg.n_subjs+1,3}.rt_a=subjdata{logg.n_subjs+1,3}.rt_a +subjdata{s,3}.rt_a;
        subjdata{logg.n_subjs+1,3}.rt_r=subjdata{logg.n_subjs+1,3}.rt_r +subjdata{s,3}.rt_r;
        subjdata{logg.n_subjs+1,3}.rt_e=  subjdata{logg.n_subjs+1,3}.rt_e +subjdata{s,3}.rt_e;
       
        % Compile pCho data too
     
        
%         % Add to overall, without ruling out cells (i.e. empty = 0
%         ws.m= subjdata{s,3}.rt_a(:);  ws.m(isnan(ws.m))=0; ws.m= reshape(ws.m,6,6);
%         subjdata{logg.n_subjs+1,3}.rt_a=subjdata{logg.n_subjs+1,3}.rt_a +ws.m;
%         ws.m= subjdata{s,3}.rt_r(:);  ws.m(isnan(ws.m))=0; ws.m= reshape(ws.m,6,6);
%         subjdata{logg.n_subjs+1,3}.rt_r=subjdata{logg.n_subjs+1,3}.rt_r + ws.m; 
%         ws.m= subjdata{s,3}.rt_e(:);  ws.m(isnan(ws.m))=0; ws.m= reshape(ws.m,6,6);
%         subjdata{logg.n_subjs+1,3}.rt_e=  subjdata{logg.n_subjs+1,3}.rt_e + ws.m; 
        
        if request.plotrt_individuals
            subtightplot(logg.n_subjs+1,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            text(0.8,0, subjdata{s,1}, 'FontSize', f.FontSize); axis off
        
            % Overall
            subtightplot(logg.n_subjs,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(subjdata{s,3}.rt, 'nancolor', [0.1 0.1 0.1]), colorbar, axis off, axis square
            if s==1; title('RT', 'FontSize', f.FontSize); end
            
            % Accept
            subtightplot(logg.n_subjs,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(subjdata{s,3}.rt_a, 'nancolor', [0.1 0.1 0.1]), colorbar, axis off, axis square
            if s==1; title('RT Accept', 'FontSize', f.FontSize); end
            
            % Reject
            subtightplot(logg.n_subjs,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(subjdata{s,3}.rt_r, 'nancolor', [0.1 0.1 0.1]), colorbar, axis off, axis square
            if s==1; title('RT Reject', 'FontSize', f.FontSize); end
            
            % Explore
            subtightplot(logg.n_subjs,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(subjdata{s,3}.rt_e, 'nancolor', [0.1 0.1 0.1]), colorbar, axis off, axis square
            if s==1; title('RT Explore', 'FontSize', f.FontSize); end
            
        end
    end 
    subjdata{logg.n_subjs+1,3}.rt=subjdata{logg.n_subjs+1,3}.rt./logg.n_subjs;
    subjdata{logg.n_subjs+1,3}.rt_a=subjdata{logg.n_subjs+1,3}.rt_a./logg.n_subjs;
    subjdata{logg.n_subjs+1,3}.rt_r=subjdata{logg.n_subjs+1,3}.rt_r./logg.n_subjs;
    subjdata{logg.n_subjs+1,3}.rt_e=subjdata{logg.n_subjs+1,3}.rt_e./logg.n_subjs;
    
    % Group level 
    figure('Color', 'w', 'Position', [100 50 700 300]);  k=1; s=logg.n_subjs; 
    subtightplot(1,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    text(0.0,0, 'Group', 'FontSize', f.FontSize); axis off; 
    %
    subtightplot(1,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagescnan(subjdata{logg.n_subjs+1,3}.rt, 'nancolor', [0.1 0.1 0.1]), colorbar,   axis square
    title('RT', 'FontSize', f.FontSize); 
    ylabel('EnvThreat','FontSize',f.FontSize, 'FontName', f.FontName);  xlabel('No. Tokens','FontSize',f.FontSize, 'FontName', f.FontName);  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
    set(gca,'FontSize',f.FontSize, 'FontName', f.FontName, 'LineWidth', 0.8);
%     colormap autumn
%         
    
    % Accept
    subtightplot(1,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagescnan(subjdata{logg.n_subjs+1,3}.rt_a, 'nancolor', [0.1 0.1 0.1]), colorbar,  axis square
    title('RT Accept', 'FontSize', f.FontSize); 
    ylabel('EnvThreat','FontSize',f.FontSize, 'FontName', f.FontName);  xlabel('No. Tokens','FontSize',f.FontSize, 'FontName', f.FontName);  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
    set(gca,'FontSize',f.FontSize, 'FontName', f.FontName, 'LineWidth', 0.8);
    % Reject
    subtightplot(1,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagescnan(subjdata{logg.n_subjs+1,3}.rt_r, 'nancolor', [0.1 0.1 0.1]), colorbar,  axis square
    title('RT Reject', 'FontSize', f.FontSize); 
    ylabel('EnvThreat','FontSize',f.FontSize, 'FontName', f.FontName);  xlabel('No. Tokens','FontSize',f.FontSize, 'FontName', f.FontName);  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
    set(gca,'FontSize',f.FontSize, 'FontName', f.FontName, 'LineWidth', 0.8);
    % Explore
    subtightplot(1,f.subplotcols,k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagescnan(subjdata{logg.n_subjs+1,3}.rt_e, 'nancolor', [0.1 0.1 0.1]), colorbar,  axis square
    title('RT Explore', 'FontSize', f.FontSize); 
    ylabel('EnvThreat','FontSize',f.FontSize, 'FontName', f.FontName);  xlabel('No. Tokens','FontSize',f.FontSize, 'FontName', f.FontName);  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
    set(gca,'FontSize',f.FontSize, 'FontName', f.FontName, 'LineWidth', 0.8);
    
%     % Plot overall, split by task 
%     figure('color','w'), f.FontSize=25; f.FontName='PT Sans';
%     barwitherr(reshape(std(d_rts)./sqrt(logg.n_subjs),4,1), reshape(mean(d_rts),4,1))
%     xlabel('Choice','FontSize', f.FontSize, 'FontName', f.FontName), ylabel('Mean RT (ms)','FontSize', f.FontSize, 'FontName', f.FontName)
%     title('Non-optimal choice, V(Chosen)<V(Best Unchosen)','FontSize', f.FontSize, 'FontName', f.FontName), legend('Exp', 'Ctrl')
%     set(gca,'FontSize', f.FontSize, 'FontName', f.FontName)


    % Plot both tasks (combine both into d_rts
    disp('Manual setup!! See code ')  
    
    openvar d_rts
    a=[]; acf=[];  act=[];  openvar acf, act
    d_rts=  [acf act];   % Manually paste cf and ct rts together 
    
     
    
    
    % RTs 
    figure('color','w'), f.FontSize=25; f.FontName='PT Sans';
    barwitherr(reshape(nanstd(d_rts)./sqrt(logg.n_subjs),3,2), reshape(nanmean(d_rts),3,2))
    xlabel('Choice','FontSize', f.FontSize, 'FontName', f.FontName), legend('Ap/Av', 'Ap/Ap')
    set(gca,'FontSize', f.FontSize, 'FontName', f.FontName), set(gca,'xticklabel',{'Accept/No bomb', 'Reject/Bomb', 'Explore'})
    ylabel('Mean RT (ms)','FontSize', f.FontSize, 'FontName', f.FontName)  % RENAME 
colormap summer
    [wf.a]=teg_repeated_measures_ANOVA(d_rts, [2 3], {'Condition', 'Choice'});  % Row=Task, Choice, TxC; Col=F, df1, df2, p
    disp('RTs differs by Task and Choice? ANOVA:')
    
    
    disp([wf.a.labels' num2cell( [cellfun(@(x)f_sig(x), num2cell(wf.a.R(:, 4))) wf.a.R(:, 4)] ) ]) 
    for c=1:4  % Stats: compare cF and ct
        [wc.h wc.p wc.ci wc.stats]= ttest(d_rts(:,c)-d_rts(:, c+4));
        disp([ ch{c} ' : t(' num2str(wc.stats.df) ')=' num2str(wc.stats.tstat,3) '  ,p=' num2str(wc.p,3)]);
    end
    
% Plot pCho 

    
    
end


%% Save ?
request.saveres=0;
if request.saveres==1
    details=logg.subjects;
    save(['res_' date '_plotbehaviour'],  'plots', 'details','subjdata', 'r_corrvals')
end