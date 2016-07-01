% Free exploration of behaviour 
clear all; close all hidden; clc



req.subjects= {};

% Inzone subjects
% req.subjects= {'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};  % Inzone6 subjects 
 

req.onlythird=0;   %  1/2/3, ohterwise ignore 
req.onlyfirsthalf=0;
req.onlylatterhalf=0; 

for o1=1:1  % Settings
    
    % Settings that don't change much 
    req.taskfile='taskfMRI'; 
    
    
    %  Anxiety scores (State, trait)    
    req.pers={'p01_GV',7,14;'p02_YY',0,6;'p04_MP',3,8;'p06_KB',9,15;'p08_SG',2,8;'p10_RC',1,4;'p13_HL',6,8;'p15_SH',4,8;'p17_SJ',1,8;'p18_MS',30,23;'p21_ES',13,26;'p23_BS',12,16;'p25_RJ',16,10;'p27_DM',0,7;'p30_KL',4,16;'p34_TB',11,9;'p35_SM',15,20;'p36_FR',15,21;'p38_MK',3,10;'p41_AC',13,20;}; 
    
    % Useful functions 
     fmat= @(mat, index)mat(index);
    fvect= @(mat)mat(:);
    f_sig=@(p)fmat(sortrows(repmat([1 0.5 zeros(1,18)]',500,1), -1), find(ceil(p*100)/100<=linspace(0.001, 1, length(repmat([1 0.5 zeros(1,18)]',500,1))),1,'first')-1);   % Convert p value to significane marker: 1 (p<=0.05), 0.5 (p<.1), 0. p values below 0.01 are treated as 0.01 for the sake of indexing
  
    % Where 
    req.inzone6= {'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};     
    w.w=pwd; if strcmp(w.w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis';  where.modelscripts=[where.where '/4 Fit computational models'];  where.where2='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI';  where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI subjdata/1 Brain data';
    else  where.where='C:\Users\e.loh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis'; where.where2='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI';  where.data_brain='G:\2 [Explore]\1 Brain data'; end
    where.scripts=[where.where fs '6 Adhoc Behavioural scripts'];  where.modelscripts=[where.where fs '4 Fit computational models'];
    where.data_beh= [where.where2 filesep '1 Behavioural data'];   path(pathdef); addpath(where.modelscripts)    
    cd(where.data_beh), w.s=dir('p*');  logg.subjects=cellstr(char(w.s.name)); logg.n_subjs=length(logg.subjects); cd(where.scripts)
    [logg.subjects logg.n_subjs w.w] = f_selectsubjects([], req.subjects,  [{'Subs' 'All'}; [logg.subjects num2cell(ones(logg.n_subjs,1))]], 'All');  % Subject selection
    
    
    
    disp(req)
    disp('---------------------------------------------------------------------------------')
end

%%  Load data for all subjects
%   subjdata: col 1=subject, col 2= conflict , col 3= control 

subjdata=cell(logg.n_subjs+1,2);   d_anxiety=nan(logg.n_subjs,2);
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
    col.NTokens=31;    
    col.pLoss = 32;
    col.Entropy= 33;
    col.EV=34;
end
for  s=1: logg.n_subjs
    for task=1:2
        ws.w=load([where.data_beh filesep logg.subjects{s} filesep logg.subjects{s} '_file_' req.taskfile]);
        switch req.taskfile
            case 'taskfMRI'   % fMRI task (cF and ct)
                switch task
                    case 1; ws.d.subjdata=ws.w.conflict; % Conflict
                    case 2; ws.d.subjdata=ws.w.control; % Control
                end
                ws.d.settings=ws.w.settings;
            case 'taskprac'  % Practice sessions for fMRI screening  (cF and ct)
                switch task
                    case 1; ws.d.subjdata=ws.w.conflict; % Conflict
                    case 2; ws.d.subjdata=ws.w.control; % Control
                end
                ws.d.settings=ws.w.settings;
                %         case  'taskconflictcontrol'  % Integrated conflict & control
                %             ws.d=ws.w.taskconflictcontrol;
                %             switch task
                %                 case 1; ws.d.subjdata=ws.d.conflictsubjdata; % Conflict
                %                 case 2; ws.d.subjdata=ws.d.controlsubjdata; % Control
                %             end
                %         case 'taskcontrol', ws.d=ws.w.taskcontrol; % Control only
                %         case 'taskconflict', eval(['ws.d=ws.w.' req.taskfile ';'])% Conflict only
            otherwise, error('Check sessions!');
        end
        if req.onlylatterhalf==1 % Latter half of session only (each task separtely)
            ws.d.subjdata= sortrows(ws.d.subjdata, col.Trialnum);
            ws.d.subjdata= ws.d.subjdata(round(size(ws.d.subjdata,1)/2)+1:end, :);
        elseif req.onlyfirsthalf==1 % First half 
            ws.d.subjdata= sortrows(ws.d.subjdata, -col.Trialnum);
            ws.d.subjdata= ws.d.subjdata(round(size(ws.d.subjdata,1)/2)+1:end, :);
        end
        
        if req.onlythird==1 % First third
            ws.d.subjdata= sortrows(ws.d.subjdata, col.Trialnum);
            ws.tn = 1: round(size(ws.d.subjdata,1)/3);
            ws.d.subjdata= ws.d.subjdata(ws.tn, :);
        elseif req.onlythird==2 % Second third
            ws.d.subjdata= sortrows(ws.d.subjdata, col.Trialnum);
            ws.tn = round(size(ws.d.subjdata,1)/3)+1:round(size(ws.d.subjdata,1)*2/3);
            ws.d.subjdata= ws.d.subjdata(ws.tn, :);
        elseif req.onlythird==3 % Last third
            ws.d.subjdata= sortrows(ws.d.subjdata, col.Trialnum);
            ws.tn = round(size(ws.d.subjdata,1)*2/3)+1: size(ws.d.subjdata,1);
            ws.d.subjdata= ws.d.subjdata(ws.tn, :);
        end
        
        % Sub-sample InZone
%         if s==1; disp('SUBSAMPLING FOR INNER TRIALS!'); end; ws.d.subjdata(ws.d.subjdata(:, col.EnvThreat)>4.5, col.TrialValid)=0;  ws.d.subjdata(ws.d.subjdata(:, col.EnvThreat)<2.5, col.TrialValid)=0;  ws.d.subjdata(ws.d.subjdata(:, col.NTokPairs)<2.5, col.TrialValid)=0; 
         
         % Fill in parameters 
        ws.d.subjdata(:, col.NTokens)= ws.d.subjdata(:, col.NTokPairs)*2; 
        switch task
            case 1; ws.d.subjdata= fpar_conflict( ws.d.subjdata, col);
            case 2; ws.d.subjdata= fpar_control( ws.d.subjdata, col);
        end
         
        % RT
        ws.d.subjdata(ws.d.subjdata(:, col.RT1)<0)=nan;
        ws.meanrt =  nanmean(ws.d.subjdata(:,col.RT1));
        ws.d.subjdata(:, col.RT1) = log(ws.d.subjdata(:, col.RT1));  if s==1 && task==1, disp('RTs logg'), end % Log RTs by subject
        % %     ws.d.subjdata(:, col.RT1) =wr.rts(wr.rts(:,1)==task,2);  if s==1, disp('RTs z scored'), end
        %     ws.d.subjdata(:, col.RT1) = ws.d.subjdata(:, col.RT1)-ws.meanrt;  if s==1, disp('RTs mean centred'), end % Mean centre RTs by subject
        
        subjdata{s,1}=logg.subjects{s};
        subjdata{s,1+ task}= ws.d.subjdata(ws.d.subjdata(:,col.TrialValid)==1, :); 
        subjdata{logg.n_subjs+1,1+ task}=[subjdata{logg.n_subjs+1,1+ task}; subjdata{s,1+ task}];
        ws=[];
    end
    
    d_anxiety(s,:)= cell2mat(req.pers(strcmp(req.pers(:,1), logg.subjects{s}), 2:3)); 
end
 
%% Manual calcs 

% Col + other settings
for o=1:1
    logg.chocomp={'AR' [1 2]; 'AE' [1 3]; 'RE' [2 3]};  % Simple effects for choice 
    col.choh=  size(subjdata{s,2+1},2 )+1; % choice entropy 
end 
        
% {1} cf means {2} ct means {3} cf > ct 
d_choh=repmat({nan(logg.n_subjs, 36)}, 1,3);   % EnvThreat x NTok in a vector (e1n1, e1n2..e2n1..)
d_pchomat=repmat({nan(logg.n_subjs, 36)}, 3,3);   % p(Choice): row= choice, col= cF, ct, cF>ct. Within each cell: EnvThreat x NTok in a vector (e1n1, e1n2..e2n1..)
d_pcho= repmat({nan(logg.n_subjs, 3)}, 1,3);  d_rt = repmat({nan(logg.n_subjs, 3)}, 1,3); 
d_chohxcho = {nan(logg.n_subjs, 6)   nan(logg.n_subjs, 6)}; % Choice entropy x Choice: cf-A, cf-R, cf-e, ct-A, ct-R, ct-e 
d_rhchoh=  nan(logg.n_subjs, 3); % r choh & h: cf, ct, cf>ct  (Correlation of choice and objective entropy)
for s=1:logg.n_subjs
    for tk=1:2
        ws.d=subjdata{s,tk+1}; 
         
        % Split by cell 
        for e=1:6
            for n=1:6
                ws.c =   ws.d(ws.d(:, col.EnvThreat)==e/6 & ws.d(:, col.NTokPairs)==n, :);
                for ch=1:3
                    ws.p(ch) =  mean( ws.c(:, col.Choice)==ch); 
                    d_choh{tk,ch}(s,(e-1)*6+n) =  ws.p(ch);  
                    
                    if ws.p(ch) ==1, ws.p(ch) =ws.p(ch) -eps; elseif ws.p(ch) ==0, ws.p(ch)=ws.p(ch) +eps; end 
                end
                 ws.d(ws.d(:, col.EnvThreat)==e/6 & ws.d(:, col.NTokPairs)==n, col.choh) =  sum ( (-ws.p).*log(ws.p) ) ; 
                
                 % Record in data structures 
                 d_choh{tk}(s,(e-1)*6+n) =  sum ((-ws.p).*log(ws.p) ); 
                  
            end
        end
        
        % Split by choice 
        for ch=1:3
            d_chohxcho{1}(s, (tk-1)*3+ch) =  mean( ws.d(ws.d(:, col.Choice)==ch, col.choh)); 
            d_pcho{tk}(s,ch)= mean( ws.d(:, col.Choice)==ch); 
            d_rt{tk}(s,ch) = mean( ws.d(ws.d(:, col.Choice)==ch, col.RT1)); 
            
        end
        
        % Compare across choice, within task 
        for cp=1:size(logg.chocomp,1)
            d_chohxcho{2}(s, (tk-1)*3+cp) =  d_chohxcho{1}(s, (tk-1)*3+logg.chocomp{cp,2}(1)) - d_chohxcho{1}(s, (tk-1)*3+logg.chocomp{cp,2}(2)) ; 
        end 
        
        % Choice entropy correlated w objective entropy? 
        d_rhchoh(s,tk) = corr( ws.d(:, col.Entropy), ws.d(:, col.choh)); 
        
        ws=[]; 
    end
end
d_choh{3}= d_choh{1} - d_choh{2}; 
d_chohxcho{2} =[-d_chohxcho{2}(:,1:3) nan(logg.n_subjs,1)  -d_chohxcho{2}(:,4:6 ) nan(logg.n_subjs,1) (-d_chohxcho{2}(:,1:3))-(-d_chohxcho{2}(:,4:6 )) nan(logg.n_subjs,1)   (d_chohxcho{1}(:,1:3)-d_chohxcho{1}(:,4:6))] ;   % comparisons are flipped relative to what they usually are because this is easier to interpret 
d_rhchoh(:, 3)= d_rhchoh(:, 1)-d_rhchoh(:, 2);
d_pcho{3}= d_pcho{1}-d_pcho{2}; 
d_rt{3}= d_rt{2}-d_rt{1}; 


%% Plot  + analyse  

for o=1:1 % Basic plots (overall pCho, RTs)
    figure('color','w','Position',[100, 0, 600 800]), 
    f.FontSize=20; f.FontName='PT Sans'; colormap default
    ch = {'Accept/No bomb', 'Reject/Bomb', 'Explore'}; 
    % RTs 
    subplot(2,1,1)
    wf.d = [d_pcho{1} d_pcho{2}]; 
    barwitherr(reshape(nanstd(wf.d)./sqrt(logg.n_subjs),3,2), reshape(nanmean(wf.d),3,2))
    xlabel('Choice','FontSize', f.FontSize, 'FontName', f.FontName), legend('Ap/Av', 'Ap/Ap')
    set(gca,'FontSize', f.FontSize, 'FontName', f.FontName), set(gca,'xticklabel',ch)
    ylabel('% Chosen','FontSize', f.FontSize, 'FontName', f.FontName)  
    [wf.a]=teg_repeated_measures_ANOVA(wf.d, [2 3], {'Condition', 'Choice'});  % Row=Task, Choice, TxC; Col=F, df1, df2, p
    disp('Choice differs by Task and Choice? ANOVA:'),  disp([wf.a.labels' num2cell(  wf.a.R(:,1:4)  ) ]) 
    for c=1:3 % Stats: compare cF and ct
        [wc.h wc.p wc.ci wc.stats]= ttest(wf.d(:,c)-wf.d(:, c+3));
        disp([ ch{c} ' : t(' num2str(wc.stats.df) ')=' num2str(wc.stats.tstat,3) '  ,p=' num2str(wc.p,3)]);
    end
    
    % RTs 
    subplot(2,1,2)
    wf.d = [d_rt{1} d_rt{2}]; 
    barwitherr(reshape(nanstd(wf.d)./sqrt(logg.n_subjs),3,2), reshape(nanmean(wf.d),3,2))
    xlabel('Choice','FontSize', f.FontSize, 'FontName', f.FontName), legend('Ap/Av', 'Ap/Ap')
    set(gca,'FontSize', f.FontSize, 'FontName', f.FontName), set(gca,'xticklabel',{'Accept/No bomb', 'Reject/Bomb', 'Explore'})
    ylabel('Mean log(RT)','FontSize', f.FontSize, 'FontName', f.FontName)  
    [wf.a]=teg_repeated_measures_ANOVA(wf.d, [2 3], {'Condition', 'Choice'});  % Row=Task, Choice, TxC; Col=F, df1, df2, p
    disp('RTs differs by Task and Choice? ANOVA:'),disp([wf.a.labels' num2cell(  wf.a.R(:,1:4)  ) ]) 
    for c=1:3  % Stats: compare cF and ct
        [wc.h wc.p wc.ci wc.stats]= ttest(wf.d(:,c)-wf.d(:, c+3));
        disp([ ch{c} ' : t(' num2str(wc.stats.df) ')=' num2str(wc.stats.tstat,3) '  ,p=' num2str(wc.p,3)]);
    end
    ylim([6 7.5])
    
    
    

end 
for o=1:1 % Choice entropy 

    f.plotcols= 3; f.plotrows=4;  f.markersize=4;
    f.figwidth= 900; f.figheight=600; f.fontsize=40;  f.subplot_VerHorz=[0.15 0.08]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05];
    figure('Name', 'Choice Entropy', 'Position', [2000 350 f.figwidth f.figheight], 'Color', 'w','number', 'off'); k=1;

    % Overall choice entropy
    wf.d= d_chohxcho{1};   wf.cluster_npercluster=3;  wf.cl_shiftx=0.23;
    wf.sc_x=  sortrows(repmat([(1:size(wf.d,2)/wf.cluster_npercluster)-wf.cl_shiftx 1:size(wf.d,2)/wf.cluster_npercluster (1:size(wf.d,2)/wf.cluster_npercluster)+wf.cl_shiftx]', logg.n_subjs,1));
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
   barwitherr(reshape( nanstd(wf.d )/sqrt(logg.n_subjs), 3,2)',  reshape(nanmean(wf.d ), 3,2)'), colormap 'summer', 
%     hold on,   scatter(wf.sc_x , wf.d(:),f.markersize)
    ylabel('Mean choice entropy','FontSize', f.fontsize), title('Choice entropy by choice','FontSize', f.fontsize); set(gca, 'xticklabel', {'Ap/Av', 'Ap/Ap'},'FontSize', f.fontsize)
    wf.x_markt = sortrows(unique(wf.sc_x));legend('Accept/No Bomb', 'Reject/Bomb', 'Explore')
    for i=1:size(wf.d,2) % Manual mark sig for clustered graph
        %             disp(mean(wf.d(:,i)-(wf.d(:,i)*0+wf.null)))
        %             [tstat(i), pvals(i)]= ttest(wf.d(:,i)-(wf.d(:,i)*0+wf.null));
        %             if pvals(i)<0.1 &&  pvals(i)>0.05, hold on, scatter(wf.x_markt(i), 0.15, '.', 'MarkerEdgeColor', 'r')
        %             elseif  pvals(i)<0.05, hold on, scatter(wf.x_markt(i), 0.15, '+', 'MarkerEdgeColor', 'r')
        %             end
    end
    [wf.a]=teg_repeated_measures_ANOVA(wf.d, [2 3], {'Condition', 'Choice'});  % Row=Task, Choice, TxC; Col=F, df1, df2, p
    disp('ChoiceH differs by Task and Choice? ANOVA:'),     disp([wf.a.labels' num2cell( [cellfun(@(x)f_sig(x), num2cell(wf.a.R(:, 4))) wf.a.R(:, 4)] ) ])
     
    % Choice entropy - simple effects 
    wf.d= d_chohxcho{2};  wf.null=0; 
    subtightplot(f.plotrows, f.plotcols , k:k+1 ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+2;
    barwitherr(nanstd(wf.d )/sqrt(logg.n_subjs),  nanmean(wf.d ), 'y'),hold on,  scatter( fvect(repmat(1:size(wf.d ,2), logg.n_subjs,1)), wf.d(:),f.markersize)
    ylabel(sprintf('Difference in choice entropy\n (Opt 1>Opt2)')), title('Choice entropy by choice - simple fx');
    set(gca, 'xticklabel', [cellfun(@(x)fliplr(x), logg.chocomp(:,1),'UniformOutput',0)'  {' ' } cellfun(@(x)fliplr(x), logg.chocomp(:,1),'UniformOutput',0)'  {' ' } cellfun(@(x)fliplr(x), logg.chocomp(:,1),'UniformOutput',0)'  {' '  'A', 'R', 'E'}] ),  xlim([0 size(wf.d,2)+1])
    [tstat, pvals]= f_markfigstat_1samt(wf.d, -0.3, 'r');  % data, null, marker y-cord, color
    xlabel('Ap/Av     Ap/Ap           Ap/Av > Ap/Ap        '),  ylim([-0.5 0.8]) 
%     [h p ci st] = ttest(wf.d); p'    
    
    % Heat plots of choice entropy    - Ap/Av
    wf.d= flipud( reshape( mean(  d_choh{1}  ),  6,6)');
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagesc( wf.d), colorbar, axis square, colormap default 
    title('[Ap/Av] Choice entropy over task space');  ylabel('EnvThreat'), xlabel('No. Tokens')
    set(gca, 'xtick', 1:6, 'xticklabel', 2:2:12, 'ytick', 1:6, 'yticklabel', {'1/6', '2/6', '3/6', '4/6', '5/6', '6/6'})
    
    % Heat plots of choice entropy    - Ap/Ap
    wf.d=  flipud( reshape( mean(  d_choh{2}  ),  6,6)');
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagesc( wf.d), colorbar, axis square 
    title('[Ap/Ap] Choice entropy over task space');  ylabel('EnvThreat'), xlabel('No. Tokens')
    set(gca, 'xtick', 1:6, 'xticklabel', 2:2:12, 'ytick', 1:6, 'yticklabel', {'1/6', '2/6', '3/6', '4/6', '5/6', '6/6'})
    
    % Heat plots of choice entropy    - Ap/Av > Ap/Av
    wf.d= flipud( reshape( mean(  d_choh{3}  ),  6,6)'); 
    [wf.h wf.p wf.ci  wf.st]= ttest(d_choh{3}) ;
     wf.d=  flipud( reshape( wf.st.tstat,  6,6)');   % t statistic  
%      wf.d=  flipud( reshape( wf.h,  6,6)');   % is t-stat significant? 
%      wf.d=  flipud( reshape( wf.p,  6,6)');   % p value  
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagesc( wf.d), colorbar, axis square 
    title('[Ap/Av  > Ap/Ap] Choice entropy over task space (see code for what is being plotted)');  ylabel('EnvThreat'), xlabel('No. Tokens')
    set(gca, 'xtick', 1:6, 'xticklabel', 2:2:12, 'ytick', 1:6, 'yticklabel', {'1/6', '2/6', '3/6', '4/6', '5/6', '6/6'})
      
    
    % NEG Heat plots of choice entropy    - Ap/Av
    wf.d= - flipud( reshape( mean(  d_choh{1}  ),  6,6)');
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagesc( wf.d), colorbar, axis square, colormap default 
    title('[Ap/Av] Negative choice entropy over task space');  ylabel('EnvThreat'), xlabel('No. Tokens')
    set(gca, 'xtick', 1:6, 'xticklabel', 2:2:12, 'ytick', 1:6, 'yticklabel', {'1/6', '2/6', '3/6', '4/6', '5/6', '6/6'})
    
    % NEG Heat plots of choice entropy    - Ap/Ap
    wf.d= - flipud( reshape( mean(  d_choh{2}  ),  6,6)');
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagesc( wf.d), colorbar, axis square 
    title('[Ap/Ap] Negative choice entropy over task space');  ylabel('EnvThreat'), xlabel('No. Tokens')
    set(gca, 'xtick', 1:6, 'xticklabel', 2:2:12, 'ytick', 1:6, 'yticklabel', {'1/6', '2/6', '3/6', '4/6', '5/6', '6/6'})
    
    % NEG Heat plots of choice entropy    - Ap/Av > Ap/Av
    wf.d=  flipud( reshape( mean(  -d_choh{1}-(-d_choh{2})   ),  6,6)');
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    imagesc( wf.d), colorbar, axis square 
    title('[Ap/Av  > Ap/Ap] Negative choice entropy over task space');  ylabel('EnvThreat'), xlabel('No. Tokens')
    set(gca, 'xtick', 1:6, 'xticklabel', 2:2:12, 'ytick', 1:6, 'yticklabel', {'1/6', '2/6', '3/6', '4/6', '5/6', '6/6'})
    
    % Do choice entropy and objective entropy correlate? 
    wf.d=  d_rhchoh ;  wf.null=0; 
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    barwitherr(nanstd(wf.d )/sqrt(logg.n_subjs),  nanmean(wf.d ), 'y'),hold on,  scatter( fvect(repmat(1:size(wf.d ,2), logg.n_subjs,1)), wf.d(:),f.markersize)
    ylabel(sprintf('Mean r statistic')), title('Correlation of entropy with choice entropy');
    set(gca, 'xticklabel', {'Ap/Av','Ap/Ap', 'Ap/Av > Ap/Ap'}),  xlim([0 size(wf.d,2)+1])
    [tstat, pvals]= f_markfigstat_1samt(wf.d, -0.3, 'r');  % data, null, marker y-cord, color
%     xlabel('Ap/Av     Ap/Ap           Ap/Av > Ap/Ap        '),  ylim([-0.5 0.8]) 
    
    
end
for o=1:1 % Misc choice analysis  
    f.plotcols= 2; f.plotrows=1;  f.markersize=4;
    f.figwidth= 900; f.figheight=600; f.fontsize=20;  f.subplot_VerHorz=[0.15 0.08]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05];
%     figure('Name', 'Misc choice analysis', 'Position', [2000 350 f.figwidth f.figheight], 'Color', 'w','number', 'off'); k=1;
 figure('color','w')
    
    % Overall % Choice proportion 
    wf.d=  [d_pcho{1} d_pcho{2}];  
%     subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    barwitherr(reshape( std(wf.d )/sqrt(logg.n_subjs), 3,2),  reshape(nanmean(wf.d ), 3,2)), colormap default
    ylabel('p(Choice)', 'FontSize', f.fontsize),legend('Ap/Av', 'Ap/Ap'),  title('Overall Choice %', 'FontSize', f.fontsize); set(gca, 'xticklabel', {'Accept/No Bomb', 'Reject/Bomb', 'Explore'}, 'FontSize', f.fontsize); 
    colormap 'summer'
%     wf.cluster_npercluster=2;  wf.cl_shiftx=0.23;    wf.sc_x=  sortrows(repmat([(1:size(wf.d,2)/wf.cluster_npercluster)-wf.cl_shiftx 1:size(wf.d,2)/wf.cluster_npercluster (1:size(wf.d,2)/wf.cluster_npercluster)+wf.cl_shiftx]', logg.n_subjs,1));
%     wf.x_markt = sortrows(unique(wf.sc_x)); hold on,   scatter(wf.sc_x , wf.d(:),f.markersize)
%     for i=1:size(wf.d,2) % Manual mark sig for clustered graph
%         %             disp(mean(wf.d(:,i)-(wf.d(:,i)*0+wf.null)))
%         %             [tstat(i), pvals(i)]= ttest(wf.d(:,i)-(wf.d(:,i)*0+wf.null));
%         %             if pvals(i)<0.1 &&  pvals(i)>0.05, hold on, scatter(wf.x_markt(i), 0.15, '.', 'MarkerEdgeColor', 'r')
%         %             elseif  pvals(i)<0.05, hold on, scatter(wf.x_markt(i), 0.15, '+', 'MarkerEdgeColor', 'r')
%         %             end
%     end
%     [wf.a]=teg_repeated_measures_ANOVA(wf.d, [2 3], {'Condition', 'Choice'});  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%     disp('ChoiceH differs by Task and Choice? ANOVA:'), disp([wf.a.labels' num2cell( [cellfun(@(x)f_sig(x), num2cell(wf.a.R(:, 4))) wf.a.R(:, 4)] ) ])
%   
end

for o=1:1 % Correlations w anxiety
    close all hidden
    
%     wf.corrcommand =  '[wf.r wf.p]=  corr(wf.d(:,1), wf.d(:,2));';   % Pearson (parametric) correlation 
    wf.corrcommand =  '[wf.r wf.p]=  corr(wf.d(:,1), wf.d(:,2),''type'', ''Spearman'')';   disp('Non par correlations!') %  Nonpar corr
    
    
    f.plotcols= 2; f.plotrows=2;  f.markersize=20;
    f.figwidth= 900; f.figheight=600; f.fontsize=20;  f.subplot_VerHorz=[0.2 0.15]; f.fig_BotTop=[0.1 0.1]; f.fig_LeftRight=[0.1 0.05];
    figure('Name', 'Correlations', 'Position', [2000 350 f.figwidth f.figheight], 'Color', 'w','number', 'off'); k=1;

    
    % Anxiety correlated w RT Anxiety cF>ct?
    wf.d=  [d_anxiety(:,2) d_rt{3}(:,1)] ; % Train anxiety
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    scatter(wf.d(:,1), wf.d(:,2),f.markersize); lsline;  eval(wf.corrcommand) 
    xlabel('Trait anxiety','FontSize', f.fontsize),  ylabel(sprintf('RT difference\n(Ap/Ap No Bomb > Ap/Av Accept)'),'FontSize', f.fontsize),
    title(['r=' num2str(wf.r) ',  p= ' num2str(wf.p)],'FontSize', f.fontsize);  set(gca, 'FontSize', f.fontsize)
    wf.d=  [d_anxiety(:,1) d_rt{3}(:,1)] ;   % State anxiety
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    scatter(wf.d(:,1), wf.d(:,2),f.markersize); lsline;  eval(wf.corrcommand) 
    xlabel('State anxiety','FontSize', f.fontsize),  ylabel(sprintf('RT difference\n(Ap/Ap No Bomb > Ap/Av Accept)'),'FontSize', f.fontsize),
    title(['r=' num2str(wf.r) ',  p= ' num2str(wf.p)],'FontSize', f.fontsize);  set(gca, 'FontSize', f.fontsize)
    
    % Choice adjustment (Accept cF>ct) correlated w RT Anxiety cF>ct?
    wf.d=  [ d_pcho{3}(:,1) d_rt{3}(:,1)] ;  
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    scatter(wf.d(:,1), wf.d(:,2),f.markersize); lsline;  eval(wf.corrcommand) 
    xlabel('% Ap/Av Accept > % Ap/Ap No Bomb)','FontSize', f.fontsize),  ylabel(sprintf('RT difference\n(Ap/Ap No Bomb > Ap/Av Accept)'),'FontSize', f.fontsize) 
    title(['r=' num2str(wf.r) ',  p= ' num2str(wf.p)],'FontSize', f.fontsize);  set(gca, 'FontSize', f.fontsize)
    
    % Choice adjustment (Reject cF>ct) correlated w RT Anxiety cF>ct?
    wf.d=  [ d_pcho{3}(:,2) d_rt{3}(:,1)] ;  
    subtightplot(f.plotrows, f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    scatter(wf.d(:,1), wf.d(:,2),f.markersize); lsline;  eval(wf.corrcommand) 
    xlabel('% Ap/Av Reject > % Ap/Ap Bomb)','FontSize', f.fontsize),  ylabel(sprintf('RT difference\n(Ap/Ap No Bomb > Ap/Av Accept)'),'FontSize', f.fontsize) 
    title(['r=' num2str(wf.r) ',  p= ' num2str(wf.p)],'FontSize', f.fontsize);  set(gca, 'FontSize', f.fontsize)
    
end


[d_pcho{1} d_pcho{2}]
openvar ans

