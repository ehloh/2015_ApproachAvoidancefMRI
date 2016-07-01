% s_checkcollinearity
clear all, close all hidden, clc

% Requested analysis
log.specificsubjects={}; 
% log.specificsubjects={'p01_CW'}; 

log.plot_ranger=[];
% log.plot_ranger=[0 0.5];  % Plot range of correlations
 

for o1=1:1 % General settings and specifications


    % Subjects
    where.where='D:\Dropbox\SCRIPPS\5 Explore fMRI';   where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data'; % where.data_beh=[where.where filesep '2 Behavioural data'];
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
%     w.onsets=load([where.data_brain filesep 'Onsetslog_m_c1_Contextall_Hit' ]); % Apply subject filter first according to validity of onsets file
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog,log.specificsubjects, [log.w.datalog  [{'OK'};  num2cell(ones(size(log.w.datalog,1)-1,1))]  ], 'OK');
    
end


%% For mixed Context + Item models (where both context mem and item mem are represented)

dothis=1;
if dothis
    
    
    req.testwhat='all'; 
    req.testwhat='choice2x3'; 
    
    
    % % WHICH MODEL?
    log.FLthread='2 First level s4Ants';
    log.FLmod='m_c13_ChoiceFull_ULPEN';

    % Fetch design matrices
    d_mat=cell(log.n_subjs,2);  % Design matrix, Correlation matrix
    for s=1:log.n_subjs
        ws.s=load([where.data_brain filesep log.subjects{s} filesep log.FLthread filesep log.FLmod '_Basic Contrasted\SPM.mat']);
%         ws.s=load([where.data_brain filesep log.subjects{s} filesep log.FLthread filesep log.FLmod ' Estimated\SPM_' log.FLmod '.mat']); disp('loadig wrong thing !')
        

        if s==1; log.RegNames= ws.s.SPM.xX.name(1:3:end)'; log.RegNames = log.RegNames(1:find(strcmp(log.RegNames, 'Sn(1) R1'))-1); log.n_regs= length(log.RegNames);  end
        d_mat{s,1} = ws.s.SPM.xX.X(:, 1:3:end);
        d_mat{s,1} = d_mat{s,1}(:, 1: log.n_regs) ;
%         log.RegNames
%         error
        
        % Which regressors to include? 
        if strcmp(req.testwhat,'choice2x3')==1, 
            d_mat{s,1}= d_mat{s,1}(:,  [find(strcmp(log.RegNames, 'Sn(1) cF_Accept*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) cF_Reject*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) cF_Explore*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) ct_NoBomb*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) ct_Bomb*bf(1)'))  find(strcmp(log.RegNames, 'Sn(1) ct_Explore*bf(1)'))  ]);
            ws.d= d_mat{s,1}; 
            ws.d=  [ws.d sum(ws.d(:,2:6),2) sum(ws.d(:,[1   3 4 5 6]),2) sum(ws.d(:,[1 2   4 5 6]),2) sum(ws.d(:,[1 2 3   5 6]),2) sum(ws.d(:,[1 2 3 4  6]),2) sum(ws.d(:,1:5),2)];
            d_mat{s,1} = ws.d; 
        end 
        
        ws=[];
    end
    
    % PLOT all
    if strcmp(req.testwhat,'choice2x3')==1,  
        cells={ 'cF_Acc' 'cF_Rej' 'cF_Exp' 'ct_Acc' 'ct_Rej' 'ct_Exp'   'xcF_Acc' 'xcF_Rej' 'xcF_Exp' 'xct_Acc' 'xct_Rej' 'xct_Exp'}; 

    elseif strcmp(req.testwhat,'all')==1, disp('cell names not specified!') 
    end
    if exist('cells','var')==0; cells = cellfun(@(x)x(7:end-6), log.RegNames, 'UniformOutput',0);  end 
    
    d_corrs=nan(log.n_subjs,5);  % Mean, min, max, absmean, absmax
    close all hidden, f.plotcols=6;  f.figwidth= 2400; f.figheight=800; f.fontsize=10; f.fontname='PT Sans Caption';  f.subplot_VerHorz=[0.05 0.03]; f.fig_BotTop=[0.05 0.15]; f.fig_LeftRight=[0.05 0.05];
%     figure('Name', 'Reg Collinearity', 'NumberTitle', 'off', 'Position', [100 50 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for s=1:log.n_subjs
        
        % Correlate
        [ws.r ws.p]=corr(d_mat{s,1});
        d_mat{s,2} =ws.r ; ws.r(ws.r==1)=nan;
        d_corrs(s,:)= [nanmean(ws.r(:)) min(ws.r(:)) max(ws.r(:)) nanmean(abs(ws.r(:))) max(abs(ws.r(:)))];
        d_outcomecorr(s,1:2)=  [min(ws.r(:)) max(ws.r(:))]; % Greatest anticorr and corr 
%         ws.r(4:end, 1:3)=nan; % Blot out cross-cell correlations (i.e. SR vs SN)
%         ws.r([1:3 7:end], 4:6)=nan;
%         ws.r([1:6 10:end], 7:9)=nan;
%         ws.r(1:9, 10:12)=nan;
%         d_outcomecorr(s,1:2)=  [min(ws.r(:)) max(ws.r(:))]; % Greatest anticorr and corr 
        
        
        % PLOT
%         subtightplot(ceil((log.n_subjs)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
%         imagescnan(d_mat{s,2}, 'nancolor',[0.1 0.1 0.1]) % , axis square   %,   title(log.subjects{s})
%         colorbar
%         set(gca, 'FontSize',10, 'XTick', 1:length(cells), 'XTickLabel', cells, 'YTick', 1:length(cells), 'YTickLabel', cells);  % xlim([0.3 2.7])

        
        
%         if (k-1)/f.plotcols==floor((k-1)/f.plotcols), colorbar, end
%         if ~isempty(log.plot_ranger), caxis(log.plot_ranger), end
    end
    % disp([num2cell(1:length(log.RegNames))' log.RegNames'] )
%     mean(d_corrs)
    
    figure('color','w'), subplot(1,2,1), barwitherr(std(d_corrs)./sqrt(log.n_subjs), mean(d_corrs), 'y'), ylabel('r correlation',  'FontSize', 20), set(gca,  'FontSize', 20, 'XTick', 1:5, 'XTickLabel', {'Mean' 'min' 'max' 'Abs mean' 'Abs Max'})
    d_mat{log.n_subjs+1,2}    =zeros(size(d_mat{1,2}));
    for s=1:log.n_subjs, d_mat{log.n_subjs+1,2}=d_mat{log.n_subjs+1,2}+ d_mat{s,2}    ; end
    d_mat{log.n_subjs+1,2}=d_mat{log.n_subjs+1,2}./log.n_subjs;
    subplot(1,2,2),
    imagescnan(d_mat{log.n_subjs+1,2}, 'nancolor',[0.1 0.1 0.1]), colorbar, axis square   %,   title(log.subjects{s})
    set(gca, 'FontSize',10, 'XTick', 1:length(cells), 'XTickLabel', cells, 'YTick', 1:length(cells), 'YTickLabel', cells);  % xlim([0.3 2.7])
    xticklabel_rotate
    
end

% title('Context event Context memory Object (recognized)', 'FontSize',20);
% % axis off
% % subplot(1,2,2)
% % % 
% a={'Context event';'Context memory';'Recognized object'}
% a= [a;a;a;a]
% 
% set(gca, 'FontSize',17, 'XTick', 1:length(cells), 'XTickLabel', cells, 'YTick', 1:length(cells), 'YTickLabel', a);  % xlim([0.3 2.7])
% 

abs(d_outcomecorr)>0.6


mean( d_outcomecorr(:,2)) 
std( d_outcomecorr(:,2)) 




%%

