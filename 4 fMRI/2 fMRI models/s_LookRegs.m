% Delete specific files/folders for all subjects 
clear all;close all hidden; clc

% Request specific
log.specificsubjects={
% 'p01_GV';'p06_KB'
}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications
   
    % Load subjects
    where.where='D:\Dropbox\SANDISK\1 Explore fMRI';  where.data_brain='G:\2 [Explore]\1 Brain data'; 
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI'; where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/1 Brain data';  
    
     
    where.data_beh=[where.where filesep '1 Behavioural data'];
    
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all'); 
    
    % Interface
    disp('=======================================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp('Requested analysis: ADHOC SCRIPT')
    disp('See direct code for exact actions to be executed !! ')
    disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location (brain): ' where.data_brain])
    disp(['Data location (behaviour): ' where.data_beh])
    disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end
 f_mat=@(A,x)A(x); 

where.data_brain = '/Users/EleanorL/Dropbox/WorkPC/Conflict'; 

mods = {
    'm_v3g_vChosenAnd_bpji08bpji11',  {'Sn(1) cF_onsetxcF_vChosen^1*bf(1)';'Sn(1) cF_onsetxcF_vBestUnchosen^1*bf(1)';'Sn(1) ct_onsetxct_vChosen^1*bf(1)';'Sn(1) ct_onsetxct_vBestUnchosen^1*bf(1)';};
    'm_v36g_RejectOrvMargChoDiff_bpji08bpji11'   {'Sn(1) cF_onsetxcF_Rej_vMargChosen^1*bf(1)';'Sn(1) cF_onsetxcF_NonRej_vMargChosen^1*bf(1)';'Sn(1) ct_onsetxct_Rej_vMargChosen^1*bf(1)';'Sn(1) ct_onsetxct_NonRej_vMargChosen^1*bf(1)';}
};

%% 

mods{1,3} ='Chosen and counterfactual value';
mods{1,4}={'Ap/Av V(Chosen)   ', 'Ap/Av V(Best Unchosen)   ', 'Ap/Ap V(Chosen)   ', 'Ap/Ap V(Best Unchosen)   '};
mods{2,3} ='Choice x Value difference'; 
mods{2,4} = {'Ap/Av Reject, V(Chosen) > V(Best Unchosen)', 'Ap/Av Accept/Explore, V(Chosen) > V(Best Unchosen)', 'Ap/Ap Bomb, V(Chosen) > V(Best Unchosen)', 'Ap/Ap No bomb/Explore, V(Chosen) > V(Best Unchosen)'}; 
 
d=cell(log.n_subjs,2);
d_corr=nan(log.n_subjs, 8*8);
for s=1: log.n_subjs  % Collect regressors 
    ws.c=clock;  disp(['Subject ' num2str(s) '   -  ' log.subjects{s} '   [' num2str(ws.c(4)) ':' num2str(ws.c(5)) ']  ------------------']);
    ws.subfol=[where.data_brain filesep log.subjects{s} filesep]; 
    ws.subfol=[ws.subfol '2 First level s4Ants' filesep];
 
    
    for m=1:size(mods,1)
        wm.s=load([ws.subfol mods{m,1} '_Basic Contrasted' fs 'SPM.mat']);
        wm.reg = [];
        for r=1:length(mods{m,2})
            wm.reg =[wm.reg   wm.s.SPM.xX.X(:,  find(strcmp(wm.s.SPM.xX.name,  mods{m,2}{r})))];
        end
        
        d{s,m} =  wm.reg; wm=[];
    end
     
    ws.d = [d{s,1} d{s,2}];
    
    
    ws.r =  corr(ws.d); 
    d_corr(s,:) = ws.r(:); 
 
ws=[];
end

 
log.nVar = 8; 

log.VarName =  [mods{1,4} mods{2,4}]';

%% Plot all 

% Req
req.r_limit= 0.4; 


close all hidden
f.plotcols=2;  f.figwidth= 800; f.figheight=800; f.fontsize=15; f.fontsize_title=30;
f.subplot_VerHorz=[0.02 0.10]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.2 0.1];
f.markersize= 4; 
for o=1:1 % Plots 
    figure('Name', 'Pairwise correlatedness', 'NumberTitle', 'off', 'Position', [200 150 f.figwidth f.figheight], 'Color', 'w');  k=1;

    % Group plot
    subtightplot(log.nVar,  f.plotcols,((1:log.nVar)-1)*f.plotcols+1,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
    %   subtightplot(log.nVar+1,  f.plotcols,k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);  k=k+1;
    imagesc( reshape(mean(d_corr), [log.nVar log.nVar])), axis square,  colorbar
    set(gca,'xtick', 1:log.nVar,'ytick', 1:log.nVar,'xticklabel', 1:log.nVar, 'yticklabel',  cellfun(@(x,n)[x '  [' num2str(n) ']   '], log.VarName, num2cell(1:log.nVar)','uniformoutput',0), 'FontSize',f.fontsize)
    title('Mean correlation r', 'FontSize', f.fontsize_title )
    
    for v=1:log.nVar   % Group plot w SEs
        subtightplot(log.nVar,  f.plotcols, (v-1)*f.plotcols+2 ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        %   subtightplot(log.nVar+1,  f.plotcols, k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);  k=k+1;
        barwitherr(f_mat(reshape(std(d_corr), [log.nVar log.nVar])/sqrt(log.n_subjs),    (v-1)*log.nVar+(1:log.nVar)),         f_mat(reshape(mean(d_corr), [log.nVar log.nVar]),    (v-1)*log.nVar+(1:log.nVar))    ,'y'   ); hold on
        scatter(f_mat(  repmat(1:log.nVar, log.n_subjs,1),   1:log.nVar*log.n_subjs),  f_mat(  d_corr(:, (1:log.nVar)+(v-1)*log.nVar),    1:log.nVar*log.n_subjs),f.markersize), hold on
        plot(0:log.nVar+1, ones(log.nVar+2,1)*req.r_limit, 'g'), hold on
        plot(0:log.nVar+1, ones(log.nVar+2,1)*-req.r_limit, 'g')
        xlim([0 length(log.VarName)+1]) 
        ylabel(log.VarName{v}, 'FontSize', f.fontsize)
        
        % Collect data in a more identifiable context  
        wv.vd=  d_corr(:, (1:log.nVar)+(v-1)*log.nVar);
        for vv=1:log.nVar   
            d_corrbeh{v,vv}= wv.vd(:, vv);
        end 
        
    end
end

er 

%% Plot for paper 

dd_mc= reshape(mean(d_corr), [log.nVar log.nVar]);
% d_mc = d_mc(5:8, 1:4); d_smc =  d_corrbeh(5:8, 1:4);

d_mc{1}  = dd_mc(5:6,  1:2);
d_mc{2}  = dd_mc(7:8,  3:4);
 



% a= d_corrbeh(5:6,  1:2);
a= d_corrbeh(7:8,  3:4);
cellfun(@(x)['mean=' num2str(mean(x)) ' se=' num2str(std(x)/sqrt(20))], a, 'uniformoutput',0)
 


cellfun(@(x) mean(x), a)


cellfun(@(x) std(x)/sqrt(20), a)

 
    
    


% close all hidden
f.plotcols=2;  f.figwidth= 1000; f.figheight=800; f.fontsize=25; f.fontsize_title=30;
f.subplot_VerHorz=[0.02 0.10]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.2 0.1];
f.markersize= 4; 
for o=1:1 % Plots 
    figure('Name', 'Pairwise correlatedness', 'NumberTitle', 'off', 'Position', [200 150 f.figwidth f.figheight], 'Color', 'w');  k=1;

    % Ap/Av condition 
    subtightplot(1,  f.plotcols, 1,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); 
    imagesc(d_mc{1}), axis square,  colorbar 
    title('Ap/Av condition', 'FontSize', f.fontsize_title )
    set(gca,'xtick', 1:2,'ytick', 1:2, ...
        'xticklabel', log.VarName(5:6), 'yticklabel', log.VarName(1:2)) 
    xticklabel_rotate 
    set(gca, 'FontSize',f.fontsize)
    
    
    % Ap/Ap condition 
    subtightplot(1,  f.plotcols, 2,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);  
    imagesc(d_mc{2}), axis square,  colorbar 
    title('Ap/Ap condition', 'FontSize', f.fontsize_title )
    set(gca,'xtick', 1:2,'ytick', 1:2, ...
        'xticklabel', log.VarName(7:8), 'yticklabel', log.VarName(3:4)) 
    xticklabel_rotate 
    set(gca, 'FontSize',f.fontsize)
end


set(gca, 'FontSize', 50)




