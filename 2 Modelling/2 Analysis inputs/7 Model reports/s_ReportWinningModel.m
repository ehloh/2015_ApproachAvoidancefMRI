
clear all; close all hidden, clc
% cd('D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical')
cd('/Users/EleanorL/Dropbox/sandisk/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/3 Hierarchical')
cf=load('res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968.mat');
ct= load('res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293.mat');
cfee=load('ETwarp conflict.mat'); ctee=load('ETwarp control.mat'); % Assumed to be the correct model
cfev=cfee.d_etwarpv; ctev=ctee.d_etwarpv;
cfe= cfee.d_etwarp; cte= ctee.d_etwarp;



%% 

% Get details
m=1; 
cf.r_res=sortrows(cf.r_res,3);  cf.spars=cf.r_res{m,2}; cf.modname=cf.r_res{m,1};
cf.modd=cf.details.models(strcmp( cf.details.models(:, 1), cf.r_res{m,1}), :);
ct.r_res=sortrows(ct.r_res,3);  ct.spars=ct.r_res{m,2}; ct.modname=ct.r_res{m,1};
ct.modd=ct.details.models(strcmp( ct.details.models(:, 1), ct.r_res{m,1}), :);
cfe=[cfe(size(cf.spars,1)+1:end, :);  cfe(1:size(cf.spars,1), :)];   cte=[cte(size(ct.spars,1)+1:end, :);  cte(1:size(ct.spars,1), :)]; 



% Parameter ranges
figure('color','w'), f.FontSize=20; f.FontName='PT Sans Caption';
for p=1: length(cf. modd{3})  % cf
    wp=cf.spars(:, 3+p);
    
    subplot(2,  6 , p); 
%     hist(wp), axis square
%     title(['[cF] ' cf.modd{3}{p}], 'FontSize',  f.FontSize, 'FontName', f.FontName)
    
    barwitherr(std(wp)/sqrt(length(cf.spars)), mean(wp), 'y')
    title([cf.modd{3}{p}], 'FontSize',  f.FontSize, 'FontName', f.FontName)
    set(gca, 'FontSize',  f.FontSize, 'FontName', f.FontName)
end
for p=1: length(ct.modd{3})  % ct 
    wp=ct.spars(:, 3+p);
    
    subplot(2, 6, p+6); 
%     hist(wp), axis square
%     title(['[ct] ' ct.modd{3}{p}], 'FontSize',  f.FontSize)

    barwitherr(std(wp)/sqrt(length(cf.spars)), mean(wp), 'y')
    title([cf.modd{3}{p}], 'FontSize',  f.FontSize, 'FontName', f.FontName)
    set(gca, 'FontSize',  f.FontSize, 'FontName', f.FontName)
end


cf. modd{3}
p=4
mean(cf.spars(:,p))
std(cf.spars(:,p))./sqrt(20)


ct. modd{3}
mean(ct.spars(:,p))
std(ct.spars(:,p))./sqrt(20)


% Shared parameters correlate between cF and ct?
sharedpar={'b';'p';'j';'i';'w'};  % HARD CODE
figure('color','w'), f.fontsize=15; f.fontsize_title=20; f.fontname='PT Sans Caption'; k=1; f.plotcols=3; 
for p=1:length(sharedpar)
    wp.cf= cf.spars(:, 3 + find(strcmp(cf.modd{3}, sharedpar{p}))    );
    wp.ct= ct.spars(:, 3 + find(strcmp(ct.modd{3}, sharedpar{p}))    );
    [wp.r wp.p]=corr(wp.cf, wp.ct);
    
    subplot(ceil(length(sharedpar)/f.plotcols),f.plotcols,  k); k=k+1; 
    scatter(wp.cf, wp.ct), lsline, axis square
    xlabel('Exp parameter value','FontSize', f.fontsize, 'FontName', f.fontname)
    ylabel('Ctrl parameter value','FontSize', f.fontsize, 'FontName', f.fontname)
    set(gca,'FontSize', f.fontsize, 'FontName', f.fontname)
    title([sharedpar{p} ': r= ' num2str(wp.r,2)    ', p= ' num2str(wp.p,3)],'FontSize', f.fontsize_title, 'FontName', f.fontname)
    %
    wp=[];
end

% Different parameters correlated across subjects?
figure('color','w'), f.fontsize=20; f.fontsize_title=25; f.fontname='PT Sans Caption'; 
[r p]= corr(cf.spars(:, 4:end));  % cf
rr=r(:); rr((p(:)>0.05))=nan; r= reshape(rr,cf.modd{2}, cf.modd{2}); 
subplot(1,2,1), imagescnan(r,'nancolor',[0 0 0]), axis square,  colorbar
ylabel('Parameter','FontSize', f.fontsize, 'FontName', f.fontname),  set(gca, 'ytick', 1:cf.modd{2}, 'yticklabel',  cf.modd{3},'FontSize', f.fontsize, 'FontName', f.fontname)
xlabel('Parameter','FontSize', f.fontsize, 'FontName', f.fontname),  set(gca, 'xtick', 1:cf.modd{2}, 'xticklabel',  cf.modd{3},'FontSize', f.fontsize, 'FontName', f.fontname)
title(['Exp  - '  cf.modd{1}],'FontSize', f.fontsize_title, 'FontName', f.fontname), set(gca, 'FontSize', f.fontsize, 'FontName', f.fontname)
[r p]= corr(ct.spars(:, 4:end));
rr=r(:); rr((p(:)>0.05))=nan; r= reshape(rr,ct.modd{2}, ct.modd{2}); 
subplot(1,2,2), imagescnan(r,'nancolor',[0 0 0]), axis square, colorbar
ylabel('Parameter','FontSize', f.fontsize, 'FontName', f.fontname),  set(gca, 'ytick', 1:ct.modd{2}, 'yticklabel',  ct.modd{3},'FontSize', f.fontsize, 'FontName', f.fontname)
xlabel('Parameter','FontSize', f.fontsize, 'FontName', f.fontname),  set(gca, 'xtick', 1:ct.modd{2}, 'xticklabel',  ct.modd{3},'FontSize', f.fontsize, 'FontName', f.fontname)
title(['Ctrl  - '  ct.modd{1}],'FontSize', f.fontsize_title, 'FontName', f.fontname), set(gca, 'FontSize', f.fontsize, 'FontName', f.fontname)


% EnvThreat distortions: is the task space correlated across cF and ct?
request.etplot={'EnvThreat';'NTok';'pLoss';'Entropy';'EntropyNTok'};
r_etcor=[request.etplot';cell(size(cf.spars,1),    length(request.etplot))]; rs_etcor=r_etcor;   rp_etcor=r_etcor; 
for p=1:length(request.etplot)
    for s=1:size(cf.spars,1)+2
        eval(['ws.cf=cfe{ s,2}.' request.etplot{p} ';'])
        eval(['ws.ct=cte{ s,2}.' request.etplot{p} ';'])
        [r_etcor{s+1, p} rp_etcor{s+1, p}]=corr(ws.cf(:), ws.ct(:));
        
        if rp_etcor{s+1, p}<0.05
            rs_etcor{s+1, p}=r_etcor{s+1, p};
        else rs_etcor{s+1, p}=nan;
        end
        ws=[];
    end
   
    % Data for plots
    r_etcorplot(:, p)=cell2mat(r_etcor(2:end, p));
    
    wp=[];
end
%
f.fontsize=15; f.fontsize_title=20; f.fontname='PT Sans Caption';
f.subplot_VerHorz=[0.025 0.1]; f.fig_BotTop=[0.03 0.01]; f.fig_LeftRight=[0.1 0.1]; k=1; 
figure('Name', 'cF vs ct task space after EnvThreat distortion', 'Color', 'w');  
ws.corrres=sortrows([cf.details.subjects  cfe(3:end,1) cte(3:end,1) num2cell(cell2mat(cfe(3:end,1))-cell2mat(cte(3:end,1))) num2cell(r_etcorplot(3:end,:))], 4);
subtightplot(2,2, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1; imagescnan(cell2mat(ws.corrres(:, 5:end)),'nancolor', [0 0 0]); axis square, colorbar
ylabel('Subject','FontSize', f.fontsize, 'FontName', f.fontname), set(gca, 'ytick', 1:size(ws.corrres,1), 'yticklabel',  ws.corrres(:,1))
xlabel('Variable','FontSize', f.fontsize, 'FontName', f.fontname), set(gca, 'xtick', 1:5, 'xticklabel', request.etplot,'FontSize', f.fontsize, 'FontName', f.fontname)
% rotateXLabels( gca(), 45)
% set(gca,'XAxisLocation', 'bottom')
k=k+1;
subtightplot(2,2, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1; imagescnan(cell2mat(ws.corrres(:, 2:3)),'nancolor', [0 0 0]); axis square, colorbar
ylabel('Subject','FontSize', f.fontsize, 'FontName', f.fontname), set(gca, 'ytick', 1:size(ws.corrres,1), 'yticklabel',  ws.corrres(:,1))
xlabel('j param','FontSize', f.fontsize, 'FontName', f.fontname), set(gca, 'xtick', 1:2, 'xticklabel', {'Exp j';'Ctrl j'},'FontSize', f.fontsize, 'FontName', f.fontname)
subtightplot(2,2, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1; imagescnan(cell2mat(ws.corrres(:, 4)),'nancolor', [0 0 0]); axis square, colorbar
ylabel('Subject','FontSize', f.fontsize, 'FontName', f.fontname), set(gca, 'ytick', 1:size(ws.corrres,1), 'yticklabel',  ws.corrres(:,1))
xlabel('j param','FontSize', f.fontsize, 'FontName', f.fontname), set(gca, 'xtick', 1:2, 'xticklabel', {'j, Exp - Ctrl'},'FontSize', f.fontsize, 'FontName', f.fontname)
%




% EnvThreat distortion: average task space experienced (cF and ct)
request.plotwhich={'EnvThreat';'NTok';'pLoss';'Entropy';'EntropyNTok';'EV'};
f.plotcols=length(request.plotwhich);   f.fontsize=15; f.fontsize_title=20; f.fontname='PT Sans Caption';
f.subplot_VerHorz=[0.02 0.02]; f.fig_BotTop=[0.01 0.03]; f.fig_LeftRight=[0.04 0.01];
figure('Name', 'Distorted details.task space', 'NumberTitle', 'off', 'Position', [130 85 900 600], 'Color', 'w');  k=1;
%
for  p=1:length(request.plotwhich)
    subtightplot(1,  f.plotcols, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
    eval(['wp.cf=cfe{end,2}.' request.plotwhich{p} ';'])
    eval(['wp.ct=cte{end,2}.' request.plotwhich{p} ';'])
    imagesc(  (wp.cf +wp.cf)./2  ); colorbar, axis square
    
%     imagesc(  wp.ct  ); colorbar, axis square
    
    
    title(request.plotwhich{p},'FontSize',f.fontsize_title,'FontName', f.fontname)
    
    
    switch request.plotwhich{p}
        case 'EnvThreat', caxis([0 1])
        case 'NTok';
            title('No. Tokens','FontSize',f.fontsize_title,'FontName', f.fontname)
        case 'pLoss', caxis([0 1])
            title('p(ActBomb)','FontSize',f.fontsize_title,'FontName', f.fontname)
        case 'Entropy'  % , caxis([0 1])
            title('Uncertainty','FontSize',f.fontsize_title,'FontName', f.fontname)
        case 'EntropyNTok', caxis([0 8])
            title('Value-scaled uncertainty','FontSize',f.fontsize_title,'FontName', f.fontname)
%         case 'EV'   % , caxis([0 1])
%             title('EV (Ctrl)','FontSize',f.fontsize_title,'FontName', f.fontname)
    end
    
    
    % Labels
    ylabel('Env Threat','FontSize',f.fontsize,'FontName', f.fontname); xlabel('No. Tokens','FontSize',f.fontsize,'FontName', f.fontname)
    set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6,'FontSize',f.fontsize,'FontName', f.fontname)
end



%% Correlations of modelling results with other measures 




error
%
f.plotcols=4;  f.subplot_VerHorz=[0.02 0.02]; f.fig_BotTop=[0.01 0.03]; f.fig_LeftRight=[0.07 0.02];
figure('Name', 'cF vs ct task space after EnvThreat distortion', 'NumberTitle', 'off', 'Position', [130 85 900 600], 'Color', 'w');  k=1;

    
    subtightplot(length(request.etplot)/f.plotcols ,  f.plotcols, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);




