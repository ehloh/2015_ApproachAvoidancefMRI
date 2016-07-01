% Look at BICs over the entire model space
clear all; close all hidden; clc

cd('3 Hierarchical');
% load('res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968')
load('res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293')


% Compare w Bayes factor?
m1=1;
m2=2;
B=exp((r_res{m1,3}-r_res{m2,3})*-0.5);
disp(['B=' num2str(B)  '  (m1=' r_res{m1,1} ', m2='  r_res{m2,1} ')     [3-10=moderate evidence, >10=strong]'])


%%
r_res=sortrows(r_res,-3);
bics=cell2mat(r_res(:,3));
% bics= bics - min(bics); disp('Plotting change in BIC')
n_params=cellfun(@(x)size(x,2)-3, r_res(:,2));
modnames=r_res(:,1); conds={'Exp';'Ctrl'};

% Sorting by # params
all=[num2cell(bics) num2cell(n_params) modnames]; all=sortrows(all, 2);
trans=[1 find(cell2mat(all(:,2))==2,1,'first') find(cell2mat(all(:,2))==3,1,'first') find(cell2mat(all(:,2))==4,1,'first') find(cell2mat(all(:,2))==5,1,'first') find(cell2mat(all(:,2))==6,1,'first') find(cell2mat(all(:,2))==7,1,'first')];

% BIC plots 
figure('color','w'); k=1; f.FontSizeTitle=30;  f.FontSize=25; f.FontName='PT Sans';
subplot(3,1,k); bar(bics-min(bics)); title('Difference in BIC from winning model', 'FontSize', f.FontSizeTitle); k=k+1; xlabel('Model', 'FontSize', f.FontSize), ylabel('Change in BIC', 'FontSize', f.FontSize)
subplot(3,1,k); bar(bics./n_params); title('BIC/No. params', 'FontSize', f.FontSizeTitle); k=k+1; xlabel('Model', 'FontSize', f.FontSize), ylabel('BIC/No. Params', 'FontSize', f.FontSize)
subplot(3,1,k); 

figure('color','w'); 

bar(cell2mat(all(:,1)) ); title('Change BIC across entire model space', 'FontSize', f.FontSizeTitle,'FontName', f.FontName), ylabel('Change in BIC', 'FontSize', f.FontSize,'FontName', f.FontName);
for i=2:7
    hx = graph2d.constantline(trans(i)-0.5, 'LineStyle', ':', 'Color','b'); changedependvar(hx,'x');
end
hold on, bar( find(cell2mat(all(:,1))== min(cell2mat(all(:,1)))), min(bics), 'r'); % highlist best 
set(gca, 'FontSize', f.FontSize,'FontName', f.FontName);
xlim([0 192])
ylim([0 30000])

% error

% BIC minus one
figure('color','w');  f.FontSizeTitle=40;  f.FontSize=30; 
allmods=details.models(:,1);
for m=1:size(allmods,1)
    allmods{m,2}= char(sortrows(details.models{m,3}))'; 
    
    if isempty(strfind(allmods{m,2}, 'w'))==0
        allmods{m,3}= allmods{m,1}(strfind(allmods{m,1}, 'w')-1);
    end
    
end
r_res=sortrows(r_res,3); pars=details.models{strcmp( details.models(:, 1), r_res{1,1}),3}; pars=char(sortrows(pars))'; winmod=r_res{1,1};
d_lessmod={}; % Name, group BIC, subject Ls, mssing param, 
for p=2:length(pars)
    wp.pars=pars;
    wp.pars(p)=[];
    
    wp.which=find(strcmp(allmods(:,2), wp.pars));
    if length(wp.which)~=1;
        wp.mods=allmods(wp.which, :);
        wp.thismod= wp.mods( strcmp(wp.mods(:,3),  winmod(strfind(winmod, 'w')-1)), :);
    else wp.thismod= allmods( wp.which, :);
    end
    
    d_lessmod{p,1} = wp.thismod{1};
    d_lessmod{p,2} = r_res{ strcmp(r_res(:,1),  d_lessmod{p,1}), 3}; 
    d_lessmod{p,3} = r_res{ strcmp(r_res(:,1),  d_lessmod{p,1}), 2}(:,2); 
    d_lessmod{p,4} = pars(p);
end
d_lessmod(1,:)=[r_res(1,1) r_res(1,3) {r_res{1,2}(:,2)} 'Winning']; 
d_lessmod=[ d_lessmod(1,:);  sortrows(d_lessmod(2:end,:) ,2)];
d_lessmod= sortrows(d_lessmod, -2); 

bar(cell2mat(d_lessmod(:, 2)), 'y'); 
xlabel('Model', 'FontSize', f.FontSize,'FontName', f.FontName), ylabel('Model BIC', 'FontSize', f.FontSize,'FontName', f.FontName)
set(gca, 'FontSize', f.FontSize,'FontName', f.FontName)
set(gca, 'xticklabel', [cellfun(@(x) ['Omit ' x],  d_lessmod(1:end-1,4), 'UniformOutput', false); {'Winning'};])
title(conds{details.tasktype} , 'FontSize', f.FontSizeTitle);
ylim([10000 13000])

