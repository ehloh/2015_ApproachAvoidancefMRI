% Plot pmods from onsets file itself
clear all, close all hidden, clc

cd('/Users/EleanorL/Dropbox')
load('p01_GV_onsets_m_v3g_vChosenAnd_bpji08bpji11.mat')
% load('p34_TB_onsets_m_v3g_vChosenAnd_bpji08bpji11.mat')

whichpmods=3:6;
% whichpmods=3:4;

%% Get pmods
onsname=cellstr(names');
pvar=[]; % nan( length(pmod(whichpmods(1)).param{1}), length(whichpmods));
pname={};
for p=1:length(whichpmods)
    if isempty(pmod(whichpmods(p)).name)==0
        pvar(1:length(pmod(whichpmods(p)).param{1}), p)=pmod(whichpmods(p)).param{1};
        pname{p,1}=[pmod(whichpmods(p)).name{1} '    ' num2str(p)];
    else pname{p,1}=[];
    end
end


% pvar=sortrows(pvar,1);


pvar=[sortrows(pvar(:,1:2),1) sortrows(pvar(:,3:4),1)];


% min(pvar(:,2))

%% Display

figure('color','w');  colors={'b';'r';'y';'c';'b';'m';'g';'p'};
%
subplot(3,1,1);  % How correlated?
[r p]= corr(pvar);
rr=r(:);   rr(p(:)>0.05)=nan; r= reshape(rr, size(pvar,2), size(pvar,2));
imagescnan(r, 'nancolor', [0.1 0.1 0.1]), colorbar, axis square
set(gca, 'ytick', 1:length(pname), 'yticklabel', pname)
subplot(3,1,2);  % Plot variable
for p=1:size(pvar,2) 
    plot(1:size(pvar,1), pvar(:,p),colors{p}), hold on
end
disp('COLOURS:'),  disp([colors(1:length(pname)) (pname)   ]);
subplot(3,1,3);  % Plot variable closeup
for p=1:size(pvar,2) 
    plot(1:size(pvar,1), pvar(:,p),colors{p}), hold on
end
xlim([0 50]);

figure, scatter(pvar(:,1), pvar(:,2))
