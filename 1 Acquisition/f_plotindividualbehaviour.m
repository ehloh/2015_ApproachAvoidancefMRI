function [b] = f_plotindividualbehaviour(p,data, subjname)
% Generate heat-plots for individual subject's behaviour
% [b] = f_plotindividualbehaviour(p,data, subjname)
%
% Plots & variables:
% 
%         b.accept
%         b.reject
%         b.explore
%         b.choicevar (in progress)
%         b.RT
%         b.outcomemean
%         b.outcomevar
%
% To plot a single subject's data, fill in the following: 
%
%         clear all; close all hidden
%         subjname='p11_SL';
%         stagenum=2;
%         stage='taskconflict'; % 'taskconflict' (#2), 'taskcontrol' (#3)
%         where='H:\5 Conflict Exploration\2 Task [Goal-conflict-exploration]';
%         %
%         load([where filesep 'Data' filesep subjname '_file_' num2str(stagenum) stage '.mat'])
%         eval(['data=' stage '.data;'])
%         eval(['p=' stage '.settings;'])
%
% -------------------------------------------------------------------------

for o1=1:1 % Documentation for par
%     Col 1:      Trial type (read # from grid, starting from bottom-left corner)
%     Col 2:      # Token-pairs offered (Trial type X, 1-6)
%     Col 3:      Set pBomb level  (Trial type Y, 1-6 in ascending order of probability)
%     Col 4:      Explored half (if applicable)
%     Col 5:      [BLANK]
%     Col 6:      Bomb present in Activated tokens?
%     Col 7:      Trial number
%     Col 8:      [1st Response] Accept, Reject, Explore?  (1=Accept, 2=Reject, 3=Explore)
%     Col 9:      [1st Response] RT
%     Col 10:     [2nd Response] Accept or Reject? 
%     Col 11:     [2nd Response] RT
%     Col 12:     Outcome (1=Gain, -1=Loss, 0=No effect) 
%     Col 13:     Trial aborted (to be repeated later) (1=Aborted)  
%     Col 14:     Did exploration reveal a bomb? (1=Yes, 0=No, 999=Not explored)
%     Col 15:     Outcome in amount (payment)
end

%% Generate data matrices
b.accept=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels);  % Resets
b.reject=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
b.explore=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
b.choicevar=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels); % how to calculate this?
b.RT=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
b.outcomemean=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
b.outcomevar=nan*ones(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
for j=1:p.task.envthreat_nlevels 
    for i=1:p.task.envthreat_nlevels
        ws=[];
        ws.choice=data((data(:,3)==j & data(:,2)==i),8);
        ws.RT=data((data(:,3)==j & data(:,2)==i),9);
        ws.outcome=data((data(:,3)==j & data(:,2)==i ),15);
        b.accept(p.task.envthreat_nlevels+1-j,i)=sum(ws.choice(:)==1)/size(ws.choice,1);
        b.reject(p.task.envthreat_nlevels+1-j,i)=sum(ws.choice(:)==2)/size(ws.choice,1);
        b.explore(p.task.envthreat_nlevels+1-j,i)=sum(ws.choice(:)==3)/size(ws.choice,1);
        b.RT(p.task.envthreat_nlevels+1-j,i)=mean(ws.RT);
        b.outcomemean(p.task.envthreat_nlevels+1-j,i)=mean(ws.outcome);
        b.outcomevar(p.task.envthreat_nlevels+1-j,i)=var(ws.outcome);
    end
end

%% Plot

figure('Name',['Behavioural performance: ' subjname],'NumberTitle','off','Position',[715,220,900,600]);
set(gcf,'Color',[1 1 1])
for i=1:p.task.envthreat_nlevels
    ylabels{p.task.envthreat_nlevels+1-i,1}=[num2str(i) '/' num2str(p.task.envthreat_nlevels)];
end
        
    subplot(2,3,1) % Accept
    imagesc(b.accept, [0 1])
    axis('square')
    title('% Accept')
    colorbar
    ylabel('p(Bomb) present');
    xlabel('# of Tokens offered');
    set(gca,'YTick',1:p.task.envthreat_nlevels)
    set(gca, 'YTickLabel', ylabels)
    set(gca,'XTick',1:p.task.envthreat_nlevels)
    set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
    
    subplot(2,3,2) % Reject
    imagesc(b.reject, [0 1])
    axis('square')
    title('% Reject')
    colorbar
    ylabel('p(Bomb) present');
    xlabel('# of Tokens offered');
    set(gca,'YTick',1:p.task.envthreat_nlevels)
    set(gca, 'YTickLabel', ylabels)
    set(gca,'XTick',1:p.task.envthreat_nlevels)
    set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
    
    subplot(2,3,3) % Explore
    imagesc(b.explore, [0 1])
    axis('square')
    title('% Explore')
    colorbar
    ylabel('p(Bomb) present');
    xlabel('# of Tokens offered');
    set(gca,'YTick',1:p.task.envthreat_nlevels)
    set(gca, 'YTickLabel', ylabels)
    set(gca,'XTick',1:p.task.envthreat_nlevels)
    set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
    
    subplot(2,3,4) % Mean RT
    imagesc(b.RT)
    axis('square')
    title('% Mean RT')
    colorbar
    ylabel('p(Bomb) present');
    xlabel('# of Tokens offered');
    set(gca,'YTick',1:p.task.envthreat_nlevels)
    set(gca, 'YTickLabel', ylabels)
    set(gca,'XTick',1:p.task.envthreat_nlevels)
    set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
    
    subplot(2,3,5) % Outcome mean
    imagesc(b.outcomemean)
    axis('square')
    title('Outcome mean')
    colorbar
    ylabel('p(Bomb) present');
    xlabel('# of Tokens offered');
    set(gca,'YTick',1:p.task.envthreat_nlevels)
    set(gca, 'YTickLabel', ylabels)
    set(gca,'XTick',1:p.task.envthreat_nlevels)
    set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
    
    subplot(2,3,5) % Outcome variance
    imagesc(b.outcomevar, [0 1])
    axis('square')
    title('Outcome variance')
    colorbar
    ylabel('p(Bomb) present');
    xlabel('# of Tokens offered');
    set(gca,'YTick',1:p.task.envthreat_nlevels)
    set(gca, 'YTickLabel', ylabels)
    set(gca,'XTick',1:p.task.envthreat_nlevels)
    set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)

end

