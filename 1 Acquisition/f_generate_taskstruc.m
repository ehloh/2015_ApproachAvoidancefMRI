function [p par Norm check]= f_generate_taskstruc(p, exec )
% [p Norm par]= f_generate_taskstruc(p, exec )
% 
% Generate structure for task, and check if requested
%
% Things to specify:
%         exec.plotfigs=1;
%         exec.checkstruc=1; 
%         % Design
%         p.task.envthreat_nlevels=6; % How many levels of p(Bomb present)
%         p.n.reps_pertrialtype=12; % If you're going to change this, CHECK the prob struc generated!
%         p.task.pBombEnv=nan*zeros(p.task.envthreat_nlevels,1);
%         for i=1:p.task.envthreat_nlevels
%             p.task.pBombEnv(p.task.envthreat_nlevels+1-i,1)=i/p.task.envthreat_nlevels; % Contingency, in descending probability
%         end
%         p.task.explorationcost=2; % # of tokens
%         p.task.fixedloss=-2*p.task.envthreat_nlevels;
%
% Updated 3/12/12 Eleanor Loh
% ---------------------------------------------------------------------
% % To use this script to generate probabilistic structures for examination,
% % run as a script with the following inputs un-commented & defined
%
%         exec.plotfigs=0;
%         exec.checkstruc=0; 
%         % Design
%         p.task.envthreat_nlevels=4; % How many levels of p(Bomb present)
%         p.n.reps_pertrialtype=12; % If you're going to change this, CHECK the prob struc generated!
%         p.task.pBombEnv=nan*zeros(p.task.envthreat_nlevels,1);
%         for i=1:p.task.envthreat_nlevels
%             p.task.pBombEnv(p.task.envthreat_nlevels+1-i,1)=i/p.task.envthreat_nlevels; % Contingency, in descending probability
%         end
%         p.task.explorationcost=2; % # of tokens
%         p.task.fixedloss=-2*p.task.envthreat_nlevels;
%
% % --------------------------------------------------------------

for o1=1:1 % Documentation for 'par'
%     Col 1:      Trial type (read # from grid, starting from bottom-left corner)
%     Col 2:      # Token-pairs offered (Trial type X, 1-6)
%     Col 3:      Set pBomb level  (Trial type Y, 1-6 in ascending order of probability)
%     Col 4:      Monitoring task trial?
%     Col 5:      Bomb present in set?
%     Col 6:      Bomb present in Activated tokens?
%     Col 7:      Trial number
%     Col 8:      [Monitor] Responded correctly? (1=Hit, 0=Miss, 999=n/a)
%     Col 9:      [Monitor] RT
%     Col 10:     [BLANK]
%     Col 11:     [BLANK]
%     Col 12:     Outcome (1=Gain, -1=Loss, 0=No effect) 
%     Col 13:     [BLANK]
%     Col 14:     [BLANK]
%     Col 15:     Position of first (leftmost) token (1-6)
%
%
% CHOICE TASK ================
%
%     Col 1: Trial number
%     Col 2: Left Colour #
%     Col 3: Right Colour #
%     Col 4: Better colour (correct choice - 1=Left, 2=Right)
%     Col 5: Choice (1=Left, 2=Right)
%     Col 6: Accuracy
end
for o1=1:1 % Task design
    
    % Generate matrices
    Norm.Reward=ones(p.task.envthreat_nlevels,1)*(2:2:p.task.envthreat_nlevels*2);
    Norm.pLoss=(p.task.pBombEnv*ones(1,p.task.envthreat_nlevels)).*(Norm.Reward/(2*p.task.envthreat_nlevels)); 
    Norm.ExpectedValue=(Norm.Reward).*(1-Norm.pLoss)+p.task.fixedloss*Norm.pLoss;
end
for o1=1:1 % Subject-specific parameters
% Assign Trial type
w.yaxis=ones(p.task.envthreat_nlevels*p.task.envthreat_nlevels, 1);
w.xaxis=ones(p.task.envthreat_nlevels*p.task.envthreat_nlevels, 1);
for i=1:p.task.envthreat_nlevels
    w.yaxis((i-1)*p.task.envthreat_nlevels+1:i*p.task.envthreat_nlevels)=i*ones(p.task.envthreat_nlevels,1);
    w.xaxis((i-1)*p.task.envthreat_nlevels+1:i*p.task.envthreat_nlevels)=1:p.task.envthreat_nlevels';
end
w.onerep(:,2)=w.xaxis;
w.onerep(:,3)=w.yaxis;
w.onerep(1:size(w.onerep,1),1)=1:size(w.onerep,1);
w.n_trialtypes=size(w.onerep,1);
% Load one_rep (one run thru of all trial types) as many times as needed
par=zeros(size(w.onerep,1)*p.n.reps_pertrialtype,p.task.envthreat_nlevels);
for i=1:p.n.reps_pertrialtype
    j=i-1;
    par(j*w.n_trialtypes+1:i*w.n_trialtypes,1:3)=w.onerep;
    par(j*w.n_trialtypes+1:i*w.n_trialtypes,4)=i; % Col 4 temporarily tracks rep #
end
% Hard-code probabilistic structure of task (pLoss)
for o2=1:1 
    % 'w.outcomes' lists the set of outcomes for all trials, arranged by trial
    % type (Col 1= all outcomes for Trial type 1). For each trial, read off the
    % outcome from the appropriate column, and strike that outcome off the list
par(:,5)=999;
w.outcomes=zeros(p.n.reps_pertrialtype,p.task.envthreat_nlevels*p.task.envthreat_nlevels); 
w.wcontin=floor(Norm.pLoss*p.n.reps_pertrialtype); % always round DOWN to make it an integer
w.contin=[];
for i=1:p.task.envthreat_nlevels
    w.contin=vertcat(w.contin, w.wcontin(p.task.envthreat_nlevels+1-i,:)');
end
for i=1:p.task.envthreat_nlevels*p.task.envthreat_nlevels % Compile outcomes for all trial types 
    w.outcomes(1:w.contin(i) ,i)=1;
end
w.outcomecounter(1:p.task.envthreat_nlevels*p.task.envthreat_nlevels)=1;w.outcomecounter=w.outcomecounter';
for i=1:size(par,1) % Read outcomes off to par
    ws.trialtype=par(i,1);
    par(i,6)=w.outcomes(w.outcomecounter(ws.trialtype), ws.trialtype);
    w.outcomecounter(ws.trialtype)=w.outcomecounter(ws.trialtype)+1; % Increment counter
    if w.outcomecounter(ws.trialtype)==p.n.reps_pertrialtype+1;
        w.outcomecounter(ws.trialtype)=1;
    end
end
end 
% --------------------------------
% Mark inactivated bombs (training only)
wb.wcontin=p.task.pBombEnv*ones(1,p.task.envthreat_nlevels)*p.task.envthreat_nlevels;
wb.contin=[];
for i=1:p.task.envthreat_nlevels
    wb.contin=vertcat(wb.contin,wb.wcontin(p.task.envthreat_nlevels+1-i,:)');
end
wb.contin=wb.contin*p.n.reps_pertrialtype/p.task.envthreat_nlevels;
wb.outcomes=zeros(p.n.reps_pertrialtype,p.task.envthreat_nlevels*p.task.envthreat_nlevels);
for i=1:size(wb.contin,1)
    wb.nbombs=wb.contin(i);
    for j=1:wb.nbombs
        wb.outcomes(j,i)=1;
    end
end
wb.counters=ones(p.task.envthreat_nlevels*p.task.envthreat_nlevels,1);
par(:,5)=0; % Load bombs to par
par=sortrows(par,[1 -6]); 
for i=1:size(par,1)
    wb.trialtype=par(i,1);
    par(i,5)=wb.outcomes(wb.counters(wb.trialtype),wb.trialtype);
    % Increment counter
    if wb.counters(wb.trialtype)==p.n.reps_pertrialtype
        wb.counters(wb.trialtype)=1;
    else
        wb.counters(wb.trialtype)=wb.counters(wb.trialtype)+1;
    end
end
end
for o1=1:1 % Check statistical structure (observed probabilities) of task
   
    % Distribution of p(Bomb present)
    for i=1:p.task.envthreat_nlevels  
       wc.wherebomb=find(and(par(:,3)==i, par(:,5)==1));
       wc.wheretrials=find(par(:,3)==i);
       check.bomb(p.task.envthreat_nlevels+1-i,1)=size(wc.wherebomb,1)/size(wc.wheretrials,1);
       clear ('wc')
   end
    check.bombdiscrepancy=check.bomb-p.task.pBombEnv;
    if sum(check.bombdiscrepancy)~=0
        disp('Error: Distribution of p(Bomb present) does not reflect the planned contingencies!')
        if sum((par(:,5)-par(:,6))<0)>0
            disp('Error: There are Activated bombs on trials where no bomb (activated or not) have been planted')
        end
    end
    
    % Structure of pLoss
    if exec.checkstruc==1;
        for j=1:p.task.envthreat_nlevels 
            for i=1:p.task.envthreat_nlevels 
                check.pLoss(p.task.envthreat_nlevels+1-j,i)=mean(par(par(:,1)==(j-1)*p.task.envthreat_nlevels+i,6));
                check.bombpresent(p.task.envthreat_nlevels+1-j,i)=mean(par(par(:,1)==(j-1)*p.task.envthreat_nlevels+i,5));
            end
        end
    end
end
for o1=1:1 % Plots if requested
      if exec.plotfigs==1
        scrsz = get(0,'ScreenSize'); % scrsz=[x offset from left edge of screen, y offset from bottom of screen, horizontal width, vertical height];
        figure('Position',[scrsz(3)/2, scrsz(4)/2, 800, 500])
        set(gcf,'Color',[1 1 1])
        for i=1:p.task.envthreat_nlevels
            ylabels{p.task.envthreat_nlevels+1-i,1}=[num2str(i) '/' num2str(p.task.envthreat_nlevels)];
        end
        
        % Reward -------------------
        subplot(2,3,1);
        imagesc(Norm.Reward);
        axis(gca,'square');
        title('Reward offered')
        colorbar
        ylabel('p(Bomb) present');
        xlabel('# of Tokens offered');
        set(gca,'YTick',1:p.task.envthreat_nlevels)
        set(gca, 'YTickLabel', ylabels)
        set(gca,'XTick',1:p.task.envthreat_nlevels)
        set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)

        % Loss -------------------
        subplot(2,3,2);
        imagesc(Norm.pLoss)
        axis(gca,'square');
        title('p(Loss)')
        colorbar
        ylabel('p(Bomb) present');
        xlabel('# of Tokens offered');
        set(gca,'YTick',1:p.task.envthreat_nlevels)
        set(gca, 'YTickLabel', ylabels)
        set(gca,'XTick',1:p.task.envthreat_nlevels)
        set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)

        % Expected value ----------------
        subplot(2,3,3);
        imagesc(Norm.ExpectedValue)
        axis(gca,'square');
        colorbar
        title('E(V)')
        ylabel('p(Bomb) present');
        xlabel('# of Tokens offered');
        set(gca,'YTick',1:p.task.envthreat_nlevels)
        set(gca, 'YTickLabel', ylabels)
        set(gca,'XTick',1:p.task.envthreat_nlevels)
        set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
        
        % GENERATED PROBABILITIES
        % Probability of Loss ----------------
        subplot(2,3,4); 
        imagesc(check.pLoss) 
        axis('square');
        title('Observed p(Loss)')
        colorbar
        set(gca,'YTick',1:p.task.envthreat_nlevels)
        set(gca, 'YTickLabel', ylabels)
        set(gca,'XTick',1:p.task.envthreat_nlevels)
        set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
        
        % p(Bomb present) ----------------
        subplot(2,3,5); 
        imagesc(check.bombpresent) 
        axis('square');
        title('Observed p(Bomb present)')
        colorbar
        set(gca,'YTick',1:p.task.envthreat_nlevels)
        set(gca, 'YTickLabel', ylabels)
        set(gca,'XTick',1:p.task.envthreat_nlevels)
        set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)
      end
end

end