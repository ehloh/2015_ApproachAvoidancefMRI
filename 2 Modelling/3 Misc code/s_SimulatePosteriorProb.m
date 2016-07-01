% Simulations for posterior probability
% clear all; close all hidden; clc
clear all; clc


% where.acquisition_script='/Volumes/SANDISK/5 Explore fMRI/2 fMRI acquisition scripts';
where.acquisition_script='I:\5 Explore fMRI\2 fMRI acquisition scripts';

addpath(where.acquisition_script)

%% Specify setup (follow acqusition scripts)

p.task.envthreat_nlevels=6;
p.n.reps_pertrialtype=1000;
p.task.pBombEnv=nan*zeros(p.task.envthreat_nlevels,1);
p.task.pBombEnv=((p.task.envthreat_nlevels:-1:1)/6)';
p.task.explorationcost=2; % # of tokens
p.task.fixedloss=-2*p.task.envthreat_nlevels;
exec.checkstruc=1;
exec.plotfigs=0;
[p par Norm check]=f_generate_taskstruc(p,exec);

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

col.NTokenPairs=2;
col.Trialtype=1;
col.EnvThreatLevel=3;
col.BombPresent=5;
col.ActBomb=6;
end

% New variables 
col.Choice=7; 
col.ExploredBomb=8;

%% Calculate probabilities
%   For all trials, choose Explore, coin-toss explored ActBomb if present

par(:, col.Choice)=3;
par(:,col.ExploredBomb)=0;
par(par(:, col.ActBomb)==1, col.ExploredBomb)=2-randi(2,sum(par(:, col.ActBomb)==1), 1);
r.p_BombExplored_if_ActBomb=mean(par(par(:, col.ActBomb)==1,col.ExploredBomb));


% Calculate observed conditional probabilities
r_pNoActBomb_NoExpBomb=nan*zeros(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
r_pActBomb_NoExpBomb=nan*zeros(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
for e=1:p.task.envthreat_nlevels
    for n=1:p.task.envthreat_nlevels
        ws.d=par(par(:,col.EnvThreatLevel)==e & par(:,col.NTokenPairs)==n, :);
        ws.d=sortrows(ws.d,[-1*col.ActBomb -1*col.ExploredBomb]);
        
        % p(NoActBomb|NoExploredBomb)
        r_pActBomb_NoExpBomb(7-e,n)=sum(ws.d(:,col.ActBomb)==1 & ws.d(:,col.ExploredBomb)==0)/sum(ws.d(:,col.ExploredBomb)==0);
        r_pNoActBomb_NoExpBomb(7-e,n)=sum(ws.d(:,col.ActBomb)==0 & ws.d(:,col.ExploredBomb)==0)/sum(ws.d(:,col.ExploredBomb)==0);
        
        %
        ws=[];
    end
end

figure; 
subplot(3,2,1)
imagesc(r_pActBomb_NoExpBomb,[0 1]); axis square; 
colorbar; title('Observed p(Act Bomb|No Explored Bomb)');ylabel('EnvThreat'); xlabel('N Tokens');
set(gca, 'YTickLabel', [num2str((p.task.envthreat_nlevels:-1:1)') repmat('/6', 6,1)])
set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)

subplot(3,2,2)
imagesc(r_pNoActBomb_NoExpBomb,[0 1]); axis square; 
colorbar; title('Observed p(No Act Bomb|No Explored Bomb)');ylabel('EnvThreat'); xlabel('N Tokens');
set(gca, 'YTickLabel', [num2str((p.task.envthreat_nlevels:-1:1)') repmat('/6', 6,1)])
set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)




%% Formulae calculation

form.postpActBomb='(ET*NTok)/(24-ET*NTok)'; % my calculation

rf_pActBomb_NoExpBomb=nan*zeros(p.task.envthreat_nlevels,p.task.envthreat_nlevels);
for e=1:p.task.envthreat_nlevels
    for n=1:p.task.envthreat_nlevels
       ET=p.task.pBombEnv(7-e);
       NTok=n*2;
        eval(['rf_pActBomb_NoExpBomb(7-e,n)=' form.postpActBomb ';'])
    end
end

subplot(3,2,3)
imagesc(rf_pActBomb_NoExpBomb,[0 1]); axis square; 
colorbar; title('Predicted p(Act Bomb|No Explored Bomb)');ylabel('EnvThreat'); xlabel('N Tokens');
set(gca, 'YTickLabel', [num2str((p.task.envthreat_nlevels:-1:1)') repmat('/6', 6,1)])
set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)


obsMpred=r_pActBomb_NoExpBomb-rf_pActBomb_NoExpBomb;
subplot(3,2,4)
imagesc(obsMpred,[-0.05 0.05]); axis square; 
colorbar; title('Observed - Predicted');ylabel('EnvThreat'); xlabel('N Tokens');
set(gca, 'YTickLabel', [num2str((p.task.envthreat_nlevels:-1:1)') repmat('/6', 6,1)])
set(gca, 'XTickLabel',2:2:2*p.task.envthreat_nlevels)

