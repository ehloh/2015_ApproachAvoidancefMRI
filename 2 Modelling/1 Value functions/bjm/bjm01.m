function [ V ] = bjm01(x, modelinput)
% 
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[]    data   fPar   col}
%        V:                   Value associated with Accept (1), Reject (2), Explore (3)
%                                   (row parampoint, col=trialnum, 3rd dimension=Accept/Reject/Explore)
%
% ---------------------------------------------------------------------------------------------------------------------------------

beta=20./(1+exp(-x(:,1)));            % 0 < beta < 20
EnvThreatDistort= 20./(1+exp(-x(:,2))) ;       % 0 < j < 20
InOptim= 20./(1+exp(-x(:,3))) ;       % 0 < m < 20

%% General setup

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};
nPP=size(x,1);  % No. of parameter points (unique combinations of parameter values for all free parameters covered by grid)

b=beta;  % Save for ct
etd=EnvThreatDistort;
io=InOptim;

% % Simplifications for debugging.
% dd=nan(36,15); k=1; for e=1:6; for n=1:6; dd(k,:)=data(find(data(:,col.EnvThreat)==e/6 & data(:,col.NTokens)==n*2, 1, 'first'),:); k=k+1; end; end
% data=dd; disp('Fake data!!'); data=sortrows(data, [ col.NTokens -col.EnvThreat]); cell=data(:,[ col.NTokens col.EnvThreat]); reshape(cell(:,1),6,6) reshape(cell(:,2),6,6)

%% Conflict task

% Set up data and parameters for cF task
data_cF=data(data(:, col.Task)==1,:); nT_cF = size(data_cF,1);  % cF
order_cF=data(data(:, col.Task)==1,col.Trialnum); 
V_cF=nan(nPP, nT_cF ,3);
FixedLoss=repmat(fPar.cF_FL,[nPP nT_cF]);  miscvar.FixedLoss=fPar.cF_FL;
ExploreCost=repmat(fPar.cF_EC,[nPP nT_cF]);   miscvar.ExploreCost=fPar.cF_EC;

% Distort EnvThreat + everything knocks on
EnvThreat=power(data_cF(:,col.EnvThreat), EnvThreatDistort);
[ov] = fcf_changeEnvThreat(EnvThreat, data_cF(:,col.NTokens), miscvar);
%
pActBomb=repmat(ov.pLoss',[nPP 1]);   % Task-space variables (knock on effects)
NTok=repmat(data_cF(:,col.NTokens)',[nPP 1]);
EnvThreat=repmat(EnvThreat',[nPP 1]); 
beta=repmat(b, [1 nT_cF]);    % Free parameters
InOptim=repmat(io, [1 nT_cF]) ;

% [1] V(Accept) & V(Reject)
V_cF(:,:,1)= pActBomb .* FixedLoss +   (1-pActBomb).*NTok;
V_cF(:,:,2)=0;

% [2a] Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
%     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
%                                      = (EnvThreat * NTokens) / (24 - EnvThreat*NTokens)   = (pActBomb/2) ./ (1- pActBomb/2) 
%     posterior p(~ActBomb)  = (1-pActBomb) ./ (1- pActBomb./2)
ppActBomb=(EnvThreat.*NTok) ./  (24-(EnvThreat.*NTok));
ppABscaled=power(ppActBomb, InOptim);


% [2b] Likelihood of Accepting, if one Explores and sees no Bomb  ***
%       Calculating EV(Explore) must weight the value by the probability of Accepting/Rejecting. Otherwise, V(Explore) will sometimes give drastically low values that would normally be rejected
posVacc=ppABscaled.* FixedLoss + (1-ppABscaled).* NTok;
posVrej=0;
pAccGivNoBombSeen=posVacc>posVrej;
EV_NoSee=pAccGivNoBombSeen .*  posVacc  + (1-pAccGivNoBombSeen) .* 0;

% [3] V(Explore) =   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
%                       = p(See)*V(See)  +    p(~See) *  [       p(Accept|~See) * { posterior p(ActBomb) * Fixed Loss    +   (1- posterior p(ActBomb))   *  NTokens   }      +     p(Reject|~See) * 0      ]   +   Exploration Cost
V_cF(:,:,3)= 0 +  (1-pActBomb.*0.5) .* EV_NoSee    + ExploreCost;

% Remove workspace variables for conflict task
clearvars EV_NoSee EnvThreat ExploreCost FixedLoss NTok V_AccGivNoSee data_cF nT_cF pAccGivNoBombSeen pActBomb posVacc posVrej ppActBomb ppABscaled miscvar EnvThreatDistort ov InOptim
if size(whos,1)>12; keyboard; disp('Workspace clear from cF to ct: more parameters in workspace than usual'); end

%% Control task

% Set up data and parameters for ct task
data_ct=data(data(:, col.Task)==2,:); nT_ct = size(data_ct,1);  % ct
order_ct=data(data(:, col.Task)==2,col.Trialnum); 
V_ct=nan(nPP, nT_ct,3);
ExploreCost=repmat(fPar.ct_EC,[nPP nT_ct]);  miscvar.ExploreCost=fPar.ct_EC;

% Distort EnvThreat + everything knocks on
EnvThreatDistort=etd;
EnvThreat=power(data_ct(:,col.EnvThreat), EnvThreatDistort);
[ov] = fct_changeEnvThreat(EnvThreat, data_ct(:,col.NTokens), miscvar);
%
pActBomb=repmat(ov.pLoss',[nPP 1]);   % Task-space variables (knock on effects)
NTok=repmat(data_ct(:,col.NTokens)',[nPP 1]);
EnvThreat=repmat(EnvThreat',[nPP 1]);
beta=repmat(b, [1 nT_ct]);    % Free parameters
InOptim=repmat(io, [1 nT_ct]) ;

% [1] V(NoBomb) & V(Bomb)
V_ct(:,:,1)= (1-pActBomb) .* NTok;
V_ct(:,:,2)= pActBomb .* NTok;

% [2a] Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
ppActBomb=(EnvThreat.*NTok) ./  (24-(EnvThreat.*NTok));
ppABscaled=power(ppActBomb, InOptim);

% [2b] Likelihood of Accepting, if one Explores and sees no Bomb
posVacc= (1-ppABscaled) .* NTok;
posVrej= ppABscaled .* NTok;
pAccGivNoBombSeen=posVacc>posVrej;
EVGivenNotSee = pAccGivNoBombSeen.* posVacc   +   (1-pAccGivNoBombSeen).* posVrej;  

% [3] V(Explore) =   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
V_ct(:,:,3) =  (pActBomb.*0.5) .* NTok  +  (1- pActBomb.*0.5) .* EVGivenNotSee   + ExploreCost;

%% Combine values for both tasks

Vv=[V_cF V_ct];
% V=Vv; 
% 
% Preserve trial order (as experienced by subject - other functions expect this)
Vv(nPP+1,:,:)=repmat([order_cF; order_ct]',[1 1 3]);
Vv(:,:,1)=sortrows(Vv(:,:,1)',nPP+1)';
Vv(:,:,2)=sortrows(Vv(:,:,2)',nPP+1)';
Vv(:,:,3)=sortrows(Vv(:,:,3)',nPP+1)';

V=nan(nPP,length(order_cF) + length(order_ct),3);
V(:,:,1)=Vv(1:nPP,:,1);
V(:,:,2)=Vv(1:nPP,:,2);
V(:,:,3)=Vv(1:nPP,:,3);

end


