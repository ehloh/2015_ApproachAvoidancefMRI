function [ outvars] = fct_changeEnvThreat(EnvThreat, NTok, miscvar)
% Given EnvThreat + NTok vectors, output all other parameter vectors
%   
%   Input variables:                     EnvThreat (E), NTokens (N)
%                                               OutcomeMagnitude (optional)
%   Additional output variables:    pLoss (P), Entropy (U), VExplore (X), EntropyNTok (K)
%                                               (if OutcomeMagnitude available) OutcomeMean,  OutcomeVariance
%   
%   Columns for EnvThreat, NTokens, pLoss & Entropy must be specified. All
%           other variables only calculated if specified in variable col.

try  FixedLoss=miscvar.FixedLoss; catch;  FixedLoss= -12; end
try  ExploreCost=miscvar.ExploreCost; catch;  ExploreCost= -2; end
EntropyCorrection=0.00001;

%%

% (1) EnvThreat and NTokens: Convert to proper variable values if necessary
if sum(NTok==12)==0; NTok=2*NTok;
end
if sum(EnvThreat<1)==0;  EnvThreat=EnvThreat/6;
end

% (2) pLoss
pLoss=(EnvThreat).*(NTok./12);

% (3) Entropy
pLossEn=pLoss;  pLossEn(pLossEn==1)=  pLossEn(pLossEn==1)- EntropyCorrection;
Entropy= -(pLossEn).*log(pLossEn) -   (1-pLossEn).*log(1-pLossEn);

% (4) EntropyNTok
EntropyNTok=Entropy.*NTok;

% (5) EV (differs from Conflict task)
% EV=(1-Entropy).*NTok;
EV= max([pLoss.*NTok (1- pLoss).*NTok], [], 2);

% (6) VExplore   ##############################

% (a) Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
%     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
%                                      = (EnvThreat * NTokens) / (24 - EnvThreat*NTokens)
ppActBomb=(EnvThreat.*NTok)  ./  (24-(EnvThreat.*NTok));

% (b) EV Accept/Reject given Explore 
EVAcc = NTok .* (1-ppActBomb);
EVRej = NTok .* ppActBomb;
EVGivenNotSee = (EVAcc >= EVRej) .* EVAcc + (EVRej > EVAcc) .* EVRej;

% (c) V(Explore) =   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
EV_Explore=  pLoss.*0.5 .* NTok  +  (1-(pLoss*0.5)) .*EVGivenNotSee  + ExploreCost;

% (d) Value of Exploring =EV(Explore) - EV(NonExplore)
%       EV(NonExplore)=EV(Accept) if >0, otherwise 0
VExplore=EV_Explore-EV;

%%  More esoteric quantites

EVGain=EV;
EVLoss=EV.*0;
EVConflict=  EVGain - EVLoss; 
PavConflict=EV.*0; 

%% Output

outvars.EnvThreat =EnvThreat ;
outvars.NTok =NTok ;
outvars.NTokens =NTok ;
outvars.pLoss =pLoss;
outvars.Entropy=Entropy;
outvars.EntropyNTok =EntropyNTok;
outvars.VExplore=VExplore;
outvars.EV=EV;
outvars.EVGain=EVGain;
outvars.EVLoss=EVLoss;
outvars.EVConflict=EVConflict;
outvars.PavConflict=PavConflict;



end

