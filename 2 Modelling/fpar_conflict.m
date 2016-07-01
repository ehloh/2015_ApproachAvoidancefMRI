function [ data ] = fpar_conflict( data, col)
% [ data ] = fpar_conflict( data, col)
% Run through CONFLICT trials and mark RL parameters for each trial
%   
%   Input variables:                     EnvThreat (E), NTokens (N)
%                                               OutcomeMagnitude (optional)
%   Additional output variables:    pLoss (P), Entropy (U), VExplore (X), EntropyNTok (K)
%                                               (if OutcomeMagnitude available) OutcomeMean,  OutcomeVariance
%   
%   Columns for EnvThreat, NTokens, pLoss & Entropy must be specified. All
%           other variables only calculated if specified in variable col.
%
% --------------------------------------------------------------------------------------------------------------

plot=0;

% Setup
FixedLoss= -12;
ExploreCost= - 2;
EntropyCorrection=0.00001;

%% Mark RL Variables (straightforward)

% (1) EnvThreat and NTokens: Convert to proper variable values if necessary
if sum(data(:, col.NTokens)==12)==0
    data(:,col.NTokens)=2*data(:,col.NTokens);
end
if sum(data(:,col.EnvThreat)<1)==0; 
    data(:,col.EnvThreat)=data(:,col.EnvThreat)/6;
end

% (2) pLoss
data(:,col.pLoss)=(data(:,col.EnvThreat)).*(data(:,col.NTokens) ./12);

% (3) Entropy
for i=1:size(data,1) % pLoss values  of 1 = corrected by -EntropyCorrection for entropy calculation
    
    % Convert if pLoss == 1 
    if data(i,col.pLoss)==1 
        data(i,col.pLoss)=data(i,col.pLoss)-EntropyCorrection;
    end
      
    % Entropy calculation:  -p(loggp)-(1-p)logg(1-p) 
    data(i,col.Entropy)= -(data(i,col.pLoss)).*log(data(i,col.pLoss)) -   (1-data(i,col.pLoss)).*log(1-data(i,col.pLoss));        
        
    % Convert pLoss back if it's been altered
    if data(i,col.pLoss)==1-EntropyCorrection;
        data(i,col.pLoss)=data(i,col.pLoss)+EntropyCorrection;
    end
    
end

% (4) Expected value
%       Value of gamble if one accepts only
data(:,col.EV)=data(:,col.pLoss).*FixedLoss + (1-data(:,col.pLoss)).*data(:,col.NTokens);


% (5) EV Gain/Loss
if isfield(col, 'EVGain')
    data(:,col.EVGain)= (1-data(:,col.pLoss)).*data(:,col.NTokens);
end
if isfield(col, 'EVLoss')
    data(:,col.EVLoss)=  data(:,col.pLoss).*FixedLoss;
end
if isfield(col, 'EVConflict')
    data(:,col.EVConflict)= data(:,col.EVGain) - data(:,col.EVLoss);
%     data(:,col.EVConflict)= abs(data(:,col.EVGain)).* abs(data(:,col.EVLoss));
%     data(:,col.EVConflict)= data(:,col.NTokens).* abs(data(:,col.EVLoss));
end


% (6) Pavlovian conflict
if isfield(col, 'PavConflict')
    data(:,col.PavConflict)= data(:,col.NTokens).* abs(data(:,col.pLoss));
end


%% Additional variables (that need thinking through)

% Value of exploration: Value gained from exploring
%       EV(AfterExplore) - EV(BeforeExplore)   [taking into account posterior probabilities after Exploring)
%       p(Accept/Reject after Explore) as binary, depending on which EV is bigger
% Note: If specification here is changed, also change functions for re-specifying parameters (e.g. after warping EnvThreat)
if isfield(col, 'VExplore')
    
    % (a) Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
    %     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
    %                                      = (EnvThreat * NTokens) / (24 - EnvThreat*NTokens)
    ppActBomb=(data(:, col.EnvThreat).*data(:, col.NTokens)) ./  (24-(data(:, col.EnvThreat).*data(:, col.NTokens)));
    
    % (b) Likelihood of Accepting, if one Explores and sees no Bomb
    %   Cannot assume heuristic Accept/Reject deterministically from Explored-Info - will produce erroneously low V(Explore) values.
    %   For now, can assume deterministic post-Exploration behaviour:  deterministically  Accept/Reject
    %       depending on which has higher value. Can also calculate p(Accept) etc on the basis of the softmax betas.
    pAccGivNoBombSeen = (1-ppActBomb)   .* data(:, col.NTokens)>=   ppActBomb .* FixedLoss;
    
    % (c) EV(Explore)=   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
    %                       = p(See)*V(See)
    %                         +    p(~See) *  [       p(Accept|~See) * { posterior p(ActBomb) * Fixed Loss    +   (1- posterior p(ActBomb))   *  NTokens   }      +     p(Reject|~See) * 0      ]
    %                         +   Exploration Cost
    %     p(See Bomb) * V(See Bomb) = 0  ( Assume Reject, thus V(See)=0)
    %     p(Accept|~See): Deterministic/binary, if V(Accept)>V(Reject) based on posterior probabilities
    %     V(~See):   pAcceptGiven~See * [   posterior pActBomb * Fixed Loss + (1- posterior pActBomb )*NTokens  ]    + pRejectGiven~See * 0
    V_AccGivNoSee=    ppActBomb  .* FixedLoss           +            (1-ppActBomb ) .* data(:,col.NTokens);
    EV_NoSee= pAccGivNoBombSeen .*  V_AccGivNoSee  + (1-pAccGivNoBombSeen) .* 0;
    EV_Explore= 0 +  (1-(data(:, col.pLoss).*0.5)) .* EV_NoSee +  ExploreCost;
    
    % (d) Value of Exploring =EV(Explore) - EV(NonExplore)
    %       EV(NonExplore)=EV(Accept) if >0, otherwise 0
    EV_NonExplore=data(:,col.EV).*(data(:, col.EV)>=0);
    data(:, col.VExplore)=EV_Explore-EV_NonExplore;
end

% EntropyNToks
if isfield(col, 'EntropyNTok')
    data(:, col.EntropyNTok)=data(:, col.Entropy).* data(:, col.NTokens);
end

% EntropyEV
if isfield(col, 'EntropyEV')
    data(:, col.EntropyEV)=data(:, col.Entropy).* data(:, col.EV);
end

% Theoretical Outcome variance/std dev and value according to mean-variance theory (D'acremont & Bossaerts 2008 CABN)
%     meanoutcome= Sum_forNoutcomes [ probability(outcome) * magnitude(outcome) ]
%     Unbiased sample variance= [1/(N-1)] .* Sum_forNoutcomes [    probability(outcome) * (outcome - meanoutcome)^2   ]
%     Std Deviation= sqrt( Unbiased sample variance )
%     * NOTE: Code here must agree with fcf_meanvar_quantities
OutcomeMean=data(:,col.EV);
StanDev=sqrt(   (1-data(:, col.pLoss)).*(data(:, col.NTokens)-OutcomeMean).^2 +  data(:, col.pLoss).*(FixedLoss -OutcomeMean).^2   );
if isfield(col,'StanDev')
    data(:, col.StanDev)=StanDev;
end
if isfield(col,'vMeanVar')
    data(:, col.vMeanVar)=  OutcomeMean - StanDev;   % V=mean - b*standev
end

% Binomial pActBomb variance
if isfield(col,'BinomVar')
    data(:, col.BinomVar)=data(:,col.pLoss).*(1-data(:,col.pLoss));
end

% Conflict: Uncertainty * (EVpositive - EVnegative)
if isfield(col, 'Conflict')
    data(:,col.Conflict)=data(:,col.Entropy).*((1-data(:,col.pLoss)).*data(:,col.NTokens) -   data(:,col.pLoss).*FixedLoss);
end


% Miscellaneous quantity (for plotting)
if isfield(col, 'Misc')
%     pls=data(:,col.pLoss);  % Rescale
%     if min(pls(:))<0, pls= pls- min(pls(:));
%     else pls=  pls - min(pls(:));
%     end
%     pls=pls./max(pls(:));
%     ens=data(:,col.Entropy);
%     if min(ens)<0, ens= ens- min(ens(:));
%     else ens=  ens- min(ens(:));
%     end
%     ens=ens./max(ens(:));
%     
% %     data(:,col.Misc)= pls-ens;
%     data(:,col.Misc)= data(:,col.pLoss) .* data(:,col.Entropy);
    data(:,col.Misc)= data(:,col.Entropy) .* data(:,col.EV);
end

%% (Experienced) outcome statistics (optional)

if isfield(col, 'OutcomeMagnitude')
    for e=1:6
        for n=2:2:12
            wo.dat=data(data(:, col.EnvThreat)==e & data(:, col.NTokens)==n,col.OutcomeMagnitude);
            data(data(:, col.EnvThreat)==e & data(:, col.NTokens)==n, col.OutcomeMean)=mean(wo.dat);
            data(data(:, col.EnvThreat)==e & data(:, col.NTokens)==n, col.OutcomeVariance)=var(wo.dat);
            wo=[];
        end
    end
end


%% Plot if requested

if plot==1
    fpar_plotvars(data, col);
end


end


