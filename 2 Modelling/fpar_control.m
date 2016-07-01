function [ data ] = fpar_control( data, col)
% [ data ] = fpar_control( data, col)
% Run through CONTROL trials and mark RL parameters for each trial
%   
%   Input variables:                     EnvThreat (E), NTokens (N)
%                                               OutcomeMagnitude (optional)
%   Additional output variables:    pLoss (P), Entropy (U), VExplore (X), Conflict (C)
%                                               (if OutcomeMagnitude available) OutcomeMean,  OutcomeVariance
%
%   Columns must be specified for all these variables (eg. col.VExplore)
%   Note: Conflict parameter is not yet specified!
%
% All parameter-specifications are the same as CONFLICT task, but the
% psychological interpretation changes. EV is the only parameter whose
% specification formula changes. 
%
%   Variables whose psychological interpretation changes:
%           - pLoss         -->         p(ActBomb)
%
% --------------------------------------------------------------------

plot=0;

% Setup 
FixedLoss=0;
ExploreCost=  -2;
EntropyCorrection=0.00001;

%% Mark RL Variables (straightforward)

% (1) EnvThreat and NTokens: Convert to proper variable values if necessary
if sum(data(:, col.NTokens)==12)==0
    data(:,col.NTokens)=2*data(:,col.NTokens);
end
if sum(data(:,col.EnvThreat)<1)==0; 
    data(:,col.EnvThreat)=data(:,col.EnvThreat)/6;
end

% (2) pLoss (actually pActBomb)
data(:,col.pLoss)=(data(:,col.EnvThreat)).*(data(:,col.NTokens)./12);

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

% (4) Expected value (differs from Conflict task)
% EV =    pLoss * NTok             if  pLoss * NTok > (1-pLoss) * NTok
%           (1-pLoss) * NTok        otherwise
% Old: EV = p(Correct)*V(Correct) + p(Wrong)*V(Wrong)
% data(:,col.EV)=(1-data(:,col.Entropy)).*data(:,col.NTokens);
data(:,col.EV)= max([data(:,col.pLoss).*data(:,col.NTokens) (1- data(:,col.pLoss)).*data(:,col.NTokens)], [],2); 

% (5) EV Gain/Loss
if isfield(col, 'EVGain')
    data(:,col.EVGain)= data(:,col.EV);
end
if isfield(col, 'EVLoss')
    data(:,col.EVLoss)=  0;
end
if isfield(col, 'EVConflict')
    data(:,col.EVConflict)= data(:,col.EVGain) - data(:,col.EVLoss);
end


% (6) Pavlovian conflict
if isfield(col, 'PavConflict')
    data(:,col.PavConflict)= data(:,col.NTokens).* 0;
end

%% Additional variables (that need thinking through)

% Value of exploration: Value gained from exploring
% Note: If specification here is changed, also change functions for re-specifying parameters (e.g. after warping EnvThreat)
if isfield(col, 'VExplore')
    
    % (a) Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
    %     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
    %                                      = (EnvThreat * NTokens) / (24 - EnvThreat*NTokens)
    ppActBomb=(data(:, col.EnvThreat).*data(:, col.NTokens))  ./  (24-(data(:, col.EnvThreat).*data(:, col.NTokens)));
    
    % (b) EV Accept/Reject given Explore 
    EVAcc = data(:, col.NTokens) .* (1-ppActBomb);
    EVRej = data(:, col.NTokens) .* ppActBomb;
    EVGivenNotSee = (EVAcc >= EVRej) .* EVAcc + (EVRej > EVAcc) .* EVRej;
    
    % (c) V(Explore) =   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
    EV_Explore=  data(:, col.pLoss).*0.5 .* data(:, col.NTokens)  +  (1-(data(:, col.pLoss)*0.5)) .*EVGivenNotSee  + ExploreCost;
    
    % (d) Value of Exploring =EV(Explore) - EV(NonExplore)
    %       EV(NonExplore)=EV(Accept) if >0, otherwise 0
    data(:,col.VExplore)=EV_Explore-data(:,col.EV);
    
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
%     * NOTE: Code here must agree with fct_meanvar_quantities
% Outcomes in ct: Correct: p=(1-Entropy),  Mag= NTokens
%                          Incorrect: p=Entropy, Mag=0
OutcomeMean=data(:,col.EV); %  EV = p(Correct)*V(Correct) + p(Wrong)*V(Wrong)
StanDev=sqrt(  (1-data(:, col.Entropy)).*(data(:, col.NTokens)-OutcomeMean).^2  +  data(:, col.Entropy).*(0 -OutcomeMean).^2);
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
    data(:,col.Conflict)=data(:, col.Entropy) .*(    (1-data(:, col.Entropy)).*data(:, col.NTokens)  - data(:, col.Entropy).*0     );
end


% Miscellaneous quantity (for plotting)
if isfield(col, 'Misc')
    pls=data(:,col.pLoss);  % Rescale
    if min(pls(:))<0, pls= pls- min(pls(:));
    else pls=  pls - min(pls(:));
    end
    pls=pls./max(pls(:));
    ens=data(:,col.Entropy);
    if min(ens)<0, ens= ens- min(ens(:));
    else ens=  ens- min(ens(:));
    end
    ens=ens./max(ens(:));
    
    data(:,col.Misc)= pls-ens;
%     data(:,col.Misc)= data(:,col.pLoss) - data(:,col.Entropy);
end

%% Outcome statistics (optional)

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