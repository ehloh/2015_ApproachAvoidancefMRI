% Changes I've made to Peter's code to make it like mine (after which, it
% does produce the same L as my bjm16)
% 
% - w applied at the same point as e.
% - beta changed to accomodate that hes working w points and im working w tokens
% - 2nd stage choice is now made deterministic
% - j and beta constraints changed - Peter's is stricter
% - Peter has a bug in his m parameter implementation in V(Accept 2nd
%   stage). m is applied to both posterior p(bomb) and p(no bomb)
% - eps distortion to p(Accept/Rej/Explore) (like the softmax lapse) turned off

function [pacc,prej,pexp,pcb,pcn,pce]=efitter(X)

% sigma=inline('1./(1+exp(-x))');
sigma=inline('5./(1+exp(-x))');  % Laxer constraint

% Email dated 23/12/14, originally named efitter.
%       Order of parameter inputs changed to match mine (bjm16_feow)
%       Original order of X: j, b, f (percentage), e, m, ow
% nlpar=sigma(X(1));   % linear (nlpar=0) to quadratic (nlpar=1) warp for env threat
% beta=exp(X(2));      % softmax temperature for 1st choice
% betaexp=1.0;         % softmax for 2nd choice - fixed, not a parameter
% cst=120*X(3);        % cost of an explosion (*10 in terms of tokens)
% expbon=10*X(4);      % excess cost/benefit of exploration
% powm=exp(X(5));      % Eleanor's power-law probability warp at the 2nd stage
% uncbon=0;
% uncbon=X(6);         % multiplier of token*entropy (my least favourite parameter)

if length(X)<6; X(length(X)+1:6)=nan; disp('control task'); end  % control task

% Order of X changed to match my code ## 
nlpar=0;
nlpar=sigma(X(2));  % j 
beta=exp(X(1));  % b
betaexp=1.0;  % beta for exploration assumed to be 1
% cst=120*X(4);   % f (percentage)
cst=  - X(4)*10;  % f is positive in dayans code
expbon=10*X(5);   % e (converted to points currency)
powm=exp(X(3));  %  m
uncbon=0;
uncbon=X(6);  % o  (converted to points currency)


beta=beta/10;  % compensate beta to be used w points (not n tokens)
disp(['[Dy] Parameter values: ']); %  num2str( [beta nlpar powm -cst expbon uncbon],3)])
disp([beta nlpar powm -cst expbon uncbon])


[pplant,ptoke]=meshgrid((6:-1:1)/6,((2:2:12))/12); %12*ptoke=ntoke
% pplant=nlp(pplant,nlpar); % #Altered #######
pplant=power(pplant,nlpar); 
pbomb=pplant.*ptoke;
hpb=-pbomb.*log(pbomb+1e-32)-(1-pbomb).*log(1-pbomb+1e-32); % Entropy
Vacc=-cst*pbomb+120*ptoke.*(1-pbomb);
Vrej=0;


% Peter's power distortion is wrong - it is applied to both posterior  p(Bomb) and posterior p(No bomb)
dacc=(1-  (((pbomb/2)./(1-pbomb/2)).^powm)).*ptoke*120+...
     -cst*((pbomb/2)./(1-pbomb/2)).^powm;  % V(accept) @ 2 if no bomb appears
% dacc=(    (     (1-pbomb)./(1-pbomb/2)   )      ).^powm.*ptoke*120+...
%      -cst*((pbomb/2)./(1-pbomb/2)).^powm;  % V(accept) @ 2 if no bomb appears
% dacc=dacc+hpb.*ptoke*120*uncbon;
%  Veacc=dacc.*sigma(betaexp*dacc);  % allow for stage 2 mistakes
% dacc=(dacc>0).*dacc;  

% Changes made (compared to Dayan's original efitter): exploration bonuses
% applied only at the end, deterministic 2nd-stage choice rather than sigma-distorted
Veacc=dacc.* (dacc>0);  % deterministic 2ndstage 
Vexp=-20+(1-pbomb/2).*Veacc+expbon; % Add fixed exploration bonus
Vexp=Vexp+hpb.*ptoke*120*uncbon;  % Add variable exploration bonus ( ow )

pacc=exp(beta*Vacc);
pexp=exp(beta*Vexp);
prej=1./(1+pexp+pacc); % reject at stage 1 gives 0
pexp=pexp.*prej;
pacc=pacc.*prej;

% pacc=(1-eps)*pacc+eps/3;
% prej=(1-eps)*prej+eps/3;
% pexp=(1-eps)*pexp+eps/3;

%%


Vcb=120*pbomb.*ptoke; % these are now for the control task [bjm01]
Vcn=120*(1-pbomb).*ptoke;
dcb=120*(pbomb/2)./(1-pbomb/2).*ptoke;
dcn=120*((1-pbomb)./(1-pbomb/2)).*ptoke;
% pcacc=sigma(betaexp*(dcb-dcn));
pcacc=dcn<dcb;   % Made deterministic + corrected bug
dcacc=pcacc.*dcb+(1-pcacc).*dcn;
Vcexp=-20+(1-pbomb/2).*dcacc+120*(pbomb/2).*ptoke;  % no exploration bonuses for ct

pcb=exp(beta*Vcb);
pcn=exp(beta*Vcn);
pce=exp(beta*Vcexp);
psum=1./(pcb+pcn+pce);
pcb=pcb.*psum;
pcn=pcn.*psum;
pce=pce.*psum;


% Output What]





% pcb=(1-eps)*pcb+eps/3;
% pcn=(1-eps)*pcn+eps/3;
% pce=(1-eps)*pce+eps/3;

%% end to dayans code

disp('Output Value rathre than probability');
pacc=Vacc; prej=Vrej; pexp=Vexp;
pcb=Vcb; pcn=Vcn; pce=Vcexp;

