function [pacc,prej,pexp,pcb,pcn,pce]=efitter(X);

nlpar=0;
nlpar=sigma(X(1));
beta=exp(X(2));
betaexp=1.0;
cst=120*X(3);
expbon=10*X(4);
powm=exp(X(5));
uncbon=0;
uncbon=X(6);

[pplant,ptoke]=meshgrid((6:-1:1)/6,((2:2:12))/12); %12*ptoke=ntoke
pplant=nlp(pplant,nlpar);
pbomb=pplant.*ptoke;
hpb=-pbomb.*log(pbomb+1e-32)-(1-pbomb).*log(1-pbomb+1e-32);
Vacc=-cst*pbomb+120*ptoke.*(1-pbomb);
Vrej=0;

dacc=(((1-pbomb)./(1-pbomb/2))).^powm.*ptoke*120+...
     -cst*((pbomb/2)./(1-pbomb/2)).^powm;  % V(accept) @ 2 if no bomb appears
dacc=dacc+hpb.*ptoke*120*uncbon;
 Veacc=dacc.*sigma(betaexp*dacc);  % allow for stage 2 mistakes
%dacc=(dacc>0).*dacc;
 
Vexp=-20+(1-pbomb/2).*Veacc+expbon; % include the cost of exploring

pexp=exp(beta*Vexp);
pacc=exp(beta*Vacc);
prej=1./(1+pexp+pacc); % reject at stage 1 gives 0
pexp=pexp.*prej;
pacc=pacc.*prej;

pacc=(1-eps)*pacc+eps/3;
prej=(1-eps)*prej+eps/3;
pexp=(1-eps)*pexp+eps/3;

Vcb=120*pbomb.*ptoke; % these are now for the control task
Vcn=120*(1-pbomb).*ptoke;
dcb=120*(pbomb/2)./(1-pbomb/2).*ptoke;
dcn=120*((1-pbomb)./(1-pbomb/2)).*ptoke;
pcacc=sigma(betaexp*(dcb-dcn));
dcacc=pcacc.*dcb+(1-pcacc).*dcn;
Vcexp=-20+(1-pbomb/2).*dcacc+120*(pbomb/2).*ptoke;

pcb=exp(beta*Vcb);
pcn=exp(beta*Vcn);
pce=exp(beta*Vcexp);
psum=1./(pcb+pcn+pce);
pcb=pcb.*psum;
pcn=pcn.*psum;
pce=pce.*psum;

pcb=(1-eps)*pcb+eps/3;
pcn=(1-eps)*pcn+eps/3;
pce=(1-eps)*pce+eps/3;

