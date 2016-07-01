
% and the function
% 


    
%%

clear all;
[pplant,ptoke]=meshgrid((6:-1:1)/6,(2:2:12)/12); %12*ptoke=ntoke


% Distort Envthreat 
nlpar=1;
pplant=nlp(pplant,nlpar);
    



pbomb=pplant.*ptoke;
sigma=inline('1./(1+exp(-x))');
beta=0.1;   % for softmax for the first choice
betaexp=1.0;% for softmax for the second choice
Vacc=-120*pbomb+120*ptoke.*(1-pbomb);
Vrej=0;

dacc=((1-pbomb)./(1-pbomb/2)).*ptoke*120+...
     -120*(pbomb/2)./(1-pbomb/2);  % V(accept) @ 2 if no bomb appears
dacc_undistorted=dacc';
 
dacc=dacc.*sigma(betaexp*dacc);  % allow for stage 2 mistakes
dacc_sig=dacc';
 %dacc=(dacc>0).*dacc;
 
Vexp=-20+(1-pbomb/2).*dacc; % include the cost of exploring


%%
pexp=exp(beta*Vexp);
pacc=exp(beta*Vacc);
prej=1./(1+pexp+pacc); % reject at stage 1 gives 0
pexp=pexp.*prej;
pacc=pacc.*prej;

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

figure(1);
clf;
subplot(2,3,1);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pplant');
colorbar;
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
title('prior bomb');
%axis square;
subplot(2,3,2);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),12*ptoke');
colorbar;
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
title('number tokens');
%axis square;
subplot(2,3,3);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pbomb');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
colorbar;
title('prob loss');
%axis square;
subplot(2,3,4);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pacc');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
colorbar;
title(['p accept; beta=' num2str(beta)]);
subplot(2,3,5);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),prej');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
%axis square;
colorbar;
title(['p reject; beta=' num2str(beta)]);
subplot(2,3,6);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pexp');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
%axis square;
colorbar;
title(['p explore; beta=' num2str(beta)]);

figure(2);
clf;
subplot(2,3,1);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pplant');
colorbar;
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
title('prior bomb');
%axis square;
subplot(2,3,2);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),12*ptoke');
colorbar;
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
title('number tokens');
%axis square;
subplot(2,3,3);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pbomb');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
colorbar;
title('prob loss');
%axis square;
subplot(2,3,4);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pcn');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
colorbar;
title(['p cno bomb; beta=' num2str(beta)]);
subplot(2,3,5);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pcb');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
%axis square;
colorbar;
title(['p cbomb; beta=' num2str(beta)]);
subplot(2,3,6);
imagesc(12*ptoke(:,1),pplant(1,6:-1:1),pce');
set(gca,'ytick',(1:6)/6,'yticklabel',1-(0:5)/6)
%axis square;
colorbar;
title(['p cexplore; beta=' num2str(beta)]);

%%

function ou=nlp(p,par);
    ou=par*(p.*p)+(1-par)*p;

