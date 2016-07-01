% Random bits of code from modelling 



%%

% Trying to understand the power law - what do values of k mean? [y=x^k]
k=1; 
x=0:0.1:1;
y=nan(length(x),1); yy=y; yyy=y;
for i=1:length(x)
    y(i)=x(i)^k;
    yy(i)=x(i)^(2*k);
    yyy(i)=x(i)^(-6*k);
end
figure, scatter(x,y); hold on; scatter( x,yy);
hold on; scatter( x(2:end),yyy(2:end));clc

% What does the power law look like? Say a range of k from 1/15 and 15 is reasonable
x=0:0.001:1; Lx = length(x);
uo = logspace(-1.2,1.2,12); 
uo =1/15:  (15-1/15)/10: 15; % evenly spaced in numerical space --> doesn't cover whole space of power law 
uo=log10(1/15) :   (log10(15)-log10(1/15))/10 :  log10(15);
uo=logspace(  log10(1/15)  ,  log10(15), 10);  % How to scale it properly? Say you want evenly spaced intervals with k=[1/5 15] as a reasonable range
Luo = length(uo);
xm = repmat(x', [1 Luo]);
uom = repmat(uo, [Lx 1]);
y = power(xm,uom);
figure, plot(x,y)
%
% set(gca, 'YScale', 'log')



% Sigmoid?
pp=0:0.1:1;
% y=1./(1+exp(-8*exp(1)*(pp) +1.5.*exp(1) ))

y=1./(1+exp(-10*(pp) +  3.*exp(1) ));
close all hidden; figure; plot(pp,y); hold on; plot(0:0.1:1,0:0.1:1);




%% Grid + map nLL 

load('D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\All data (18-Mar-2014).mat')

beta=logspace(-6,1,12);
epsilon = 0.05:0.02:0.332;

beta = log(beta);
epsilon=-log( 1./(3.*epsilon) -1 );

% Quick gridding + look at nLL surface. Highest nLL value should be nTrials*log(1/3)  [ 1/3 =chance]
probbias = -3:0.2:2;
bll = nan([length(beta) length(epsilon) length(probbias)]);
for iB=1:length(beta)
    disp(beta(iB))
    for iE=1:length(epsilon)
        for iP = 1:length(probbias)
            bll(iB,iE,iP) = f_nllsoftmax_lapse([beta(iB) epsilon(iE) probbias(iP)] ,  {'bpa01' subjdata{1,2}  details.fixedpar  details.col});
        end
    end
end
figure, surf(exp(beta), 1./(3+3.*exp(-epsilon)), squeeze(bll(:,:,10))'), xlabel('beta'), ylabel('epsilon')


% Trinomial softmax?
L= length(-5:0.1:5);
V1=repmat((-5:0.1:5)', [1 L L]);
V2=repmat(-5:0.1:5, [L 1 L]);
V3=repmat(shiftdim(-5:0.1:5,-1), [L L 1]);
sm = exp(V1)./(exp(V1)+exp(V2)+exp(V3));
figure, surf(-5:0.1:5, -5:0.1:5, squeeze(sm(:,:,find((-5:0.1:5)==5))))  % Change V(3), see how the choices change

