% Frame to compare L from my vs dayans value functions
clear all;  clc;

task=1; model='bjm16_feow'; modelpars={'b';'j';'m';'f';'e';'w'};
% task=2; model='bj01'; modelpars={'b';'j';};

for o=1:1 % Set up
%     where.where='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
    where.where='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
    where.inputs=[where.where filesep '2 Analysis inputs'];    
    path(pathdef); addpath(where.where), addpath([where.where filesep '1 Value functions Det Jpower']);
    addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bjm']);
    addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bj']);
    addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bpm']);
    addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bm']);    
    
    % Setup on my side
    d_ts=load([where.inputs filesep 'All data (09-May-2014).mat']);
    dcol=d_ts.details.col; d_ts=d_ts.subjdata;
    details.fixedpar.cF_FL=  -12;  % Task parameters
    details.fixedpar.cF_EC=  -2;
    details.fixedpar.ct_FL=0;
    details.fixedpar.ct_EC=  -2;
end
for o=1:1 % Hard coded somethings, 'hc'
%     where.res='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\Det Jpower';
    where.res='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/Det Jpower';
    if task==1; rr=load([where.res filesep 'res_fitmodels_cF (02-Dec-2014) j.mat']);
    else rr=load([where.res filesep 'res_fitmodels_ct (02-Dec-2014) j.mat']); end
    
    hc.dayans_nll_bjm16={ % wrong anyway
        'p01_GV'  306.4927
        'p02_YY'  216.4332
        'p04_MP'  131.4657
        'p06_KB'  252.2438
        'p08_SG' 285.8671
        'p10_RC'  275.2847
        'p13_HL'  217.2757
        'p17_SJ'  180.3606
        'p18_MS'  196.2476
        'p21_ES'  373.4595
        'p23_BS'  152.9958
        'p25_RJ'  291.4573
        'p27_DM'  420.4694
        'p30_KL'  326.2263
        'p34_TB'  254.9597
        'p35_SM'  112.4145
        'p36_FR'  331.8405
        'p38_MK'  228.7145
        'p41_AC'  249.9471
        };
end

%%

hc.mymodres=rr.r_res{strcmp(rr.r_res(:,1), model),2};
% hc.mymodres=rr.r_res{strcmp(rr.r_res(:,1), model),2}(find(1-strcmp(rr.details.subjects, 'p15_SH')),:); 
hc.modL=hc.mymodres(:,2); hc.subpars=hc.mymodres(:,4:end);

Lok=nan(size(hc.subpars,1),3);
for s=1:size(hc.subpars,1)
    disp(s)

    % Artificial inputs
    hc.subpars(s,:)=[hc.subpars(s,1) 1-eps 1-eps -12 0 0];   disp('All parameters except beta turned off')
    hc.subpars(s,:)=[10 1 1 -12 0 0];   disp('All parameters  turned off')
%         hc.subpars(s,:)=[hc.subpars(s,1) 1-eps 1-eps];   disp('All parameters except beta turned off')
    % d_ts{s,task+1}=d_ts{s,task+1}(1:100,:);
%     hc.subpars(s,1:3)=[1.2 1 1];
    % hc.subpars(s,2)=1-eps;
%     hc.subpars(s,4)=-12;
%     hc.subpars(s,5:6)=0;

    % Feed thru my code
    ws.origparvals=hc.subpars(s,:);
    ws.parvals=f_transpar(modelpars, ws.origparvals, 'from');
    [ws.mynll pch]=f_nllsoftmax(ws.parvals, {model d_ts{s,task+1} details.fixedpar dcol}); pch=pch';
    disp(['[My] Parameter values: ' num2str(ws.origparvals,3)])
    disp(['          ws.mynll:' num2str(ws.mynll )])

    % Feed thru Dayan's code
    ws.parvals_dayan=ws.origparvals; 
    ws.parvals_dayan(strcmp(modelpars,'j')) =   - log(  (5/  ws.origparvals(strcmp(modelpars,'j'))   )      -1  );  % Inverse parameter transforms 
    ws.parvals_dayan(strcmp(modelpars,'b'))=log(ws.origparvals(strcmp(modelpars,'b')));
    ws.parvals_dayan(strcmp(modelpars,'m'))=log(ws.origparvals(strcmp(modelpars,'m')));
    [pacc,prej,pexp,pcb,pcn,pce]=efitter_makelikemine(ws.parvals_dayan); % p(Each Choice)
    ws.pchoice{1}{1}=pacc'; ws.pchoice{1}{2}=prej'; ws.pchoice{1}{3}=pexp';  
    ws.pchoice{2}{1}=pcn'; ws.pchoice{2}{2}=pcb'; ws.pchoice{2}{3}=pce'; % return to my frame of ref (env x ntok framing)
    ws.dd=d_ts{s,task+1}(:, [dcol.Choice dcol.EnvThreat dcol.NTokens]); % Eval choice according to dayans calcs    
     ws.dd(:, 2)=ws.dd(:, 2)*6; ws.dd(:, 3)=ws.dd(:, 3)/2;  % choice, e, n, p(observed choice)
    for t=1:size(ws.dd,1)
        ws.dd(t,4)=ws.pchoice{task}{ws.dd(t,1)}(  7-ws.dd(t,2) ,   ws.dd(t,3));
    end
    ws.dynnll=sum(-log(ws.dd(:,4))); disp(ws.dynnll)

    % Compare my ws.mynll to Dayans 
    Lok(s,1:2)=[ws.mynll ws.dynnll];
    if abs(ws.mynll  -ws.dynnll)<0.001, Lok(s,3)=1; disp('            L ok'); else disp('            nLLs don''t match!'); Lok(s,3)=0; end 
    disp([ws.mynll    ws.dynnll])
    
    
    for o=1:1 % What's the best choice?
        for t=1:2
            
            for e=1:6
                for n=1:6
                    ws.pc=[ws.pchoice{t}{1}(e,n)  ws.pchoice{t}{2}(e,n) ws.pchoice{t}{3}(e,n)];
                    if sum(max(ws.pc)==ws.pc)==1
                        ws.bestchoice{t}(e,n)=find(max(ws.pc)==ws.pc);
                    elseif  sum(find(max(ws.pc)==ws.pc))==3 % A + R
                        ws.bestchoice{t}(e,n)=1.5;
                    elseif  sum(find(max(ws.pc)==ws.pc))==4 % A + E
                        ws.bestchoice{t}(e,n)=1.8;
                    elseif  sum(find(max(ws.pc)==ws.pc))==5 % E + R
                        ws.bestchoice{t}(e,n)=2.5;
                    end
                end
            end
        end
        bestcf= ws.bestchoice{1};
        bestct= ws.bestchoice{2};
%         figure('color','w')
%         subplot(1,2,1)
%         imagesc(ws.bestchoice{1}); axis square, colorbar
%         title('cF Best choice'); ylabel('EnvThreat');  xlabel('No. Tokens');  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
%         subplot(1,2,2)
%         imagesc(ws.bestchoice{2}); axis square, colorbar
%         title('ct Best choice');ylabel('EnvThreat');  xlabel('No. Tokens');  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
    
    
    end
    
    error
    
    ws=[];

end

%%


% Output values to compare
dopt=[
{pacc/10}
{prej/10}
{pexp/10}
{pcn/10}
{pcb/10}
{pce/10}];  

dp=d_plotmatrix(1,2:end);



dopt{5} - dp{4}'