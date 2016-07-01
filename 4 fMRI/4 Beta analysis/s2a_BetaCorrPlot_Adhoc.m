
% % AD HOC PLOTS ------------------------------
% beh='per.cF_Reject';
beh='Trait anxiety scores';
% beh='State anxiety scores';
% beh='m.cF_i';
% beh='Probability of following null information\n from Exploring in experimental task';
% beh='BIS scores';
% beh='infoexnull.all';
% beh='per.NotAccept_cF-ct';
% beh='per.Accept_ct-cF';
% beh='per.cF_Accept';
% beh='infoexnull.ct'; 
% beh='RT for Rejecting (Exp condition)';
% beh='Difference in % Rejecting (Exp>Ctrl)'; 

% rrrr={'per.Accept_ct-cF';}
% rrrr={'RT for Rejecting (Exp condition)';}
% rrrr={'Probability of following null information\n from Exploring in experimental task'};

rrrr={'HPC_aL_tv001-NonRej_vMargCho_cF-ct'}
% rrrr={'HPC_aR_tv001-cF_vMargCho_NonRej-Rej'}
    

roicon=rrrr{1}; wr.rc= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),roicon)) ));  wr.b= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),beh)) ));
f.scatter_dotsize=100;   f.scatter_linewidth=4;   f.FontSize=25; f.fontname='PT Sans Caption';   
% figure('color', 'w')
scatter(wr.rc, wr.b, f.scatter_dotsize,'LineWidth', 3); h=lsline; set(h,'LineWidth', f.scatter_linewidth); xlabel(sprintf(roicon),'FontSize',20), ylabel(sprintf(beh),'FontSize',20)
if      request.kendalls_correlation==0, [r p]= corr(wr.rc, wr.b); title(['r='  num2str(r) ', p=' num2str(p)],'FontSize',f.FontSize)
else      [r p]= corr(wr.rc, wr.b,'type', 'Kendall'); title(['tau='  num2str(r) ', p=' num2str(p)],'FontSize',f.FontSize)
end; set(gca,'FontSize',f.FontSize)



xlabel(sprintf('Difference in left hippocampal parameter estimates\n  on Non-Reject trials (Exp>Ctrl)'), 'FontSize',f.FontSize,'FontName', 'PT Sans Caption')
xlabel(sprintf('Difference in right hippocampal parameter \n  estimates (Exp, Non-Reject>Reject)'), 'FontSize',f.FontSize,'FontName', 'PT Sans Caption')
ylabel(sprintf(beh),'FontSize',f.FontSize)

request.kendalls_correlation=0; 
request.kendalls_correlation=1; 


%% MEDIAN split 


% split_iv={'Probability of following null information\n from Exploring in experimental task'}; 
split_iv={'State anxiety scores'};

wm.ivbeh=cell2mat(d_betas(2:end, strcmp(d_betas(1,:), split_iv))) ; 
sum(wm.ivbeh >median(wm.ivbeh))
wm.subs_high= find(wm.ivbeh > median(wm.ivbeh));   % HIGH score group 
wm.subs_low= find(1- (wm.ivbeh >median(wm.ivbeh)));   % HIGH score group 


wm.subs_low



