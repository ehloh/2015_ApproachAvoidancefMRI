% Get results from fminunc model fitting, for manual manipulation
clear all; close all hidden; clc

for o1=1:1 % Documentation r_res
% variables (results columns in 'rc')
% 'r_res'
%       Col 1: Model name
%       Col 2: Subject fit parameters
%             Col i: BIC
%             Col ii: nLL
%             Col iv onwards: parameters (beta first)
%       Col 3: Model BIC (summed across subjects)
%       Col 4: Hessians
%       Col 5: 
%       Col 6 onwards: mean parameter values
path(pathdef)
addpath('2 Analysis inputs');  sdata=load(['2 Analysis inputs' filesep 'All Data (09-May-2014).mat']);
end

% Fit settings
Valfxn_type=' Det'; addpath(['2 Analysis inputs' filesep 'Det'])

% Which fits?
rrr=load(['res_fitmodels_cF (18-Jul-2014) all20iter.mat']); mod='bpmi14_feuw';  % cF
% mod='bpmi16_feow';
mod='bpm16_feow';
% mod='bpmi08_fuw';
% rrr=load(['res_fitmodels_ct (18-Jul-2014) all20iter.mat']);  mod='bpmi11_euw'; % ct

% Get details & model pars
disp([ rrr.details.task  '  ' mod]);
d=rrr.details.models(strcmp(rrr.details.models(:,1),mod),:);    % model details: d
rr=rrr.r_res(strcmp(rrr.r_res(:,1),mod),:);         % group results: rr
r=rr{1,2};                                                      % sub results: r
dcol=rrr.details.col;

% Alter parameters
% r(r(:,3+4)<-20, 3+4)= -20;  rr(6:end)=num2cell(mean(r(:,4:end))); disp('Cut out extreme values of a')  
% rr{5+ find(strcmp(d{3}, 'm'))}=1; disp('Set a to 1');
% % rr{5+ find(strcmp(d{3}, 'a'))}=0; disp('Set a to 0');
% rr{5+ find(strcmp(d{3}, 'e'))}=0; disp('Set e to 0');
% rr{5+ find(strcmp(d{3}, 'w'))}=0; disp('Set w to 0');

for o1=1:1 % Set up
    
    path(pathdef)
    addpath(['1 Value functions'  Valfxn_type]);   
    addpath(['1 Value functions'  Valfxn_type filesep 'b']);   
    addpath(['1 Value functions'  Valfxn_type filesep 'bp']);
    addpath(['1 Value functions'  Valfxn_type filesep 'bm']);
    addpath(['1 Value functions'  Valfxn_type filesep 'bi']);   
    addpath(['1 Value functions'  Valfxn_type filesep 'bpi']); 
    addpath(['1 Value functions'  Valfxn_type filesep 'bmi']); 
    addpath(['1 Value functions'  Valfxn_type filesep 'bpm']);
    addpath(['1 Value functions' Valfxn_type filesep 'bpmi']);  
    
    % Set up fake trialstats
    d6x6(:, [dcol.EnvThreat dcol.NTokens dcol.Choice])=[sortrows(repmat((1:6)', 6,1))./6 repmat((1:6)', 6,1).*2  randi(3, 36,1)];
    d6x6(:, [dcol.Task])=rrr.details.tasktype;
    switch rrr.details.tasktype
        case  1; d6x6=fpar_conflict(d6x6, dcol);
        case 2; d6x6=fpar_control(d6x6, dcol);
    end
end

%% Parameter correlations/tradeoffs

% Across-subjects correlation matrix
% pars={'m';'i';'w'};
    
d_corrmat=nan(length(pars), length(pars));
f.figheight= 100; f.figheight=800; f.subplotcols=3; f.subplot_VerHorz=[0.1 0.07]; f.fig_BotTop=[0.01 0.05]; f.fig_LeftRight=[0.05 0.1];
figure('Name', 'Parameters correlated across subjects?', 'Position',[200, 200,f.figheight,f.figheight]); set(gcf,'Color',[1 1 1]);  k=1;
for p1=1:length(pars)
    for p2=1:length(pars)
        
        
        subtightplot(length(pars), length(pars),  k ,  f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        if p1~=p2
            scatter(r(:, 3+ find(strcmp(d{3}, pars{p1}))), r(:, 3+ find(strcmp(d{3}, pars{p2}))))
            xlabel(pars{p1}); ylabel(pars{p2})
            [rr p]=corr(r(:, 3+ find(strcmp(d{3}, pars{p1}))), r(:, 3+ find(strcmp(d{3}, pars{p2}))));
            if p>0.1;  % title([pars{p1}  ' &  ' pars{p2}  ': correlation nsf']);
            else
                title([pars{p1}  ' &  ' pars{p2}  ': r='  num2str(rr) ',   p= ' num2str(p, 2)]);
                d_corrmat(p1,p2)=rr;
            end
        else
            title(pars{p1}, 'FontSize',20); axis off
        end
        
        k=k+1;
    end
end

figure, imagescnan(d_corrmat ,'NanColor', [0 0 0]); colorbar
set(gca,'YTick',1:length(pars), 'YTickLabel', pars,'XTick',1:length(pars), 'XTickLabel',pars)


%% How does parameter act across plausible values?
% compar(Par, Valstep, Choice, EnvThreat, NTok)

request.alterpar_steps=10;
request.alterpar={
% % %     'b'     linspace(1.3, 20, request.alterpar_steps);      % cF parameter limits (min-max)
% % %     'p'     linspace(0.001, 0.039, request.alterpar_steps);
%     'm'    logspace(log10(0.49),  log10(17), request.alterpar_steps);
%     'a'     linspace(-1.54, 13.5, request.alterpar_steps);
% % %     'f'      linspace(-16, -6, request.alterpar_steps);
%     'e'      linspace(-10, 10, request.alterpar_steps);
% %     'w'     linspace(-0.53, 0.64, request.alterpar_steps);
    %
%     'm'   logspace(log10(0.01),  log10(1.5), request.alterpar_steps);  % Mean +/- 2*SD
    'i'     linspace(-10, 15, request.alterpar_steps);
%     'f'      linspace(-16, -6, request.alterpar_steps);
%     'e'      linspace(-14, 7, request.alterpar_steps);
%     'w'     linspace(-2,2, request.alterpar_steps);
% %     
    
% % %     'b'     linspace(0.56, 3.36, request.alterpar_steps);      % ct parameter limits (min-max)
% % %     'p'     linspace(0.001, 0.02, request.alterpar_steps);
%     'm'    logspace(log10(0.02),  log10(1.8), request.alterpar_steps);
%     'a'     linspace(-20, 20, request.alterpar_steps);
% %     'e'      linspace(-26, 3.53, request.alterpar_steps);
% %     'w'     linspace(-0.92, 17.33, request.alterpar_steps);
%     %
%     'a'     linspace(-20, 20, request.alterpar_steps);
%     'e'      linspace(-13, 2, request.alterpar_steps);
%     'w'     linspace(1.2, 11.43, request.alterpar_steps);
};

mean(r(:,4:end))-1.*std(r(:,4:end))
mean(r(:,4:end))+1.*std(r(:,4:end))



% Compile DV: To change the quantity/DV that you're looking at, alter contents of variable 'compar'
%       Other data variables:  t_ = trialstats format, m_ = matrix, otherwise d_
compar=nan(size(request.alterpar,1), request.alterpar_steps, 3,6,6);  % compar(Par, Valstep, Choice, EnvThreat, NTok)
scompar=nan(20, size(request.alterpar,1), request.alterpar_steps, 3,6,6);  % scompar(Subject, Par, Valstep, Choice, EnvThreat, NTok)
for iPar=1:size(request.alterpar,1)
    for iVal=1: request.alterpar_steps
        
        sfitpars=r(:,4:end);
        for s=1:size(sfitpars,1);
%             ws.d=sdata.subjdata{s, rrr.details.tasktype+1} ;
%             ws.pars=sfitpars(s,:);
%             ws.pars(find(strcmp(d{3}, request.alterpar{iPar, 1})))=request.alterpar{iPar, 2}(iVal);
%             
%             % Cell-specific nLLs
%             for e=1:6
%                 for n=1:6
%                     wc.d=ws.d(ws.d(:, dcol.EnvThreat)==e/6 & ws.d(:, dcol.NTokens)==n*2,:);
%             
%                     % Fit parameters ####################################
%                     if isempty(strfind(mod, 'p'))==0
%                         [wc.nll wc.pch]=f_nllsoftmax_lapse(f_transpar(d{3}, ws.pars, 'from'), {mod wc.d rrr.details.fixedpar dcol});
%                     else [wc.nll wc.pch]=f_nllsoftmax(f_transpar(d{3}, ws.pars, 'from'), {mod wc.d rrr.details.fixedpar dcol});
%                     end
%                     
%                     % DV = nLL (cell - by - cell measured)
%                      scompar(s, iPar, iVal, 1:3, 7-e, n) = wc.nll;
%                      if isnan(scompar(s, iPar, iVal, 1, 7-e,   n))==1; error('NAN'); end
%                      
%                      wc=[];
%                 end
%             end
%             wd=[];
        end
%         compar(iPar, iVal, 1, :, :)=squeeze(mean(scompar(:, iPar, iVal, 1, :, :)));         % Re-assemble group from subjects
%         compar(iPar, iVal, 2, :, :)=compar(iPar, iVal, 1, :, :);
%         compar(iPar, iVal, 3, :, :)=compar(iPar, iVal, 1, :, :);
        
        % Group level ONLY ##############################################
        % Fitted params
        fitpars=[rr{6:end}];
        fitpars(find(strcmp(d{3}, request.alterpar{iPar,1}))) =    request.alterpar{iPar,2}(iVal);
        
        % Set up data that you want to collect  (t_ = trialstats format, m_ = matrix, otherwise d_)
        t_v=nan(36,1);
        m_v=repmat({nan(6, 6)}, 1,3);
        
        % [nll pch]=f_nllsoftmax(f_transpar(d{3}, fitpars, 'from'), {mod  d6x6 rrr.details.fixedpar dcol});
        eval(['t_v= squeeze('  mod '(f_transpar(d{3}, fitpars, ''from''), {[] d6x6 rrr.details.fixedpar, dcol}));'])
        for e=1:6
            for n=1:6
                m_v{1}(7-e, n)= t_v(d6x6(:, dcol.EnvThreat)==e/6 & d6x6(:, dcol.NTokens)==n*2, 1);
                m_v{2}(7-e, n)= t_v(d6x6(:, dcol.EnvThreat)==e/6 & d6x6(:, dcol.NTokens)==n*2, 2);
                m_v{3}(7-e, n)= t_v(d6x6(:, dcol.EnvThreat)==e/6 & d6x6(:, dcol.NTokens)==n*2, 3);
            end
        end
        % Keep values!
        compar(iPar, iVal, 1,:,:)=m_v{1};
        compar(iPar, iVal, 2,:,:)=m_v{2};
        compar(iPar, iVal, 3,:,:)=m_v{3};
        if isreal(sum(m_v{1}(:)+m_v{1}(:)+m_v{3}(:))) ~=1; error('Imaginary values!!') ; end
        
        
    end
end
% k=compar(:, :, :, :, :); kk=k(:);    input('nLL: correct for anti-chance?');
% kk(kk>mean(cellfun(@(x)size(x,1)/(6*6),  sdata.subjdata(:, 1+rrr.details.tasktype))))=mean(cellfun(@(x)size(x,1)/(6*6),  sdata.subjdata(:, 1+rrr.details.tasktype)));
% compar=reshape(kk,size(compar)); disp('Correcting for anti-chance nLLs!')
        

% Plot DV for each choice across all param ranges (compar)
f.figheight= 100; f.figheight=800; f.subplotcols=3; f.subplot_VerHorz=[0.05 0.07]; f.fig_BotTop=[0.01 0.05]; f.fig_LeftRight=[0.05 0.1];
% for iPar=1:size(request.alterpar,1)
%     figure('Name', ['[' request.alterpar{iPar,1} '] Values across par steps'], 'NumberTitle', 'off', 'Position',[200+iPar*15,00+iPar*15,f.figheight,f.figheight]); set(gcf,'Color',[1 1 1]); k=1; 
    figure('Name', ['Values across par steps'], 'NumberTitle', 'off', 'Position',[200+iPar*15,00+iPar*15,f.figheight,f.figheight]); set(gcf,'Color',[1 1 1]); k=1; 
for iPar=1:size(request.alterpar,1); f.subplotcols=size(request.alterpar,1);
    for iVal=1: request.alterpar_steps
        wp.par=squeeze(compar(iPar, :, :,:,:));  wp.par=wp.par(:);
        wp.parval=squeeze(compar(iPar, iVal, :,:,:));   wp.parval= wp.parval(:);
        wp.parchoice{1}=squeeze(compar(iPar, :, 1,:,:));    wp.parchoice{1}=wp.parchoice{1}(:);
        wp.parchoice{2}=squeeze(compar(iPar, :, 1,:,:));    wp.parchoice{2}=wp.parchoice{2}(:);
        wp.parchoice{3}=squeeze(compar(iPar, :, 1,:,:));    wp.parchoice{3}=wp.parchoice{3}(:);
        
        for c=3:3
            
            k=(iVal-1)*f.subplotcols+iPar;
            
            subtightplot(request.alterpar_steps, f.subplotcols,  k ,  f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
            imagesc(squeeze(compar(iPar, iVal, c,:,:))); axis square; colorbar; k=k+1; axis off
            title(['(' num2str(iVal) ') ' request.alterpar{iPar,1} ' = '  num2str(request.alterpar{iPar,2}(iVal), 2) ])
            
            
            % Standardize ranges 
            caxis([min(wp.par) max(wp.par)])   % within par (everything)
%             caxis([min(wp.parval) max(wp.parval)])   % within par (everything)
%             caxis([min(wp.parchoice{1}) max(wp.parchoice{1})])   % within par (everything)
%             input(' ');    
        end
    end
end

% Difference/subtraction plots (compdiff) 
request.parsteps_diff=[3 1]; g=1;
figure('Name', ['[' request.alterpar{iPar,1} '] Values across par steps'], 'NumberTitle', 'off', 'Position',[200+iPar*15,00+iPar*15,f.figheight,f.figheight]); set(gcf,'Color',[1 1 1]); k=1;
for iPar=1:size(request.alterpar,1); f.subplotcols=size(request.alterpar,1);
    
    for iDiff=1:request.alterpar_steps-request.parsteps_diff(1) % Compile differences
        g=iDiff+1; % if change here, change in plot too!
        compdiff(iPar, iDiff, 1,:,:)=squeeze(compar(iPar, request.parsteps_diff(1)+iDiff, 1, :,:) -  compar(iPar, request.parsteps_diff(2)+g, 1, :,:));
        compdiff(iPar, iDiff, 2,:,:)=squeeze(compar(iPar, request.parsteps_diff(1)+iDiff, 2, :,:) -  compar(iPar, request.parsteps_diff(2)+g, 2, :,:));
        compdiff(iPar, iDiff, 3,:,:)=squeeze(compar(iPar, request.parsteps_diff(1)+iDiff, 3, :,:) -  compar(iPar, request.parsteps_diff(2)+g, 3, :,:));
        
    end
    wp.pardiff{3}=squeeze(compdiff(iPar,:, 3,:,:)); % For ranges
    
    for iDiff=1:request.alterpar_steps-request.parsteps_diff(1)
        g=iDiff+1; % if change here, change in compile too!
        k=(iDiff-1)*f.subplotcols+iPar;
        
        
        subtightplot(request.alterpar_steps-request.parsteps_diff(1), f.subplotcols,  k ,  f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        imagesc(squeeze(compdiff(iPar, iDiff, 3,:,:))); axis square, axis off; colorbar       
        title([request.alterpar{iPar,1} ': step ' num2str(request.parsteps_diff(1)+iDiff) ' - ' num2str(request.parsteps_diff(2)+g)])
        
     
        caxis([min(wp.pardiff{3}(:)) max(wp.pardiff{3}(:))])
        

        
    end
    
end

 
% L=f(par) plots - how does DV change with par, in each ExN cell?
for iPar=1:size(request.alterpar,1); f.subplotcols=6; k=1; 
    figure('Name', ['[' request.alterpar{iPar,1} '] V=f(par)'], 'NumberTitle', 'off', 'Position',[400+20*iPar, 100+20*iPar,f.figheight,f.figheight]); set(gcf,'Color',[1 1 1]); k=1;
    for e=1:6
        for n=1:6
            subplot(6,6, k);
            plot(1:request.alterpar_steps, compar(iPar, :, 3, e, n)); axis square; 
           ylim([min(min(squeeze(min(compar(iPar, :, 3, :, :)))))  max(max(squeeze(max(compar(iPar, :, 3, :, :)))))])
            k=k+1;
        end
    end
    
end


%% BIC-wise, we do need every parameter?
%   This bit of code is very useful. For each parameter, compare the BICs
%   to the corresponding model without that particular parameter. Is a good
%   way of looking at whether each parameter is meaningfully contributing
%   or not, ceterus paribus

a=r_res(:,[1 3]);  % Get from results model names and BICs
a=sortrows(a,1);
a(:,3)=num2cell(sortrows(repmat((1:8)',16,1)));  % mark families (16 in each family, 8 families)
p1=~cellfun(@isempty, strfind(a(:,1), 'p'));    % Parameter indices: 
m1=~cellfun(@isempty, strfind(a(:,1), 'm'));
i1=~cellfun(@isempty, strfind(a(:,1), 'i'));
f1=~cellfun(@isempty, strfind(a(:,1), 'f'));
e1=~cellfun(@isempty, strfind(a(:,1), 'e'));
uw1=~cellfun(@isempty, strfind(a(:,1), 'uw'));
vw1=~cellfun(@isempty, strfind(a(:,1), 'vw'));
ow1=~cellfun(@isempty, strfind(a(:,1), 'ow'));

bic=[a{:,2}];
inFamily=find(e1); outFamily=find(~e1);
family1=find(uw1); family2 = find(vw1); family3 = find(ow1); family4 = find(~(vw1|uw1|ow1));
% figure, bar(inFamily,bic(inFamily), 'r', 'BarWidth', 0.4), hold all, bar(outFamily, bic(outFamily), 'b', 'BarWidth', 0.4)
% figure, bar(1:2:128,bic(inFamily), 'r', 'BarWidth', 0.4), hold all, bar(2:2:128, bic(outFamily), 'b', 'BarWidth', 0.4)


figure, bar(1:4:128, bic(family1), 'r', 'BarWidth', 0.4), hold all,     % for w param: map all 4
bar(2:4:128, bic(family2), 'g', 'BarWidth', 0.4)
bar(3:4:128, bic(family3), 'b', 'BarWidth', 0.4)
bar(4:4:128, bic(family4), 'y', 'BarWidth', 0.4)
legend({'uw1', 'vw1','ow1', 'none'})


dec2bin(0:127)  % In binary, finding every 2ndth
mod(round((0:127)/2),2)


%% Plot all BICs (for manuscript figures)

error('Manually load BICs ############# ')
a=sortrows(a,-2);  % col 1= model, col 2=BIC


% Figure settings
fontsize=25;
fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms
yslack=1000*0.5;

% Plot
bar(1:size(a,1), cell2mat(a(:,2)))
set(gcf,'Color','w')
set(gca,'FontSize',fontsize, 'FontName', fontname);
ylabel('Model BIC'); xlabel('Model')
% title('BICs across model space in Experimental task'); yrange=[-12897 * log(1/3)     max(cell2mat(a(:,2)))];
title('BICs across model space in Control task'); yrange=[-12884 * log(1/3)     max(cell2mat(a(:,2)))];
yrange=[min(cell2mat(a(:,2)))-yslack   max(cell2mat(a(:,2)))+yslack];
axis([0 size(a,1)+2  yrange(1) yrange(2)])


% Bayes factor?
bic_best=cell2mat(a(end,2));
bic_2ndbest=cell2mat(a(end-1,2));



B=(bic_best-bic_2ndbest)*-0.5; B


% 
% 12897 * log(1/3)  % Chance for cF
% 12884 * log(1/3)  % Chance for ct


%%

% load('D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical\res_hierarfitmodels_ct (21-Jul-2014) top10.mat')
% 
% 
% for s=1: size(r_res{1,2},1)
% details.subj_ntrials
% 
% 
% exp(  - r_res{1,2}(s,2)/details.subj_ntrials(s)  )
% 
% 
% 
% exp(   sum(r_res{1,2}(:,2))/12897   )
