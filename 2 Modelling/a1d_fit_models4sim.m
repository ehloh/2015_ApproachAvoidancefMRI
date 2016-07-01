% Fit models to the choice data, calculate BICs & parameter values
%   To run: all subjects, try/catch
clear all; close all hidden; clc

% Requested analysis
details.tasktype=2; % 1=cF, 2=ct
%
details.n_iterations=10;    

for o1=1:1 % Model families
    details.modelfams.b={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow';'b17_yw';'b18_fyw';'b19_eyw';'b20_feyw'; 'b21_kw';'b22_fkw';'b23_ekw';'b24_fekw'};
    details.modelfams.bp={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow';'bp17_yw';'bp18_fyw';'bp19_eyw';'bp20_feyw'; 'bp21_kw';'bp22_fkw';'bp23_ekw';'bp24_fekw'};
    details.modelfams.bm={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow';'bm17_yw';'bm18_fyw';'bm19_eyw';'bm20_feyw'; 'bm21_kw';'bm22_fkw';'bm23_ekw';'bm24_fekw'};
    details.modelfams.bpm ={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow';'bpm17_yw';'bpm18_fyw';'bpm19_eyw';'bpm20_feyw'; 'bpm21_kw';'bpm22_fkw';'bpm23_ekw';'bpm24_fekw'};
    %
    details.modelfams.bi={'bi01';'bi02_f';'bi03_e';'bi04_uw';'bi05_vw';'bi06_ow';'bi07_fe';'bi08_fuw';'bi09_fvw';'bi10_fow';'bi11_euw';'bi12_evw';'bi13_eow';'bi14_feuw';'bi15_fevw';'bi16_feow';'bi17_yw';'bi18_fyw';'bi19_eyw';'bi20_feyw'; 'bi21_kw';'bi22_fkw';'bi23_ekw';'bi24_fekw'};
    details.modelfams.bpi ={'bpi01';'bpi02_f';'bpi03_e';'bpi04_uw';'bpi05_vw';'bpi06_ow';'bpi07_fe';'bpi08_fuw';'bpi09_fvw';'bpi10_fow';'bpi11_euw';'bpi12_evw';'bpi13_eow';'bpi14_feuw';'bpi15_fevw';'bpi16_feow';'bpi17_yw';'bpi18_fyw';'bpi19_eyw';'bpi20_feyw'; 'bpi21_kw';'bpi22_fkw';'bpi23_ekw';'bpi24_fekw'};
    details.modelfams.bmi ={'bmi01';'bmi02_f';'bmi03_e';'bmi04_uw';'bmi05_vw';'bmi06_ow';'bmi07_fe';'bmi08_fuw';'bmi09_fvw';'bmi10_fow';'bmi11_euw';'bmi12_evw';'bmi13_eow';'bmi14_feuw';'bmi15_fevw';'bmi16_feow';'bmi17_yw';'bmi18_fyw';'bmi19_eyw';'bmi20_feyw'; 'bmi21_kw';'bmi22_fkw';'bmi23_ekw';'bmi24_fekw'};
    details.modelfams.bpmi ={'bpmi01';'bpmi02_f';'bpmi03_e';'bpmi04_uw';'bpmi05_vw';'bpmi06_ow';'bpmi07_fe';'bpmi08_fuw';'bpmi09_fvw';'bpmi10_fow';'bpmi11_euw';'bpmi12_evw';'bpmi13_eow';'bpmi14_feuw';'bpmi15_fevw';'bpmi16_feow';'bpmi17_yw';'bpmi18_fyw';'bpmi19_eyw';'bpmi20_feyw'; 'bpmi21_kw';'bpmi22_fkw';'bpmi23_ekw';'bpmi24_fekw'};
    details.ct_modnums=[1 3 4 5 6 11 12 13 17 19]; % true models for ct
    %
    details.modelfams.bj={'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow';'bj17_yw';'bj18_fyw';'bj19_eyw';'bj20_feyw'; 'bj21_kw';'bj22_fkw';'bj23_ekw';'bj24_fekw'};
    details.modelfams.bpj={'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow';'bpj17_yw';'bpj18_fyw';'bpj19_eyw';'bpj20_feyw'; 'bpj21_kw';'bpj22_fkw';'bpj23_ekw';'bpj24_fekw'};
    details.modelfams.bjm={'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow';'bjm17_yw';'bjm18_fyw';'bjm19_eyw';'bjm20_feyw'; 'bjm21_kw';'bjm22_fkw';'bjm23_ekw';'bjm24_fekw'};
    details.modelfams.bpjm={'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm09_fvw';'bpjm10_fow';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow';'bpjm17_yw';'bpjm18_fyw';'bpjm19_eyw';'bpjm20_feyw'; 'bpjm21_kw';'bpjm22_fkw';'bpjm23_ekw';'bpjm24_fekw'};
end
details.whichmodels=details.modelfams.bpjm(19);

for o1=1:1 %% Set up
    
details.folder4saving=['5 Simulation' filesep '3 Paramfits refit'];
    details.dataset_date=details.whichmodels{1};

%     details.dataset_date='(09-May-2014)'; % Actual
    Valfxn_type=' Det Jpower';   % details.folder4saving=Valfxn_type(2:end);
    details.Vanfxn_type=Valfxn_type;
    w=pwd;  if strcmp(w(2), ':')==1; where.where='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models'; else where.where='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models'; end
    path(pathdef); addpath([where.where filesep '1 Value functions' Valfxn_type]);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'b']);   
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bp']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bm']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpm'])
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bi']);   
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpi']); 
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bmi']); 
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpmi']);  
    %
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bj']);   
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpj']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bjm']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpjm']);
    if details.tasktype==1; details.task='cF'; elseif details.tasktype==2; details.task='ct'; else error('Invalid details.tasktype requested (1=cF, 2=ct)'); end    
%     w=load([where.where filesep '2 Analysis inputs' filesep 'All Data ' details.dataset_date '.mat']);
    

    % Load simulated data
    w=load([where.where filesep '2 Analysis inputs' filesep '5 Simulation' filesep '2 Simulated data' filesep 'Simdata Fit ' details.task ' - ' details.dataset_date '.mat']);

    col=w.details.col; subjdata=w.subjdata; details.n_subjs=w.details.n_subjs;  details.subjects=w.details.subjects;  details.col=col;
    details.dataset=['All Data ' details.dataset_date '.mat']; % Selection of specific subjects not implemented
    
    % Task parameters
    details.fixedpar.cF_FL=  -12;
    details.fixedpar.cF_EC=  -2;
    details.fixedpar.ct_FL=0;
    details.fixedpar.ct_EC=  -2;
    
    diary([where.where filesep '2 Analysis inputs' filesep '1 Fit logs' filesep 'diary_fit_' details.task ' ' Valfxn_type ' (' date ')'])
    
    % Model settings + fetch requested model details
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels,details.n_iterations);
    details.n_models=length(details.whichmodels); errorlog={}; e=1;
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
    
    disp('=============================================================='); w.c=clock; w.c1=w.c;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ');
    disp(['Requested: ' num2str(details.n_subjs) ' subjects, ' num2str(details.n_iterations) ' iterations per model'])
    disp(['details.task type:  ' details.task ]); disp(' ');
    disp([num2str(length(details.whichmodels))  ' Models:']); disp(details.whichmodels); disp(' ');
    disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
    input('Hit enter to start                                   ');
    disp('==============================================================')
end

%% Run all models

startm=1;

% 'r_res'
%       Col 1: Model name
%       Col 2: Subject fit parameters
%             Col i: BIC
%             Col ii: nLL
%             Col iii: fminunc exceeded default iterations?
%             Col iv onwards: parameters (beta first)
%       Col 3: Model BIC (summed across subjects)
%       Col 4: Hessians
%       Col 5: 
%       Col 6 onwards: mean parameter values
for o1=1:1 % Columns 
    rc.modname=1;
    rc.subpars=2;
    rc.sp.bic=1;
    rc.sp.nll=2;
    rc.sp.p1=4;
    rc.modelbic=3;
    rc.hessians=4;
    rc.mean_p1=6;
    %
    if exist('r_res', 'var')==0; if startm~=1; input('Requested start is NOT from model #1. Proceed?    '); end; r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2)))); r_iterations=cell(details.n_models,4);  end
end
for m=startm:details.n_models
    w.c=clock; disp(['Model ' num2str(m) ' - ' details.models{m,1} '       ['  num2str(w.c(4))  ':'  num2str(w.c(5))  '] ##############' ])
    r_res{m,1}=details.models{m,1};
    r_iterations{m,1}=details.models{m,1};
    r_iterations{m,2}=zeros(details.n_subjs*details.n_iterations, details.models{m,2}+2);
    r_iterations{m,3}=cell(details.n_subjs*details.n_iterations, 2);
    if details.n_iterations>1
        for p=1:details.models{m, 2} % Randomize order of par steps
            details.models{m, 6}(p) =cellfun(@(x)x(:,1),  {sortrows([details.models{m, 6}{p}(:) rand(details.n_iterations,1)],2)}, 'UniformOutput', 0);
        end
    end
    % 'r_iterations' - model fits for each subject x iteration
    %       Col 1=Model name
    %       Col 2=Iteration results (nLL, parameters)
    %             Col 1: Subject
    %             Col 2: nLL
	%             Col 3: Iterations exceeded fminunc?
    %             Col 4 onwards: model parameters (1st parameter is beta/inverse temperature)
    %       Col 3=Iteration hessians
    %             Col 2: Hessian for this iteration
    %       Col 4=nLL histogram data
    %       Col 5= Subject hessians for best fit (col 1), ok? (non-singular; col 2)
    for s= 1:  size(subjdata,1)
        disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
        ws.data=subjdata{s, details.tasktype+1};
%         ws.data=sortrows(ws.data, [col.EnvThreat col.NTokens]); disp('Sorting by trial type!!');
%         tn=1:18:size(ws.data,1);  ws.data=ws.data(tn, :);  ws.data(:, col.Trialnum)=1:36;; disp('Data simplifications: 1 per trial type');
        
        for i=1:details.n_iterations
%             try 
                %  Starting parameters (read seed + apply inverse constraint if applicable)
                wi.startpar=nan(1,details.models{m,2});
                for p=1:details.models{m,2}
                    x=details.models{m,6}{p}(i); eval(['wi.startpar(p)=' details.models{m,5}{p} ';']);
                    if strcmp(details.models{m,3}{p}, 'm')~=1;  wi.startpar(p)=round(wi.startpar(p)*10)/10; end
                end
                
                % Fit parameters ####################################
                if strcmp(details.models{m,1}(2), 'p')==1         
                    [wi.par, wi.L, wi.exit, wi.output, wi.grad, wi.hessians]=fminunc(@(x)f_nllsoftmax_lapse(x, {details.models{m,1} ws.data details.fixedpar col},0,0), wi.startpar, w.options);
                else
                    [wi.par, wi.L, wi.exit, wi.output, wi.grad, wi.hessians]=fminunc(@(x)f_nllsoftmax(x, {details.models{m,1} ws.data details.fixedpar col},0,0), wi.startpar, w.options);
%                     [nll pch]=f_nllsoftmax(wi.startpar, {details.models{m,1} ws.data details.fixedpar col});
                end
                
%                 error
                
                % Convert parameters to parameter space
                for p=1:details.models{m,2}
                    x=wi.par(p); eval(['wi.par(p)=' details.models{m,4}{p} ';']);
                end
                
                % Write to array
                wi.rownum=(s-1)*details.n_iterations+i;
                r_iterations{m,2}(wi.rownum,1)=s;
                r_iterations{m,2}(wi.rownum,2)=wi.L;
                r_iterations{m,2}(wi.rownum,3)= wi.exit==0;
                r_iterations{m,2}(wi.rownum,4:3+length(wi.par))=wi.par;
                r_iterations{m,3}{wi.rownum,1}=s;
                r_iterations{m,3}{wi.rownum,2}=wi.hessians;
                wi=[];
                
%                 catch MExc
%                     errorlog{e,1}=['Failed: ' details.subjects{s} '  -  ' details.models{m,1} '    iteration ' num2str(i)];
%                     errorlog{e,2}=wi.startpar;
%                     errorlog{e,3}=MExc;
%                     disp(errorlog{e,1}); e=e+1;
%                     
%                     % Fake inputs
%                     wi.rownum=(s-1)*details.n_iterations+i;
%                     r_iterations{m,2}(wi.rownum,1)=s;
%                     r_iterations{m,2}(wi.rownum,2)=nan;
%                     r_iterations{m,2}(wi.rownum,4:3+length(wi.startpar))=wi.startpar;
%                     r_iterations{m,3}{wi.rownum,1}=s;
%                     r_iterations{m,3}{wi.rownum,2}=[];
%             end
        end
        
        % Choose best iteration for this subject + record details
        ws.iters=r_iterations{m,2}(  (s-1)*details.n_iterations+1:s*details.n_iterations, :);
        ws.iterhess=r_iterations{m,3}((s-1)*details.n_iterations+1:s*details.n_iterations, :);
        ws.bestiter= ws.iters(find(ws.iters(:,2)== min(ws.iters(:,2)), 1, 'first'),:);
        r_res{m,2}(s,2)=ws.bestiter(2);
        r_res{m,2}(s,4:3+length(ws.bestiter)-3)=ws.bestiter(4:end);
        if ws.bestiter(3)==1;  r_res{m,2}(s, 3)=1; end
        r_res{m,4}{s,1}=ws.iterhess{find(ws.iters(:,2)== min(ws.iters(:,2)), 1, 'first'),2};
        r_iterations{m,5}{s,1}=ws.iterhess{find(ws.iters(:,2)== min(ws.iters(:,2)), 1, 'first'),2};
        r_iterations{m,5}{s,2}= sum(isnan(r_iterations{m,5}{s,1}(:)))>0 + rank(r_iterations{m,5}{s,1})<details.models{m,2};
        
        % Calculate BIC: 2*nll+ K*ln(n_trials)
        r_res{m,2}(s,1)= 2*r_res{m,2}(s,2) +  details.models{m,2} * log(size(ws.data,1));
        
        %
        ws=[];
    end
    
    % Calculate overall BIC & mean parameters for all subjects 
    r_res{m,3}=sum(r_res{m,2}(:,1));
    r_res(m, 6:5+details.models{m,2})=num2cell(mean(r_res{m,2}(:,4:end),1));
    r_iterations{m,4}= mean(reshape(r_iterations{m,2}(:,2), details.n_iterations, details.n_subjs),2);  % nLL histograms
    
    % Save in partial fits in-between models
    wi.c=clock; save(['2 Analysis inputs' filesep 'PartialFits' filesep 'Partial fit ' details.task ' (' date ' - ' num2str(wi.c(4)) ' hrs)']);
end
r_res=sortrows(r_res,3);

%% END 

disp('===================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp('INSTRUCTION: Results variables (prefix r):r_iterations, r_res, r_models')
disp('                     [r_res: Col 1= Model, Col 2=Subject fits, 3=Model BIC, Col 4=N valid subjects, Col 6 onwards=Mean parameter values]'); disp(' ')
disp('Check command-window log for poor fits; See ''further instructions'' at end of script for Bayes factor & model weights'); disp(' ')
disp('====================================')

% Save
resfilewhere=[where.where filesep '2 Analysis inputs' filesep details.folder4saving  filesep];
resfilename=['Fitsim ' details.task ' - ' details.models{1}];
save([resfilewhere resfilename], 'details', 'r_iterations','r_res', 'rc','errorlog'); diary off



%% Further instructions (Copy to command window to calculate)

openvar r_res; diary off;  details.whichmodels
error('Done :)')


try  f_sendemail('kurzlich', ['Modelling fitting (using fminunc) done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end; char(errorlog{:,1})

% nLL histograms: can we trust the minimizations? Iterations should reliably find the lowest nLL fit
f.subplot_rowcols=[4 4];
% f.subplot_rowcols=[4 2];
f.figwidth= 500; f.figheight=1000; ff=1; mm=1;
f.subplot_VerHorz=[0.1 0.03]; f.fig_BotTop=[0.05 0.03]; f.fig_LeftRight=[0.05 0.1];
figure('Name', ['nLL Histograms Plot ' num2str(ff)], 'NumberTitle', 'off', 'Position',[200,00,f.figheight,f.figwidth], 'Color',[1 1 1]);
for m=1:details.n_models
    subtightplot(f.subplot_rowcols(1), f.subplot_rowcols(2),   mm,  f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
    hist(round(r_iterations{m,4})); title(r_iterations{m,1}, 'FontSize', 14)
%     axis off
    if m~=details.n_models && m/(f.subplot_rowcols(1)*f.subplot_rowcols(2)) == round(m/(f.subplot_rowcols(1)*f.subplot_rowcols(2)))
        ff=ff+1; mm=1; 
        figure('Name', ['nLL Histograms Plot ' num2str(ff)], 'NumberTitle', 'off', 'Position',[200,00,f.figheight,f.figheight], 'Color',[1 1 1]);
    else mm=mm+1; 
    end
end

