% Fit models to the choice data, calculate BICs & parameter values
clear all; close all hidden; clc

% Requested analysis
details.tasktype=1; % 1=cF, 2=ct
details.dataset_date='(09-May-2014)';
details.par_steps=5;
% details.folder4saving='2 High iter';
details.folder4saving=[];

for o1=1:1 % Request models 
    details.modelfams{1}={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow'};
    details.modelfams{2}={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow'};
    details.modelfams{3}={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow'};
    details.modelfams{4} ={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow'};
    %
    details.modelfams{5}={'bi01';'bi02_f';'bi03_e';'bi04_uw';'bi05_vw';'bi06_ow';'bi07_fe';'bi08_fuw';'bi09_fvw';'bi10_fow';'bi11_euw';'bi12_evw';'bi13_eow';'bi14_feuw';'bi15_fevw';'bi16_feow'};
    details.modelfams{6} ={'bpi01';'bpi02_f';'bpi03_e';'bpi04_uw';'bpi05_vw';'bpi06_ow';'bpi07_fe';'bpi08_fuw';'bpi09_fvw';'bpi10_fow';'bpi11_euw';'bpi12_evw';'bpi13_eow';'bpi14_feuw';'bpi15_fevw';'bpi16_feow'};
    details.modelfams{7} ={'bmi01';'bmi02_f';'bmi03_e';'bmi04_uw';'bmi05_vw';'bmi06_ow';'bmi07_fe';'bmi08_fuw';'bmi09_fvw';'bmi10_fow';'bmi11_euw';'bmi12_evw';'bmi13_eow';'bmi14_feuw';'bmi15_fevw';'bmi16_feow'};
    details.modelfams{8} ={'bpmi01';'bpmi02_f';'bpmi03_e';'bpmi04_uw';'bpmi05_vw';'bpmi06_ow';'bpmi07_fe';'bpmi08_fuw';'bpmi09_fvw';'bpmi10_fow';'bpmi11_euw';'bpmi12_evw';'bpmi13_eow';'bpmi14_feuw';'bpmi15_fevw';'bpmi16_feow'};
    details.ct_modnums=[1 3 4 5 6 11 12 13]; % true models for ct
    % NOTE: k param models cannot be grided (21-24)
    % NOTE: j parameter models cannot be gridded
    details.modelfams{9}={'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow'};
    details.modelfams{10}={'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow'};
    details.modelfams{11}={'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow'};
    details.modelfams{12}={'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm09_fvw';'bpjm10_fow';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow'};
end
details.whichmodels=[details.modelfams{1}; details.modelfams{2}; details.modelfams{4}; details.modelfams{5}; ];
% details.whichmodels=details.modelfams{5};
% which=[1]; details.whichmodels=[ details.modelfams{1}(which); details.modelfams{2}(which); details.modelfams{4}(which);  details.modelfams{5}(which) ];
details.whichmodels=details.modelfams{3}(1);
% details.whichmodels=[details.modelfams{1}; details.modelfams{4}];
details.whichmodels={ 
% %     'bpm11_euw'
    'bpm16_feow'
    };
% details.whichmodels=details.modelfams{9}(1);
% details.whichmodels={ 
% 'bio01'; 'bpio01'; 'bpio04_uw';'bpio08_fuw'

%     'bpi04_uw_o.m';'bpi08_fuw_o.m';'bpji08_fuw_o.m'  % Additional models from CerCortex Reviews 1
%     };   % Additional models from CerCortex Reviews 1



for o1=1:1 %% Set up
     
    Valfxn_type=[];   details.folder4saving=Valfxn_type;
    details.Vanfxn_type=Valfxn_type;
      w=pwd;  if strcmp(w(2), ':')==1;  where.where='C:\Users\e.loh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis\4 Fit computational models'; else  where.where='/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis/4 Fit computational models'; end
    path(pathdef); addpath([where.where filesep '1 Value functions' Valfxn_type]);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'misc']);   
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'misc' fs 'Base']);   

    
%     
%     w.w=pwd;  if strcmp(w.w(2), ':')==1; where.where='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models'; else where.where='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models'; end
%     path(pathdef); addpath([where.where filesep '1 Value functions' Valfxn_type]);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'b']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bp']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bm']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpm']);
    %
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bj']);   
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bpj']);
    addpath([where.where filesep '1 Value functions' Valfxn_type filesep 'bjm']);
    
    w=load([where.where filesep '2 Analysis inputs' filesep 'All Data ' details.dataset_date '.mat']);
    col=w.details.col; subjdata=w.subjdata; details.n_subjs=w.details.n_subjs;  details.subjects=w.details.subjects;  details.col=col;
    details.dataset=['All Data ' details.dataset_date '.mat']; % Selection of specific subjects not implemented
    
    % Task parameters
    details.fixedpar.cF_FL=-12;
    details.fixedpar.cF_EC=-2;
    details.fixedpar.ct_FL=0;
    details.fixedpar.ct_EC=-2;
    if details.tasktype==1; details.task='cF'; elseif details.tasktype==2; details.task='ct'; else error('Invalid details.tasktype requested (1=cF, 2=ct)'); end    
    rand('state',sum(100*clock));
    
    % Model settings + fetch requested model details
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels,details.par_steps);
    details.n_models=length(details.whichmodels);
    
    disp('===================================='); w.c=clock; w.c1=w.c;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ');
    disp(['Requested: ' num2str(details.n_subjs) ' subjects'])
    disp(['Task type:  ' details.task ]); disp(' ');
    disp('Models:'); disp(details.whichmodels); disp(' ');
    disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
    input('Hit enter to start                                   ');
    disp('====================================')
 end
 
%% Grid search all parameter combinations

details.n_subjs=2; input('Just testing out low n subjects!'); 


startm=1;

% Fit grid: 'r_grid' - model grids for each subject
%       Col 1=Model name
%       Col 2=Grid search results (nLL, parameters)
%             Col 1: Subject
%             Col 2: nLL
%             Col 3: Full grid of nLL (all param-points, p-dimensions)
%             Col 4: Indices of best fit (calculated later)
%             Col 6 onwards: best fit model parameters (1st parameter is beta/inverse temperature)
if exist('r_grid', 'var')==0; r_grid=cell(details.n_models,1); end
for m=startm:details.n_models
    disp(['Model ' num2str(m) ' - ' details.models{m,1} '  ##############' ])
    r_grid{m,1}=details.models{m,1};
        
    % Apply inverse constraints to parameter ranges (correct for constraints
    %       being applied within the softmax/value functions). Col 7=Parameters
    %       with inverse constraints
    wi.startpar=cell(1,details.models{m,2});
    for p=1:details.models{m,2}
        x=details.models{m,6}{p};
        eval(['wi.startpar{p}=' details.models{m,5}{p} ';']);
    end
    details.models{m,7}=wi.startpar;
    
    % Fit all subjects
    for s=1: details.n_subjs
        disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
        r_grid{m,2}{s,1}=details.subjects{s};
        
        if strcmp(details.models{m,1}(2), 'p')==1
            eval(['[ ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax_lapse, {'''  details.models{m,1} ''' subjdata{s, details.tasktype+1} details.fixedpar col},  details.models{m,7});'])
        else
            eval(['[ ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax, {'''  details.models{m,1} ''' subjdata{s, details.tasktype+1} details.fixedpar col},  details.models{m,7});'])            
%              [ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax, {'b01' subjdata{s, details.tasktype+1} details.fixedpar col},  details.models{m,7});
        end
        r_grid{m,2}{s,3}=ws.nLL;
    end
    
end

%% For each subject, find best fit (+ calculate model vals)

% 'r_res'
%       Col 1: Model name
%       Col 2: Subject fit parameters
%             Col i: BIC
%             Col ii: nLL
%             Col iv onwards: parameters (beta first)
%       Col 3: Model BIC (summed across subjects)
%       Col 4:
%       Col 5:
%       Col 6 onwards: mean parameter values
r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2))));
for o1=1:1 % Columns 
    rc.modname=1;
    rc.subpars=2;
    rc.sp.bic=1;
    rc.sp.nll=2;
    rc.sp.p1=4;
    rc.modelbic=3;
    rc.mean_p1=6;
    %
    r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2)))); r_iterations=cell(details.n_models,1); 
end
for m= 1:details.n_models
    disp(['Model ' num2str(m) ' - ' details.models{m,1} '  ##############' ])
    r_res{m,1}=details.models{m,1};
    r_res{m,2}=nan(details.n_subjs, 3+details.models{m,2});
    
    wm.indices=[]; % Indices for best fit
    for i=1: details.models{m,2}; wm.indices= [wm.indices ' ii' num2str(i)]; end
    
    for s=1:details.n_subjs
        disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
        r_grid{m,2}{s,1}=details.subjects{s};
        ws.gridnLL=r_grid{m,2}{s,3};
        
        % Locate best match
        ws.nmatch=length(find(ws.gridnLL==min(ws.gridnLL(:))));
        eval(['[' wm.indices ']=ind2sub(size(ws.gridnLL), find(ws.gridnLL==min(ws.gridnLL(:))));'])
        if ws.nmatch~=1; disp('Multiple matches! Assume 1st (arbitrary)'); for i=1: details.models{m,2}; eval(['ii' num2str(i) '=ii' num2str(i) '(1);']); end; end
        eval(['ws.index=num2cell([' wm.indices ']);']);
        ws.ii=  [cellfun(@(x)[num2str(x) ','], ws.index(1:length(ws.index)-1), 'UniformOutput',0) num2str(ws.index{length(ws.index)}) ];
        r_grid{m,2}{s,4}=ws.index;
        
        % Record details of best fit (r_grid)
        eval(['r_grid{m,2}{s,2}=ws.gridnLL(' strcat(ws.ii{:}) ');'])
        for p=1:details.models{m,2}
            r_grid{m,2}{s,5+p}=details.models{m,6}{p}(ws.index{p});
        end
        
        % Record details of best fit (r_res{m,2})
        r_res{m,2}(s,2)=r_grid{m,2}{s,2};
        r_res{m,2}(s,4:3+details.models{m,2})=cell2mat(r_grid{m,2}(s,6:end));
        r_res{m,2}(s,1)= 2*r_res{m,2}(s,2) +  details.models{m,2} * log(size(subjdata{s,details.tasktype+1},1)); % Calculate BIC: 2*nll+ K*ln(n_trials)
    end
    
    % Calculate overall BIC & mean parameters for all subjects
    r_res{m,3}=sum(r_res{m,2}(:,1));
    r_res(m, 6:5+details.models{m,2})=num2cell(mean(r_res{m,2}(:,4:end),1));
    
end
r_res=sortrows(r_res,3);

% Save
resfilename=['res_gridmodels_' details.task ' (' date ')'];
resfilewhere=[where.where filesep '2 Analysis inputs' filesep details.folder4saving  filesep ];
resfilesthere=cellstr(spm_select('List', resfilewhere, ['res_gridmodels_' details.task '*.*' date]));
if isempty(resfilesthere)~=1 && strcmp(resfilesthere{end}(strfind(resfilesthere{end}, ')')+1), '.')==0 % number re existing fit from this date
    k=1; kk=0;
    while kk==0 
        if sum(strcmp(resfilesthere,[resfilename num2str(k) '.mat']));  k=k+1;
        else resfilename=[resfilename num2str(k)]; kk=1;
        end
    end
end
save([resfilewhere resfilename],  'r_grid','r_res', 'details')

%% END 

disp('===================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp('INSTRUCTION: Results variables (prefix r):r_grid, r_res, r_models')
disp('                     [r_res: Col 1= Model, Col 2=Subject fits, 3=Model BIC, Col 4=N valid subjects, Col 6 onwards=Mean parameter values]'); disp(' ')
disp('Check command-window log for poor fits; See ''further instructions'' at end of script for Bayes factor & model weights'); disp(' ')
disp('====================================')


%% Further instructions (Copy to command window to calculate

% Calculate Bayes factor ----------------------
m1=1;  m2=2;
B=(r_res{m1,3}-r_res{m2,3})*-0.5;

disp(['B=' num2str(B)  '  (m1=' r_res{m1,1} ', m2='  r_res{m2,1} ')'])
openvar r_res;
% Interpreting B (conventions): 3-10=moderate evidence, >10=strong evidence (in favour of m1)

try f_sendemail('kurzlich', ['Modelling fitting (using grid) done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end
%
% % ----------------------------------------------


% Fits pars?
m=1; disp(['Model:   '   r_grid{m,1} '     ' details.models{m,1}     '       '  r_res{m,1} ' (if mismatch-DELETE!)'])
parnum=2; pars=r_res{m,2}(:, 3+parnum);
hist(pars); details.models{m,6}{parnum}

%% Look directly at nLL surface
%       Chance nLL: nTrials x log(1/3)

m=1; s=4; 

if strcmp(r_grid{m,1}, details.models{m,1})+strcmp(r_grid{m,1}, r_res{m,1})~=2; error('Check order of models!'); end
disp(['Model name:   '  details.models{m,1}]);
gg=r_grid{m,2}; g=gg{s,3};
for p=4: size(r_res{m,2},2);  % Which grid point for best fit?
    disp([details.models{1,3}{p-3} ':  '  num2str(find(details.models{1,6}{p-3}==r_res{m,2}(s,p))) ])
end

% Plot
surf(g(:, :, 7,2,3,10))
ylabel('Par 1'); xlabel('Par 2');  
title([r_grid{m,1} '  (' r_grid{m,2}{s,1} ', min='   num2str(r_res{m,2}(s,2), 4) ')'])


