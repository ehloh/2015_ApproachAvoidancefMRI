function [ matlabbatch contrastfiles] = RL_linearFX(logg,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch contrastfiles] = RL_linearFX(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% (1) Trial-type contrasts (Task x EnvThreat x NTokens) are fed into 2x6x6
%       factorial
% (2) Specify new contrasts using RL weights (linear effect)
%       (For all RL variables)
%
% ------------------------------------------------------------------------------------------------------------

% Execute following to use as script:   where.data_brain; subjectlist=log.subjects; firstlevelmodel=log.onsetsmodel; secondlevelmodel=log.secondlevelmodel; logg=log; clear 'log';

%% (1) Details for this model

% Run contrasts only?
disp('------------------------------------------------------------------------------')
disp('Which 2nd level steps to run:')
disp('(1)   Specify & Estimate model, then add Contrasts')
disp('(2)   JUST run Contrasts (assuming model is spec & est already)')
contrastsonly=input('Enter option (1 or 2)   ');

% % Set up folder 
secondlevelfolder=[where.resultsfolder filesep secondlevelmodel];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); mkdir([secondlevelfolder filesep 'ROI']); end

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

% Available contrasts
w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
spm=w.SPM;  
contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% Construct design first (RL variables)
RLvariables={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'EV';}; 
d.Env=[((1:6)'*ones(1,6))/6]'; % No longer in visualization space
d.N=[2*ones(6,1)*(1:6)]';
d.pLoss=[d.Env.*(d.N/12)]';
d.Unc=nan*zeros(6,6); d.Conf=nan*zeros(6,6);
d.EV_Conflict=(1-d.pLoss).*d.N + d.pLoss*(-12); 
for e=1:6 % Uncertainty, Conflict
    for n=1:6 
        if d.pLoss(e,n)==1 % Correction for infinite values
            d.pLoss(e,n)=d.pLoss(e,n)-0.00001;            
        end
        
        % Uncertainty/Entropy=-p (logp) - (1-p) logg(1-p)
        d.Unc(e,n)= - d.pLoss(e,n).* log(d.pLoss(e,n)) - (1-d.pLoss(e,n)) *(log(1-d.pLoss(e,n)));
        
        % Conflict= Uncertainty * NTokens
        d.Conf(e,n)=d.Unc(e,n)*d.N(e,n);
        
        if d.pLoss(e,n)==1-0.000001 % Correction for infinite values
            d.pLoss(e,n)=d.pLoss(e,n)+0.000001;
        end
    end
end
d.EV_Conflict=(1-d.pLoss).*d.N + d.pLoss*(-12); 
d.EV_Control= (1-d.Unc).*d.N;
d.cF_design=zeros(6*6,7); 

% Mean-centre all variables (except d.Env & d.N)
d.cF_design(:,1)=1; 
d.cF_design(:,2)=d.Env(:);
d.cF_design(:,3)=d.N(:);
d.cF_design(:,4)=d.pLoss(:)-mean(d.pLoss(:));
d.cF_design(:,5)=d.Unc(:)-mean(d.Unc(:));
d.cF_design(:,6)=d.Conf(:)-mean(d.Conf(:));
d.cF_design(:,7)=d.EV_Conflict(:)-mean(d.EV_Conflict(:));
d.ct_design=d.cF_design; % Control task: Follows Conflict parameter specification, except for EV
d.ct_design(:,1)=2;
d.ct_design(:,7)=d.EV_Control(:)-mean(d.EV_Control(:));

% Factor assignments + RL value for each contrast (Task x EnvThreat x NTokens contrasts only):
%   'design': Col 1=Task, Col 2=EnvThreat, Col 3=NTokens, Col 4=pLoss, Col 5=Entropy, Col 6= Conflict, Col 7= EV
%   'design_cell': Col 1=Trialtype, Col 2= con image, Col 3=Task ... (same as above)
contrasts=contrasts([7:78, 85:2*6*6],:);
design_cell=[contrasts num2cell(vertcat(d.cF_design, d.ct_design))];
col.b4design=3;

%% (3) Specify model for Factorial analysis

% General settings
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% 2x6*6 ANOVA: Task x EnvThreat x NTokens
disp('NOTE: ANOVA SETTINGS ARE DEFAULT')
matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0; % non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'EnvThreat';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 6;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'NTokens';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 6;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).ancova = 0;

% Include subjects as covariates
for s=1:length(subjectlist)
    sub=zeros(1,length(subjectlist)); sub(s)=1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 2*6*6])';
   matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
end
% matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});

% Assign for each cell: Scan (con image) & Factor levels
for c=1:size(design_cell,1)
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(c).levels =[design_cell{c,3} design_cell{c,4}*6  design_cell{c,5}/2]; % EnvThreat + N Tokens, categorical coding
    for s=1:logg.n_subjs
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(c).scans{s,1}= [where.data_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep design_cell{c,2} ',1'];
    end
end

%% RUN

% [EXECUTE] Specify  ##############################
disp('Specifying model (fake factorial) --------------------------')
if contrastsonly==1
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
end
matlabbatch=[];

%  [EXECUTE] Estimate model ##############################
disp('Estimating model (fake factorial) ------------------------------')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
if contrastsonly==1
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
end
matlabbatch=[];

%% (4) Re-specify sensible contrasts

% Set up requested contrasts in Matlab
disp('Running contrasts to specify linear effects ####################')
matlabbatch{1}.spm.stats.con.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 1; % Delete existing contrasts?
if matlabbatch{1}.spm.stats.con.delete==1; disp('(Deleting fake-factorial''s contrasts)'); end

% (1) Contrasts for Task --------------
j=1;
matlabbatch{1}.spm.stats.con.consess{j}.tcon.name ='Conflict Task - Pos';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = [ones(1,36) zeros(1,36)]';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
j=j+1;
matlabbatch{1}.spm.stats.con.consess{j}.tcon.name ='Conflict Task - Neg';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = -1*[ones(1,36) zeros(1,36)]';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
j=j+1;
matlabbatch{1}.spm.stats.con.consess{j}.tcon.name ='Control Task - Pos';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec =  [zeros(1,36) ones(1,36)]';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
j=j+1;
matlabbatch{1}.spm.stats.con.consess{j}.tcon.name ='Control Task - Neg';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = -1*[zeros(1,36) ones(1,36)]';
matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
j=j+1;

% (2) Contrasts for RL variables (ignoring task) ------------
for i=1:length(RLvariables)
    % Specify contrasts according to instructions 
    
    % Positive
    matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = [RLvariables{i} ' - Pos'];
    matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = cell2mat(design_cell(:, col.b4design+i)); % weights
    matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
    j=j+1;
    % Negative
    matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = [RLvariables{i} ' - Neg'];
    matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = -1* cell2mat(design_cell(:, col.b4design+i)); % weights
    matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
    j=j+1;
end


% Contrasts for Task x RL variables
for t=1:2
    if t==1; prefix='cF_'; sel=[ones(1,36) zeros(1,36)]';
    elseif t==2; prefix='ct_'; sel=[zeros(1,36) ones(1,36) ]';
    end
    
    for i=1:length(RLvariables)
        % Specify contrasts according to instructions
        
        % Positive
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = [prefix RLvariables{i} ' - Pos'];
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = cell2mat(design_cell(:, col.b4design+i)).*sel; % weights
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        j=j+1;
        
        % Negative
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = [prefix RLvariables{i} ' - Neg'];
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = (-1* cell2mat(design_cell(:, col.b4design+i))).*sel; % weights
        matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        j=j+1;
    end
end


disp('DONE ------------------------------------- ')

%%

% Execute
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];
contrastfiles=design_cell; % For the function output

% Save details of this model + contrasts
details.design_cell=design_cell;
details.design_parameters=d;
details.RLvariables=RLvariables;
save([secondlevelfolder filesep 'details_2ndlevel.mat'], 'details', 'logg'); % Save details in 2nd level folder

end

