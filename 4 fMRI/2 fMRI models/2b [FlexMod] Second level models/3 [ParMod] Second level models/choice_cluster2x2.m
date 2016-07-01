function [ matlabbatch contrastfiles] = choice_cluster2x2(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch contrastfiles] = choicecluster_2x2(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)

% Execute following to use as script:   subjectlist=log.subjects; firstlevelmodel=log.onsetsmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% % Set up folder 
secondlevelfolder=[where.resultsfolder filesep secondlevelmodel];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); mkdir([secondlevelfolder filesep 'ROI']); end

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

% Available contrasts
w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
spm=w.SPM;  
contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% Assign Choices (Reject/Bomb & Explore only)
choices={'in_cF_Reject' [1 1];'in_cF_Explore'  [1 2];'in_ct_Bomb'  [2 1];'in_ct_Explore'  [2 2]};
for i=1:size(choices,1)
       a.which=strcmp(contrasts(:,1), choices{i});
        if sum(a.which)>1
            error(['ERROR: More than 1 contrast file found for RL variable ' choices{i}])
        elseif sum(a.which)==0
            error(['ERROR: Could not find contrast file for RL variable ' choices{i}])
        end
        choices{i,3}=contrasts{find(a.which(:,1)),2};
end

%% (2) Specify model for Factorial analysis

% 2x3 ANOVA: Task x Choice (Accept/Reject/Explore, NoBomb/Bomb/Explore)
matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Choice';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Specify Contrast files + Factorial cell 
contrastfiles=cell(size(choices,1),1);
for i=1:size(choices,1)
    % Assign level
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels = choices{i,2};
    
    % Assign contrast files
    contrastfiles{i}=cell(length(subjectlist),1);
    for s=1:length(subjectlist)
        contrastfiles{i}{s}=[where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread  filesep firstlevelmodel ' Contrasted' filesep choices{i,3} ',1'];
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans =contrastfiles{i};
end

% Include subjects as covariates
for s=1:length(subjectlist)
    sub=zeros(1,length(subjectlist)); sub(s)=1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(choices,1)])';
   matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
end

%% Run model

% Execute (Specify only)
disp('Specifying model --------------------------')
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (2) Estimate model ##############################
disp('Estimating model ------------------------------')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (3) Add sensible contrasts
disp('Adding sensible contrasts ------------------------------')
matlabbatch{1}.spm.stats.con.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 0; c=1;

% Identity matrix
matlabbatch{1}.spm.stats.con.consess{c}.fcon.name='Identity'; % Identity matrix (F Contrast)
matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec ={[1 0 0 0]; [0 1 0 0];[0 0 1 0]; [0 0 0 1];};
matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none'; c=c+1; 

% Sensible contrasts
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'cF-ct'; % Compare tasks
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 1 -1 -1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'ct-cF';
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [ -1 -1 1 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'Rej-Explore'; % Compare Reject & Explore
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 -1 1 -1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'Explore-Reject';
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 1 -1 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'cF_Rej-Exp'; % Comparing cells (Ex vs VEx)
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 -1 0 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'cF_Exp-Rej'; 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 1 0 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'ct_Bomb-Exp'; 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 0 1 -1 ];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'ct_Exp-Bomb'; 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 0 -1 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'Rej_cF-ct';  % Comparing cells (Task)
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'Rej_ct-cF';  
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 0 1 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'Exp_cF-ct';
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 1 0 -1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'Exp_ct-cF';  
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 -1 0 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;

%
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (4) Save details to 2nd level folder
details.RLvariables=RLvariables;
save([secondlevelfolder filesep 'details_2ndlevel.mat'], 'details', 'log','contrastfiles'); % Save details in 2nd level folder


end

