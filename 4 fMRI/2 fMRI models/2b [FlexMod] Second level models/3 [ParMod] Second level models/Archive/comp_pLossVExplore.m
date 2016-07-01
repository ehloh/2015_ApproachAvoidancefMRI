function [ matlabbatch contrastfiles] = comp_pLossVExplore(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch contrastfiles] = comp_pLossVExplore(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% Second level models:
%     one sample ttests             pLoss-VExplore
%                                             VExplore-pLoss
%     ANOVA (Task x Comp)     pLoss-VExplore_cF    vs    pLoss-VExplore_ct
%                                            VExplore-pLoss_cF    vs    VExplore-pLoss_ct
%                                            (simpleFX & one-sample ttest for each cell)

% Execute following to use as script:   subjectlist=log.subjects; firstlevelmodel=log.onsetsmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% % Set up folder 
secondlevelfolder=[where.resultsfolder filesep secondlevelmodel];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); end

% Correct 1st level model?
OKmods={'m_c4_CompeteFull_XUPEN';};
if strcmp(OKmods, firstlevelmodel)==0; error('Invalid first-level model selected. Change first level model, or amend settings in script.'); end
    
% Empty output
contrastfiles=[];

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

% Available contrasts
w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']); spm=w.SPM;  
contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% [SampleT] Identify target contrasts
ConInstruc.onesamT={'pLoss-VExplore' []    [secondlevelfolder filesep 'pLoss-VExplore both'];
                                    'VExplore-pLoss' []    [secondlevelfolder filesep 'VExplore-pLoss both']};
for i=1:size(ConInstruc.onesamT,1)
    ConInstruc.onesamT{i,2}=contrasts{find(strcmp(ConInstruc.onesamT{i,1}, contrasts(:,1))),2};
end

% [ANOVA] Identify target contrasts
ConInstruc.anova(1,1:2)={{'pLoss-VExplore_cF' [];'pLoss-VExplore_ct' []}   [secondlevelfolder filesep 'pLoss-VExplore Anova']};
ConInstruc.anova(2,1:2)={{'VExplore-pLoss_cF' [];'VExplore-pLoss_ct' []}   [secondlevelfolder filesep 'VExplore-pLoss Anova']};
for a=1:length(ConInstruc.anova)
    for i=1:size(ConInstruc.anova{a},1)
        ConInstruc.anova{a,1}{i,2}=contrasts{find(strcmp(ConInstruc.anova{a,1}{i,1}, contrasts(:,1))),2};
    end
end

%% (3a) Implement model for one-sample ttests
% Specify & estimate 2nd level model, implement 2nd-level contrasts
% (positive and negative)

for t=1:size(ConInstruc.onesamT,1)
    disp(['Running model for one sample T #'  num2str(t) '  :  ' ConInstruc.onesamT{t,1} '---------------------------'])
    mkdir(ConInstruc.onesamT{t,3}); mkdir([ConInstruc.onesamT{t,3} filesep 'ROI']); 
    
    % (i) Specify model  ##############################
    matlabbatch{1}.spm.stats.factorial_design.dir = ConInstruc.onesamT(t,3);
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    for s=1:log.n_subjs % Load contrast files  
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{s}=[where.data_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep ConInstruc.onesamT{t,2} ',1'];
    end
    %
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
    % (ii) Estimate model  ##############################
    disp('Estimating model ------- ')
    matlabbatch{1}.spm.stats.fmri_est.spmmat ={[ConInstruc.onesamT{t,3} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
    % (iii) Add 2nd level contrast (Positive, negative) ##############################
    disp('Adding 2nd level contrasts ------ ')
    matlabbatch{1}.spm.stats.con.spmmat = {[ConInstruc.onesamT{t,3} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name =  'Pos';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name =  'Neg';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = -1;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    %
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
    % Save
    details.RLvariables=RLvariables;
    details.FLcons=ConInstruc.onesamT(t,:);
    save([ConInstruc.onesamT{t,3} filesep 'details_2ndlevel.mat'], 'details', 'log'); details=[];
end

%% (3b) Implement model for ANOVAs

for a=1:size(ConInstruc.anova,1)
    disp(['Running ONE-WAY anova #'  num2str(a) '  :  ' ConInstruc.anova{a,2} '---------------------------'])
    mkdir(ConInstruc.anova{a,2}); mkdir([ConInstruc.anova{a,2} filesep 'ROI']); 
    
    % (i) Specify One way ANOVA ##############################
    matlabbatch{1}.spm.stats.factorial_design.dir = ConInstruc.anova(a,2);
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    for i=1:size(ConInstruc.anova{a,1},1) % Allocate to cells: ONE-way ANOVA only
        
        % Cell allocation
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels =i;
        
        % Contrast image allocation
        for s=1:log.n_subjs 
            matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans{s}=[where.data_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep ConInstruc.anova{a,1}{i,2} ',1'];
        end
        
    end
    %
    for s=1:length(subjectlist) % Include subjects as covariates
        sub=zeros(1,length(subjectlist)); sub(s)=1;
        matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(ConInstruc.anova{a,1},1)])';
        matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
        matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
    end
    %
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
       
    % (ii) Estimate model ##############################
    disp('Estimating model ------------------------------')
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {[ConInstruc.anova{a,2} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
    % (iii) Add sensible contrasts ##############################
    disp('Adding sensible contrasts ------------------------------')
    matlabbatch{1}.spm.stats.con.spmmat = {[ConInstruc.anova{a,2} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0; c=1;
    % Identity
    matlabbatch{1}.spm.stats.con.consess{c}.fcon.name='Identity'; % Identity matrix (F Contrast)
    matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec = {[1 0]; [0 1]}';
    matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none'; c=c+1;
    % Sensible contrasts
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'cF'; % Identity, each task
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 0];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'ct'; % Identity, each task
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'cF-ct'; % Compare tasks
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  'ct-cF';
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [ -1 1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
    %
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
    % Save
    details.RLvariables=RLvariables;
    details.FLcons=ConInstruc.anova(a,:);
    save([ConInstruc.anova{a,2} filesep 'details_2ndlevel.mat'], 'details', 'log'); details=[];
end

end

