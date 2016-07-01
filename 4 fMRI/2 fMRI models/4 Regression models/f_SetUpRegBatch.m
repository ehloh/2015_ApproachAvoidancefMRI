function [ matlabbatch] =f_SetUpRegBatch(log,where, behvar, TargetConName)
% [ matlabbatch] =f_TrialTypeVariables(log,where, behvar, RLvar)
%   Set up regression model for a contrast from the t1_TrialType model
% ------------------------------------------------------------------------------------

%%

% Set up batch
matlabbatch{1}.spm.stats.factorial_design.dir = {where.secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Load behavioural variables as covariates
for i=1:size(log.BehaviouralVariables,1)
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).c = behvar(:,i);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(i).cname = log.BehaviouralVariables{i};
end


% Identify contrast image (sample 1st subject)
p=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep log.FLfolname filesep 'SPM.mat']);
SubCons=cell(size(p.SPM.xCon,2),2);
for i=1:size(p.SPM.xCon,2);
    SubCons{i,1}=p.SPM.xCon(i).name;
    SubCons{i,2}=p.SPM.xCon(i).Vcon.fname;
end
ConImage =char(SubCons(find(strcmp(SubCons(:,1),TargetConName)),2));

% Load scans 
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans=cell(log.n_subjs,1);
for s=1:log.n_subjs
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level'  filesep log.FLfolname filesep ConImage ',1'];
end

end

