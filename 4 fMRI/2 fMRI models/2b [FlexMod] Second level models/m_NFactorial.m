function [ matlabbatch ] = m_NFactorial(where, log, Des)
% [ matlabbatch ] = m_NFactorial(where, log, Des)
%  Set up 2nd-level N-way factorial ANOVA in batch system (Specify &
%  Estimate). Additional sensible contrasts may be added in the pipeline,
%  but they are NOT added here. 
%
% Des.Directory:    Results folder(first level / second level)
% Des.Design:     Cell, with each cell holding details for each factor in the factorial 
%                       (.name, .levels, .dept, .variance, .gmsca, .ancova)
% Des.Cells:        Cell matrix, allocating contrast images to factorial cells
%                       (Col 1=Contrast image, Col 2=Factor 1 level, Col 2=Factor 2 level ...)
% Des.Covar(n):    Covariate details (subject-specific)
%                             .Ncells=No. of contrast images per subject
%                             .CoVector=Array containing subject-specific covariate number
%                             .CovNames=Cell vector containing subject-specific covariate names
%
% ------------------------------------------------------------------------------------

% Execute to debug: m=1; Des=Con.SecondLevelAnalysis{m,3};


for o1=1:1 % Fake inputs 
% Des.Directory='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\testfolder';
% Des.Design{1}.name='Task';  % Factorial design
% Des.Design{1}.levels=2;
% Des.Design{1}.dept=0; Des.Design{1}.variance=1; Des.Design{1}.gmsca=0;  Des.Design{1}.ancova=0;
% Des.Design{2}.name='Choice';
% Des.Design{2}.levels=2;
% Des.Design{2}.dept=0; Des.Design{2}.variance=1; Des.Design{2}.gmsca=0;  Des.Design{2}.ancova=0;
% Des.Cells={'con0006.img' 1 1;                 'con0007.img' 1 2;                 'con0008.img' 1 3;                 'con0009.img' 2 1;                 'con0010.img' 2 2;                'con0011.img' 2 3; };
% Des.Covar=[];
end

if isdir(Des.Directory)==0; mkdir(Des.Directory); end

%% Set up Factorial analysis in batch

% (1) Load details for the factorial
matlabbatch{1}.spm.stats.factorial_design.dir = {Des.Directory};
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
for f=1:length(Des.Design) % Details about factors 
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).name=Des.Design{f}.name;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).levels=Des.Design{f}.levels;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).dept=Des.Design{f}.dept;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).variance=Des.Design{f}.variance;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).gmsca=Des.Design{f}.gmsca;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).ancova=Des.Design{f}.ancova;
end

% (2) Covariates?
if isempty(Des.Covar)==0
    c=1;
    for j=1:length(Des.Covar)
        switch Des.Covar(j).Type
            
            % Split Covariates ----------------------------------
            case 'split'
                
                for s=1:log.n_subjs
                    cov=zeros(1,log.n_subjs); cov(s)=Des.Covar(j).CoVector(s);
                    %
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).c = repmat(cov,[1 Des.Covar(j).Ncells]);
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = Des.Covar(j).CovNames{s};
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1;
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;
                    %
                    c=c+1;
                end
                
            % Single-vector Covariates ----------------------------------
            case 'single'                
                matlabbatch{1}.spm.stats.factorial_design.cov(c).c =Des.Covar(j).CovSingle_Vector;
                matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = Des.Covar(j).CovSingle_Name;
                matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1;
                matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;     
                c=c+1;
        end
    end
else
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end

% (3) Factorial cells: Assign contrast images to factorial cells
for c=1:size(Des.Cells,1)
    
    % Cell assignment (for N factors)
    for f=1: size(Des.Cells,2)-1
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(c).levels(f)= Des.Cells{c,f+1};
    end

    % Load contrast images for all subjects for each cell
    for s=1:log.n_subjs
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(c).scans{s,1}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep Des.Cells{c,1} ',1'];
    end
end

% Save details in 2nd level folder
save([Des.Directory filesep 'details_2ndlevel.mat'], 'log','Des'); 

%% Estimate model

matlabbatch{2}.spm.stats.fmri_est.spmmat = {[Des.Directory filesep 'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


end

