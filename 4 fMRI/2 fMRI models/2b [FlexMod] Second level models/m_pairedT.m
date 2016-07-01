function [ matlabbatch] = m_pairedT(where, log, Des)
% [ matlabbatch] = m_pairedT(where, log, Des)
%
%   Des.Directory:        Folder for this 2nd level model result
%   Des.ConImage1:      Con image 1
%   Des.ConImage2:      Con image 2
%   Des.Covar(n):          Covariate details (subject-specific)
%                                         .Ncells=No. of contrast images per subject
%                                         .CoVector=Array containing subject-specific covariate number
%                                         .CovNames=Cell vector containing subject-specific covariate names
%
%
% ---------------------------------------------------------------------------------

if isdir(Des.Directory)==0; mkdir(Des.Directory); end

%%

% Set up 
matlabbatch{1}.spm.stats.factorial_design.dir = {Des.Directory};
matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

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

% Load each subject's con images
for s=1:log.n_subjs 
    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{1}=[where.data_brain filesep log.subjects{s} '2 First level' filesep  'm_' log.firstlevel_contraststype ' Contrasted' filesep Des.ConImage1 ',1'];
    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{2}=[where.data_brain filesep log.subjects{s} '2 First level' filesep  'm_'  log.firstlevel_contraststype ' Contrasted' filesep Des.ConImage2 ',1'];
end

% Save details in 2nd level folder
save([Des.Directory filesep 'details_2ndlevel.mat'], 'log','Des'); 

end

