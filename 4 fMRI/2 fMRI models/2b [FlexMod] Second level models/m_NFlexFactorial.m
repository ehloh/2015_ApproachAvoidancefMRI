function [ matlabbatch] = m_NFlexFactorial(where, log, Des)
% [ matlabbatch ] = m_NFactorial(where, log, Des)
%  Set up 2nd-level Flexible factorial ANOVA in batch system (Specify & Estimate). 
%   Additional sensible contrasts may be added in the pipeline, but they are NOT added here. 
% 
% Des.Directory:                        Results folder(first level / second level)
% 
% Des.Design:                            Cell, with each cell holding details for each factor in the factorial 
%                                                     (.name, .dept, .variance, .gmsca, .ancova)
%                                                     
% Des.Cells:                              Cell matrix, allocating contrast images to factorial cells
%                                                     (Col 1=Contrast image, Col 2=Factor 1 level, Col 2=Factor 2 level ...)
% 
% Des.Cell_FacAssignment:       Cell vector, with factor assignment for each cell 
%                                                 (correponding to con images in Des.Cells, i.e. which factor assignment for each con image
% 
% Des.Covar(n):                       Covariate details (subject-specific)
% 
%                                                 .Ncells=No. of contrast images per subject
%                                                 .CoVector=Array containing subject-specific covariate number
%                                                 .CovNames=Cell vector containing subject-specific covariate names
%
% --------------------------------------------------------------------------------------------------------------------------------------------------------

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

if isdir(Des.Directory)==0; mkdir(Des.Directory); mkdir([Des.Directory filesep 'ROI']); end

%% Set up Factorial analysis in batch

% General settings
matlabbatch{1}.spm.stats.factorial_design.dir = {Des.Directory};
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Specify Factorial design
%         Note: These can be entered as the true factorial design. BUT, can also be
%         entered as a 1 x NCells design, and then the comparisons corresponding to
%         a factorial design added in the 2nd-level Contrasts. Always have
%         'Subject' as a separate factor, however.
for f=1:length(Des.Design)
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).name=Des.Design{f}.name;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).dept = Des.Design{f}.dept;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).variance =Des.Design{f}.variance;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).gmsca = Des.Design{f}.gmsca;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).ancova = Des.Design{f}.ancova;
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
                matlabbatch{1}.spm.stats.factorial_design.cov(c).c =repmat(Des.Covar(j).CovSingle_Vector, Des.NFactorialCells,1) ;
                matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = Des.Covar(j).CovSingle_Name;
                matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1;
                matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;     
                c=c+1;
        end
    end
else
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end



% Factors to include in the model (not necessarily factors of interest, i.e. all factors here)
for f=1:length(Des.ModelFac)
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{f}.fmain.fnum = Des.ModelFac(f);
end

% Load Scans + Conditions
for s=1:log.n_subjs
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).conds=Des.Cells_FacAssignment;
    for c=1:size(Des.Cells,1)
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(s).scans{c,1}= [where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep Des.Cells{c,3}  ',1' ];       
    end
end

% Save details in 2nd level folder
save([Des.Directory filesep 'details_2ndlevel.mat'], 'log','Des'); 

%% Estimate model

matlabbatch{2}.spm.stats.fmri_est.spmmat = {[Des.Directory filesep 'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Archive

for o1=1:1 % ORIGINAL BATCH CODE (archive) 
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Cell';
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 1;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 1;
% % matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
% % matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Subject';
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(log.n_subjs).conds = 1:NFactorialCells; 
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(1).scans = {
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p01_GV\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0001.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p01_GV\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0002.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p01_GV\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0003.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p01_GV\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0004.img,1'
%                                                                                   };
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(1).conds = [1 2 3 4];
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(2).scans = {
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p06_KB\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0001.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p06_KB\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0002.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p06_KB\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0003.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p06_KB\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0004.img,1'
%                                                                                   };
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(2).conds = [1 2 3 4];
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(3).scans = {
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p08_SG\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0001.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p08_SG\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0002.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p08_SG\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0003.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p08_SG\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0004.img,1'
%                                                                                   };
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(3).conds = [1 2 3 4];
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(4).scans = {
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p10_RC\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0001.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p10_RC\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0002.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p10_RC\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0003.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p10_RC\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0004.img,1'
%                                                                                   };
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(4).conds = [1 2 3 4];
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(5).scans = {
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p13_HL\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0001.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p13_HL\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0002.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p13_HL\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0003.img,1'
%                                                                                   'C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p13_HL\2 First level\m_f1_ChoicexTrialtype Contrasted flex4_ExploreOr_FixWind_4Most\con_0004.img,1'
%                                                                                   };
end

end



