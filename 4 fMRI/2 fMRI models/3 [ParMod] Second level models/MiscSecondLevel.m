function [  matlabbatch details ] = MiscSecondLevel(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% function [ matlabbatch details ] = MiscSecondLevel( input_args )
%   Manually specify a second level!!  


%% (1) Details for this model

% % Set up folder 
secondlevelfolder=[where.resultsfolder filesep secondlevelmodel];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); mkdir([secondlevelfolder filesep 'ROI']); end
FLfols = cellfun(@(x)[where.data_brain filesep x filesep '2 First level' log.FirstLevelThread  filesep firstlevelmodel ' Contrasted' filesep ], subjectlist, 'UniformOutput',0); 
details.designcells={};

%% (2) Choose contrasts to include (sampled from first subject)


% Available contrasts
w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
spm=w.SPM;   contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% Choose contrasts ######
% FLcons={'cF_vChosen'; 'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen';  'ct_vChosen'; 'ct_vBestUnchosen';};      % <--- Specify123 -------
% % FLcons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg';  'ct_vChosen'; 'ct_vBestUnchosen_pos';};
% FLcons={'cF_vGamblePos'; 'cF_vGambleNeg'; 'ct_vGamble'};
% FLcons={'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vBestUnchosen_pos'};
FLcons={'cF_vBUpos_vMargChosen'; 'cF_vBUneg_vMargChosen'; 'ct_vBUpos_vMargChosen'};
% FLcons={'cF_EVGain'; 'cF_EVLoss'; 'ct_EVGain'};


disp('Ingoing contrasts:'), disp(FLcons),  input('Continue?');
for i=1:size(FLcons,1)
       a.which=strcmp(contrasts(:,1), FLcons{i});
        if sum(a.which)>1, error(['ERROR: More than 1 contrast file found for RL variable ' FLcons{i}])
        elseif sum(a.which)==0,  error(['ERROR: Could not find contrast file for RL variable ' FLcons{i}])
        end
        FLcons{i,3}=contrasts{find(a.which(:,1)),2};
end

%% SECOND LEVEL options
% All lines that need to be manually specified are marked: <--- Specify123 ------- 


onewayanova=1; 
for o=1:1 % ONE WAY ANOVA

    if onewayanova
        % Manually specify 1-way ANOVA contrasts                                                                        % <--- Specify123 -------
        cellspex={1 [1] 'cF_vBUpos_vMargChosen'};  % cell #, cell level in the ANOVA, Contrast name
        details.designcells= [details.designcells; cellspex];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellspex{1}).levels = cellspex{2};   % [1 2]
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellspex{1}).scans =cellfun(@(x)[x FLcons{strcmp(FLcons(:,1), cellspex{3}),3} ',1'], FLfols , 'UniformOutput', 0);
        %
        cellspex={2 [2] 'cF_vBUneg_vMargChosen'};  % cell #, cell level in the ANOVA, Contrast name
        details.designcells= [details.designcells; cellspex];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellspex{1}).levels = cellspex{2};   % [1 2]
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellspex{1}).scans =cellfun(@(x)[x FLcons{strcmp(FLcons(:,1), cellspex{3}),3} ',1'], FLfols , 'UniformOutput', 0);
        %
        cellspex={3 [3] 'ct_vBUpos_vMargChosen'};  % cell #, cell level in the ANOVA, Contrast name
        details.designcells= [details.designcells; cellspex];
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellspex{1}).levels = cellspex{2};   % [1 2]
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellspex{1}).scans =cellfun(@(x)[x FLcons{strcmp(FLcons(:,1), cellspex{3}),3} ',1'], FLfols , 'UniformOutput', 0);
        
        
        % Batch
        matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'vBUtype';     % <--- Specify123 -------
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = size(details.designcells,1);    
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;     % <----------
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;       % <----------
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % Include subjects as covariates
        for s=1:length(subjectlist)
            sub=zeros(1,length(subjectlist)); sub(s)=1;
            matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(details.designcells,1)])';
            matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
            matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
        end
        
        % SPM: Specify, Estimate, Identity SL contrast
        disp('Specifying + estimating model --------------------------')
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
        matlabbatch{3}.spm.stats.con.delete = 0; c=1;
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.name='Identity';
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.convec = {eye(size(details.designcells,1))};   % {[1 0 0 0 0 0 ]; [0 1 0 0 0 0 ];[0 0 1 0 0 0]; [0 0 0 1 0 0]; [0 0 0 0 1 0]; [0 0 0 0 0 1]};
        matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'none'; c=c+1;
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        
        
        % RECORD DETAILS
        details.descrip='One-way ANOVA w vMargCho, split by vBUPosNeg and task';   % <--- Specify123 -------
        details.name='vMargCho x vBUpn';                           % <--- Specify123 -------
        details.matlabbatch=matlabbatch;
        save([secondlevelfolder filesep 'details2ndlevel.mat'], 'details')
        disp(['DONE! Details in SL folder. Manually rename SL folder to: ' details.name]);
        
    end
end

end

