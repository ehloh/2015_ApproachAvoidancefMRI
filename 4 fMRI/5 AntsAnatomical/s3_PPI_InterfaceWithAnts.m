% Interface with ANTs, transferring files between SPM data folder and Ants (Neurodebian shared) folder
clear all;close all hidden; clc


where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.exp_folder='C:\Users\eloh\Desktop\2 [Explore]'; 
where.data_brain=[where.exp_folder filesep '1 Brain data'];  where.antsfolder='D:\host\cF_study\2b_AdjustCons';
where.data_brainppi='F:\2 Explore fMRI\1 Brain data PPI';

% Steps?
log.specificsubjects={}; 
request.TransformedSeedImgs_ToPC=1;
request.ConsToAnts=0;
request.AdjustedConsToPC=0;

% PPI details
request.seed='sph_HPC_aL_cmSRvSN';
request.psych='DisN';
request.ppi_conditiontype=1; % 1=Context onsets, 2=Context pmod, 3= ContextOnsetMem

%
request.WhichAnalysisThread='s4Ants';
request.WhichFLmodel='m_ci3_ContextItemNomotor_Hit';
request.FL_shortname='s4FullCardiacWithDeriv_ci3Hit';
request.AdjustConMethod='Basic';


for o1=1:1 % General settings and specifications
    
    % Load subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
%     
%     % PPI details
%     switch request.ppi_conditiontype
%         case 1; request.ppi_condition_prefix='co';
%         otherwise; error('not specified yet')
%     end
%     log.ppiname=['PPI ' request.ppi_condition_prefix 'Fam_' request.seed '_psy_' request.psych];
%     log.pc_FLfol=['2 First level ' request.WhichAnalysisThread  filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep log.ppiname filesep];
%     
%     % Ants
%     log.AntsTypes={'Basic';'Landmarks';'Template';'Landmarks2';'Landmarks4';'Landmarks5';'Landmarks6'};
%     if isempty(strfind(request.FL_shortname, request.WhichAnalysisThread(1:end-4)))==1; error('Analysis Thread not in FL short name (for neurodebian folder)'); end
%     if isempty(strfind(request.FL_shortname, request.WhichFLmodel(3:5)))==1;     error('Requested FL folder not in FL short name (for neurodebian folder)'); end
%     if isempty(strfind(request.FL_shortname, request.WhichFLmodel(max(strfind(request.WhichFLmodel, '_'))+1:end)))==1;     error('Requested FL memtype not in FL short name (for neurodebian folder)'); end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ');     
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain]);
    disp(' ');  input('Hit Enter to start      ')
    disp('=======================================================')

end

%% (1) Bring ants-adjusted Seed imgs to subject PPI folders
%   Details of VOIs are in 2nd level PPI folder

if request.TransformedSeedImgs_ToPC
  
    request.Imgs2Bring={
        'sph_HPC_aL_STC';  'sph_HPC_aR_sTC';  % c3 TxC
        'sph_BA10_C'; 'sph_Precuneus'; 'sph_Striatum_R_C'; 'sph_SupMidFG'; % Choice Explore
        'sph_HPC_aL_sC'; 'sph_HPC_aR_saC'; % Choice HPC
        };
    
    % -----------------------------------------------------------------
    disp('VOIs to move over (file name, new voi name):'); disp(request.Imgs2Bring); input('OK?  ');
    disp('Bring VOIs (in subject space) to PC ################')
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.from=[where.antsfolder filesep log.subjects{s} filesep request.AdjustConMethod filesep 'InverseTransform' filesep];
        ws.to=[where.data_brainppi filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep 'VOIs' filesep];
        
        % Transfer all vois for this subject
        for v=1:size(request.Imgs2Bring,1)
            copyfile([ws.from log.subjects{s} '_' request.Imgs2Bring{v} '.nii'], [ws.to request.Imgs2Bring{v} '.nii'])
        end
        
        % Save details
        if exist([ws.to 'voi_details.mat'],'file')==0; vois=[];
        else;  ws.v=load([ws.to 'voi_details.mat']); vois=ws.v.vois;
        end
        vois=[vois; [{date} {request.FL_shortname} {request.AdjustConMethod} {request.Imgs2Bring}]];
        save([ws.to 'voi_details.mat'],'vois')
        
        
        ws=[]; vois=[];
    end
end

%%  (2) Unadjusted cons --> Ants 

if request.ConsToAnts
    request.saveoriginal=1;
    disp('Bring UnAdjusted contrasts from PC to Ants ################')
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.ants_fol=[where.antsfolder filesep log.subjects{s} filesep request.AdjustConMethod filesep];
        ws.ants_ppifol=[ws.ants_fol request.FL_shortname filesep log.ppiname(5:end) filesep]; mkdir(ws.ants_ppifol);
        ws.pc_ppifol=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep log.ppiname filesep];
        ws.pc_ppifol_cleancons=[ws.pc_ppifol 'PreAnts cons' filesep];
        
        
        % Save original contrasts imgs
        if request.saveoriginal
            if isdir([ws.pc_ppifol 'PreAnts cons'])==0
                mkdir([ws.pc_ppifol 'PreAnts cons'])
                movefile([ws.pc_ppifol 'con_0001.hdr'], [ws.pc_ppifol_cleancons 'con_0001.hdr'])
                movefile([ws.pc_ppifol 'con_0001.img'], [ws.pc_ppifol_cleancons 'con_0001.img'])
                movefile([ws.pc_ppifol 'con_0002.hdr'], [ws.pc_ppifol_cleancons 'con_0002.hdr'])
                movefile([ws.pc_ppifol 'con_0002.img'], [ws.pc_ppifol_cleancons 'con_0002.img'])
                % 
            else
                disp('PreAnts folder found. Not saving orignal cons, assumed already saved. OK?')
                if s==1; input('OK?   '); end
            end
        end
        
        % Convert
        v=spm_vol([ws.pc_ppifol_cleancons 'con_0001.img']) ;ima=spm_read_vols(v);
        v.fname=[ws.pc_ppifol_cleancons 'con_0001.nii']; spm_write_vol(v,ima);
%         v=spm_vol([ws.pc_ppifol_cleancons 'con_0002.img']) ;ima=spm_read_vols(v);
%         v.fname=[ws.pc_ppifol_cleancons 'con_0002.nii']; spm_write_vol(v,ima);
        
        % Copy over
        copyfile([ws.pc_ppifol_cleancons  'con_0001.nii'],  [ws.ants_ppifol 'con_0001.nii'])
%         copyfile([ws.pc_ppifol_cleancons  'con_0002.nii'],  [ws.ants_ppifol 'con_0002.nii'])
        
        %
        ws=[];
    end
end

%% (2) Adjust cons, from Ants --> SPM first-level models on PC

if request.AdjustedConsToPC
    
    disp('Bring Adjusted contrasts from Ants to PC ################')
    request.ConvertFirst=1;
    
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.ants_fol=[where.antsfolder filesep log.subjects{s} filesep request.AdjustConMethod filesep];
        ws.ants_ppifol=[ws.ants_fol request.FL_shortname filesep log.ppiname(5:end) filesep]; 
        ws.pc_ppifol=[where.data_brain filesep log.subjects{s} filesep '2 First level ' request.WhichAnalysisThread filesep request.WhichFLmodel '_' request.AdjustConMethod ' Contrasted' filesep log.ppiname filesep];
        ws.pc_ppifol_cleancons=[ws.pc_ppifol 'PreAnts cons' filesep];
        
        % Convert
            if request.ConvertFirst
                matlabbatch{1}.spm.util.imcalc.input = {[ws.ants_ppifol 'r_con_0001.nii,1']};
                matlabbatch{1}.spm.util.imcalc.output = 'r_con_0001.img';
                matlabbatch{1}.spm.util.imcalc.outdir = {ws.ants_ppifol};
                matlabbatch{1}.spm.util.imcalc.expression = 'i1';
                matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
%                 matlabbatch{2}.spm.util.imcalc.input = {[ws.ants_ppifol 'r_con_0002.nii,1']};
%                 matlabbatch{2}.spm.util.imcalc.output = 'r_con_0002.img';
%                 matlabbatch{2}.spm.util.imcalc.outdir = {ws.ants_ppifol};
%                 matlabbatch{2}.spm.util.imcalc.expression = 'i1';
%                 matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
%                 matlabbatch{2}.spm.util.imcalc.options.mask = 0;
%                 matlabbatch{2}.spm.util.imcalc.options.interp = 1;
%                 matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
                spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                matlabbatch=[];
            end
            
            % Copy to FL folder on PC
            copyfile([ws.ants_ppifol 'r_con_0001.hdr'], [ws.pc_ppifol 'con_0001.hdr']);
            copyfile([ws.ants_ppifol 'r_con_0001.img'], [ws.pc_ppifol 'con_0001.img']);
%             copyfile([ws.ants_ppifol 'r_con_0002.hdr'], [ws.pc_ppifol 'con_0002.hdr']);
%             copyfile([ws.ants_ppifol 'r_con_0002.img'], [ws.pc_ppifol 'con_0002.img']);
            
            ws=[];
    end
end





