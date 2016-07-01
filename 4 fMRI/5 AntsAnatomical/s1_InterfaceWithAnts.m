% Interface with ANTs, transferring files between SPM data folder and Ants (Neurodebian shared) folder
clear all;close all hidden; clc

% Steps?
% request.AntsMethod='Template';
request.AntsMethod='Basic';
% request.AntsMethod='LM3';
for o=1:1 
% log.onsetsmodel='m_c3_ChoiceFull_OULPEN';  
% log.onsetsmodel='m_c4_Choice_OUPEN'; 
% log.onsetsmodel='m_c6_ChCluster4Full_OULPEN'; request.WhichContrasts=1:4;
% log.onsetsmodel='m_c7_ChCluster6Full_OULPEN'; request.WhichContrasts=1:4;
% log.onsetsmodel='m_v1c_vChoice_bpm16bpmi11'; request.WhichContrasts=1:5; 
% log.onsetsmodel='m_v3c_vChosenAnd_bpm16bpmi11'; request.WhichContrasts=1:4; 
% log.onsetsmodel='m_v4c_ChoicevChosenAnd_bpm16bpmi11'; request.WhichContrasts=1:10; 
% log.onsetsmodel='m_v4d_ChoicevChosenAnd_bpm16bpm11'; request.WhichContrasts=1:10; 
% log.onsetsmodel='m_v5c_ChoiceXvChosenAnd_bpm16bpmi11';  request.WhichContrasts=20:23;
% log.onsetsmodel='m_v6_ChoicevChosenAndposneg2_bpmi16bpmi11'; request.WhichContrasts=1:11; 
% log.onsetsmodel='m_v6c_ChoicevChosenAndposneg2_bpm16bpmi11'; request.WhichContrasts=7:11; 
% log.onsetsmodel='m_v6e_ChoicevChosenAndposneg2_b01b01'; request.WhichContrasts=1:5;
% log.onsetsmodel='m_v8c_ChoicesubEVposneg_bpm16bpmi11'; request.WhichContrasts=1:3; 
% log.onsetsmodel='m_t1_1_Trialtype';   request.WhichContrasts=1:78; % NOTE: all t1 models in same folder. Not differentiated between t1_1 & t1_2 etc.
% log.onsetsmodel='m_v9c_vChosenAndposneg2_bpm16bpmi11';   request.WhichContrasts=[2 3 5]; % 6:8;
% log.onsetsmodel='m_v11c_vGambleOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v11c_vGambleOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v12c_vModalchoiceOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v12c_vModalchoiceOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v13c_OutcomevGambleOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v14c_OutcomevModalchoiceOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoPE_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoVal_bpm16bpmi11';
% log.onsetsmodel='m_v16c_vChoicePE_bpm16bpmi11';
% log.onsetsmodel='m_v17c_vChoiceOutcomeMag_bpm16bpmi11';
% log.onsetsmodel='m_v18c_vChoiceAtOutcomeMag_bpm16bpmi11';
% log.onsetsmodel='m_v19c_VExploreInfoPEOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v20c_pExplore_bpm16bpmi11';
% log.onsetsmodel='m_v21c_ExploreGamInfoOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v1e_vChoice_b01b01';
% log.onsetsmodel='m_v22f_OutcomevChosenOutcomeMagnitude_b02b01';
% log.onsetsmodel='m_c1_Choice_ENU';
% log.onsetsmodel='m_v10c_RejectOrvChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_t2_1_TrialtypeNc';
% log.onsetsmodel='m_v1g_vChoice_bpji08bpji11';  
% log.onsetsmodel='m_c13_ChoiceFull_ULPEN';
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
% log.onsetsmodel='m_v23g_RejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v6g_ChoicevChosenAndposneg2_bpji08bpji11';
% log.onsetsmodel='m_v25g_RejectOrvGamble_bpji08bpji11';
% log.onsetsmodel='m_c14_Choice';
% log.onsetsmodel='m_v26e_pvBestChoice_b01b01';
% log.onsetsmodel='m_v24g_predChoice_bpji08bpji11';
% log.onsetsmodel='m_c15g_NextChoice_ULPEN';
% log.onsetsmodel='m_v4g_ChoicevChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v27g_ChoiceRejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v28g_ChoicePredChoice_bpji08bpji11';
% log.onsetsmodel='m_v29g_ChoiceRejectOrvGamble_bpji08bpji11';
% log.onsetsmodel='m_v30g_cFRejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v31g_ChoicecFRejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_c16_ChCluster4';
% log.onsetsmodel='m_c17_ChCluster6';
% log.onsetsmodel='m_v32g_vGamble_bpji08bpji11';
% log.onsetsmodel='m_v33g_ChoicevGamble_bpji08bpji11';
% log.onsetsmodel='m_v8g_vGamblePosNeg_bpji08bpji11'; 
% log.onsetsmodel='m_v34g_ChoiceXvGamble_bpji08bpji11';
% log.onsetsmodel='m_v35g_ChoiceChoiceXvGamble_bpji08bpji11';
% log.onsetsmodel='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11';
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v36g_RejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v38g_EVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v40g_vBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v42g_pLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v39g_ChoiceEVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v41g_ChoicevBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v43g_ChoicepLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v45g_ChoicevBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11'; 
% log.onsetsmodel='m_c13_ChoiceFull_ULPEN';
% log.onsetsmodel='m_c6_ChCluster4Full_ULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_ULPEN'; 
% log.onsetsmodel='m_c21_RejExpFull_ULPEN'; 
% log.onsetsmodel='m_c22_ChoiceFull_HLPEN'; 
end
log.onsetsmodel='m_c23_ChoiceFullRT_ULPEN'; 
% log.onsetsmodel='m_c24_ChCluster6RT_ULPEN'; 


request.FLthread=' s4Ants';
% request.WhichContrasts=7:16; % [10 11 15];  % Blank to process all
request.WhichContrasts=[];
% request.WhichContrasts=5:8;
%
request.ConvertLandmarks_nii=0;
request.TransferCons_Pc2Ants=0; 
request.TransferCons_Ants2Pc=1;

% Which subjects
log.specificsubjects={};
% log.specificsubjects={'p01_GV'}; 
% log.specificsubjects={
% %     'p01_GV';
% 'p02_YY';'p04_MP';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';
%     'p35_SM';'p36_FR';'p38_MK';'p41_AC'};

for o1=1:1 % General settings and specifications
    % Settings
    request.ConvertLandmarks_nii=0;
    log.AntsTypes={'Template';'Basic';'LM1'; 'LM2'; 'LM3';'LM4'};
    if sum(strcmp(log.AntsTypes, request.AntsMethod))~=1; error(['Invalid contrast-adjustment method selected: '  request.AntsMethod]); end
    w=strfind(log.onsetsmodel, '_'); request.FLmodShortname=log.onsetsmodel(1:w(2)-1);
    if isempty(strfind(request.FLthread, 'Ants')); error('Chose a FL model without Ants!!!'); end
    if isempty(strfind(log.onsetsmodel, 'Outcome'))==0
        disp('Altering shortname for outcome!')
        request.OutcomeType=log.onsetsmodel(strfind(log.onsetsmodel, 'Outcome')+7:18+strfind(log.onsetsmodel(20:end), '_b'));
        request.FLmodShortname =[request.FLmodShortname '_' request.OutcomeType];
    elseif isempty(strfind(log.onsetsmodel, 'ExploreInfo'))==0
        request.InfoType=log.onsetsmodel(strfind(log.onsetsmodel, 'ExploreInfo')+7:18+strfind(log.onsetsmodel(20:end), '_b'));
        request.FLmodShortname =[request.FLmodShortname '_' request.InfoType];
    end
    
    % Where
    where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI';  where.exp_folder='G:\2 [Explore]'; where.data_brain=[where.exp_folder filesep '1 Brain data']; 
    where.antsfolder='C:\Users\e.loh\Documents\Neurodeb\1_cF';
    if  isempty(strfind(log.onsetsmodel, 'Trialtype')) ==0; where.data_brain='I:\1 Explore fMRI'; end
    where.ants_AdjustCons=[where.antsfolder '\2b_AdjustCons'];  % In Host Neurodebian folder, where to? All subjects' T1s & cons go here (no subject folders)
    addpath(where.where)    
    
    % Load subjects (no selection of specific subjects!)
    log.exploreinfo_oksubs={'p02_YY';'p04_MP';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p35_SM';'p36_FR';'p38_MK';'p41_AC'};
    w.modelsneedingsubselect={'m_c6_';'m_c7_';'m_c8_';'m_c9_';'m_c10_';'m_v7_';'m_c17';'m_c16';'m_c24'};      % Apply further subject selection for some models
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    if sum(strcmp(log.onsetsmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep '4 Set up models' filesep 'i_Subjectdataok_SpecificModels.xlsx']); 
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.onsetsmodel);
    else [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ');      disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); 
    disp(' '); disp(['Ants adjusted by: ' request.AntsMethod]); disp(' '); 
    disp(['First level model:   ' log.onsetsmodel])
    disp(' ');  input('Hit Enter to start      ')
    disp('=======================================================')

end


%% Archived steps 1-2 (Template construction, adjust functional scans)

request.TransferT1sForTemplateConstruc=0;
request.Transfer_Func2Ants=0;
for o1=1:1 
% Step 1: Template constcruction
if request.TransferT1sForTemplateConstruc
    where.to_T1w=[where.antsfolder '\CreateTemplate'];  % In Host Neurodebian folder, where to?
%     where.to_T1w=[where.antsfolder '\3a_Data_CreateTemplatePartial']; disp('Partial volume template!!')
    
    where.T1wFol='1 Preprocessed\MPM\Antst1w';
    disp('Transferring T1ws for Template construction -----------------------------------')
    for s=1: log.n_subjs
        disp([log.subjects{s} ' ------------------']);
        ws.subfol_pc=[where.data_brain filesep log.subjects{s} filesep ];
        
        
        % Mean functiona (ubf), is coregistered to all functionals (non-normalized,
        % prefix ubf). Convert to nii format. - used to create a partial
        % volume template!        
%         f=spm_select('List', [ws.subfol_pc '1 Preprocessed\Func_b1\Preproc_b1'], '^meanubf.*.img'); if size(f,1)~=1; error('More than 1 mean functional found!'); end
%         v=spm_vol([ws.subfol_pc '1 Preprocessed\Func_b1\Preproc_b1\' f(1,:)]);
%         ima=spm_read_vols(v);
%         v.fname=[ws.subfol_pc where.T1wFol filesep f(1,1:length(f(1,:))-4) '.nii'];
%         spm_write_vol(v,ima);
        
        try
            f=spm_select('List', [ws.subfol_pc filesep where.T1wFol], '^bsMQ.*_T1w.nii'); % Bias corrected T1w structurals
%             f=spm_select('List', [ws.subfol_pc filesep where.T1wFol], '^meanubf.*.nii'); disp('partial volume template!')
            if size(f,1)~=1; error('No. found files ~=1'); end
            copyfile([ws.subfol_pc where.T1wFol filesep f], [where.to_T1w filesep log.subjects{s} '_' f]);
        catch
            disp('Did not')
        end
        
        ws=[];
        
    end
end

% Step 2: Adjust functional data directly
%   All scans must be coregistered with the T1s used for template construction. Perform this step just before smoothing
if request.Transfer_Func2Ants
    log.preants_prefix='ubfM';
    log.wherefuncs=[where.antsfolder filesep '2a_AdjustFuncs' filesep];
    
    disp('Transferring functional scans to be adjusted -----------------------------------')
    for s=1: log.n_subjs
        disp([log.subjects{s} ' ------------------']);
        ws.subfol_pc_func=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_r'];
        ws.subfol_ants=[log.wherefuncs log.subjects{s} filesep];
        ws.funclist={}; ff=1;
        
        for r=1:6
            f=spm_select('List', [ws.subfol_pc_func num2str(r)], ['^' log.preants_prefix '.*.img']);
            disp(['   run ' num2str(r) '  -  '  num2str(size(f,1)) ' files'])
            
            % Convert to nii format + Write to ants folder 
            for i=1:size(f,1)
                v=spm_vol([ws.subfol_pc_func num2str(r) filesep f(i,:)]);
                ima=spm_read_vols(v);
                v.fname=[ws.subfol_ants   f(i,1:length(f(i,:))-4) '.nii'];
                spm_write_vol(v,ima);
                
                % Create list 
                ws.funclist{ff,1}=[f(i,1:length(f(i,:))-4) '.nii']; ff=ff+1; 
            end
        end
        
        ws=[];
    end
    
    % Create list 
    
end

end


%% Landmarks files: hdr/img --> nii

if request.ConvertLandmarks_nii
    request.LMname=request.AntsMethod;
    for s=1:log.n_subjs
        disp(log.subjects{s}); cd([where.ants_AdjustCons filesep log.subjects{s}])
        
        % Make sure Landmarks are in same space as the subject structural
        spm_get_space([request.AntsMethod '.hdr'],  spm_get_space([log.subjects{s} '_T1w_coreg.nii']));
       
        % Convert
        v=spm_vol([request.LMname '.img']);
        ima=spm_read_vols(v);
        v.fname=[request.LMname '.nii'];
        spm_write_vol(v,ima);
        
        % Archive old file
        if isdir('LM archive')==0;  mkdir('LM archive'); end
        movefile([request.LMname '.img'], ['LM archive' filesep request.LMname '.img']);
        movefile([request.LMname '.hdr'], ['LM archive' filesep request.LMname '.hdr']);
       
        % Set up LM folder if not present, and move LM mask file there
        if isdir(request.LMname)==0;  mkdir(request.LMname); end
        movefile([request.LMname '.nii'], [request.LMname  filesep request.LMname '.nii']);
    end
    
    
    % Group template: make sure in correct space and format
    cd([where.antsfolder filesep '1a_CreateTemplate'])
    spm_get_space([where.antsfolder filesep '1a_CreateTemplate' filesep 'Landmarks' filesep request.AntsMethod '.hdr'],  spm_get_space('hc_template.nii')); 
    cd('Landmarks');
    v=spm_vol([request.LMname  '.img']); % Convert to .nii
    ima=spm_read_vols(v);
    v.fname=[request.LMname '.nii'];
    spm_write_vol(v,ima);
    movefile([request.LMname '.img'], ['LM archive' filesep request.LMname '.img']);  % Move
    movefile([request.LMname '.hdr'], ['LM archive' filesep request.LMname '.hdr']);
end


%% Step 3: Adjust first-level contrasts themselves
%   Must be coregistered with the T1s used for template construction

% ##### Transfer un-adjusted contrasts from PC to Ants ##############################
if request.TransferCons_Pc2Ants
    alreadytransferedsome=0;
    
    %
    if alreadytransferedsome==0;
        ws=load([where.data_brain filesep log.subjects{1} '\2 First level' request.FLthread '\' log.onsetsmodel ' Contrasted\SPM.mat']);
    else ws=load([where.data_brain filesep log.subjects{1} '\2 First level' request.FLthread '\' log.onsetsmodel '_' request.AntsMethod ' Contrasted\SPM.mat']);
    end
    
    if isempty(request.WhichContrasts); request.WhichContrasts=1:size(ws.SPM.xCon,2); end
    disp('Requested contrasts (to be moved to Ants for adjustment):   '); disp([repmat('       ', length(request.WhichContrasts),1) char(ws.SPM.xCon(request.WhichContrasts).name)]); input('Proceed?    ');
    
    for s=1:  log.n_subjs
        disp([log.subjects{s} ' ------------------']);
        if alreadytransferedsome==0;
            ws.pc_FLfol=[where.data_brain filesep log.subjects{s} '\2 First level' request.FLthread '\' log.onsetsmodel ' Contrasted\'];
        else ws.pc_FLfol=[where.data_brain filesep log.subjects{s} '\2 First level' request.FLthread '\' log.onsetsmodel '_' request.AntsMethod ' Contrasted\'];  % If this isnt the first ant stransfer
        end
        %
        ws.ants_Adjustfol=[where.ants_AdjustCons filesep log.subjects{s} filesep request.AntsMethod filesep];
        ws.ants_AdjustModfol=[ws.ants_Adjustfol        strtrim(request.FLthread(1:strfind(request.FLthread, 'Ants')-1)) request.FLmodShortname '\'];
        if isdir(ws.ants_Adjustfol)==0; mkdir(ws.ants_Adjustfol); end
        if isdir(ws.ants_AdjustModfol)==0; mkdir(ws.ants_AdjustModfol); end
        if isdir([ws.ants_Adjustfol 'InverseTransform'])==0; mkdir([ws.ants_Adjustfol 'InverseTransform']); end
        
        % Save original contrasts (in .nii format) if not done already (un-normalized, un-adjusted)
        if isdir([ws.pc_FLfol 'Original contrasts\'])==0
            mkdir([ws.pc_FLfol 'Original contrasts\'])            
            f=spm_select('List', ws.pc_FLfol, '^con.*.img');
            for i=1:size(f,1)
                
                movefile([ws.pc_FLfol f(i,:)], [ws.pc_FLfol 'Original contrasts\' f(i,:)])
                movefile([ws.pc_FLfol f(i,1:end-4) '.hdr'], [ws.pc_FLfol 'Original contrasts\' f(i,1:end-4) '.hdr'])
                %
                v=spm_vol([ws.pc_FLfol 'Original contrasts\' f(i,:)]);
                ima=spm_read_vols(v);
                v.fname=[ws.pc_FLfol 'Original contrasts\' f(i,1:end-4) '.nii'];
                spm_write_vol(v,ima);
                %
                delete([ws.pc_FLfol 'Original contrasts\' f(i,:)])
                delete([ws.pc_FLfol 'Original contrasts\' f(i,1:end-4) '.hdr'])
            end    
        elseif alreadytransferedsome 
            f=cellfun(@(x)['con_' x(length(x)-3:end) '.img'], cellfun(@(x)num2str(100000+x), num2cell(request.WhichContrasts), 'UniformOutput',0),'UniformOutput',0)';
            f=char(f); 
            for i=1:size(f,1)
                movefile([ws.pc_FLfol f(i,:)], [ws.pc_FLfol 'Original contrasts\' f(i,:)])
                movefile([ws.pc_FLfol f(i,1:end-4) '.hdr'], [ws.pc_FLfol 'Original contrasts\' f(i,1:end-4) '.hdr'])
                %
                v=spm_vol([ws.pc_FLfol 'Original contrasts\' f(i,:)]);
                ima=spm_read_vols(v);
                v.fname=[ws.pc_FLfol 'Original contrasts\' f(i,1:end-4) '.nii'];
                spm_write_vol(v,ima);
                %
                delete([ws.pc_FLfol 'Original contrasts\' f(i,:)])
                delete([ws.pc_FLfol 'Original contrasts\' f(i,1:end-4) '.hdr'])
            end  
        end
        
        % Move over!
        for c=1:length(request.WhichContrasts)
            ws.cname=num2str(10000+request.WhichContrasts(c));
            copyfile([ws.pc_FLfol 'Original contrasts\con_' ws.cname(2:end) '.nii'], [ws.ants_AdjustModfol 'con_' ws.cname(2:end) '.nii'])
        end
        ws=[];
    end
end

% ##### Transfer successfully-adjusted contrasts from Ants to PC ##############################
if request.TransferCons_Ants2Pc
    request.Convert2HdrFirst=1; 

    alreadytransferedsome=0;
    
    if alreadytransferedsome
        ws=load([where.data_brain filesep log.subjects{1} '\2 First level' request.FLthread '\' log.onsetsmodel '_' request.AntsMethod ' Contrasted\SPM.mat']);
    else ws=load([where.data_brain filesep log.subjects{1} '\2 First level' request.FLthread '\' log.onsetsmodel ' Contrasted\SPM.mat']);
    end
    
    if isempty(request.WhichContrasts);  f=cellstr(spm_select('List', [where.ants_AdjustCons filesep log.subjects{1} filesep request.AntsMethod filesep  request.FLthread(2:3) request.FLmodShortname  filesep], '^r.*.nii'));  request.WhichContrasts=cell2mat(cellfun(@(x)str2double(x(8:10)), f, 'UniformOutput',0)); end
    disp('Requested adjusted contrasts (to be moved back to PC):   '); disp([repmat('       ', length(request.WhichContrasts),1) char(ws.SPM.xCon(request.WhichContrasts).name)]); input('Proceed?    ');
    for s=1: log.n_subjs
        disp(['Subject ' num2str(s) ' :  ' log.subjects{s} ' ------------------']);
        ws.pc_FLfol=[where.data_brain filesep log.subjects{s} '\2 First level' request.FLthread '\' log.onsetsmodel ' Contrasted\'];
        ws.pc_FLfoladjusted=[where.data_brain filesep log.subjects{s} '\2 First level' request.FLthread '\' log.onsetsmodel '_' request.AntsMethod ' Contrasted\'];
        ws.ants_Adjustfol=[where.ants_AdjustCons filesep log.subjects{s} filesep request.AntsMethod filesep];
        ws.ants_AdjustModfol=[ws.ants_Adjustfol  strtrim(request.FLthread(1:strfind(request.FLthread, 'Ants')-1)) request.FLmodShortname '\'];
        
        % Convert files .nii --> .hdr/.img  (stays in Ants folder)
        if request.Convert2HdrFirst
            if isdir([ws.ants_AdjustModfol 'Formatted for SL'])==0; mkdir([ws.ants_AdjustModfol 'Formatted for SL']); end
            for c=1:length(request.WhichContrasts)
                ws.c=num2str(10000+request.WhichContrasts(c));
                matlabbatch{1}.spm.util.imcalc.input = {[ws.ants_AdjustModfol 'r_con_' ws.c(2:end) '.nii']};
                matlabbatch{1}.spm.util.imcalc.output = ['con_' ws.c(2:end) '.img'];
                matlabbatch{1}.spm.util.imcalc.outdir = {[ws.ants_AdjustModfol 'Formatted for SL']};
                matlabbatch{1}.spm.util.imcalc.expression = 'i1';
                matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                matlabbatch=[];
            end
        end
        
        % Create new first-level folder for this Ants adjustment x FL model !!! ####
        if isdir(ws.pc_FLfoladjusted)==0; copyfile(ws.pc_FLfol, ws.pc_FLfoladjusted);
        elseif s==1 && isempty(spm_select('List', ws.pc_FLfoladjusted, '^con.*.img'))==0; input('There are some con files in the ants-adjusted FL folder already. Proceed?'); 
            % Are they existing cons in the new FL folder?
        end
        
        % Transfer over!!
        for c=1:length(request.WhichContrasts)
            ws.c=num2str(10000+request.WhichContrasts(c));
            copyfile([ws.ants_AdjustModfol 'Formatted for SL\con_' ws.c(2:end) '.hdr'], [ws.pc_FLfoladjusted 'con_' ws.c(2:end) '.hdr'])
            copyfile([ws.ants_AdjustModfol 'Formatted for SL\con_' ws.c(2:end) '.img'], [ws.pc_FLfoladjusted 'con_' ws.c(2:end) '.img'])
        end
        
        ws=[];
    end
    

    % Notify!!
    disp('======================================================='); w.c=clock;
    disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(' '); disp(['Data location (brain): ' where.data_brain]); disp(' ')
    disp(['First level model thread:    ' request.FLthread])
    disp(['Ants adjustment type:   '  request.AntsMethod])
    disp(['First level model:  ' log.onsetsmodel])
    disp('FL contrasts adjusted:  '); disp(request.WhichContrasts)
    disp('=======================================================')
    
    try % Notify researcher
        f_sendemail('kurzlich', strcat('Analysis script is complete (s1InteraceWithAnts)'), ' ',1);
    end
end


