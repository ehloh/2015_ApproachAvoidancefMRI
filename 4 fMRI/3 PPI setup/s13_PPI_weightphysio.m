% PPI - single (in SPM)
clear all;close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_brain='/Volumes/SANDISK/2 EXPLORE Brain data';
where.where='I:\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
% where.where='C:\Users\eloh\Desktop\2 [Explore]\5 Explore fMRI';

% Requested analysis ###############
log.specificsubjects={};
%
request.weightphysio_firstlevel=1;
request.weightphysio_secondlevel=1;
%
log.firstlevelmodel='m_c3_CompeteFull_XUVPEN';
log.ppimodel='compFam_HPC_L_strictanatTcf_psy_cF_Rej-Exp';


for o1=1:1 % General settings and specifications 
   
    % Add paths
    addpath(where.where); 
    log.func_prefix='s';
    
   % Load subjects
   log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
   [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
   %
   w.modelsneedingsubselect={'m_c6_Cluster4CompeteFull_XUVPEN';'m_c7_Cluster6CompeteFull_XUVPEN';'m_c8_Cluster4MovCompeteFull_XUVPEN';'m_c9_Cluster6MovCompeteFull_XUVPEN';'m_c10_Cluster6CompeteRT_XUVPEN'};
    if sum(strcmp(log.firstlevelmodel, w.modelsneedingsubselect))==1
        disp('Executing additional subject selections for this univariate model');
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep '4 Set up models' filesep 'i_Subjectdataok_SpecificModels.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    end
    where.subFLfol=cell(log.n_subjs,1);    
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.firstlevelmodel ' (' date ')' ])
    errorlog=cell(1,1); e=1;

    for s=1:log.n_subjs  % Locations of 1st level models  
        where.subFLfol{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep ];
%         if strcmp(log.firstlevelmodel(1), 'f')==1;
%             where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_f1 Contrasted   ' log.firstlevelmodel filesep];
%             log.univmod_type='flex';
%         elseif strcmp(log.firstlevelmodel(1), 't')==1;
%             where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_t1 Contrasted   ' log.firstlevelmodel filesep];
%             log.univmod_type='flex';
%         else
%             where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.firstlevelmodel ' Contrasted' filesep];
%             log.univmod_type='par';
%         end
    end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.experiment_folder]);
    disp(' '); disp('Requested analysis:'); disp(request); disp(' ');
    disp(' '); disp('             -------------------  CHECK HERE  -------------------'); disp(' ')
    disp(['Univariate model: '  log.firstlevelmodel ]);
    disp(['PPI model:   ' log.ppimodel])
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
end

%%

% Physiological variable = 3rd regressor
if request.weightphysio_firstlevel
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) '  - ' log.subjects{s} ' -----------'])
    ws.whereppi=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.firstlevelmodel ' PPIs' filesep log.ppimodel];
    matlabbatch{1}.spm.stats.con.spmmat = {[ws.whereppi filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0;
    %
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Physio';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [0 0 1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    %
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];ws=[];    
end

end



if request.weightphysio_secondlevel
where.secondlevelres=[where.experiment_folder filesep '2 Second level results' filesep log.firstlevelmodel filesep log.ppimodel filesep];
log.PPImodelname=log.ppimodel;

    if isdir(where.secondlevelres)==0
        mkdir(where.secondlevelres);
        mkdir([where.secondlevelres  'ROI' filesep]);
    end

    disp(['Running 2nd level model: ' log.ppimodel ' ################################'])
    [batch]=onesampleT(where,log,'Physio');
end


%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
% disp(['GLM model: ' log.AnalysisType '      ' log.firstlevelmodel '     ' log.secondlevelmodel]); disp(' ')
disp('Analysis completed:'); disp(request);
disp('Errors:'); disp(errorlog)
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', ['Analysis batchscript is complete (' mfilename ')'], ' ',1);
end
