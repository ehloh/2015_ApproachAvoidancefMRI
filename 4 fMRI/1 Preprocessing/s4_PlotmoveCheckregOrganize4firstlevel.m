% Plot movement, Organize, Checkreg
clear all; close all hidden; clc

% Requested analysis
request.plotmovement=1;
request.checkreg=0;
request.organize4firstlevel=0;
log.specificsubjects={}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications    
    
    % Load subjects
% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_brain='/Volumes/SANDISK/2 EXPLORE Brain data';
where.where='D:\Dropbox\SANDISK\'; where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data';  log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;

    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    errorlog=cell(1,1); e=1;
    
    % Settings that don't really change
    request.organize_preprocessedfiles='swu';
    request.organize_zipunused=0;
    request.organize=0;
    
    % Interface
    disp('===================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested:'); disp(request)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.specificsubjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); disp(' ')
    input('Hit Enter to start      ')
    disp('====================================')
    
end

%% (1) Plot movement regressors

if request.plotmovement==1
    disp('----------------- Collecting movement parameters ----------------- ')
    movement=cell(log.n_subjs,7);
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '  (' log.subjects{s} ')'])
        for b=1:6
            f=spm_select('List', [where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep 'Func_r' num2str(b)], '.txt$');
            if isempty(f)==1
                movement{s,b}=[];
            else
                movement{s,b}=load([where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep 'Func_r' num2str(b) filesep f]);
            end
        end
        movement{s,7}=log.subjects{s};
    end
    
    % Plot
    figure('Name', 'Movement parameters', 'Position', [400 100 1200 1000])
    for s=1:log.n_subjs
        subplot(ceil(log.n_subjs/5), round(log.n_subjs/ceil(log.n_subjs/5)), s);
        for i=1:6
            plot(movement{s,i}); hold on;
        end
        text(0,1, log.subjects{s}); text(0,1, log.subjects{s}); text(0,1, log.subjects{s})
        axis tight
    end
    
end

%% (2) CheckReg 

if request.checkreg==1
    disp('--------------- Collecting scans for CheckReg --------------------')
    w.nScansPerFigure=10;
    %
    images=cell(log.n_subjs*2,1); i=1;
    for s=1:log.n_subjs
        for b=1:6
            wb.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_r' num2str(b) filesep];
            % Identify scans
            f=spm_select('List', wb.where, ['^' request.organize_preprocessedfiles '.*img$']);
            if isempty(f)==1
                disp(['Error: No available scans for ' log.subjects{s} '  block ' num2str(b) ' ---'])
            else
                images{i,1}=[wb.where f(randi(size(f,1)),:) ',1'];
                images{i,2}=log.subjects{s};  i=i+1;
            end
        end
    end
    
    % Collate into figures
    f=1; i=1; figs=cell(ceil(size(images,1)/w.nScansPerFigure),1);
    for d=1:size(images,1)
        figs{f}{i,1}=char(images{d,1});
        images{d,3}=f;
        if i==w.nScansPerFigure
            i=1; f=f+1;
        else
            i=i+1; 
        end
    end
    
    % Display / Instructions for display
    if w.nScansPerFigure>=size(images,1)
        matlabbatch{1}.spm.util.checkreg.data=figs{f};
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    else
        disp(' --------------------- INSTRUCTIONS FOR CHECKREG -------------------------------------')
        disp(['Multiple figures to display (' num2str(size(figs,1)) ' figures)'])
        disp('Note: To change no. of images per display, specify in script (at beginning of module)')
        disp('See variable ''images'' to identify subject scans (Col 3)')
        disp(' ')
        disp('To display each batch, specify value of ''f'' and execute following command:')
        disp(' ')
        disp('f = # ; eval(checkregcommand) ')
        disp(' ')        
        disp(' -----------------------------------------------------------------------------------------------')
        spm_jobman('initcfg');  checkregcommand='matlabbatch{1}.spm.util.checkreg.data=figs{f}; spm_jobman(''run'' , matlabbatch);';
    end
end

%% (3) Prep for First level

if request.organize4firstlevel==1
    disp('--------------- Preparing folder for 1st level --------------------')
    for s=1:log.n_subjs
        disp(['Subject  ' num2str(s) '  (' log.subjects{s} ') '])
        wb.where=[where.data_brain filesep log.subjects{s} filesep];
        wb.wherefrom=[wb.where '1 Preprocessed' filesep];
        wb.whereto=[wb.where '2 First level' filesep];
        if isdir(wb.whereto)==0; mkdir(wb.whereto); end;
        if isdir([wb.whereto 'Preproc functionals'])==0; mkdir([wb.whereto 'Preproc functionals']); end
        % Collect functionals
        for r=1:6
            cd([wb.wherefrom 'Func_r' num2str(r)])
            wb.files=dir([request.organize_preprocessedfiles '*']);
            if size(wb.files,1)==9
                errorlog{e}=['Error in prep for 1st level: Can''t find preproc functionals   ('  log.subjects{s} ' -  run  ' num2str(r) ')']; disp(errorlog{e}); e=e+1;
            end
            for i=1:size(wb.files,1)
                copyfile([wb.wherefrom 'Func_r' num2str(r) filesep wb.files(i).name],[wb.whereto 'Preproc functionals' filesep wb.files(i).name]);
            end
        end
    end
end

%%  

disp(' ###################################################')
disp(' DONE ' )
if request.checkreg==1
    disp('CHECKREG: See command window (earlier command) for display instructions')
end
disp(' ###################################################')



