% Preprocessing: realign & unwarp, coregister, segment, normalze, smooth
clear all;close all hidden; clc

% Requested analysis
process.realignunwarp=0;
process.coregister=0;
process.segment=0;
process.normalize=0;
process.smooth=0;
process.setup_newfol_4FL=1;
%
process.smoothingsize=3; % Smoothing size?
log.specificsubjects={}; % BLANK to process all subjects

where.data='F:\2 Explore fMRI\1 Brain data Preproc PPIold';
% where.data='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data';

for o1=1:1 % General settings and specifications
    
    % Load subjects
    where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; 
    where.data_beh=[where.where filesep 'Behavioural data'];
    addpath(where.where)
    log.w=load([where.data filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
  
    % Log
    diary([where.data filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:'); disp(process)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.specificsubjects); end
    disp(' '); disp(['Data location: ' where.data]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
end

%% Step 1: Realign & unwarp (prefix u)

if process.realignunwarp==1
    disp(' ########## (1) REALIGN & UNWARP: Create VDM files ##########')
    default.fieldmap=[where.spm filesep 'toolbox' filesep 'FieldMap' filesep 'pm_defaults_Trio_eFoV.m'];
    for s=1:length(log.subjects)
        disp(['Subject ' num2str(s)   '   (' log.subjects{s} ')  --------------- '])
        ws.where=[where.data filesep log.subjects{s} filesep '1 Preprocessed' filesep];
        % Which Fieldmap for which Run?
        ws.fm_run=[log.datalog{s+1,11} log.datalog{s+1,12} log.datalog{s+1,13} log.datalog{s+1,14} log.datalog{s+1,15} log.datalog{s+1,16}];
        ws.fm_nums=unique(ws.fm_run);
        ws.fm=cell(length(ws.fm_nums),1);
        for i=1:length(ws.fm_nums)
            ws.fm{i,1}=ws.fm_nums(i);
            ws.fm{i,2}=find(ws.fm_run(:)==ws.fm{i,1});
        end
        % Create VDM files
        for i=1:size(ws.fm,1)
%             try
                wb.n_mag=ws.fm{i,1};
                wb.n_phase=wb.n_mag+1;
                wb.whereFM=[ws.where 'Fieldmap_' num2str(ws.fm{i,1}) filesep];
                wb.wherefirstEPI=[ws.where 'Func_r' num2str(ws.fm{i,2}(1)) filesep];
                % Choose files
                f=spm_select('List', wb.whereFM, ['^sM.*.00' num2str(wb.n_phase) '-00001.*\.img$']); % Phase
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase=cellstr([wb.whereFM f ',1']);
                f=spm_select('List', wb.whereFM, ['^sM.*.00' num2str(wb.n_mag) '-00001-000001-01.img']); % Mag
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude=cellstr([wb.whereFM f ',1']);
                for r=1:length(ws.fm{i,2}) % pick first volume of each run
                    f=spm_select('list',[ws.where 'Func_r' num2str(ws.fm{i,2}(r)) filesep], '^bfM.*img$'); % first epi corresponding to that field map
                    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(r).epi=cellstr([ws.where 'Func_r' num2str(ws.fm{i,2}(r)) filesep f(1,:) ',1']);
                end
                %
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsfile = {default.fieldmap};
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
                matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;
                % Run !
                spm_jobman('initcfg'); spm_jobman('run',matlabbatch);
                matlabbatch=[];
                wb=[];
%             catch
%                 errorlog{e,1}=['Failed: Realign & Unwarp: Create VDM file  --- ' log.datalog{s+1,1} ',  Fieldmap no. ' num2str(ws.fm{i,1})]; disp(errorlog{e,1}); e=e+1;
%             end
        end
        ws=[];
    end
    
    % Execute unwarping ------------------
    disp(' ########## (1) REALIGN & UNWARP: Execute realigning & unwarping ##########')
    for o1 =1:1 % Realign & unwarp settings
        settings.realignunwarp.eoptions.quality = 0.9;
        settings.realignunwarp.eoptions.sep = 4;
        settings.realignunwarp.eoptions.fwhm = 5;
        settings.realignunwarp.eoptions.rtm = 0;
        settings.realignunwarp.eoptions.einterp = 2;
        settings.realignunwarp.eoptions.ewrap = [0 0 0];
        settings.realignunwarp.eoptions.weight = '';
        settings.realignunwarp.uweoptions.basfcn = [12 12];
        settings.realignunwarp.uweoptions.regorder = 1;
        settings.realignunwarp.uweoptions.lambda = 100000;
        settings.realignunwarp.uweoptions.jm = 0;
        settings.realignunwarp.uweoptions.fot = [4 5];
        settings.realignunwarp.uweoptions.sot = [];
        settings.realignunwarp.uweoptions.uwfwhm = 4;
        settings.realignunwarp.uweoptions.rem = 1;
        settings.realignunwarp.uweoptions.noi = 5;
        settings.realignunwarp.uweoptions.expround = 'Average';
        settings.realignunwarp.uwroptions.uwwhich = [2 1];
        settings.realignunwarp.uwroptions.rinterp = 4;
        settings.realignunwarp.uwroptions.wrap = [0 0 0];
        settings.realignunwarp.uwroptions.mask = 1;
        settings.realignunwarp.uwroptions.uwroptions.prefix = 'u';
    end
    for s=1:length(log.subjects)
%         try
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            ws.where=[where.data filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            
            % Identify fieldmaps-run matches
            ws.fm_run=[log.datalog{s+1,11} log.datalog{s+1,12} log.datalog{s+1,13} log.datalog{s+1,14} log.datalog{s+1,15} log.datalog{s+1,16}];
            ws.fm_nums=unique(ws.fm_run); ws.fm=cell(length(ws.fm_nums),1);
            for i=1:length(ws.fm_nums); ws.fm{i,1}=ws.fm_nums(i); ws.fm{i,2}=find(ws.fm_run(:)==ws.fm{i,1});end
            
            % For each run, select EPIs & vdm file (from corresponding fieldmap folder)
            matlabbatch{1}.spm.spatial.realignunwarp=settings.realignunwarp;
            for r=1:6
                wb.fieldmapnum=log.datalog{s+1,10+r}; % Identify fieldmap + fieldmap-session
                wb.sessionnum=find(ws.fm{ws.fm_nums(:)==wb.fieldmapnum,2}==r);
                f=spm_select('List', [ws.where 'Fieldmap_' num2str(wb.fieldmapnum) filesep], ['^vdm5_.*.session' num2str(wb.sessionnum) '*.img$']); % VDM file
                matlabbatch{1}.spm.spatial.realignunwarp.data(r).pmscan = { [ws.where 'Fieldmap_' num2str(wb.fieldmapnum) filesep f ',1']};
                f=spm_select('List', [ws.where 'Func_r' num2str(r)], '^bfM.*.img'); % Choose EPIs
                for i=1:size(f,1)
                    wb.r{i,1}= [ws.where 'Func_r' num2str(r) filesep f(i,:) ',1'];
                end
                matlabbatch{1}.spm.spatial.realignunwarp.data(r).scans =wb.r;
                %
                wb=[];
            end
            %
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            matlabbatch=[];
%         catch
%             errorlog{e,1}=['Failed: Realign & Unwarp  --- ' log.datalog{s+1,1}]; disp(errorlog{e,1}); e=e+1;
%         end
        ws=[];
    end
end

%% Step 2: Coregister (no prefix)

if process.coregister==1
    disp(' ############### (2) COREGISTRATION ############ ##########')
    for o2=1:1 % Settings for Coregistration
        settings.coreg.other = {''};
        settings.coreg.eoptions.cost_fun = 'nmi';
        settings.coreg.eoptions.sep = [4 2];
        settings.coreg.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        settings.coreg.eoptions.fwhm = [7 7];
    end
    for s=1:length(log.subjects)
        try
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.coreg.estimate=settings.coreg;  % Specifications 
            wb.where=[where.data filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            % Choose files
            f=spm_select('List', [wb.where 'MPM'], '^sM*.*T1w.nii'); % Reference: Structural
            matlabbatch{1}.spm.spatial.coreg.estimate.ref={[wb.where 'MPM' filesep f ',1']};
            f=spm_select('List', [wb.where 'Func_r1'], '^meanubfM*.*.img'); % Source image: meanubfM file (or, 1st volume from 1st run) 
            matlabbatch{1}.spm.spatial.coreg.estimate.source={[wb.where 'Func_r1' filesep f ',1']};
            % Collect functionals from all runs
            for r=1:6 
                f=spm_select('List', [wb.where 'Func_r' num2str(r)], '^ubfM*.*img');
                wb.r{r}=[repmat([wb.where 'Func_r' num2str(r) filesep],size(f,1),1) f repmat(',1',size(f,1),1)];
            end
            matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(vertcat(wb.r{1},wb.r{2},wb.r{3},wb.r{4},wb.r{5},wb.r{6}));
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Coregister & Reslice --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% Step 3: Segment (no prefix)

if process.segment==1 
    disp(' ############### (3) Normalization: Segmentation  ############ ##########')
    for o1=1:1 % Settings for Segmentation (performed only on T1)
        settings.segment.output.GM = [0 0 1];
        settings.segment.output.WM = [0 0 1];
        settings.segment.output.CSF = [0 0 0];
        settings.segment.output.biascor = 1;
        settings.segment.output.cleanup = 0;
        settings.segment.opts.tpm = {[where.spm filesep 'tpm' filesep 'grey.nii'];[where.spm filesep 'tpm' filesep 'white.nii'];[where.spm filesep 'tpm' filesep 'csf.nii']};
        settings.segment.opts.ngaus = [2;2;2; 4];
        settings.segment.opts.regtype = 'mni';
        settings.segment.opts.warpreg = 1;
        settings.segment.opts.warpco = 25;
        settings.segment.opts.biasreg = 0.0001;
        settings.segment.opts.biasfwhm = 60;
        settings.segment.opts.samp = 3;
        settings.segment.opts.msk = {''};
    end
    for s=1:length(log.subjects)
        try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.preproc=settings.segment;
            f=spm_select('List',  [where.data filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep 'MPM' ], '^sM*.*_T1w.nii');
            matlabbatch{1}.spm.spatial.preproc.data = {[where.data filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep 'MPM' filesep f ',1' ]};
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Segmentation of T1 --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%% Step 4: Normalization (prefix 'w')

if process.normalize==1  
    disp(' ############### (3) Normalization: Execute normalization ############ ##########')
    for o1=1:1 % Settings for Normalization 
        settings.normalization.preserve = 0;
        settings.normalization.bb = [-78 -112 -50;78 76 85];
        settings.normalization.vox = [2 2 2];
        settings.normalization.interp = 1;
        settings.normalization.wrap = [0 0 0];
        settings.normalization.prefix = 'w';
    end
    for s=1:length(log.subjects)
        try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.normalise.write.roptions=settings.normalization;
            wb.where=[where.data filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            f=spm_select('List', [wb.where 'MPM'], '^sM*.*_T1w_seg_sn.mat'); % Parameter file
            matlabbatch{1}.spm.spatial.normalise.write.subj.matname = {[wb.where 'MPM' filesep f]}; 
            % Collect functionals from all runs
            for r=1:6 
                f=spm_select('List', [wb.where 'Func_r' num2str(r)], '^ubfM*.*img');
                wb.r{r}=[repmat([wb.where 'Func_r' num2str(r) filesep],size(f,1),1) f repmat(',1',size(f,1),1)];
            end
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample=cellstr(vertcat(wb.r{1},wb.r{2},wb.r{3},wb.r{4},wb.r{5},wb.r{6}));
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Normalization  --- ' log.datalog{s+1,1}]; e=e+1;
        end
    end
end

%% Step 5: Smoothing
% Always saved in different folders



if process.smooth==1  
    disp([' ############### (4) SMOOTHING (size: ' num2str(process.smoothingsize) 'x' num2str(process.smoothingsize) 'x' num2str(process.smoothingsize) ')   ############ ##########'])
    input('Check code! Are we smoothing normalized or un-normalized functional data?');
    for o1=1:1 % Settings for Smoothing
        settings.smooth.fwhm = [process.smoothingsize process.smoothingsize process.smoothingsize]; % Smoothing size
        settings.smooth.dtype = 0;
        settings.smooth.im = 0;
%         settings.smooth.prefix='s';
        settings.smooth.prefix = ['s' num2str(process.smoothingsize)];
    end
    for s=1:length(log.subjects)
        try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.smooth=settings.smooth;
            wb.where=[where.data filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            % Collect functionals from all runs
            for r=1:6 
%                 f=spm_select('List', [wb.where 'Func_r' num2str(r)], '^wubfM*.*img');
                f=spm_select('List', [wb.where 'Func_r' num2str(r)], '^ubfM*.*img'); % for ants
                wb.r{r}=[repmat([wb.where 'Func_r' num2str(r) filesep],size(f,1),1) f repmat(',1',size(f,1),1)];
            end
            matlabbatch{1}.spm.spatial.smooth.data=cellstr(vertcat(wb.r{1},wb.r{2},wb.r{3},wb.r{4},wb.r{5},wb.r{6}));
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            %
            matlabbatch=[];
            wb=[];
        catch
            errorlog{e,1}=['Failed: Smoothing (size: ' num2str(process.smoothingsize) ' --- ' log.datalog{s+1,1}];
            e=e+1;
        end
    end
end

%%

if process.setup_newfol_4FL
    request.newfol.func_prefix='s3ubf';
    request.newfol.FL_type=' s3Ants';
    
    disp('Setting up new FL folder for:'); disp(request.newfol); input('Proceed?   ');
    for s=1:log.n_subjs
        disp(log.subjects{s})
        ws.wherefunc=['F:\2 Explore fMRI\1 Brain data Preproc PPIold\' log.subjects{s} '\1 Preprocessed\Func_r'];
        ws.toFL=['C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\'    log.subjects{s} filesep '2 First level' request.newfol.FL_type filesep];
        ws.tofunc=[ws.toFL 'Preproc functionals' filesep];
        mkdir(ws.toFL); mkdir(ws.tofunc)
        
        % Fetch functional scans
        for r=1:6
            ws.from=[ws.wherefunc num2str(r) filesep];
            f=spm_select('List', ws.from, ['^' request.newfol.func_prefix '.*']);
            for i=1:size(f,1)
                copyfile([ws.from f(i,:)], [ws.tofunc  f(i,:)])
            end
        end
        
        % Fetch onsets and other data
        ws.origFL=['C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\'    log.subjects{s} filesep '2 First level s6'  filesep];
        f=spm_select('List', ws.origFL, ['^' log.subjects{s} '.*']);
        for i=1:size(f,1)
            copyfile([ws.origFL filesep f(i,:)], [ws.toFL filesep f(i,:)])
        end
    end
end


%% END

disp('=======================================================')
w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp('Analysis completed:')
disp(process); disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' '); disp(errorlog); disp(' ')
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s3_Preprocessing)'), ' ',1);
catch
end