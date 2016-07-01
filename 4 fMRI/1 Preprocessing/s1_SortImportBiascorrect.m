% SortImportDeletedummiesBiascorrect
clear all;close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_brain='/Volumes/SANDISK/2 EXPLORE Brain data';
where.where='I:\5 Explore fMRI'; where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data'; % where.data_beh=[where.where filesep 'Behavioural data'];

% Requested analysis
process.sortimport=1;
process.deletedummyvolumes=1;
process.biascorrect=1;
%   
log.specificsubjects={'p39_TW'}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications
    
    % Load subjects
    log.w=load([where.data filesep 'datalogexplore_allsubsALL.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    scan=load([where.where filesep 'i_scanningdetails.mat']);
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('===================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:'); disp(process); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data]); disp(' ')
    input('Hit Enter to start      ')
    disp('====================================')
    
end

%% STEP 1: Sort and import together

if process.sortimport==1
    w.unzip_funcs=1; w.unzip_mpms=1;
    disp('############### STEP 1: SORT & IMPORT ###############')
    for s=1:log.n_subjs
        try
            disp(['Subject ' num2str(s) ' (' log.datalog{s+1} ')  --------- '])
            ws.from=[where.data filesep log.datalog{s+1} filesep '0 Original' filesep];
            ws.to=[where.data filesep log.datalog{s+1} filesep '1 Preprocessed' filesep];
            if isdir(ws.to)==0; mkdir(ws.to); end
            % Functional runs  -----
            if w.unzip_funcs==1
                for i=1:6
                    f=spm_select('List', ws.from, [log.datalog{s+1,4} '.' num2str(log.datalog{s+1,4+i}) '.tar$']);
                    Import_Archive([ws.from f(1,:)],ws.to);
                    ws.rename='java.io.File([ws.to log.datalog{s+1,4} ''.'' num2str(log.datalog{s+1,4+i})]).renameTo(java.io.File([ws.to ''Func_r'' num2str(i)]));';
                    eval(ws.rename)
                end
                % Field maps (indexed by scan #) -------
                ws.fms=unique([log.datalog{s+1,11} log.datalog{s+1,12} log.datalog{s+1,13} log.datalog{s+1,14} log.datalog{s+1,15} log.datalog{s+1,16}]);
                for i=1:length(ws.fms)
                    f=spm_select('List', ws.from, [log.datalog{s+1,4} '.' num2str(ws.fms(i)) '.tar']);
                    Import_Archive([ws.from f(1,:)],ws.to);
                    eval('java.io.File([ws.to log.datalog{s+1,4} ''.'' num2str(ws.fms(i))]).renameTo(java.io.File([ws.to ''Fieldmap_'' num2str(ws.fms(i))]));')
                    f=spm_select('List', ws.from, [log.datalog{s+1,4} '.' num2str(ws.fms(i)+1) '.tar']);  % 2nd scan
                    Import_Archive([ws.from f(1,:)],ws.to);
                    cd([ws.to log.datalog{s+1,4} '.' num2str(ws.fms(i)+1)])
                    ws.fmfiles=dir('sM*');
                    for j=1:size(ws.fmfiles,1)
                        movefile([ws.to log.datalog{s+1,4} '.' num2str(ws.fms(i)+1) filesep ws.fmfiles(j).name],[ws.to 'Fieldmap_' num2str(ws.fms(i)) filesep ws.fmfiles(j).name]);
                    end
                    cd(ws.to); rmdir([ws.to log.datalog{s+1,4} '.' num2str(ws.fms(i)+1)],'s')
                end
            end
            % MPMs ---------
            if w.unzip_mpms==1
                for i=1:5
                    for j=1:length(scan.MPM_scans{i,2})
                        f=spm_select('List', ws.from, [log.datalog{s+1,17} '.' num2str(log.datalog{s+1,18}-1+scan.MPM_scans{i,2}(j)) '.tar$']);
                        Import_Archive([ws.from f(1,:)],ws.to);
                    end
                    eval('java.io.File([ws.to log.datalog{s+1,17} ''.'' num2str(log.datalog{s+1,18}-1+ scan.MPM_scans{i,2}(1))]).renameTo(java.io.File([ws.to ''MPM_'' scan.MPM_scans{i,1}]));')
                    if i==2 % B1, 2 scan files
                        ws.b1_2=[ws.to log.datalog{s+1,17} '.' num2str(log.datalog{s+1,18}-1+ scan.MPM_scans{i,2}(2)) filesep];
                        cd(ws.b1_2)
                        ws.b1_2files=dir('sM*');
                        for j=1:size(ws.b1_2files,1) % Combine B0 files
                            movefile([ws.b1_2 ws.b1_2files(j).name],[ws.to 'MPM_b0' filesep ws.b1_2files(j).name]);
                        end
                        cd(ws.to); rmdir(ws.b1_2,'s')
                    end
                end
            end
            %
            ws=[];
        catch
            errorlog{e,1}=['Failed: Sort & Import --- ' log.datalog{s+1,1}];  disp(errorlog{e,1}); e=e+1;
        end
    end
end

%% STEP 2: Delete dummy volumes + adjust for Special cases

if process.deletedummyvolumes==1
    disp('############### STEP 2: DELETING DUMMIES ###############')
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) ' (' log.subjects{s} ')  --------- '])
        for b=1:6
            wb.where=[where.data filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_r' num2str(b) filesep];
            for j=1:scan.nDummyVols
                f=spm_select('List', wb.where, ['^fM.*00000' num2str(j) '-01']); if isempty(f)==1; erorrlog{e,1}=['Error. Could not find dummy to delete  - ' log.subjects{s} ' -  b' num2str(b) '  dummy # ' num2str(j)]; end; e=e+1;
                for i=1:size(f,1)
                    delete([wb.where f(i,:)]);
                end
            end
            wb=[];
        end
    end
end

% SPECIAL CASES?
% if sum(strcmp(log.subjects, 'p09_CN'))==1
%     errorlog{e}='NOTE: p09_CN ends at volume 224'; disp(errorlog{e}); e=e+1;
%     wb.where=[where.data filesep 'p09_CN' filesep '1 Preprocessed' filesep];
%     delete([wb.where 'Func_r1' filesep 'fMQ01089-0006-00225-000225-01.hdr'])
%     delete([wb.where 'Func_r1' filesep 'fMQ01089-0006-00225-000225-01.img'])
%     wb=[];
% end

%% STEP 3: Bias correction (Functionals only)

if process.biascorrect==1
    spm_jobman('initcfg')
    disp('############### STEP 3: BIAS CORRECTION ###############')
    for s=1:log.n_subjs
        try
            disp(['Subject ' num2str(s) ' (' log.subjects{s} ')  --------- '])
            ws.to=[where.data filesep log.datalog{s+1} filesep '1 Preprocessed' filesep];
            for i=1:6
                spm_biascorrect([ws.to filesep 'Func_r' num2str(i)]);
            end
        catch
            errorlog{e,1}=['Failed: Bias correction --- ' log.datalog{s+1,1}]; disp(errorlog{e,1}); e=e+1;
        end
    end
end

%% END

disp('====================================')
w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ')
disp('Analysis completed:')
disp(process)
disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' ')
disp('Error log?')
disp(errorlog)
disp(' ')
disp('====================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s1_SortImportBiascorrect)'), ' ',1);
catch
end
