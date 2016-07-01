% Get spike regressors + combine with movement parameters + Prep for 1st level
clear all;close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_brain='/Volumes/SANDISK/2 EXPLORE Brain data';
where.where='I:\5 Explore fMRI'; where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data'; where.data_beh=[where.where filesep 'Behavioural data'];

% Requested analysis
log.specificsubjects={'p01_GV'}; % BLANK to process all subjects 
request.getspike=0;
request.format_movementNphysio_regressors=1;
%
request.spikelogtoload=[]; % '(16-May-2013) Spikelog';

for o1=1:1 % General settings and specifications
    
    % Load subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');

    % Add paths
    [spm_path,name]=fileparts(which('spm'));
    physiopath=sprintf('%s%s%s',spm_path,filesep,'toolbox',filesep,'physio');
    addpath(physiopath);
    sonpath=sprintf('%s%s%s%s%s',spm_path,filesep,'toolbox',filesep,'physio',filesep,'son');
    addpath(sonpath);
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(['Data location (behaviour): ' where.data_beh]); disp(' ')
    disp('Requested analysis:'); disp(request); disp(' ')
    input('Hit enter to start                 ');
    disp('=======================================================')
    
end

%% STEP 1: Get Spike ------------------

if request.getspike==1
    
    disp('############## GET SPIKE ############## ')
    spikelog=cell(log.n_subjs,7);
    for o1=1:1 % Scanning details 
        w.scan=load([where.where filesep 'i_scanningdetails.mat']);
        nslices=w.scan.nSlicesPerVol;  % Number of slices in volume
        ndummies=w.scan.nDummyVols;  % Number of scans excluded from analysis
        TRms=w.scan.TRms;     % TR in ms
        TR=TRms/1000;       % Slice TR in secs
        nsessions=1; % Number of scanning sessions in the file
        slicenum=nslices/2;   % Slice number to time-lock regressors to
        % The above slice number can be determined from
        % data converted to nifti format. By default, slices
        % will be numbered from bottom to top but the acquisition
        % order can be ascending, descending or interleaved.
        % If slice order is descending or interleaved, the slice number
        % must be adjusted to represent the time at which the slice of
        % interest was acquired:
        %
        % For 3D acquisition sequences: Choose centre slice
        sliceorder='ascending'; % Ascending by default (verify with physics that this applies for this sequence: DONE)
        slicenum=get_slicenum(slicenum,nslices,sliceorder);
        
        % Channels
        channel.scanner=1;
        channel.cardialTTL=2;
        channel.cardiacQRS=[];
        channel.respiration=4;
        %
    end
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------------------'])
        wb.where=[where.data_beh filesep log.subjects{s} filesep 'Spike' filesep];
        spikelog{s,1}=log.subjects{s};
        
        if log.datalog{s+1,19}==1
            for b=1:6
                
                % Identify physio files
                f=spm_select('List', wb.where, [ '^' log.subjects{s} '.*.b' num2str(b) '.smr']);
                physiofile=[wb.where f];
                slicenum=get_slicenum(slicenum,nslices,sliceorder); % Need slice information
                if ~isempty(f)
                    try
                        disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')  -   b' num2str(b) '  ---------------------------------------------------'])
                        spikelog{s,b+1}=1;
                        % #######################################################
                        %
                        % From here, script content is entirely pasted from Physics' script 'physio_script.m'. See physics wiki for details
                        %
                        % The channel numbers must be assigned as they have been in spike.
                        % Unused channels should be set to empty using [];
                        % The channel numbers can be checked using the routines
                        % show_channels and check_channels as demonstrated below.
                        % Once the channels have been set correctly, they should
                        % stay the same when using spike with the same set-up and
                        % configuration file.
                        show_channels(physiofile);
                        scanner_channel=channel.scanner;
                        cardiacTTL_channel=channel.cardialTTL;
                        cardiacQRS_channel=channel.cardiacQRS;
                        resp_channel=channel.respiration;
                        check_channels(physiofile,scanner_channel,cardiacTTL_channel,cardiacQRS_channel,resp_channel);
                        
                        % Call the main routine for calculating physio regressors
                        % NB - currently the cardiacqrs calculation is disabled.
                        [cardiac,cardiacqrs,respire,rvt]=make_physio_regressors(physiofile,nslices,ndummies,TR,...
                            slicenum,nsessions,scanner_channel,cardiacTTL_channel,cardiacQRS_channel,resp_channel);
                        
                        % Save a record of parameters used for the regressors
                        %eval(['save ' spm_str_manip(physiofile,'r') '_physioparams physiofile nslices ndummies TRms slicenum nsessions']);
                        filename=[spm_str_manip(physiofile,'r'), '_physioparams'];
                        save(filename, 'physiofile', 'nslices', 'ndummies', 'TRms','slicenum','nsessions','sliceorder');
                        
                        % Save a record of parameters used for the regressors
                        %eval(['save ' spm_str_manip(physiofile,'r') '_physioparams physiofile nslices ndummies TRms slicenum nsessions']);
                        filename=[spm_str_manip(physiofile,'r'), '_physioparams'];
                        save(filename, 'physiofile', 'nslices', 'ndummies', 'TRms','slicenum','nsessions','sliceorder');
                        
                        % For each session, put regressors in a matrix called R.
                        % Each individual set of regressors are saved and also all regressors are saved with the name 'physiofile_R_session%d'.
                        % These files can be loaded into an SPM design matrix using the 'Multiple Regressors' option.
                        % NB motion parameters can also be concatenated with the physio regressors
                        % and saved as a set of regressors called R (see below for example)
                        
                        for sessnum=1:nsessions
                            R=[];
                            if ~isempty(cardiac{sessnum}) & ~isempty(cardiac{sessnum})
                                cardiac_sess = cardiac{sessnum};
                                filename = sprintf('%s_cardiac_session%d',spm_str_manip(physiofile,'r'),sessnum);
                                save(filename, 'cardiac_sess');
                                R=cat(2,R,cardiac{sessnum}(:,1:6));
                            end
                            if ~isempty(cardiacqrs{sessnum}) & ~isempty(cardiacqrs{sessnum})
                                cardiacqrs_sess = cardiacqrs{sessnum};
                                filename = sprintf('%s_cardiacqrs_session%d',spm_str_manip(physiofile,'r'),sessnum);
                                save(filename, 'cardiacqrs_sess');
                                R=cat(2,R,cardiacqrs{sessnum}(:,1:6));
                            end
                            if ~isempty(respire) & ~isempty(respire{sessnum})
                                respire_sess = respire{sessnum};
                                filename = sprintf('%s_respire_session%d',spm_str_manip(physiofile,'r'),sessnum);
                                save(filename, 'respire_sess');
                                R=cat(2,R,respire{sessnum}(:,1:6));
                            end
                            if ~isempty(rvt) & ~isempty(rvt{sessnum})
                                rvt_sess = rvt{sessnum};
                                filename = sprintf('%s_rvt_session%d',spm_str_manip(physiofile,'r'),sessnum);
                                save(filename,'rvt_sess');
                                R=cat(2,R,rvt{sessnum}(:,1:size(rvt{sessnum},2)));
                            end
                            nfiles=size(R,1);
                            % Save R for all physio only
                            if nfiles>0
                                oR=R;
                                Rname = sprintf('%s_R_session%d',spm_str_manip(physiofile,'r'),sessnum);
                                R=R-repmat(mean(R),nfiles,1);
                                save(Rname, 'R');
                            end
                            %    If required, also concatenate with motion params, e.g.
                            %  load rp_example.txt
                            %   RP=rp_example;
                            %   R=cat(2,oR,RP);
                            %   R=R-repmat(mean(R),nfiles,1);
                            %   Rname = sprintf('%s_R_session%d',spm_str_manip(physiofile,'r'),sessnum);
                            %   save(Rname, 'R');
                        end
                    catch
                        spikelog{s,b+1}=0;
                        errorlog{e}=['SPIKE FAILED: Could not process ' log.subjects{s} '  block ' num2str(b)]; disp(errorlog{e}); e=e+1;
                    end
                    
                else
                    disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')        -   b' num2str(b) '     [missing]   ------------------------------------------'])
                    spikelog{s,b+1}=0;
                end
                % --------------------- END OF PASTED SCRIPT --------------------------------
                % ############################################
            end
            spikelog{s,8}=floor((spikelog{s,2}+spikelog{s,3}+spikelog{s,4}+spikelog{s,5}+spikelog{s,6}+spikelog{s,7})/6);
        else
            % Exclude SPIKE for subjects that I know something's up
            spikelog{s,8}=0;
        end
    end
    % Save datalog to avoid having to run Spike again
%     save([where.data_beh filesep '(' date ') Spikelog'],'spikelog')
end

% Exclude Spike regressors (Typically disused!)
request.excludespikeregressors=1;
if  request.excludespikeregressors==1
    disp('############## EXCLUDING SPIKE REGRESSORS ############## ')
    for s=1:log.n_subjs
        spikelog{s,1}=log.subjects{s};
        for i=2:8
            spikelog{s,i}=0;
        end
    end
else
    load([where.data_beh filesep request.spikelogtoload]);
end


%% STEP 2: Format Spike & Movement parameters

if request.format_movementNphysio_regressors==1
    
    disp(' ############## FORMAT SPIKE & MOVEMENT ##############')
    physiomovereg=cell(log.n_subjs,   4); % Col 1: Subject, Col 2= Movement regressors + block regressors, Col 3=Spike regressors, Col 4: Concatenated everything
    printregressorsok=cell(log.n_subjs,2);
    if request.getspike~=1 && request.excludespikeregressors==0
        load([where.data_beh filesep request.spikelogtoload]);
    end
    for i=1:log.n_subjs
        s=find(strcmp(spikelog, log.subjects{i})==1);
        disp(['Subject ' num2str(s) '   (' spikelog{s,1} ') -------------------------------'])
        physiomovereg{s,1}=spikelog{s,1}; 
        
        % Load Movement parameters (from preprocessing) + Block ---------------------------------------------------------------
        for b=1:6 
             wb.where=[where.data_brain filesep spikelog{s,1} filesep '1 Preprocessed' filesep 'Func_r' num2str(b) filesep];
             f=spm_select('List', wb.where, '^rp.*txt$');
            wb.motion=load([wb.where f(1,:)]); if isempty(f)==1; errorlog{e}=(['ERROR: Could not find motion parameters for   ' spikelog{s,1} '  block ' num2str(b)]); disp(errorlog{e}); e=e+1; end
            wb.blocks=zeros(size(wb.motion,1),5); if b~=6; wb.blocks(:,b)=1; end % Append run regressor
            physiomovereg{s,2}{b}=[wb.motion wb.blocks];
        end
        
        % Load Respiration from SPIKE  ---------------------------------------------------------------------------------------
        if spikelog{s,8}==1 && request.excludespikeregressors ==0
            ws.includespike=1;
            wb.wbreathe=[where.data_beh filesep spikelog{s,1} filesep 'Spike' filesep];
            for b=1:6
                f=spm_select('List', wb.wbreathe, [ num2str(b) '_R_session1.mat$']);
                if ~isempty(f)
                    physiomovereg{s,3}{b}=load([wb.wbreathe f]);
                else
                    ws.includespike=0;
                    physiomovereg{s,3}{b}=[];
                    errorlog{e}=['NOTE: Respiration regressors not found   (' spikelog{s,1} '  run no. ' num2str(b) ')']; disp(errorlog{e}); e=e+1;
                end
            end
        else
            ws.includespike=0;
            for b=1:6
                physiomovereg{s,3}{b}=[];
            end
        end
        
        % Combine SPIKE & Motion, Combine blocks
        wb.allblocks_motion=vertcat(physiomovereg{s,2}{1},physiomovereg{s,2}{2},physiomovereg{s,2}{3},physiomovereg{s,2}{4},physiomovereg{s,2}{5},physiomovereg{s,2}{6});
        wb.allblocks_spike=vertcat(physiomovereg{s,3}{1},physiomovereg{s,3}{2},physiomovereg{s,3}{3},physiomovereg{s,3}{4},physiomovereg{s,3}{5},physiomovereg{s,3}{6});
        physiomovereg{s,4}=horzcat(wb.allblocks_spike, wb.allblocks_motion);
        
        % Check no. of regressors (14 respiration, 6 Motion, 5 blocks)
        if spikelog{s,8}==1 && size(physiomovereg{s,4},2)~=25 || spikelog{s,8}==0 && size(physiomovereg{s,4},2)~=11
            errorlog{e}=['ERROR: No. columns for regressor txt file is wrong     (' spikelog{s,1} ')']; disp(errorlog{e}); e=e+1; 
        end
        
        % Output to each subject's 1st level folder
        if isdir([where.data_brain filesep spikelog{s,1} filesep '2 First level'])==0; mkdir([where.data_brain filesep spikelog{s,1} filesep '2 First level']); end
        [w.printok]=print2txt([where.data_brain filesep spikelog{s,1} filesep '2 First level'], [spikelog{s,1} '_reg_physiomovement'], physiomovereg{s,4});
        printregressorsok{s,1}=w.printok.fileprint; printregressorsok{s,2}=w.printok.filemoved;

%             catch
%                 errorlog{e}=['ERROR: Could not format Spike+Motion regressors for  ' spikelog{s,1} '  block  ' num2str(b)]; disp(errorlog{e}); e=e+1;
%             end
    end
    
    % Save overall? To avoid re-processing SPIKE
    request.save_regressors=0;
    if request.save_regressors==1
        save([where.data_beh filesep '(' date ') PhysioMovement regressors.mat'], 'physiomovereg');
    end
    
end

%%  End

disp('------------------------------------------------------------------------')
disp('SPIKE regressors successfully constructed [block 1-6; Col 8: overall to include spike/not]')
if request.getspike==1; disp('SPIKELOG:'); disp(spikelog); end
disp('Error log:');  disp(errorlog)
disp('Note: Regressor files are saved with brain data')
disp('------------------------------------------------------------------------')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s5_GetspikeandMovementregressors)'), ' ',1);
end
