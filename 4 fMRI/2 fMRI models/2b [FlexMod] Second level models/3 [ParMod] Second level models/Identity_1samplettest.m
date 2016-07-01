function [ matlabbatch contrastfiles] = Identity_1samplettest(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch secondlevelfolder contrastfiles] = Identity_1samplettest(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)

%  % Execute following to use as script:   where.data_brain=where.data_brain; subjectlist=log.subjects; firstlevelmodel=log.firstlevelmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% % Set up folder 
% cd(where.data_brain); cd .. ; here=pwd; 
% secondlevelfolder=[here filesep '2 Second level results' filesep firstlevelmodel];
% if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); end

% Details of regressors/contrasts, for this particular first-level model
% prefix={'cF_' 'ct_'};
prefix={'cF_'}; disp('Running cF only!!!!')
% prefix={'ct_'}; disp('Running ct only!!!!')

% % (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
spm=w.SPM;  

contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% RLvariables: Col 1=RL variable, Col 2=Contrast file for Conflict task, Col 3=Contrast file for Control task
RLvector=cell(size(RLvariables,1)*length(prefix),2);
for i=1:length(RLvariables) 
    for k=1:length(prefix)
        
        a.which=strfind(contrasts,[prefix{k} RLvariables{i}]);
        
        for j=1:size(a.which,1)
            if isempty(a.which{j,1})==0 && strcmp(contrasts{j,1}, [prefix{k} RLvariables{i}])==1
                a.namewhich(j,1)=1;
            else
                a.namewhich(j,1)=0;
            end
        end
        if sum(a.namewhich)>1
            error(['ERROR: More than 1 contrast file found for RL variable ' prefix{k} RLvariables{i}])
        elseif sum(a.namewhich)<1
            error(['ERROR: Could not find contrast file for RL variable ' prefix{k} RLvariables{i}])
        end
        RLvariables{i,k+1}=contrasts{find(a.namewhich),2};
        
        RLvector{(i-1)*length(prefix)+k,1}=[prefix{k} RLvariables{i}];
        RLvector{(i-1)*length(prefix)+k,2}=contrasts{find(a.namewhich),2};
    end
end

% Choices
if isempty(choices)==1; choices=cell(0,2); end
for i=1:length(choices)
        a.which=strfind(contrasts, choices{i});
        
        for j=1:size(a.which,1)
            if isempty(a.which{j,1})==0
                a.namewhich(j,1)=1;
            else
                a.namewhich(j,1)=0;
            end
        end
        if sum(a.namewhich)>1
            error(['ERROR: More than 1 contrast file found for RL variable ' choices{i}])
        elseif sum(a.namewhich)<1
            error(['ERROR: Could not find contrast file for RL variable ' choices{i}])
        end
        choices{i,2}=contrasts{find(a.namewhich),2};
end

% Variables to read from this model
Variables=vertcat(choices,RLvector);

%% (3) Specify and estimate models for all Ttests

contrastfiles=cell(length(Variables),2);
for a=1:size(Variables,1)
        
        % (1) Specify model ##############################
        disp(['Specifying model no. '  num2str(a) '  :  ' Variables{a,1} '---------------------------'])
        testfolder=[where.resultsfolder filesep 'i ' Variables{a,1}];
        if isdir(testfolder)==0; mkdir(testfolder); mkdir([testfolder filesep 'ROI']); end

        matlabbatch{1}.spm.stats.factorial_design.dir = {testfolder};
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % Fetch subjects' scans/contrast files
        contrastfiles{a,1}=Variables{a,1};
        contrastfiles{a,2}=cell(length(subjectlist),1);
        for s=1:length(subjectlist)
            contrastfiles{a,2}{s}= [where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread  filesep firstlevelmodel ' Contrasted' filesep Variables{a,2} ',1'];
        end
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans=contrastfiles{a,2};

        % Execute (Specify only)
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];

        % (2) Estimate model ##############################
        disp(['Estimating model no. '  num2str(a) '  :  ' Variables{a,1} '---------------------------'])
        
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[testfolder filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];
        
        
        % Save
        details.RLvariables=RLvariables;
        save([testfolder filesep 'details_2ndlevel.mat'], 'details', 'log','contrastfiles'); % Save details in 2nd level folder

end

% In 2nd level folder as well
details.RLvariables=RLvariables;
save([where.resultsfolder  filesep 'details_2ndlevel.mat'], 'details', 'log'); % Save details in 2nd level folder

%% (4) Specify the contrast within each model 

disp('Specifying contrasts within each model (only 1 available) ###############')

for  a=1:size(Variables,1)
       
    % Specify the contrast (only + and - available)
    matlabbatch{1}.spm.stats.con.spmmat = {[where.resultsfolder filesep 'i ' Variables{a,1} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0; 
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name =  ['P_' Variables{a,1}]; % Positive contrast
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name =  ['N_' Variables{a,1}]; % Negative contrast
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = -1;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    % Execute
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
end

end

