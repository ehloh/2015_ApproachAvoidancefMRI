function [ matlabbatch contrastfiles] = Identity_pairedttest(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch secondlevelfolder contrastfiles] = Identity_1samplettest(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)

%  % Execute following to use as script:   subjectlist=log.subjects; firstlevelmodel=log.onsetsmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% Instruc: Col 1=Name of comparison, Col 2&3: This minus That (names)
%               Col 4 & 5: This minus That (contrasts)
InstrucChoice={'cF-ct_Accept'     'cF_Accept'     'ct_NoBomb' [] [];
                       'cF-ct_Reject'     'cF_Reject'     'ct_Bomb' [] [];
                       'cF-ct_Explore'     'cF_Explore'     'ct_Explore' [] [];
                       'ct-cF_Accept'     'ct_NoBomb'     'cF_Accept' [] [];
                       'ct-cF_Reject'     'ct_Bomb'     'cF_Reject' [] [];
                       'ct-cF_Explore'     'ct_Explore'     'cF_Explore' [] [];};
c1=cell(size(RLvariables,1),5);c2=cell(size(RLvariables,1),5); % Compile instructions for RL variables
for i=1:size(RLvariables,1)
    c1{i,1}=['cF-ct_' RLvariables{i,1}];  % Conflict-Control
    c1{i,2}=['cF_' RLvariables{i,1}];
    c1{i,3}=['ct_' RLvariables{i,1}];
    c2{i,1}=['ct-cF_' RLvariables{i,1}];   % Control-Conflict
    c2{i,2}=['ct_' RLvariables{i,1}];
    c2{i,3}=['cF_' RLvariables{i,1}];
end
Instruc=vertcat(InstrucChoice, c1,c2);


% Other details of regressors/contrasts, for this particular first-level model
prefix{1}='cF_'; prefix{2}='ct_';

% % Set up folder 
% cd(where.data_brain); cd .. ; here=pwd; 
% secondlevelfolder=[here filesep '2 Second level results' filesep firstlevelmodel];
% if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); end % Should already be made

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

% What're we working with?
w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);spm=w.SPM;  
contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% Fill out details of instructions 
for i=1:size(Instruc,1)
    Instruc{i,4}=contrasts{find(strcmp(Instruc{i,2},contrasts)),2};
    Instruc{i,5}=contrasts{find(strcmp(Instruc{i,3},contrasts)),2};
end

%% (3) Specify and estimate models for all Ttests

contrastfiles=cell(size(Instruc,1),3);
for a=1:size(Instruc,1)
        
        % (1) Specify model ##############################
        disp(['Specifying model no. '  num2str(a) '  :  ' Instruc{a,1} '---------------------------'])
        testfolder=[where.resultsfolder filesep 'ip ' Instruc{a,1}];
        if isdir(testfolder)==0; mkdir(testfolder); mkdir([testfolder filesep 'ROI']);end

        matlabbatch{1}.spm.stats.factorial_design.dir = {testfolder};
        matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % Fetch subjects' scans/contrast files
        contrastfiles{a,1}=Instruc{a,1};
        contrastfiles{a,2}=cell(length(subjectlist),1);
        contrastfiles{a,3}=cell(length(subjectlist),1);
        for s=1:length(subjectlist)
            contrastfiles{a,2}{s}= [where.data_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep Instruc{a,4} ',1'];
            contrastfiles{a,3}{s}= [where.data_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep Instruc{a,5} ',1'];
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{1}=contrastfiles{a,2}{s};
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{2}=contrastfiles{a,3}{s};
        end
        
        % Execute (Specify only)
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];

        % (2) Estimate model ##############################
        disp(['Estimating model no. '  num2str(a) '  :  ' Instruc{a,1} '---------------------------'])
        
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

disp('Specifying contrasts witin each model (only 1 available) ###############')

for  a=1:size(Instruc,1)
       
    % Specify the contrast (only + and - available)
    matlabbatch{1}.spm.stats.con.spmmat = {[where.resultsfolder filesep 'ip ' Instruc{a,1} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0; 
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name =  ['P_' Instruc{a,1}]; % Positive contrast
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name =  ['N_' Instruc{a,1}]; % Negative contrast
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%     
%     % Overall, positive & negative
%     matlabbatch{1}.spm.stats.con.consess{2}.tcon.name =  [Instruc{a,1}]; % Negative contrast
%     matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
%     matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%     
%     
    
    
    % Execute
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
end

end

