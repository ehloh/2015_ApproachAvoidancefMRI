function [ matlabbatch contrastfiles] = RL_2xN_fakefactorial(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch secondlevelfolder contrastfiles] = RL_2xN_fakefactorial(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% (1) Feed in all contrasts via a fake factorial (Task x RLvariable)
% (2) Remove/Conserve fake factorial's contrasts
% (3) Specify new contrasts that make sense
%
% Sensible contrasts to specify: (i) identity matrix (ii) Combined across tasks 
%                                              (ii) > in Conflict (iv) > in Control
% ------------------------------------------------------------------------------------------------------------

% Execute following to use as script:   where_brain=where.data_brain; subjectlist=log.subjects; firstlevelmodel=log.firstlevelmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% % Set up folder 
secondlevelfolder=[where.resultsfolder filesep secondlevelmodel];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); mkdir([secondlevelfolder filesep 'ROI']); end

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

% Available contrasts
w=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
spm=w.SPM;  
contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% Assign RLvariables
%       Col 1=RL variable, Col 2=Contrast file for Conflict task, Col 3=Contrast file for Control task
N=length(RLvariables);
prefix{1}='cF_'; prefix{2}='ct_'; RLvector=cell(size(RLvariables,1)*2,3);
for k=1:2
    for i=1:length(RLvariables)
        a.which=strfind(contrasts,[prefix{k} RLvariables{i,1}]);
        
        for j=1:size(a.which,1)
            if isempty(a.which{j,1})==0
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
        
        RLvector{(k-1)*length(RLvariables)+i,1}=[prefix{k} RLvariables{i}];
        RLvector{(k-1)*length(RLvariables)+i,2}=[k i];
        RLvector{(k-1)*length(RLvariables)+i,3}=contrasts{find(a.namewhich),2};
    end
end

% Assign Choices
for i=1:size(choices,1)
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
        choices{i,3}=contrasts{find(a.namewhich),2};
end

% Variables to read from this model
% Variables=vertcat(choices,RLvector);
                                                                                 
%% (3) Specify model for Factorial analysis

% 2x3 ANOVA: Task x Choice (Accept/Reject/Explore, NoBomb/Bomb/Explore)
matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'RLvar';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = N;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Specify Contrast files + Factorial cell
contrastfiles=cell(size(RLvector,1),1);
for i=1:size(RLvector,1)
        % Assign level
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels = RLvector{i,2};
        
        % Assign contrast files
        contrastfiles{i}=cell(length(subjectlist),1);
        for s=1:length(subjectlist)
            contrastfiles{i}{s}=[where.data_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep RLvector{i,3} ',1'];
        end
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans =contrastfiles{i};

end
% 
% % Include subjects as covariates
% for s=1:length(subjectlist)
%     sub=zeros(1,length(subjectlist)); sub(s)=1;
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(choices,1)])';
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
% end

% [EXECUTE] Specify  ##############################
disp('Specifying model (fake factorial) --------------------------')
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

%  [EXECUTE] Estimate model ##############################
disp('Estimating model (fake factorial) ------------------------------')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];


%% (4) Re-specify sensible contrasts
% Sensible contrasts to specify: (i) identity matrix (ii) Combined within tasks
%                                              (iii) RL variables, acoss tasks
%                                              (ii) > in Conflict (iv) > in Control

% Read details about this model (available conditions for weighting, RLvariables included)
ConNames=cell(size(RLvector,1),1);cRLvector=cell(size(RLvector,1),1);  cRLvariables=cell(size(RLvariables,1),1);  nCon=size(ConNames,1);
for i=1:size(RLvector,1); ConNames{i}=RLvector{i,1}; cRLvector{i}=RLvector{i,1}; end
for i=1:size(RLvariables,1); cRLvariables{i}=RLvariables{i,1}; end

for o1=1:1 % Compile contrast weights (cInstruc, 2nd col) for these contrasts  
    disp('Re-specifying contrasts that make sense ####################')
    c=1; cInstruc=cell(size(cRLvector,1),2); z=zeros(1,nCon);
    
    % Identity matrix
    for i=1:size(cRLvector,1) 
        cInstruc{c,1}=cRLvector{i};
        cInstruc{c,2}=z; cInstruc{c,2}(i)=1;
        c=c+1;
    end
    
    % Task, ignoring RL variables
    cInstruc{c,1}='Conflict_task';
    cInstruc{c,2}=z; cInstruc{c,2}(1:length(cRLvariables))=1; c=c+1;
    cInstruc{c,1}='Control_task';
    cInstruc{c,2}=z; cInstruc{c,2}(length(cRLvariables)+1:end)=1; c=c+1;
    
    % For each RL variable: Combined Conflict & Control
    for i=1:length(cRLvariables)
        cInstruc{c,1}=cRLvariables{i};
        cInstruc{c,2}=z;
        cInstruc{c,2}([i i+length(cRLvariables)])=1;
        c=c+1;
    end
    
    % For each RL variable: Conflict - Control
    for i=1:length(cRLvariables)
        cInstruc{c,1}=['cF-ct_' cRLvariables{i}];
        cInstruc{c,2}=z;
        cInstruc{c,2}(i)=1;
        cInstruc{c,2}(i+length(cRLvariables))=-1;
        c=c+1;
    end
    
    % For each RL variable: Control - Conflict
    for i=1:length(cRLvariables)
        cInstruc{c,1}=['ct-cF_' cRLvariables{i}];
        cInstruc{c,2}=z;
        cInstruc{c,2}(i)=-1;
        cInstruc{c,2}(i+length(cRLvariables))=1;
        c=c+1;
    end
    
    % Display
    disp('Conditions for weighting: '); disp(ConNames); disp('+ subject covariates');
    for i=1:size(cInstruc,1); disp(cInstruc{i,1});disp(cInstruc{i,2}); end
    
end

% Set up requested contrasts in Matlab
disp('Running contrasts that make sense ####################')
matlabbatch{1}.spm.stats.con.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 0; % Delete contrasts from the fake factorial?
if matlabbatch{1}.spm.stats.con.delete==1; disp('(Deleting fake-factorial''s contrasts)'); end
for i=1:size(cInstruc,1) % Specify contrasts according to instructions 
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.name =  cInstruc{i,1};
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.convec = cInstruc{i,2}; % weights
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
end

% Execute
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% Save details of this model + contrasts
details.RLvariables=cRLvariables;
details.DesignRegCols=ConNames;
details.ReqContrasts=cInstruc;
save([secondlevelfolder filesep 'details_2ndlevel.mat'], 'details', 'log', 'contrastfiles'); % Save details in 2nd level folder

end

