function [ matlabbatch secondlevelfolder contrastfiles] = rl_2xN(where_brain,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
%  [ matlabbatch secondlevelfolder contrastfiles] = rl_2xN(where_brain,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)

% Execute following to use as script:   where_brain=where.data_brain; subjectlist=log.subjects; firstlevelmodel=log.firstlevelmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% Set up folder 
cd(where_brain); cd .. ; here=pwd; 
secondlevelfolder=[here filesep '2 Second level analysis' filesep ['results     ' firstlevelmodel '     ' secondlevelmodel] ];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); end

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

% Available contrasts
w=load([where_brain filesep subjectlist{1} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
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
                                                                                 
%% (2) Specify model for Factorial analysis

% 2x3 ANOVA: Task x Choice (Accept/Reject/Explore, NoBomb/Bomb/Explore)
matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'RLvar';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = N;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;% non-default
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
            contrastfiles{i}{s}=[where_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep RLvector{i,3} ',1'];
        end
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans =contrastfiles{i};

end

% Include subjects as covariates
for s=1:length(subjectlist)
    sub=zeros(1,length(subjectlist)); sub(s)=1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(choices,1)])';
   matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
end

%% Run model

% Execute (Specify only)
disp('Specifying model --------------------------')
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (2) Estimate model ##############################
disp('Estimating model ------------------------------')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[secondlevelfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];


end

