function [ matlabbatch secondlevelfolder contrastfiles] = c_2xN(where_brain,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)
% [ matlabbatch secondlevelfolder contrastfiles] = c_2xN(where_brain,  subjectlist, firstlevelmodel, secondlevelmodel, RLvariables,choices)

% Execute following to use as script:   where_brain=where.data_brain; subjectlist=log.subjects; firstlevelmodel=log.firstlevelmodel; secondlevelmodel=log.secondlevelmodel; 

%% (1) Details for this model

% Set up folder 
cd(where_brain); cd .. ; here=pwd; 
secondlevelfolder=[here filesep '2 Second level analysis' filesep ['results     ' firstlevelmodel '     ' secondlevelmodel] ];
if isdir(secondlevelfolder)==0; mkdir(secondlevelfolder); end

% Details of regressors/contrasts, for this particular first-level model
prefix{1}='cF_'; prefix{2}='ct_';

%% (2) Which contrasts no. for the requested RLvariables? (Sample 1st subject

w=load([where_brain filesep subjectlist{1} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
spm=w.SPM;  

contrasts=cell(size(spm.xCon,2),2);
for i=1:size(spm.xCon,2); 
    contrasts{i,1}=spm.xCon(i).name;
    contrasts{i,2}=spm.xCon(i).Vcon.fname;
end

% RLvariables: Col 1=RL variable, Col 2=Contrast file for Conflict task, Col 3=Contrast file for Control task
RLvector=cell(size(RLvariables,1)*2,2);
for i=1:length(RLvariables) 
    for k=1:2
        
        a.which=strfind(contrasts,[prefix{k} RLvariables{i}]);
        
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
        
        RLvector{(i-1)*2+k,1}=[prefix{k} RLvariables{i}];
        RLvector{(i-1)*2+k,2}=contrasts{find(a.namewhich),2};
    end
end

% Choices
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

                                                                                 
%% (2) Specify model for Factorial analysis

matlabbatch{1}.spm.stats.factorial_design.dir = {secondlevelfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'RLvariable';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = size(RLvariables,1);
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1; 1;1];
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Specify contrast files + factorial cell

a=1; contrastfiles=[];
for t=1:2
    for v=1:size(RLvariables,1)
        
        % Specify which factorial cell
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels=[t; v];
        
        % Specify contrast files
        contrastfiles{a,1}=RLvariables{v,1};
        contrastfiles{a,2}=[t v];
        contrastfiles{a,3}=cell(length(subjectlist),1);
        for s=1:length(subjectlist)
            contrastfiles{a,3}{s}= [where_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep RLvariables{v,t+1} ',1'];
        end
        
        
        matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(a).scans=contrastfiles{a,3};
        a=a+1;
    end
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




% 
% 
% % Specify cells to Factorial Design
% for i=1:size(factorialcells,1)
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels=factorialcells{i,2};
% end

%% (3) Choose contrast files for all subjects
% 
% % Identify contrast files for each subjects    
% disp('Specifying contrast files for 2nd level ############# ')
% contrastfiles=cell(length(subjectlist),1);
% for s=1:length(subjectlist) 
%     disp(['Subject ' num2str(s) '  -  ' subjectlist{s}])
%     ws.where1st=[where_brain filesep subjectlist{s} filesep '2 First level' filesep firstlevelmodel ' Contrasted' filesep];
%     ws.s=load([ws.where1st 'SPM.mat']); ws.SPM=ws.s.SPM;% Load contrasted SPM file
%     
%     % Identify & Asisgn the correct contrasts for this cell
%     ws.cellcontrastnums=factorialcells; contrastfiles{s}=cell(4,1);
%     for i=1:length(factorialcells)
%         disp(['Identifying contrast no. ' num2str(i) '  -  ' factorialcells{i,1}])
%         if strcmp(memtype,'Roc')==1
%             [ws.contrastnames] = f_findcontrastname(ws.SPM, {[factorialcells{i,1} 'xCMem_' memtype] });
%         else
%             [ws.contrastnames] = f_findcontrastname(ws.SPM, {[factorialcells{i,1} 'xCMem_' memtype{1}] });
%         end
%         if isempty(ws.contrastnames)==1
%             error(['Error: Could not find the correct contrast file for specified cell   -- ' subjectlist{s} '  cell: ' factorialcells{i,1}])
%         end
%         contrastfiles{s}{i}= [ws.where1st ws.contrastnames{1} ,',1']; % Assign
%     end
%     %
%     ws=[];
% end
% disp('Specifying contrast files for 2nd level: DONE ############# ')
% 
% % Compile contrasts, between-subjects
% for i=1:4
%     for s=1:length(subjectlist)
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans{s,1}=contrastfiles{s}{i};
%     end
% end

end

