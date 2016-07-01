% Delete specific files/folders for all subjects 
clear all;close all hidden; clc

% Request specific
log.specificsubjects={}; % BLANK to process all subjects

for o1=1:1 % General settings and specifications
   
    % Load subjects
    where.where='D:\Dropbox\SANDISK\1 Explore fMRI'; 
    where.data_brain='G:\2 [Explore]\1 Brain data'; where.data_beh=[where.where filesep '1 Behavioural data'];
    addpath(where.where)
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all'); 
    
    % Interface
    disp('=======================================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp('Requested analysis: ADHOC SCRIPT')
    disp('See direct code for exact actions to be executed !! ')
    disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location (brain): ' where.data_brain])
    disp(['Data location (behaviour): ' where.data_beh])
    disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

whereto='C:\Users\e.loh\Documents\Neurodeb\1_cF\Backup BasicConFiles\';
wherefrom='C:\Users\e.loh\Documents\Neurodeb\1_cF\2b_AdjustCons\';

files=  {'_AF_Affine.txt';'_PRAF_Affine.txt';'_PRAF_InverseWarpxvec.nii.gz';'_PRAF_InverseWarpyvec.nii.gz';'_PRAF_InverseWarpzvec.nii.gz';'_PRAF_Warpxvec.nii.gz';'_PRAF_Warpyvec.nii.gz';'_PRAF_Warpzvec.nii.gz'}; 


%%
for s=1: log.n_subjs
    ws.c=clock;  disp(['Subject ' num2str(s) '   -  ' log.subjects{s} '   [' num2str(ws.c(4)) ':' num2str(ws.c(5)) ']  ------------------']);
%     ws.subfol=[where.data_brain filesep log.subjects{s} filesep]; 
%     ws.subfol=[ws.subfol '2 First level s4Ants' filesep];

%%



for f=1:size(files,1) 
    copyfile( [wherefrom log.subjects{s}  fs 'Basic' fs log.subjects{s}  files{f}], [whereto log.subjects{s} files{f}])
end





%%
%         eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
ws=[];
end



try % Notify researcher
    f_sendemail('kurzlich', strcat('DONE with ad hoc script'), ' ',1);
end


