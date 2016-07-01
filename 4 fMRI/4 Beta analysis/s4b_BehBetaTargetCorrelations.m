% Perform TARGETTED correlations between betas and behaviour.
clear all; close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data']; where.parameter_scripts='/Volumes/PENNYDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
where.where='I:\5 Explore fMRI'; where.data_beh=[where.where filesep '1 Behavioural data'];  where.parameter_scripts='I:\4 Explore experiment\3 Analysis\4 Fit computational models'; where.secondlevelresults='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results';

%
log.specificsubjects={};
log.behmodelparam_file='Behavioural stats modelpars (01-Nov-2013)'; % txt
log.behdecisionboundary_file='LD Decision boundary parameters (04-Nov-2013)'; % txt
log.personality_file='Personality scores'; % xlsx
log.behmisc_file='Misc behavioural scores'; % xlsx

% Beta file *** Edit here ***   
% compare='cF_Rej-Exp';
% % compare='Rej-Exp';
% seed='BA11';
% name=[compare  ' PPI seeded by ' seed ];

log.betatextfile_where='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c7_Cluster6CompeteFull_XUVPEN\choice_cluster2x2\ROI\2 Selected c7 rois';
log.betatextfile_name='(03-Nov-2013) Extracted betas';


behprofile=1; % 1 = Beh Inhibition, 2=Explore
name=[] % 'Behavioural inhibition';
% name='Explore';
log.ppi=0;

for o1=1:1 % General setup
    
    % Load subjects
    addpath(where.where);
    log.w=load([where.data_beh filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects_all log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
end


%% Load beta file & behaviour files (data in variable d, other details in variable log)
%             Data held in d_beh_table and d_beta_table. Row=subject, Col=Variable. 
%             Variable names included as headers in first row

for o1=1:1 % Load all them thangs. 
% (1) Load beta file
w.bet=importdata([log.betatextfile_where filesep log.betatextfile_name '.txt']);
d_beta_table=vertcat([{'Subject'} strtrim(w.bet.textdata(1, 2:end))], sortrows([strtrim(w.bet.textdata(2:end,1)) num2cell(w.bet.data)],1));

% (2) Load all behavioural files (append here if there are more datasets)
%         (a) Behavioural model parameters
%         (b) Personality scores
%         (c) Misc behavioural scores 
w.behmodpar=importdata([where.data_beh filesep 'Group behaviour files' filesep log.behmodelparam_file '.txt']);
d_rawbeh.behmodpar= [strtrim(w.behmodpar.textdata(1,:));  sortrows([strtrim(w.behmodpar.textdata(2:end,1)) num2cell(w.behmodpar.data)],1)];
[w.p1 w.p2 w.personality]=xlsread([where.data_beh filesep 'Group behaviour files' filesep log.personality_file '.xlsx']);
d_rawbeh.personality=[strtrim(w.personality(1,:));  sortrows([strtrim(w.personality(2:end,1)) w.personality(2:end, 2:end)],1)];
[w.m1 w.m2 w.miscbeh]=xlsread([where.data_beh filesep 'Group behaviour files' filesep log.behmisc_file '.xlsx']);
d_rawbeh.miscbeh=[strtrim(w.miscbeh(1,:));  sortrows([strtrim(w.miscbeh(2:end,1)) w.miscbeh(2:end, 2:end)],  1)]; % assume 1 dummy row
w.behdecpar=importdata([where.data_beh filesep 'Group behaviour files' filesep log.behdecisionboundary_file '.txt']);

d_rawbeh.behdecpar= [strtrim(w.behdecpar.textdata(1,:));  sortrows([strtrim(w.behdecpar.textdata(2:end,1)) num2cell(w.behdecpar.data)],1)];
for i=2:size(d_rawbeh.behmodpar,1) % Check subjects's order matching for all behavioural files. If error, change in the behaviour files.
    if  strcmp(d_rawbeh.behmodpar(i,1),  d_rawbeh.personality(i,1))~=1 || strcmp(d_rawbeh.behmodpar(i,1) , d_rawbeh.miscbeh(i,1)) ~=1 || strcmp(d_rawbeh.behmodpar(i,1) , d_rawbeh.behdecpar(i,1)) ~=1
        error('Subjects not consisten across the 3 behaviour profiles! Append raw data files.')
    end
end
d_beh_table=[d_rawbeh.personality d_rawbeh.behmodpar(:,2:end) d_rawbeh.miscbeh(:,2:end) d_rawbeh.behdecpar(:,2:end)];

% (3) Apply subject selection form beta table to behaviour
[log.subjects log.n_subjs d_beh_table]=f_selectsubjects(d_beh_table,  d_beta_table(2:end,1), [d_beh_table(:,1) [{'all'}; num2cell(ones(size(d_beh_table,1)-1,1))]], 'all');
log.subjects=d_beh_table(2:end,1);
d_beh_table(:,1)=[];
d_beta_table(:,1)=[]; % From here on, no more subject names!
disp('Original beh list:'); disp(d_beh_table(1,:)');
disp('Original bet list:'); disp(d_beta_table(1,:)');
end

%% Ad-hoc Transformations: New scores (for beta and behaviour
%       Note: another option is to change scores in excel/texts tables

% Which ROIs to subject to further analysis at all
instruc.ROIs={  % These are used to generate new scores & correlations themselves 
% 'BA11';'BA46';'Caudate_L_anat';'Caudate_R_anat';'Putamen_L_anat';'HPC_saC_L';'HPC_saC_R';'HPC_saTc_L';'HPC_saTc_R';  % Roi battery (PPIs)
'HPC_L_tc';'BA10_c';'BA46_EmR';'SFG_med_RmE';


% 'SFG_BA11';'SFG_BA46';'MFG_BA9';'OFC_Mid';'Med_FG'; 'Cingulate_Post1';'Cingulate_Post2'; 'Precuneus1';'Precuneus2';'PostcentralGyr';'Occipit';'Cereb'    % c3 ME Choice
%
% 'Caudate_L';%  'Caudate_R';
%     'HPC_L';'HPC_R'

};
instruc.roi_suffix=[]; % ['__' seed '-sd' ];
instruc.roi_prefix= 'in_';
instruc.ROIs=cellfun(@(x)[x instruc.roi_suffix], instruc.ROIs, 'UniformOutput',0);

disp('################# Calculating new scores ###########')
for o1=1:1 % Ad hoc BEHAVIOURAL scores  
    
    
    instruc.newbehscores={
        
        % [Behavioural inhibition scores] 
        
        'per.cF_NotAccept'                      'per.cF_Reject'                 '+'     'per.cF_Explore';
        'per.ct_NotAccept'                       'per.ct_Bomb'                 '+'     'per.ct_Explore';
        'per.NotAccept_cF-ct'                   'per.cF_NotAccept'           '-'     'per.ct_NotAccept';
        'per_i4.Reject_cF-ct'                   'per_i4.cF_Reject'              '-'     'per_i4.ct_Bomb';
        'per_i6.Reject_cF-ct'                   'per_i6.cF_Reject'              '-'     'per_i6.ct_Bomb';
        
        'dbao.ao_const_cF-ct'          'dbao.cF_ao_const'     '-'     'dbao.ct_ao_const';
        'dbao.ao_slope_cF-ct'               'dbao.cF_ao_slope'            '-'     'dbao.ct_ao_slope';
        'dbao.ao_bprod_cF-ct'       'dbao.cF_ao_bprod'     '-'     'dbao.ct_ao_bprod';
        'dbf.ar_bprod_cF-ct'                   'dbf.cF_ar_bprod'                '-'     'dbf.ct_ar_bprod';
        
        'rt.Reject_cF-ct'                        'rt.cF_Reject'                      '-'     'rt.ct_Bomb';
        'rt_i4.Reject_cF-ct'                   'rt_i4.cF_Reject'                  '-'     'rt_i4.ct_Bomb';
        'rt_i6.Reject_cF-ct'                   'rt_i6.cF_Reject'                  '-'     'rt_i6.ct_Bomb';
        
        };
    
    % ########## Execute ##########################
    k=size(d_beh_table,2)+1;
    for i=1:size(instruc.newbehscores,1) % By instruction
        
        % Find inputs
        if sum(strcmp(d_beh_table(1,:), instruc.newbehscores{i,2}))~=1; error(['[ad hoc beh score] could not find beh variable ' instruc.newbehscores{i,2}]); end
        if sum(strcmp(d_beh_table(1,:), instruc.newbehscores{i,4}))~=1; error(['[ad hoc beh score] could not find beh variable ' instruc.newbehscores{i,4}]); end
        wb.beh1=cell2mat(d_beh_table(2:end, find(strcmp(d_beh_table(1,:), instruc.newbehscores{i,2})   )));
        wb.beh2=cell2mat(d_beh_table(2:end, find(strcmp(d_beh_table(1,:), instruc.newbehscores{i,4})   )));
        
        % Apply mathematical transformation
        eval(['d_beh_table(2:end,k)= num2cell(  wb.beh1 '  instruc.newbehscores{i,3} ' wb.beh2       );'])
        %
        d_beh_table{1,k}=instruc.newbehscores{i,1};
        disp(['New beh variable: ' instruc.newbehscores{i,1} '     (' instruc.newbehscores{i,2} '  ' instruc.newbehscores{i,3}  '  ' instruc.newbehscores{i,4} ')'])
        wb=[]; k=k+1;
    end
    
    for o2=1:1 % Manually specified new beh scores
        
        % (a) Mean transformations ------------------------------------
        instruc.newbeh_mean={
%             'per.Explore'                 {'per.cF_Explore';'per.ct_Explore'};
            };
        
        k=size(d_beh_table,2)+1;
        for i=1:size(instruc.newbeh_mean,1) % via instruction
            
            % Find inputs
            wb.behs=nan*zeros(log.n_subjs,length(instruc.newbeh_mean{i,2}));
            for b=1:length(instruc.newbeh_mean{i,2})
                if sum(strcmp(d_beh_table(1,:), [instruc.newbeh_mean{i,2}{b}]) )~=1; error(['[ad hoc beh score (mean)] could not find beh variable ' instruc.newbeh_mean{i,2}{b}]); end
                wb.behs(:,b)=cell2mat(d_beh_table(2:end,find(strcmp(d_beh_table(1,:), [instruc.newbeh_mean{i,2}{b}]))));
            end
            
            % Write to beh table
            d_beh_table{1,k}=[ instruc.newbeh_mean{i,1}];
            d_beh_table(2:end,k)=num2cell(mean(wb.behs,2));
            
            
            % Display to interface
            disp(['New beh variable: '  instruc.newbeh_mean{i,1}  '    (mean of inputs, following)'])
            disp(instruc.newbeh_mean{i,2})
            
            %
            wb=[]; k=k+1;
        end
        
        % (b) Other manual transformations  ------------------------------------
        instruc.newbeh_2ndtransform={
            %             'Rej-Exp'        'Reject'     '-'     'Explore';       % Explore scores
            };
        
        k=size(d_beh_table,2)+1;
        for i=1:size(instruc.newbeh_2ndtransform,1) % via instruction
            % Find inputs
            if sum(strcmp(d_beh_table(1,:), [instruc.newbeh_2ndtransform{i,2}]    ))~=1; error(['[ad hoc beh score 2nd transform] could not find beh variable '  instruc.newbeh_2ndtransform{i,2}]   ); end
            if sum(strcmp(d_beh_table(1,:),  [instruc.newbeh_2ndtransform{i,4}]))~=1; error(['[ad hoc beh score 2nd transform] could not find beh variable ' instruc.newbeh_2ndtransform{i,4}]); end
            wb.beh1=cell2mat(d_beh_table(2:end, find(strcmp(d_beh_table(1,:), [instruc.ROIs{r} '-' instruc.newbeh_2ndtransform{i,2}])   )));
            wb.beh2=cell2mat(d_beh_table(2:end, find(strcmp(d_beh_table(1,:), [instruc.ROIs{r} '-' instruc.newbeh_2ndtransform{i,4}])   )));
            
            % Apply mathematical transformation
            eval(['d_beh_table(2:end,k)= num2cell(  wb.beh1 '  instruc.newbeh_2ndtransform{i,3} ' wb.beh2       );'])
            %
            d_beh_table{1,k}=instruc.newbeh_2ndtransform{i,1};
            disp(['New beh variable: ' instruc.newbeh_2ndtransform{i,1} '     (' instruc.newbeh_2ndtransform{i,2} '  ' instruc.newbeh_2ndtransform{i,3}  '  ' instruc.newbeh_2ndtransform{i,4} ')'])
            wb=[]; k=k+1;
        end
        
    end
    
    
    
end
for o1=1:1 % Ad hoc BETA scores  
    
    % Meta-instructions: to generate instructions (specify rois & contrasts)
   if log.ppi==1
       instruc.newbeta_transformperROI={};
   else
    instruc.newbeta_transformperROI={
        
        % [New betas for CHOICE models]
        'cF_Rej-Exp'        [instruc.roi_prefix 'cF_Reject']     '-'     [instruc.roi_prefix 'cF_Explore'];       % Behavioural inhibition scores
        'Reject_cF-ct'       [instruc.roi_prefix 'cF_Reject']     '-'     [instruc.roi_prefix 'ct_Bomb'];
        'ct_Bmb-Exp'        [instruc.roi_prefix 'ct_Bomb']     '-'     [instruc.roi_prefix 'ct_Explore'];
        %
        };
   end
    
    % Generate instructions from meta instructions. 
    % Eventual form: Col 1= new name, Col 2=Var 1, Col 3= mathematical transform, Col 4=Var 2
    instruc.newbetascores=cell(length(instruc.ROIs)*size(instruc.newbeta_transformperROI,1),4); k=1;
    for r=1:length(instruc.ROIs)
        for t=1:size(instruc.newbeta_transformperROI,1)
            instruc.newbetascores(k,1:4)=strtrim( cellstr([strtrim(char({instruc.ROIs{r}; instruc.ROIs{r}; ' '; instruc.ROIs{r}})) ['-';'-';' ';'-'] char(instruc.newbeta_transformperROI{t,:})])');
            k=k+1;
        end
    end
    
    % ########## Execute to calc new beta scores ##########################
    k=size(d_beta_table,2)+1;
    for i=1:size(instruc.newbetascores,1) % via instruction
        
        % Find inputs
        if sum(strcmp(d_beta_table(1,:), instruc.newbetascores{i,2}))~=1; error(['[ad hoc beta score] could not find beta variable ' instruc.newbetascores{i,2}]); end
        if sum(strcmp(d_beta_table(1,:), instruc.newbetascores{i,4}))~=1; error(['[ad hoc beta score] could not find beta variable ' instruc.newbetascores{i,4}]); end
        wb.beta1=cell2mat(d_beta_table(2:end, find(strcmp(d_beta_table(1,:), instruc.newbetascores{i,2})   )));
        wb.beta2=cell2mat(d_beta_table(2:end, find(strcmp(d_beta_table(1,:), instruc.newbetascores{i,4})   )));
        
        % Apply mathematical transformation
        eval(['d_beta_table(2:end,k)= num2cell(  wb.beta1 '  instruc.newbetascores{i,3} ' wb.beta2       );'])
        %
        d_beta_table{1,k}=instruc.newbetascores{i,1};
        disp(['New beta variable: ' instruc.newbetascores{i,1} '     (' instruc.newbetascores{i,2} '  ' instruc.newbetascores{i,3}  '  ' instruc.newbetascores{i,4} ')'])
        wb=[]; k=k+1;
    end
    
    
    
    for o2=1:1 % Manually specified new beta scores
        
        % (a) Mean transformations ------------------------------------
        if log.ppi==1
            instruc.newbeta_mean={};
        else
        instruc.newbeta_mean={
            'Reject'                  {[instruc.roi_prefix 'cF_Reject'];[instruc.roi_prefix 'ct_Bomb']};
            'Explore'                 {[instruc.roi_prefix 'cF_Explore'];[instruc.roi_prefix 'ct_Explore']};
            };
        end
        k=size(d_beta_table,2)+1;
        for r=1:length(instruc.ROIs)
            for i=1:size(instruc.newbeta_mean,1) % via instruction
                
                % Find inputs
                wb.betas=nan*zeros(log.n_subjs,length(instruc.newbeta_mean{i,2}));
                for b=1:length(instruc.newbeta_mean{i,2})
                    if sum(strcmp(d_beta_table(1,:), [instruc.ROIs{r} '-'  instruc.newbeta_mean{i,2}{b}]) )~=1; error(['[ad hoc beta score (mean)] could not find beta variable ' instruc.newbeta_mean{i,2}{b}]); end
                    wb.betas(:,b)=cell2mat(d_beta_table(2:end,find(strcmp(d_beta_table(1,:), [instruc.ROIs{r} '-'  instruc.newbeta_mean{i,2}{b}]))));
                end
                
                % Write to beta table
                d_beta_table{1,k}=[instruc.ROIs{r} '-' instruc.newbeta_mean{i,1}];
                d_beta_table(2:end,k)=num2cell(mean(wb.betas,2));
                
                
                % Display to interface
                disp(['New beta variable: '  instruc.ROIs{r} '-'  instruc.newbeta_mean{i,1}  '    (mean of inputs, following)'])
                disp(instruc.newbeta_mean{i,2})
                
                %
                wb=[]; k=k+1;
            end
        end
        
        % (b) Other manual transformations  ------------------------------------
      if log.ppi==1
          instruc.newbeta_2ndtransform={};
      else
      
        instruc.newbeta_2ndtransform={
            'Rej-Exp'        'Reject'     '-'     'Explore';       % Explore scores
            };
      end 
        k=size(d_beta_table,2)+1;
        for r=1:length(instruc.ROIs)
            for i=1:size(instruc.newbeta_2ndtransform,1) % via instruction
                % Find inputs
                if sum(strcmp(d_beta_table(1,:), [instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,2}]    ))~=1; error(['[ad hoc beta score 2nd transform] could not find beta variable ' instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,2}]   ); end
                if sum(strcmp(d_beta_table(1,:),  [instruc.ROIs{r} '-'  instruc.newbeta_2ndtransform{i,4}]))~=1; error(['[ad hoc beta score 2nd transform] could not find beta variable ' instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,4}]); end
                wb.beta1=cell2mat(d_beta_table(2:end, find(strcmp(d_beta_table(1,:), [instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,2}])   )));
                wb.beta2=cell2mat(d_beta_table(2:end, find(strcmp(d_beta_table(1,:), [instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,4}])   )));
                
                % Apply mathematical transformation
                eval(['d_beta_table(2:end,k)= num2cell(  wb.beta1 '  instruc.newbeta_2ndtransform{i,3} ' wb.beta2       );'])
                %
                d_beta_table{1,k}=[instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,1}];
                disp(['New beta variable: ' instruc.ROIs{r} '-' instruc.newbeta_2ndtransform{i,1} '     (' instruc.newbeta_2ndtransform{i,2} '  ' instruc.newbeta_2ndtransform{i,3}  '  ' instruc.newbeta_2ndtransform{i,4} ')'])
                wb=[]; k=k+1;
            end
            
        end
        
    end
    
    
    
    
    
end

% Consolidate data tables  (beh and beta). Retain headers with variable names
log.beh_list=d_beh_table(1,:)'; log.n_beh=length(log.beh_list);
log.beta_list=d_beta_table(1,:)'; log.n_beta=length(log.beta_list);
disp('################# Final list of scores (beh and beta) ###########')
disp('Beh list: '); disp(log.beh_list);
disp('Beta list: '); disp(log.beta_list);

%% Specify targetted correlations
% Instructions are compiled both manually and by meta-instruction (specifying ROI, contrast, behavioural score)

% (1) [Meta] Beta & ROI
request.metainstruc_rois=instruc.ROIs;
for o1=1:1
    if log.ppi==1
        request.metainstruc_contrasts={''};
        
    elseif behprofile==1 % Which contrasts?
        
        % [Behavioural inhibition]
        request.metainstruc_contrasts={
            [instruc.roi_prefix 'cF_Reject']; 'cF_Rej-Exp'; 'Reject_cF-ct'; [instruc.roi_prefix 'ct_Bomb']; 'ct_Bmb-Exp';
            };
        
    elseif behprofile==2
        
        % [Explore]
        request.metainstruc_contrasts={
            'Explore'; [instruc.roi_prefix 'cF_Explore']; [instruc.roi_prefix 'ct_Explore']; 'Rej-Exp';'cF_Rej-Exp'; 'ct_Bmb-Exp';
            };
    end
end

% (2) [Meta] Behaviour scores
if behprofile==1
    
    % [Behavioural Inhibition]
    request.metainstruc_beh={
        'st.State';'st.Trait';'B.Bis';'B.Bas';'B.Drive';'B.FS';'B.RR'; % Personality
        'per.cF_Reject';'per.cF_Explore'; 'per.cF_NotAccept'; 'per.NotAccept_cF-ct';   % Choice probability
        'per_i4.cF_Reject';'per_i6.cF_Reject'; 'per_i4.Reject_cF-ct';'per_i6.Reject_cF-ct';
        'dbao.cF_ao_const';'dbao.cF_ao_slope';'dbao.cF_ao_bprod';        % Linear discrim decision boundary
        'dbao.ao_const_cF-ct';'dbao.ao_slope_cF-ct';'dbao.ao_bprod_cF-ct'; 
        'dbf.cF_ar_bprod'; 'dbf.ar_bprod_cF-ct'; 
        'mcF_2fevx.f'; 'mcF_2feux.f';        % RL model
        'rt.cF_Reject';'rt_i4.cF_Reject';'rt_i6.cF_Reject';  % RTs
        'rt.Reject_cF-ct';'rt_i4.Reject_cF-ct'; 'rt_i6.Reject_cF-ct';
        };
    
elseif behprofile ==2
    
    % [Explore profile ]
    request.metainstruc_beh={
        'B.Bas';'B.Drive';'B.FS';'B.RR';          % Personality
        'per.cF_Explore';'per.ct_Explore';  'per.Explore';  % Choice probability
        'mcF_2fevx.e';'mcF_2fevx.x';  'mct_2evx.e'; 'mct_2evx.x';  'mct_2eux.e'; 'mct_2eux.x';          % RL model
        'dbf.cF_ae_bprod'; 'dbf.ct_ae_bprod'; 'dbf.cF_re_bprod'; 'dbf.ct_re_bprod';  % Linear discrim decision boundary
        };
    
end

% (3) Manual instructed correlations (Col 1= beta, col 2= Beh)
request.manual_correlations={
%     'Cereb_L-cF_Explore'             'Dummy';
%     ' '                                           'st.Trait';
};

% ###############################################

for o1=1:1 % Compile instructions for executing correlations
    request.compiled_correlations=[];
    request.compiled_corr_byroi=cell(length(request.metainstruc_rois)*length(request.metainstruc_contrasts),1);
    for r=1:length(request.metainstruc_rois)
        for c=1:length(request.metainstruc_contrasts)
            if log.ppi
                request.compiled_corr_byroi{r,c}=[cellstr(repmat([request.metainstruc_rois{r} request.metainstruc_contrasts{c}], length(request.metainstruc_beh),1) )  request.metainstruc_beh];
            else
            request.compiled_corr_byroi{r,c}=[cellstr(repmat([request.metainstruc_rois{r} '-' request.metainstruc_contrasts{c}], length(request.metainstruc_beh),1) )  request.metainstruc_beh];
            end
            request.compiled_correlations=vertcat(request.compiled_correlations, request.compiled_corr_byroi{r,c});
        end
    end
    
end

% Final list of requested correlations:
request.correlations=vertcat(request.compiled_correlations, request.manual_correlations);
% openvar r_corr

% openvar request.correlations

%% Execute correlations; results in r_corr 
%     Col 1=Beta, Col 2=Beh, Col 3= corr stats, Col 4 =mark significance

disp('################# Executing requested correlations ###########')
r_corr=cell(size(request.correlations,1),3); r_corrT=cell(size(request.correlations,1),1);
d_targetedbeh=cell(log.n_subjs+1, size(request.correlations,1)); 
d_targetedbeta=cell(log.n_subjs+1, size(request.correlations,1)); 
k=1; j=1; 
for c=1:size(request.correlations,1)
    
    % Add a line between roi x contrasts
    if c>1; if strcmp(request.correlations{c,1},request.correlations{c-1,1})==0; k=k+1; end; end
    
    r_corr{k,1}=request.correlations{c,1};
    r_corr{k,2}=request.correlations{c,2};
    if strcmp(request.correlations{c,1}, ' ')==1; r_corr{k,1}=r_corr{k-1,1}; end
   
    
    
    % Load beta variable
    if sum(strcmp(log.beta_list, r_corr{k,1}))==1
        wc.beta=cell2mat(d_beta_table(2:end, find(strcmp(log.beta_list, r_corr{k,1}))));
    else error(['Correlation ' num2str(c) ': cannot find beta variable ' r_corr{k,1}]);
    end
    
    % Load behavioural variable
    if sum(strcmp(log.beh_list, r_corr{k,2}))==1
        wc.beh=cell2mat(d_beh_table(2:end, find(strcmp(log.beh_list, r_corr{k,2}))  ));
    else error(['Correlation ' num2str(c) ': cannot find beh variable ' r_corr{k,2}]);
    end
    
    % Run correlation & record results
    [wc.r wc.p]=corr(wc.beh, wc.beta);
    if wc.p<0.001
        r_corr{k,3}=strtrim(['r= ' num2str(wc.r,3) ', p=' num2str(wc.p, 4) ' ***']);
    elseif wc.p<0.01
        r_corr{k,3}=strtrim(['r= ' num2str(wc.r,3) ', p=' num2str(wc.p, 3) '  **']);
    elseif wc.p<0.05
        r_corr{k,3}=strtrim(['r= ' num2str(wc.r,3) ', p=' num2str(wc.p, 2) '   *']);
    elseif wc.p<0.1
        r_corr{k,3}=strtrim(['r= ' num2str(wc.r,3) ', p=' num2str(wc.p, 2)]);
    else
        r_corr{k,3}=[];
    end
    
    
    % 1-cell text
    if isempty(r_corr{k,3})==1
        r_corrT{k,1}=[r_corr{k,1} ' - ' r_corr{k,2}];
    else
        r_corrT{k,1}=[r_corr{k,1} ' - ' r_corr{k,2} '    ' r_corr{k,3}];
    end
    
    
%     Record Behavioural variable
    d_targetedbeh{1,j}=request.correlations{c,2};
    d_targetedbeh(2:end,j)=num2cell(wc.beh);
    d_targetedbeta{i,j}=request.correlations{c,1};
    d_targetedbeta(2:end,j)=num2cell(wc.beta);
    
    
    openvar d_targetedbeh
    
    k=k+1;j=j+1; 
    
    
    wc=[];
end   

if behprofile==1
    corrtype='BEHAVIOURAL INHIBITION';
else
    corrtype='EXPLORE';
end

% Export (append details of beta file)
% dat=[{['Beta file name: ' log.betatextfile_name] ' ' ' ';['Beta file loc: ' log.betatextfile_where] ' ' ' '; ' ' ' ' ' '}; r_corr];
dat=[{['Beta file name: ' log.betatextfile_name] ;['Beta file loc: ' log.betatextfile_where] ;' ';[corrtype ' correlations']; ' ';}; r_corrT];
openvar dat
% openvar d_targetedbeh


print=0;
if print
    input('Continue to print? ');
    printok=print2txt(log.betatextfile_where, ['(' date ') Targeted beta correlations - ' name], dat);
    disp(printok)
end

%% Plot?

log.targetedbeh=d_targetedbeh(1,:)';
log.targetedbeta=d_targetedbeta(1,:)';
% openvar log.targetedbeh; openvar log.targetedbeta


% instruc.plot={'per.Explore'         'Caudate_R-Explore'}; % Behaviour - Beta
instruc.plot={'per.cF_Reject'   'SFG_BA46-ct_Explore'};
% instruc.plot={'mcF_2fevx.e'         'SFG_BA46-cF_Explore'};
newbeh={'dbao.cF_ao_const'};


% Execute followng
close all hidden
p.beh=cell2mat(d_targetedbeh(2:end, find(strcmp(log.targetedbeh,  instruc.plot{1}), 1,'first')));
p.beh2=cell2mat(d_targetedbeh(2:end, find(strcmp(log.targetedbeh,  newbeh{1}), 1,'first')));
p.beta=cell2mat(d_targetedbeta(2:end, find(strcmp(log.targetedbeta,  instruc.plot{2}), 1,'first')));
%


% 
% scatter(p.beh,  p.beh2); hold on;
% title(['x= ' instruc.plot{1} ',  y= beta '  instruc.plot{2}])


% p.mdl=fitlm(p.beh,p.beta);
% plot(p.mdl)


