function [ con] = f_firstlevelcontrasts_standard(where,log)
%  [ con] = f_firstlevelcontrasts_standard(where,log)
%  Bog-standard first-level contrasts - read details off excel table 'i_firstlevel_contrasts'
%   Different excel sheet depending on which model (specified in this  function)
%
% ----------------------------------------------------------------------------------------


%% Which contrast table (from the excel sheet) to use?

for o1=1:1 % Load default specifications for requested model
    
    % Which contrast table to use?
    if strcmp(log.onsetsmodel(1:3), 'm_c')==1  % Competing models
        contrasts.contrasttable='Compete';
    elseif strcmp(log.onsetsmodel(1:3), 'm_o')==1 %Orthogonalized models
        contrasts.contrasttable='Orthog';
    elseif sum(strcmp(log.onsetsmodel,{'m_t2_Trialtype';}))==1 % TrialType models
        contrasts.contrasttable='Trialtype';
    elseif sum(strcmp(log.onsetsmodel,{'m_t1_ChoicexTrialtype';}))==1 % ChoicexTrialType models
        error('Error in Contrasts set up. Wrong contrast type selected!')
    else
        error('Error in Contrasts setup: Could not find requested onsets model. Which contrast table (i.e. excel sheet) to use?')
    end
    
    % Details for requested contrast table (excel file)
    switch contrasts.contrasttable
        case 'Compete'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=20;
            row.conditionName=1;
            row.conditionType=2; % Condition or Pmod
        case 'Orthog'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=20;
            row.conditionName=1;
            row.conditionType=2; % Condition or Pmod
        case 'Trialtype'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=78;
            row.conditionName=1;
            row.conditionType=2; % Condition or Pmod
        case 'ChoicexTrialtype'
            col.contrastname=2;
            col.weightstart=3;
            col.num_conds=191;
            row.conditionName=1;
            row.conditionType=2; % Condition or Pmod
        otherwise
            error('Error in Contrasts setup: Requested contrast table (i.e. excel sheet) not found')
    end
    col.weightend=col.weightstart+col.num_conds-1;
    col.contrastnum=col.weightend+1;
    col.requested=col.contrastnum+1;
    
end

% input('Continue to compile requested Contrasts?    ')

%% Compile requested Contrasts (Read from excel table)

% Load default details regarding available contrasts
[w.a1, w.a, w.req]=xlsread([where.where filesep 'i_firstlevel_contrasts'], contrasts.contrasttable);
for i=col.contrastname+1:col.contrastnum-1 % Read details
    contrasts.condnames{i-col.contrastname}=w.req{row.conditionName,i};
    contrasts.condtypes{i-col.contrastname}=w.req{row.conditionType,i};
end

% Construct details for requested contrasts
i=1;
for r=row.conditionType+1:size(w.req,1);
    if isnan(w.req{r,col.requested})==0
        contrasts.contrastnames{i,1}=w.req{r,col.contrastname};
        contrasts.requested(i,1)=w.req{r,col.requested};
        
        % Load condition weights
        for c=col.contrastname+1:col.contrastnum-1 % Load all weights
            contrasts.weights(i,c-col.contrastname)=w.req{r,c};
        end
        %
        i=i+1;
    end
end

% Un-request non-applicable contrasts (e.g. for models without all RL variables)
for o2=1:1
    if strcmp(log.onsetsmodel(1:3), 'm_t')==0
        % Trialtype models do not exclude RL variables
        
        % Which RL variables to exclude?
        allvar={'EnvThreat' 0;'NTokens' 0; 'pLoss' 0; 'Entropy' 0; 'Conflict' 0; 'OutcomeMean' 0; 'OutcomeVariance' 0};
        w.spm=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']);
        log.RLvariables=w.spm.SPM.RLvariables; j=1;
        for i=1:size(allvar,1)
            if isempty(cell2mat(strfind(log.RLvariables, allvar{i})))==1
                varremove{j,1}=allvar{i}; j=j+1;
            end
        end
        w.isthere=cell(length(varremove),1);
        for i=1:length(varremove)
            w.whichtoremove=strfind(contrasts.condnames',varremove{i});
            w.isthere{i}=zeros(size(w.whichtoremove,1),1);
            for j=1:size(w.whichtoremove,1)
                if isempty(w.whichtoremove{j})==0
                    w.isthere{i}(j)=1;
                end
            end
        end
        
        % Un-request
        for i=1:length(varremove)
            contrasts.requested(find(w.isthere{i}))=0;
        end
        
        disp('#################################################')
        disp('[RL VARIABLES REMOVED]       Default contrasts included all RL variables. ')
        disp('The following RL variables are not included in this model:')
        disp(' '); disp(varremove); disp(' ');
        disp('Contrasts including these variables are un-requested')
        disp('#################################################')
    end
end

% input('Continue to execute Contrasts?    ')

%% (2) Execute Contrasts

disp('Running contrasts #######################################')
for s=1:log.n_subjs
    
    disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
    
    % Organize folder for this model
    wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep];
    wb.wheremodel=[wb.where log.onsetsmodel ' Contrasted' filesep];
    if isdir(wb.wheremodel)==0;
        disp('Copying folder')
        copyfile([wb.where filesep log.onsetsmodel ' Estimated'], [wb.wheremodel]);
        disp('Done copying')
    end
    
    % Grab details about this model
    f   = spm_select('List', wb.wheremodel, 'SPM.mat');
    wb.spm  =load([wb.wheremodel filesep f]);
    wb.regnamelist=wb.spm.SPM.xX.name';
    
    % Specify requested contrasts
    disp('Specifying Contrasts ---------------------------');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr([wb.wheremodel filesep f]); c=1; 
    con.weights=cell(length(contrasts.contrastnames),1); con.chosenregs=cell(length(contrasts.contrastnames),1); con.names=cell(length(contrasts.contrastnames),1);
    for k= 1:length(contrasts.contrastnames)
        if contrasts.requested(k)==1 && mean(contrasts.weights(k,:))~=0
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = contrasts.contrastnames{k};
            
            % Which regressors to actively weight?
            wc.condnums=find(contrasts.weights(k,:)~=0);
            for i=1:length(wc.condnums)
                wc.regnames{i,1}=contrasts.condnames{wc.condnums(i)};
                wc.regnames{i,2}=contrasts.condtypes{wc.condnums(i)};
            end
            [ wc.regdetails ] = f_GetRegressorNums(wb.regnamelist, wc.regnames );
            
            % Compile weights for this contrast
            wc.weights=zeros(1,length(wb.regnamelist));
            for i=1:size(wc.regdetails ,1)
                wc.weights(wc.regdetails{i,2}(1))=contrasts.weights(k,wc.condnums(i));
            end
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=wc.weights;
            
            % Record
            con.weights{k}=wc.weights; con.chosenregs{k}=wc.regdetails; con.names{k}=contrasts.contrastnames{k};
            
            wc=[]; c=c+1;
        end
    end
    
    disp('Running contrasts -----------------')
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];wb=[];
    
end

end

