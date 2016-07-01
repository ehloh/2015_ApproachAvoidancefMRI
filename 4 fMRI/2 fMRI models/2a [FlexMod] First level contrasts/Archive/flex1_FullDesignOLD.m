function [con] = flex1_FullDesign(where, log)
% [con] = f_firstlevelcontrasts_TaskEnvNtokens( where, log)
% From TaskxChoicexTrialType 1st level model, apply contrasts to derive full task design
%
% This script is NOT for generating flexible contrasts between certain
% cells (e.g. All Explore vs all Non-explore, controlling for Design cells)
%
%   con.subj_regweights: (Subject-specific) Which regressors were weighted for each contrast?
%                 Col 1: Name of contrast
%                 Col 2: Search criteria (applied to regressor names)
%                 Col 3: Regressor nos. that were positively weighted (+1)
%                 Col 4: Regressor names that were weighted in this contrast
%
% ---------------------------------------------------------------------------------------

% Check: Is this the correct Contrasts table for the chosen onsets model?
onsetsmodels4thiscontrasttype={'m_f1_ChoicexTrialtype'};
if sum(strcmp(log.onsetsmodel,onsetsmodels4thiscontrasttype))~=1 
    error('Wrong contrasts type chosen for this onsets model!')
end

%% (1) Establish general instructions: Which regressors to weight (via search criteria) for each contrast?

con.taskchoices={'ConflictTrial'           'cF_';                               % ME Task
                            'ControlTrial'           {'ct_NoBomb'; 'ct_Bomb'; 'ct_Explore'};
                            'Accept'                  {'Accept'; 'NoBomb'};      % ME Choice
                            'Reject'                   {'Reject'; 'ct_Bomb'};
                            'Explore'                 'Explore';
                            'cF_Accept'            'cF_Accept';                      % Task x Choice
                            'cF_Reject'             'cF_Reject';
                            'cF_Explore'           'cF_Explore';
                            'ct_NoBomb'          'ct_NoBomb';
                            'ct_Bomb'              'ct_Bomb';
                            'ct_Explore'            'ct_Explore'};
                        
% Variable levels (Env, NTokens)
for i=1:6
    w.e{i,1}=['Env' num2str(i)]; w.e{i,2}=['_t' num2str(i)];
    w.n{i,1}=['N' num2str(i)]; w.n{i,2}=['-' num2str(i)];
end
con.varlevels=vertcat(w.e,w.n);

% Variable cells (Env x N Tokens)
k=1;con.varcell=cell(36,2);
for e=1:6 
    for n=1:6
        con.varcell{k,1}=['e' num2str(e) '_n' num2str(n)];
        con.varcell{k,2}=['_t' num2str(e) '-' num2str(n)]; k=k+1;
    end
end

% Task Variable cells (Task x Env x N Tokens)
disp('Search criteria for Task x Env x N Tokens contrasts #########################')
k=1;con.taskvarcell=cell(36*2,2); 
tasks={'cF_' {'Accept'; 'Reject'; 'Explore'};
            'ct_' {'NoBomb'; 'Bomb'; 'Explore'}};
for t=1:2 
    for e=1:6
        for n=1:6
            con.taskvarcell{k,1}=[tasks{t,1} 'e' num2str(e) '_n' num2str(n)];
            for i=1:length(tasks{t,2})
                con.taskvarcell{k,2}{i}=[tasks{t} tasks{t,2}{i} '_t' num2str(e) '-' num2str(n)];
            end
            
            % Display
            disp(['Criteria for       ' con.taskvarcell{k,1} ' ------------------------'])
            disp(con.taskvarcell{k,2}')
            disp(' ')
            %
            k=k+1;
        end
    end
end

%% (2) Set up and Execute Subject-specifict Contrasts

for s=1:log.n_subjs
    
    disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------'])
    ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Contrasted ' log.firstlevelcontrasts_type filesep];
    
    % Load (edited) list of regressors
    ws.s=load([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']); ws.SPM=ws.s.SPM;
    f_whatderiv=@(x)x(length(x)-2:length(x));                     % Remove derivatives
    w.a=cellfun(f_whatderiv, ws.SPM.xX.name','UniformOutput',0);
    ws.SPM.xX.name(strcmp(w.a,'(2)'))={'null'}; ws.SPM.xX.name(strcmp(w.a,'(3)'))={'null'};
    f_readrealname=@(x)x(7:length(x)-6);                          % Simplify names
    ws.regnames=cellfun(f_readrealname,ws.SPM.xX.name,'UniformOutput',0)';
    
    % Determine subject-specific weights, to derive entire design (con.instruc, 3rd col)
    disp('Specifying Contrasts (determining subject-specific regressor weights to derive entire design) ---------------------------');
    con.instruc=vertcat(con.taskchoices,con.varlevels, con.varcell, con.taskvarcell);
    ws.con_regweights=cell(size(con.instruc,1),1); 
    for c=1:size(con.instruc,1); % Compile weights (con.instruc, 3rd col)
        if iscell(con.instruc{c,2})==0
            wc.w=strfind(ws.regnames,con.instruc{c,2});
            wc.w= double(~cellfun(@isempty,wc.w));  % Positive weights applied to all regressors that find a string match
            con.instruc{c,3}=wc.w;
        else
            wc.a=zeros(length(ws.regnames),1);
            for i=1:length(con.instruc{c,2})
                wc.w{i}=strfind(ws.regnames,con.instruc{c,2}{i});
                wc.w{i}= double(~cellfun(@isempty,wc.w{i}));  % Positive weights applied to all regressors that find a string match
                wc.a=wc.a+wc.w{i};
            end
            con.instruc{c,3}=wc.a;
        end
        wc.w=[];
        
        % Record details (for checking)
        ws.con_regweights{c,1}=con.instruc{c,1};
        ws.con_regweights{c,2}=con.instruc{c,2};
        ws.con_regweights{c,3}=find(con.instruc{c,3});
        ws.con_regweights{c,4}=ws.regnames(find(con.instruc{c,3}));
        disp(['contrast ' num2str(c) '   ' con.instruc{c,1} '   -    No. of regressors weighted +1: ' num2str(sum(con.instruc{c,3}))]);
        disp(ws.con_regweights{c,4})
    end
    con.subj_regweights{s,1}=ws.con_regweights;
    
    % [Check] How many regressors in each condition: for i=1:59; disp([con.instruc{i,1} '    ' num2str(sum(abs(con.instruc{i,3})))]); end
    
    disp('ELEANOR! You really should script in some checks here.')
    input('Continue?')
    
    % Set up batch (all contrasts together)
    if isdir(ws.where_contrasts)==0; 
        disp('Copying folder for contrasts . .'); copyfile([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Estimated'],ws.where_contrasts); disp('Done')
    else
        input('Folder with the target name (onsets + contrasts type) already exists. Assume correct?')
    end
    matlabbatch{1}.spm.stats.con.spmmat = cellstr([ws.where_contrasts 'SPM.mat' ]);
    for c=1:size(con.instruc,1)
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =con.instruc{c,1};
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=con.instruc{c,3};
    end
    disp('Running contrasts -----------------')
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[]; ws=[]; con.instruc=[];
    
    % Set up batch (each contrasts on its own)
    for c=1:1 % size(con.instruc,1) % DISUSED 
        %     matlabbatch{1}.spm.stats.con.spmmat = cellstr([ws.where_contrasts 'SPM.mat' ]);
        %     matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        %     matlabbatch{1}.spm.stats.con.consess{1}.tcon.name =con.instruc{c,1};
        %     matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec=con.instruc{c,3};
        %     %
        %     disp(['Running contrast ' num2str(c) '  out of  ' num2str(size(con.instruc,1)) ' -----------------'])
        %     spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        %     matlabbatch=[];
    end
end

end

