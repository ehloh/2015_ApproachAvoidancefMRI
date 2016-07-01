function [Con log] = f1_2_Trialtype(where, log)
% % [Con log] =  f1_TaskChoiceTrialtype(where, log)
% From flexible 1st level model, apply contrasts to derive full task design 
%                                   (Task x Choice x Trial type, but ignoring ME choice)
%
% Second level models included:
%   (1) 2 x 6 x 6 Factorial, Task x EnvThreat x NTokens
%
% Output 'ConInstruc' follows specifications of RegWeights output from the
% function xfx_2GetRegweightsFromSearchCriteria, details the (subject-specific) 
% regerssor weights for each contrast for each subject
%
% ---------------------------------------------------------------------------------------


%% (1) Establish general instructions: Which regressors to weight (via search criteria) for each First-Level contrast?

if  log.execute_contrasts==1
    con.taskchoices={'ConflictTrial'           {1 {'cF_'}};                               % ME Task
                                'ControlTrial'           {1 {'ct_NoBomb'; 'ct_Bomb'; 'ct_Explore'}};
                                'Accept'                  {1 {'Accept'; 'NoBomb'}};      % ME Choice
                                'Reject'                   {1 {'Reject'; 'ct_Bomb'}};
                                'Explore'                 {1 {'Explore'}};
                                'cF_Accept'            {1 {'cF_Accept'}};                      % Task x Choice
                                'cF_Reject'             {1 {'cF_Reject'}};
                                'cF_Explore'           {1 {'cF_Explore'}};
                                'ct_NoBomb'          {1 {'ct_NoBomb'}};
                                'ct_Bomb'              {1 {'ct_Bomb'}};
                                'ct_Explore'            {1 {'ct_Explore'}}   };
    
    % Variable levels (Env, NTokens)
    for i=1:6
        w.e{i,1}=['Env' num2str(i)];
        w.e{i,2}={1 {['_t' num2str(i)]}};
        w.n{i,1}=['N' num2str(i)];
        w.n{i,2}={1 {['-' num2str(i)]}};
    end
    con.varlevels=vertcat(w.e,w.n);
    
    % Variable cells (Env x N Tokens)
    k=1;con.varcell=cell(36,2);
    for e=1:6
        for n=1:6
            con.varcell{k,1}=['e' num2str(e) '_n' num2str(n)];
            con.varcell{k,2}={1 {['_t' num2str(e) '-' num2str(n)]}}; k=k+1;
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
                con.taskvarcell{k,1}=[tasks{t,1} 'e' num2str(e) '_n' num2str(n)]; % Name of contrast
                con.taskvarcell{k,2}{1}=1; % Weight value to apply
                for i=1:length(tasks{t,2}) % Search criteria
                    con.taskvarcell{k,2}{1,2}{i}=[tasks{t} tasks{t,2}{i} '_t' num2str(e) '-' num2str(n)];
                end
                
                % Display
                disp(['Criteria for       ' con.taskvarcell{k,1} ' ------------------------'])
                disp(con.taskvarcell{k,2}{1,2}')
                disp(' ')
                %
                k=k+1;
            end
        end
    end
end

%% (2) Using general search criteria (applied to regressor names), derive subject-specific regressor weights for each contrast x each subject
% For this contrast type, search criteria is not subject-specific, because  the ultimate model
% ignores choice (the only thing to vary across subjects).
%
% But because the 1st level model is split according to Task x Env x NTokens x Choice, 
% the specific regressors to be included to derive  each design cell differs from subject to subject 
% (e.g. if for a cell a subject never chose Explore, only the other choice regressors for that
% cell will be included, since a Explore regressor for that cell does not exist).

if  log.execute_contrasts==1
    RegLists=cell(log.n_subjs,1); SearchCriterias=cell(log.n_subjs,1); ConInstruc=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------'])
        [RegLists{s}] = xfx_1LoadSimpleReglist(where, log,s); % Load (edited) list of regressors
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep];
        
        % Determine subject-specific weights, to derive entire design
        disp('####### Applying search criteria to derive regressor weights for contrasts (deriving entire design)  ##########'); disp(' ');
        SearchCriterias{s}=vertcat(con.taskchoices, con.varlevels, con.varcell, con.taskvarcell);
        [ConInstruc{s}] = xfx_2GetRegweightsFromSearchCriteria(SearchCriterias{s}, RegLists{s});
        ws=[];
    end
    %
    disp('ELEANOR! You really should script in some checks here.')
    
    % Record/output details
    Con.ConInstruc=ConInstruc;
    Con.RegLists=RegLists;
    Con.SearchCriterias=SearchCriterias;
end

%% (3) Execute contrasts in batch if requested - or Fetch already-contrasted model details

if  log.execute_contrasts==1
    disp('############### Executing contrasts #####################')
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep];
        [ matlabbatch] = xfx_3Setup1stlevelContrastbatch(where,log, ws, s,ConInstruc{s});
        ws=[];
    end
    Con.matlabbatch=matlabbatch;
end

%% (4) Set up the 2nd level models (executed next) that go with this model
% Establish instructions for 2nd level analyses to run ################
% Col 1=Name of 2nd-level model, Col 2=Which analysis function, Col 3=Inputs to function
% Instructions are executed in main script to run all possible 2nd level models

% Fetch details from the contrasted first level model
disp(' ################ Prepping for 2nd level ##################')
disp('Fetching details from (already-contrasted) first level model - sampling 1st subject')
c=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep  'SPM.mat']);
for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
disp(Con.ConNames);

% Models to run? 
%   1: [Factorial] Task x Env x NTokens (Full design)
%   2: [Paired ttest] Comparing tasks
Nmodels=1;
Analysis=cell(Nmodels,3);

% ################################################################

% [Model #1] Task x EnvThreat x NTokens ################
for o1=1:1 
    m=1; Analysis{m,1}='TaskxEnvxNtokens';
    Analysis{m,2}='NFactorial';
    Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
    %
    Analysis{m,3}.Covar=[];
    Analysis{m,3}.Design{1}.name='Task';      % Factor 1: Task
    Analysis{m,3}.Design{1}.levels=2;
    Analysis{m,3}.Design{1}.dept=0;
    Analysis{m,3}.Design{1}.variance=1;
    Analysis{m,3}.Design{1}.gmsca=0;
    Analysis{m,3}.Design{1}.ancova=0;
    Analysis{m,3}.Design{2}.name='EnvThreat';      % Factor 2: EnvThreat
    Analysis{m,3}.Design{2}.levels=6;
    Analysis{m,3}.Design{2}.dept=0;
    Analysis{m,3}.Design{2}.variance=1;
    Analysis{m,3}.Design{2}.gmsca=0;
    Analysis{m,3}.Design{2}.ancova=0;
    Analysis{m,3}.Design{3}.name='NTokens';      % Factor 3: N Tokens
    Analysis{m,3}.Design{3}.levels=6;
    Analysis{m,3}.Design{3}.dept=0;
    Analysis{m,3}.Design{3}.variance=1;
    Analysis{m,3}.Design{3}.gmsca=0;
    Analysis{m,3}.Design{3}.ancova=0;
    
    % Assign contrast images to correct design cells (via search instruction)
    Analysis{m,3}.Cells=cell(2*6*6, 4); c=1;
    for t=1:2 
        for e=1:6
            for n=1:6
                switch t
                    case 1
                        Analysis{m,3}.Cells{c,1}=['cF_e' num2str(e) '_n' num2str(n)];
                    case 2
                        Analysis{m,3}.Cells{c,1}=['ct_e' num2str(e) '_n' num2str(n)];
                end
                Analysis{m,3}.Cells{c,2}=t;
                Analysis{m,3}.Cells{c,3}=e;
                Analysis{m,3}.Cells{c,4}=n;
                c=c+1;
            end
        end
    end
    for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
        if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
            Analysis{m,3}.Cells{i,1}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
        else
            error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
        end
    end
end

% DISUSED - [Models #2 and #3] Main effects of Conflict & Control task trials
for o1=1:1
%     
%     Conflict ------------------------------------------
%     m=2; Analysis{m,1}='ME_ConflictTrial';
%     Analysis{m,2}='onesampleT';
%     Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%     
%     Analysis{m,3}.Covar=[];
%     if length(find(strcmp(Con.ConNames(:,1), 'ConflictTrial')))==1 % Fetch correct con image
%         Analysis{m,3}.ConImage=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'ConflictTrial')),2};
%     else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%     end
%     
%     Control ------------------------------------------
%     m=3; Analysis{m,1}='ME_ControlTrial';
%     Analysis{m,2}='onesampleT';
%     Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%     
%     Analysis{m,3}.Covar=[];
%     if length(find(strcmp(Con.ConNames(:,1), 'ControlTrial')))==1 % Fetch correct con image
%         Analysis{m,3}.ConImage=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'ControlTrial')),2};
%     else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%     end
%    
end

% ################################################################

% Output to main script
Con.SecondLevelAnalysis=Analysis;

end

