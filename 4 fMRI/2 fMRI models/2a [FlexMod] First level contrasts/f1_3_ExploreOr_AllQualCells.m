function [Con log] = f2_ExploreOr_AllQualCells(where, log)
% [Con log] = f2_ExploreOr_AllQualCells(where, log)
% Compare Explore vs NonExplore, limited only to cells where both nExplore 
% and nOr have enough events (minimum is specified in function). 
% All qualifying cells are included. Weights do NOT adjust for 
% no. of events in each condition. 
%
% Second level models included:
%   (1) Task x Explore (2x2 Factorial)
%   (2) Add sensible contrasts to the 2x2 (Task x Explore) Factorial
%
% Output 'Con.ConInstruc' follows specifications of RegWeights output from the
% function xfx_2GetRegweightsFromSearchCriteria, details the (subject-specific) 
% regerssor weights for each contrast for each subject
%
% ---------------------------------------------------------------------------------------

%% (1) General changeable settings

settings.min_num_events=1;

%% (2) Generate subject-specific search criteria (First level contrasts)
% Search criteria is fed into function xfx_2GetRegweightsFromSearchCriteria,
% to derive regressor weights for each contrast. Format of search criteria
% must therefore follow format specified in this search function.

if  log.execute_contrasts==1
    disp('########### Specifying subject-specific search criteria (to find regressors for weighting) #####################')
    RegLists=cell(log.n_subjs,1); SearchCriterias=cell(log.n_subjs,1); Cells2include=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        
        disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------'])
        [RegLists{s}] = xfx_1LoadSimpleReglist(where, log,s); % Load (edited) list of regressors
        ws.b=load([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.subjects{s} '_behstats.mat']);
        
        for o1=1:1 % Which cells to include? (i.e. > Min. trials for both Explore AND Non-explore)
            % cells2include 'index': Col 1=EnvThreat, Col 2=NTokens,
            %                                 Col 3=nExp, Col 4=nOr
            cells2include.cF=zeros(6,6);
            cells2include.cF_nExp=zeros(6,6);
            cells2include.cF_nOr=zeros(6,6);
            cells2include.ct=zeros(6,6); cells2include.ct_nExp=zeros(6,6); cells2include.ct_nOr=zeros(6,6);
            f=1; t=1;
            for e=1:6 % Note: EnvThreat always kept in matrix format for visualization (descending threat)
                for n=1:6
                    % Conflict task
                    if ws.b.cF_stats{7-e,n}.n_Explore>settings.min_num_events && (ws.b.cF_stats{7-e,n}.n_Accept+ws.b.cF_stats{7-e,n}.n_Reject)>settings.min_num_events
                        cells2include.cF(7-e,n)=1;
                        cells2include.cF_nExp(7-e,n)=ws.b.cF_stats{7-e,n}.n_Explore;
                        cells2include.cF_nOr(7-e,n)=ws.b.cF_stats{7-e,n}.n_Accept+ws.b.cF_stats{7-e,n}.n_Reject;
                        cells2include.cF_index(f,1)=e; % Instructions list
                        cells2include.cF_index(f,2)=n;
                        cells2include.cF_index(f,3)=cells2include.cF_nExp(7-e,n);
                        cells2include.cF_index(f,4)=cells2include.cF_nOr(7-e,n); f=f+1;
                    end
                    
                    if ws.b.ct_stats{7-e,n}.n_Explore>settings.min_num_events && (ws.b.ct_stats{7-e,n}.n_NoBomb+ws.b.ct_stats{7-e,n}.n_Bomb)>settings.min_num_events
                        cells2include.ct(7-e,n)=1; % Visualization matrices
                        cells2include.ct_nExp(7-e,n)=ws.b.ct_stats{7-e,n}.n_Explore;
                        cells2include.ct_nOr(7-e,n)=ws.b.ct_stats{7-e,n}.n_NoBomb+ws.b.ct_stats{7-e,n}.n_Bomb;
                        cells2include.ct_index(t,1)=e; % Instructions list
                        cells2include.ct_index(t,2)=n;
                        cells2include.ct_index(t,3)=cells2include.ct_nExp(7-e,n);
                        cells2include.ct_index(t,4)=cells2include.ct_nOr(7-e,n); t=t+1;
                    end
                end
            end
            
        end
        
        % ######### Format for the search function #############
        ws.Crit=cell(7,2); c=1;
        
        % [Contrasts 1-6] Choices ################
%         choices={'cF_Accept'; 'cF_Reject'; 'cF_Explore'; 'ct_NoBomb'; 'ct_Bomb';'ct_Explore'};
%         for i=1:6
%             ws.Crit{c,1}=choices{i};
%             ws.Crit{c,2}{1,1}=1;
%             ws.Crit{c,2}{1,2}{1}=choices{i};
%             c=c+1;
%         end
        
        % -------------- CELLS -----------------------------------------------------------------------------------
        % [Contrast 6+1] Conflict task, Explore ################
        ws.Crit{c,1}='cF_Explore'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.cF_index,1) ,2);
        for i=1:size(cells2include.cF_index,1)
            tt=['_t' num2str(cells2include.cF_index(i,1)) '-' num2str(cells2include.cF_index(i,2))];
            %
            ws.Crit{c,2}{v,1}=1;  % cells2include.cF_index(i,4);
            ws.Crit{c,2}{v,2}={['cF_Explore' tt];};
            v=v+1;            
        end
        c=c+1;
        
        % [Contrast 6+2] Conflict task, NonExplore ################
        ws.Crit{c,1}='cF_NonExplore'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.cF_index,1) ,2);
        for i=1:size(cells2include.cF_index,1)
            tt=['_t' num2str(cells2include.cF_index(i,1)) '-' num2str(cells2include.cF_index(i,2))];
            %
            ws.Crit{c,2}{v,1}=1; % cells2include.cF_index(i,3) ;
            ws.Crit{c,2}{v,2}={['cF_Accept' tt];['cF_Reject' tt];};
            v=v+1;
        end
        c=c+1;
        
        % [Contrast 6+3] Control task, Explore ################
        ws.Crit{c,1}='ct_Explore'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.ct_index,1) ,2);
        for i=1:size(cells2include.ct_index,1)
            tt=['_t' num2str(cells2include.ct_index(i,1)) '-' num2str(cells2include.ct_index(i,2))];
            %
            ws.Crit{c,2}{v,1}=1;  % cells2include.cF_index(i,4);
            ws.Crit{c,2}{v,2}={['ct_Explore' tt];};
            v=v+1;            
        end
        c=c+1;
        
        % [Contrast 6+4] Control task, NonExplore ################
        ws.Crit{c,1}='ct_NonExplore'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.ct_index,1) ,2);
        for i=1:size(cells2include.ct_index,1)
            tt=['_t' num2str(cells2include.ct_index(i,1)) '-' num2str(cells2include.ct_index(i,2))];
            %
            ws.Crit{c,2}{v,1}=1; % cells2include.cF_index(i,3) ;
            ws.Crit{c,2}{v,2}={['ct_NoBomb' tt];['ct_Bomb' tt];};
            v=v+1;
        end
        c=c+1;
        
        % -------------- Explore vs Not  -----------------------------------------------------------------------------------
        % [Contrast 6+5] Conflict task, Explore vs Not (weighted by n counter-events) ##########
        ws.Crit{c,1}='cF_ExploreVsNon'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.cF_index,1) * 2 ,2);
        for i=1:size(cells2include.cF_index,1)
            tt=['_t' num2str(cells2include.cF_index(i,1)) '-' num2str(cells2include.cF_index(i,2))];
            
            % Positively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.cF_index(i,4);              % + weight=n NotExplore
            ws.Crit{c,2}{v,2}={['cF_Explore' tt];};
            v=v+1;
            
            % Negatively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.cF_index(i,3)  *  -1 ;      % - weight=nExplore
            ws.Crit{c,2}{v,2}={['cF_Accept' tt];['cF_Reject' tt];};
            v=v+1;
            
        end
        c=c+1;
        
        % [Contrast 6+6] Control task, Explore vs Not (weighted by n counter-events) ##########
        ws.Crit{c,1}='ct_ExploreVsNon'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.ct_index,1) * 2 ,2);
        for i=1:size(cells2include.ct_index,1)
            tt=['_t' num2str(cells2include.ct_index(i,1)) '-' num2str(cells2include.ct_index(i,2))];
            
            % Positively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.ct_index(i,4);              % + weight=n NotExplore
            ws.Crit{c,2}{v,2}={['ct_Explore' tt];};
            v=v+1;
            
            % Negatively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.ct_index(i,3)  *  -1 ;      % - weight=nExplore
            ws.Crit{c,2}{v,2}={['ct_NoBomb' tt];['ct_Bomb' tt];};
            v=v+1;
            
        end
        c=c+1;
        
        % [Contrast 6+7] Both tasks, Explore vs Not (weighted by n counter-events) ################
        ws.Crit{c,1}='ExploreVsNon'; v=1;
        ws.Crit{c,2}=cell(size(cells2include.cF_index,1) * 2  + size(cells2include.ct_index,1) * 2 ,2);
        for i=1:size(cells2include.cF_index,1) % Conflict task
            tt=['_t' num2str(cells2include.cF_index(i,1)) '-' num2str(cells2include.cF_index(i,2))];
            
            % Positively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.cF_index(i,4);              % + weight=n NotExplore
            ws.Crit{c,2}{v,2}={['cF_Explore' tt];};
            v=v+1;
            
            % Negatively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.cF_index(i,3)  *  -1 ;      % - weight=nExplore
            ws.Crit{c,2}{v,2}={['cF_Accept' tt];['cF_Reject' tt];};
            v=v+1;
            
        end
        for i=1:size(cells2include.ct_index,1)  % Control task
            tt=['_t' num2str(cells2include.ct_index(i,1)) '-' num2str(cells2include.ct_index(i,2))];
            
            % Positively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.ct_index(i,4);              % + weight=n NotExplore
            ws.Crit{c,2}{v,2}={['ct_Explore' tt];};
            v=v+1;
            
            % Negatively weighted (Explore)
            ws.Crit{c,2}{v,1}=1; % cells2include.ct_index(i,3)  *  -1 ;      % - weight=nExplore
            ws.Crit{c,2}{v,2}={['ct_NoBomb' tt];['ct_Bomb' tt];};
            v=v+1;
            
        end
        c=c+1;
        
        %
        SearchCriterias{s}=ws.Crit;
        Cells2include{s}=cells2include;
        ws=[]; cells2include=[];
    end
    
end

%% (3) Get regressor weights for First-level contrast

% error('Pause')

if  log.execute_contrasts==1
    disp('####### Applying search criteria to derive regressor weights for contrasts  ##########'); disp(' ');
    ConInstruc=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##################'])
        [ConInstruc{s}] = xfx_2GetRegweightsFromSearchCriteria(SearchCriterias{s}, RegLists{s});
    end
    
end

%% (4) Generate batch for First-level contrasts

if  log.execute_contrasts==1
    
    input('Continue to execute contrasts?    ')
    disp('############### Executing contrasts #####################')
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep];
        [matlabbatch] = xfx_3Setup1stlevelContrastbatch(where,log, ws, s,ConInstruc{s});
    end
    
    % Record/output details
    Con.RegLists=RegLists;
    Con.Cells2include=Cells2include;
    Con.SearchCriterias=SearchCriterias;
    Con.ConInstruc=ConInstruc;
    Con.matlabbatch=matlabbatch;
end

%% (5) Set up the 2nd level models (executed next) that go with this model
% Establish instructions for 2nd level analyses to run ################
% Col 1=Name of 2nd-level model, Col 2=Which analysis function, Col 3=Inputs to function
% Instructions are executed in main script to run all possible 2nd level models

% Fetch details from the contrasted first level model
disp(' ################ Prepping for 2nd level ##################')
disp('Fetching details from (already-contrasted) first level model - sampling 1st subject')
c=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep 'SPM.mat']);
for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
disp(Con.ConNames);

% Models to run? 
%   1: Task x Explore factorial
%   2: Add sensible contrasts to the 2x2 (Task x Explore) Factorial
Nmodels=2;
Analysis=cell(Nmodels,3);

% ################################################################

% [Model #1] Task x Explore ################
for o1=1:1 
    m=1; Analysis{m,1}='TaskxExplore';
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
    Analysis{m,3}.Design{2}.name='Explore';      % Factor 2: Explore
    Analysis{m,3}.Design{2}.levels=2;
    Analysis{m,3}.Design{2}.dept=0;
    Analysis{m,3}.Design{2}.variance=1;
    Analysis{m,3}.Design{2}.gmsca=0;
    Analysis{m,3}.Design{2}.ancova=0;
    % Assign contrast images to correct design cells (via search instruction)
    Analysis{m,3}.Cells=cell(2*2, 3); 
    Analysis{m,3}.Cells={'cF_Explore'             1 1;
                                      'cF_NonExplore'       1 2;
                                      'ct_Explore'              2 1;
                                      'ct_NonExplore'        2 2;};
    for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
        if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
            Analysis{m,3}.Cells{i,1}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
        else
            error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
        end
    end
end

% [Model #2] Add sensible contrasts to the Task x Explore factorial
for o1=1:1 
    m=m+1; Analysis{m,1}='Adding Factorial contrasts to TaskxExplore';
    Analysis{m,2}='AddContrasts';
    %
    Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep 'TaskxExplore'];
    Analysis{m,3}.DeleteOldCon=0;
    Analysis{m,3}.ConInstruc={'ConflictTask'                [1 1 0 0]         1;     % Main factor means
                                                'ControlTask'                [0 0 1 1]         1;
                                                'Explore'                       [1 0 1 0]         1;
                                                'NonExplore'                [0 1 0 1]         1;
                                                'cF_Explore'                [1 0 0 0]         1;        % Cell means
                                                'ct_Explore'                [0 0 1 0]         1;
                                                'cF_NonExplore'           [0 1 0 0]         1;
                                                'ct_NonExplore'             [0 0 0 1]         1;
                                                'Task_cF-ct'             [1 1 -1 -1]         1;        % Factor effects
                                                'Task_ct-cF'             [-1 -1 1 1]         1;
                                                'Explore-Non'             [1 -1 1 -1]         1;
                                                'Non-Explore'             [-1 1 -1 1]         1;
                                                'cF_Explore-Non'        [1 -1 0 0]         1;   % Cell effects
                                                'cF_Non-Explore'        [-1 1 0 0]         1;
                                                'ct_Explore-Non'        [0 0 1 -1]         1;
                                                'ct_Non-Explore'        [0 0 -1 1]         1;
                                                'cF_Explore_vs'         [3 -1 -1 -1]         1;   % Cell versus
                                                'cF_NonExplore_vs'    [-1 3 -1 -1]         1;
                                                'ct_Explore_vs'          [-1 -1 3 -1]         1;
                                                'ct_NonExplore_vs'     [-1 -1 -1 3]         1;
                                                'TaskxExp_Inter1'       [1 -1 -1 1]         1; % Interaction terms
                                                'TaskxExp_Inter2'       [-1 1 1 -1]         1;};
end

for o1=1:1  % DISUSED MODELS 
%     % [Model #2-4] Comparisons: cF_ExporeOr, cf_ExploreOr, ExploreOr ################
%     for o1=1:1 
% 
%         % cF_Explore_Vs_Or
%         m=m+1; Analysis{m,1}='cF_Explore_Vs_Or';
%         Analysis{m,2}='onesampleT';
%         Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%         %
%         Analysis{m,3}.Covar=[];
%         if length(find(strcmp(Con.ConNames(:,1), 'cF_ExploreOr')))==1 % Fetch correct con image
%             Analysis{m,3}.ConImage=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'cF_ExploreOr')),2};
%         else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%         end
% 
%         % ct_Explore_Vs_Or
%         m=m+1; Analysis{m,1}='ct_Explore_Vs_Or';
%         Analysis{m,2}='onesampleT';
%         Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%         %
%         Analysis{m,3}.Covar=[];
%         if length(find(strcmp(Con.ConNames(:,1), 'ct_ExploreOr')))==1 % Fetch correct con image
%             Analysis{m,3}.ConImage=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'ct_ExploreOr')),2};
%         else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%         end
% 
%         % Explore_Vs_Or
%         m=m+1; Analysis{m,1}='Explore_Vs_Or';
%         Analysis{m,2}='onesampleT';
%         Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%         %
%         Analysis{m,3}.Covar=[];
%         if length(find(strcmp(Con.ConNames(:,1), 'ExploreOr')))==1 % Fetch correct con image
%             Analysis{m,3}.ConImage=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'ExploreOr')),2};
%         else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%         end
%     end
% 
%     % [Model #5] Comparing comparisons across task (cF_ExporeOr vs cf_ExploreOr) ################
%     m=m+1; Analysis{m,1}='cF_vs_ct_ExploreOr';
%     Analysis{m,2}='pairedT';
%     Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%     %
%     Analysis{m,3}.Covar=[];
%     if length(find(strcmp(Con.ConNames(:,1), 'cF_ExploreOr')))==1 % Fetch correct con image
%         Analysis{m,3}.ConImage1=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'cF_ExploreOr')),2};
%     else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%     end
%     if length(find(strcmp(Con.ConNames(:,1), 'ct_ExploreOr')))==1 % Fetch correct con image
%         Analysis{m,3}.ConImage2=Con.ConNames{find(strcmp(Con.ConNames(:,1), 'ct_ExploreOr')),2};
%     else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
%     end
end

% ################################################################

% Output to main script
Con.SecondLevelAnalysis=Analysis;

end

