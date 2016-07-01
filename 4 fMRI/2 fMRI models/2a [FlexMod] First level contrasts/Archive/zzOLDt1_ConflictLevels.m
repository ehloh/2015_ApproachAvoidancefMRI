function [Con logg] = t1_ConflictLevels(where, logg)
% [Con logg] = t1_ConflictLevels(where, logg)
% Generate contrasts for levels of conflict (16 semi-arbitrary), ignoring
% all other variables
%
% Second level models included:
%   (1) 2x6 (Task x Conflict) levels Factorial
%
% Output 'Con.ConInstruc' follows specifications of RegWeights output from the
% function xfx_2GetRegweightsFromSearchCriteria, details the (subject-specific) 
% regerssor weights for each contrast for each subject
%         
% ---------------------------------------------------------------------------------------

% Execute (First time only): logg=log; clear('log')

% Check: Is this the correct Contrasts table for the chosen onsets model?
onsetsmodels4thiscontrasttype={'m_t1_Trialtype'};
if sum(strcmp(logg.onsetsmodel,onsetsmodels4thiscontrasttype))~=1
    error('Wrong contrasts type chosen for this onsets model! Or wrong onsets model for this contrast.')
end

%% Set up details
% From true values of conflict, derive levels (via bins, e.g. x decimal places)
% Formulation: Conflict=NTokens x Entropy; Entropy= -p(log p) - (1-p) (log (1-p) )

% Establish design first (RL variables)
d.Env=((6:-1:1)'*ones(1,6))/6;
d.N=2*ones(6,1)*(1:6);
d.pLoss=d.Env.*(d.N/12);
d.Unc=nan*zeros(6,6); d.Conf=nan*zeros(6,6);
for e=1:6
    for n=1:6
        if d.pLoss(e,n)==1 % Correction for infinite values
            d.pLoss(e,n)=d.pLoss(e,n)-0.000001;
        end
        
        % Uncertainty/Entropy=-p (logp) - (1-p) log(1-p)
        d.Unc(e,n)= - d.pLoss(e,n).* log(d.pLoss(e,n)) - (1-d.pLoss(e,n)) *(log(1-d.pLoss(e,n)));
        
        % Conflict= Uncertainty * NTokens
        d.Conf(e,n)=d.Unc(e,n)*d.N(e,n);
    end
end

% Bin values of conflict
d.Conf_2dp=round(d.Conf*100)/100; % 32 unique values
d.Conf_1dp=round(d.Conf*10)/10;  % 30 unique values)
d.Conf_integ=round(d.Conf); % 9 unique values
d.Conf_halfpoint=round(2*d.Conf)/2; % 16 unique values (arbitrary 0.5 value steps)

% Which type of binning to use?
d.Conf_score=d.Conf_halfpoint;
d.Conf_levels=num2cell(unique(d.Conf_score));
disp('Raw conflict values (by trial type):'); disp(d.Conf)
disp('Conflict levels (by trial type):'); disp(d.Conf_score)
figure('NumberTitle', 'off', 'Name', 'Raw conflict values and Conflict levels (see comm window for values)');
subplot(2,1,1); imagesc(d.Conf); colorbar; axis square;
subplot(2,1,2); imagesc(d.Conf_score); colorbar; axis square;

%% Compile instructions and search criteria

% Compile instructions for each Conflict level (d.Conf_levels)
% Col 1= Conflict level, Col 2: Cells with this value (EnvThreat, NTokens),
% Col 3=Weight for each cell in this conflict level (if each level should
% be equivalent, contrast wise)

if  logg.execute_contrasts==1
    
    % Which cells for each level/bin, or which bin is each cell in?
    for i=1:length(d.Conf_levels)
        [d.Conf_levels{i,2}(:,1)  d.Conf_levels{i,2}(:,2) ]=find(d.Conf_score==d.Conf_levels{i});
        d.Conf_levels{i,2}(:,1) =7-d.Conf_levels{i,2}(:,1) ; % Convert EnvThreat from visualization space to trialtype
        w.ncells(i,1)=size(d.Conf_levels{i,2},1);
    end
    for i=1:size(d.Conf_levels,1)
        d.Conf_levels{i,3}=lcms(w.ncells)/size(d.Conf_levels{i,2},1);
    end
    
    % Format for search function (generalized search criteria, since choice doesn't matter)
    SearchCriteria=cell(6+size(d.Conf_levels,1), 2); c=1;
    
    % Choice regressors
    choices={'cF_Accept'; 'cF_Reject'; 'cF_Explore'; 'ct_NoBomb'; 'ct_Bomb';'ct_Explore'};
    for i=1:6 
        SearchCriteria{c,1}=choices{i};
        SearchCriteria{c,2}{1,1}=1;
        SearchCriteria{c,2}{1,2}{1}=choices{i};
        c=c+1;
    end
    
    for i=1:size(d.Conf_levels,1) % Conflict task
        SearchCriteria{c,1}=['cF_ConflictLevel-' num2str(i)];
        SearchCriteria{c,2}{1,1}=1; % d.Conf_levels{i,3};
        for s=1:size(d.Conf_levels{i,2},1) % Search criteria
            SearchCriteria{c,2}{1,2}{s,1}=['cF_t' num2str(d.Conf_levels{i,2}(s,1)) '-' num2str(d.Conf_levels{i,2}(s,2))];
        end
        c=c+1;
    end
    for i=1:size(d.Conf_levels,1) % Control task
        SearchCriteria{c,1}=['ct_ConflictLevel-' num2str(i)];
        SearchCriteria{c,2}{1,1}=1; % d.Conf_levels{i,3};
        for s=1:size(d.Conf_levels{i,2},1)
            SearchCriteria{c,2}{1,2}{s,1}=['ct_t' num2str(d.Conf_levels{i,2}(s,1)) '-' num2str(d.Conf_levels{i,2}(s,2))];
        end
        c=c+1;
    end
    
    % Display
    disp('######### Search criteria for this model (Which cells are included in what?) #########')
    disp(['NOTE: Each Conflict level''s weights should add up to ' num2str(lcms(w.ncells))]); disp(' ')
    for i=1:size(SearchCriteria,1)
        disp(['Contrast ' num2str(i) ':    ' SearchCriteria{i,1} ' ------------------------'])
        disp('Trial types/regressors weighted for this contrast (i.e this Conflict level)'); disp(SearchCriteria{i,2}{1,2})
        disp(['Weight value for each cell/trialtype at this level: ' num2str(SearchCriteria{i,2}{1,1})])
        disp(' ')
    end
end

%% Generate regressor weights using search critera + search function

if  logg.execute_contrasts==1
    % Get names of regressors
    [RegList] = xfx_1LoadSimpleReglist(where, logg,1); % Load (edited) list of regressors
    disp('List of regressors in this model (sampled from 1st subject)');  disp(RegList)
    
    % Feed into search function
    disp('####### Applying search criteria to derive regressor weights for contrasts  ##########'); disp(' ');
    [ConInstruc] = xfx_2GetRegweightsFromSearchCriteria(SearchCriteria,RegList);
end

%% Execute contrasts (batch)

if  logg.execute_contrasts==1
    input('Continue to execute contrasts?    ')
    disp('############### Executing contrasts #####################')
    
    for s=1:logg.n_subjs
        disp(['Subject ' num2str(s) '   -  '  logg.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep logg.subjects{s} filesep '2 First level' filesep logg.onsetsmodel(1:4) ' Contrasted   ' logg.firstlevel_contraststype filesep];
        [matlabbatch] = xfx_3Setup1stlevelContrastbatch(where,logg, ws, s,ConInstruc);
    end
    
    % Record/output details
    Con.RegList=RegList;
    Con.SearchCriterias=SearchCriteria;
    Con.ConInstruc=ConInstruc;
    Con.matlabbatch=matlabbatch;
    Con.Conflict_levels=d.Conf_levels;
    Con.Design=d;
end


%% (5) Set up the 2nd level models (executed next) that go with this model
% Establish instructions for 2nd level analyses to run ################
% Col 1=Name of 2nd-level model, Col 2=Which analysis function, Col 3=Inputs to function
% Instructions are executed in main script to run all possible 2nd level models

% Fetch details from the contrasted first level model
disp(' ################ Prepping for 2nd level ##################')
disp('Fetching details from (already-contrasted) first level model - sampling 1st subject')
c=load([where.data_brain filesep logg.subjects{1} filesep '2 First level' filesep logg.onsetsmodel(1:4) ' Contrasted   ' logg.firstlevel_contraststype filesep 'SPM.mat']);
for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
disp(Con.ConNames);

% Models to run? 
%   1: Task x Conflict level factorial
%   2-4: One-sample t-test on 3 comparison contrasts
%   5: Paired t-test comparing comparisons (ExploreOr) across task
Nmodels=1;
Analysis=cell(Nmodels,3);

% ################################################################

% error('Pause')

% [Model #1] Task x Conflict level ################
for o1=1:1 
    m=1; Analysis{m,1}='TaskxConflictLevel';
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
    Analysis{m,3}.Design{2}.name='ConflictLevel';      % Factor 2: Conflict Level
    Analysis{m,3}.Design{2}.levels=length(d.Conf_levels);
    Analysis{m,3}.Design{2}.dept=0;
    Analysis{m,3}.Design{2}.variance=1;
    Analysis{m,3}.Design{2}.gmsca=0;
    Analysis{m,3}.Design{2}.ancova=0;
    % Assign contrast images to correct design cells (via search instruction)
    Analysis{m,3}.Cells=cell(2*length(d.Conf_levels), 3); c=1;
    for t=1:2
        for f=1:length(d.Conf_levels)
            switch t
                case 1
                    Analysis{m,3}.Cells{c,1}=['cF_ConflictLevel-' num2str(f)];
                case 2
                    Analysis{m,3}.Cells{c,1}=['ct_ConflictLevel-' num2str(f)];
            end
            Analysis{m,3}.Cells{c,2}=t;
            Analysis{m,3}.Cells{c,3}=f;
            c=c+1;
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

% ################################################################

% Output to main script
Con.SecondLevelAnalysis=Analysis;

end

