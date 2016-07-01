function [variables c] = m_c8_ChCluster4MovFull_OULPEN(datastruc, col, c)
% Full model, but with Choice split as a function of exploration cluster or not
% Exploration cluster is determined on a subject-by-subject manner
%
% Competing: Each variable has its own onset + pmod
%                   Duplicate onsets are ignored in first-level contrasts
%
% Variables:  EntropyNTok (O), Entropy (U), EV (L), pLoss (P), EnvThreat (E), NTokens (N)
%
%       datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
%
% Evaluate to debug: datastruc=ws.d;

RLvariables={'EntropyNTok'; 'Entropy';  'EV'; 'pLoss'; 'EnvThreat'; 'NTokens';}; 

Min_percent=0.0;% Minimum % of responses of Reject/Bomb & Explore for entire cluster
N_cells=4; % Size of inner cluster

%% Define cluster (subject-specific basis)

% Col 1/6=EnvThreat, Col 2/7=NTokens, Col 3=%Accept, Col 4=%Reject, 
% Col 5=%Explore, Col 8=%NoBomb, Col 9=%Bomb, Col 10=%Explore
stat=zeros(36,10);
for t=1:2
    wt.d=datastruc{t}(datastruc{t}(:,col.TrialValid)==1,:);
    wt.count=zeros(36,5); k=1;
    for e=1:6
        for n=2:2:12
            wd=wt.d(wt.d(:,col.EnvThreat)==e & wt.d(:,col.NTokens)==n,:);
            %
            wt.count(k,1)=e;
            wt.count(k,2)=n;
            wt.count(k,3)=sum(wd(:,col.Resp1)==1)/size(wd,1);
            wt.count(k,4)=sum(wd(:,col.Resp1)==2)/size(wd,1);
            wt.count(k,5)=sum(wd(:,col.Resp1)==3)/size(wd,1);
            %
            wd=[]; k=k+1;
        end
    end
    stat(:,(t-1)*5+1:5*t)=wt.count;
end
stat=sortrows(stat,-5);

% Check: Good no. of Reject/Explore & Bomb/Explore?
w.a=sum(stat(1:N_cells,:))/N_cells;
if w.a([3 4 5 8 9 10])>Min_percent; 
elseif w.a([4 5 9 10])>Min_percent; 
    disp('Inner cluster ok - outer cluster not'); 
else
    disp('Problematic subject - % events ');
end

% Define cluster
stat(:,[2 7])=stat(:,[2 7])/2; % Cluster's NTokens is in levels (i.e. Token pairs)
InnerCluster=stat(1:N_cells,1:2);

%% 

% CHOICE regressors
[variables c] =  i_choice_conditions_cluster(datastruc{1}, [], 'cF_', {'Accept';'Reject';'Explore'}, col, c,InnerCluster ); % Conflict
[variables c] =  i_choice_conditions_cluster(datastruc{2}, variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c,InnerCluster ); % Control

% RL VARIABLE regressors
[variables c] = i_RLvariables_Competing(datastruc{1}, 'cF_', RLvariables, variables, col, c ); % Conflict
[variables c] = i_RLvariables_Competing(datastruc{2}, 'ct_', RLvariables, variables,col,  c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';'EntropyNTok'}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );


end


