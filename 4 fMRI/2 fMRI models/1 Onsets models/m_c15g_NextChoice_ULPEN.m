function [variables c] = m_c3_ChoiceFull_ULPEN(datastruc, col, c)
% Full model: Entropy (U), EV (L), pLoss (P), EnvThreat (E), NTokens (N)
% Competing: Each variable has its own onset + pmod
%                   Duplicate onsets are ignored in first-level contrasts
%
%       datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
% Evaluate to debug: datastruc=ws.d;

RLvariables={'EnvThreat'; 'NTokens';'pLoss'; 'Entropy'; 'EV'; };

% Exclude errors?
%     datastruc{1}=datastruc{1}(datastruc{1}(:,col.TrialValid)==1,:);
%     datastruc{2}=datastruc{2}(datastruc{2}(:,col.TrialValid)==1,:); 

%% 

% Shift choice: label each trial by the next trial's choice
datanext=datastruc;
datanext{1}(:, col.Resp1) =  [datanext{1}(2:end, col.Resp1); nan];
datanext{2}(:, col.Resp1) =  [datanext{2}(2:end, col.Resp1); nan];
for b=1:6   % Nan-out the last trial of each block (it's not predicting anything on the next trial)
    datanext{1}(  find(datanext{1}(:, col.Block)==b, 1, 'last'), col.Resp1)=4; % cf
    datanext{1}(  find(datanext{1}(:, col.Block)==b, 1, 'last'), col.TrialValid)=0;
    datanext{2}(  find(datanext{2}(:, col.Block)==b, 1, 'last'), col.Resp1)=4;  % ct
    datanext{2}(  find(datanext{2}(:, col.Block)==b, 1, 'last'), col.TrialValid)=0;
end
datanext{1}(end, col.TrialValid)=0;
datanext{2}(end, col.TrialValid)=0;

% CHOICE regressors
[variables c] =  i_choice_conditions(datanext{1}(datanext{1}(:, col.TrialValid)==1, :), [], 'cF_', {'Accept';'Reject';'Explore'}, col, c ); % Conflict
[variables c] =  i_choice_conditions(datanext{2}(datanext{2}(:, col.TrialValid)==1, :), variables,  'ct_', {'NoBomb';'Bomb';'Explore'}, col, c ); % Control

% RL VARIABLE regressors
[variables c] = i_RLvariables_Competing(datastruc{1}, 'cF_', RLvariables, variables, col, c ); % Conflict
[variables c] = i_RLvariables_Competing(datastruc{2}, 'ct_', RLvariables, variables,col,  c ); % Control

% EXPLORE information regressors
[variables c] = i_ExploredInfo(datastruc{3}, {'Entropy';}, variables, col, c );

% NO-INTEREST regressors
[variables c] = i_nointerest(datastruc{3}, {'Error'; 'Motor'; 'OutcomePresented'}, variables, col, c );

end

