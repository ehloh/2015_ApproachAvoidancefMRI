% Check that the new parameter functions match up with the old
clear all; clc
where='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
cd(where); addpath(where); addpath([where '1 Value functions Det Jpower']);



% Base setup
et=(repmat(1:6,1, 6)/6)';
ntok=2*6*sortrows(et);
%
et=power(et, 0.2);


%% Conflict task

% New function: fcf_changeEnvThreat
% [ ov] = fcf_changeEnvThreat(et, ntok, []);
[ ov] = fct_changeEnvThreat(et, ntok, []);

% ORIGINAL: fpar_conflict
col.EnvThreat=6;
col.NTokens=2;
col.pLoss=1;
col.Entropy=4;
col.VExplore=5;
col.EV=10;
col.EntropyNTok=7;
data(:, [col.EnvThreat  col.NTokens])= [et ntok];
[ data ] = fpar_conflict( data, col);
% [ data ] = fpar_control( data, col);


% Compare
empty=nan(size(data,1),1); 
newd=[ov.pLoss  ov.NTok empty ov.Entropy ov.VExplore ov.EnvThreat  ov.EntropyNTok empty empty  ov.EV];
sum(data-newd)


%%

