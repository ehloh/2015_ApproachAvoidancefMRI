function [RegList] = xfx_1LoadSimpleReglist(where, log,s)
% [RegList] = f_loadeditedreglist(where, log,s)
%   Load list of included regressors associated with this model
%   Edits to list: (1) Remove derivatives (blank)
%                      (2) Simplify names of regressors
%
% ------------------------------------------------------------------------

ws=load([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']); 

% Remove derivatives
f_whatderiv=@(x)x(length(x)-2:length(x));
w.a=cellfun(f_whatderiv, ws.SPM.xX.name','UniformOutput',0);
ws.SPM.xX.name(strcmp(w.a,'(2)'))={'null'}; ws.SPM.xX.name(strcmp(w.a,'(3)'))={'null'};

% Simplify names
f_readrealname=@(x)x(7:length(x)-6);
RegList=cellfun(f_readrealname,ws.SPM.xX.name,'UniformOutput',0)';

end

