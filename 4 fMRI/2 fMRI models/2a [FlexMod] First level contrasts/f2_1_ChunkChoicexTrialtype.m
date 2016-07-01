function [Con log] = f2_1_ChunkChoicexTrialtype(where, log)
% [Con log] = f2_1_ChunkChoicexTrialtype(where, log)
% From flexible 1st level model, apply contrasts to derive full task design 
%                                   (Task x Choice x Trial type)
%
% Unlike other flexmod contrasts, this model just weights every
% non-nuisance regressor in the first-level model (1st derivative only).
% Nuisance regressors assumed prefixed n_;
%
% ---------------------------------------------------------------------------------------

log.FLest=['2 First level' filesep log.onsetsmodel ' Estimated' filesep];
log.FLcon=['2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep];

RegLists=cell(log.n_subjs,1); Con.ConNames=cell(log.n_subjs,1);
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')  ---------'])
    ws.FLcon=[where.data_brain filesep log.subjects{s} filesep log.FLcon]; c=1;
    if isdir(ws.FLcon)==0 ;  copyfile([where.data_brain filesep log.subjects{s} filesep log.FLest],  ws.FLcon); end
    [RegLists{s}] = xfx_1LoadSimpleReglist(where, log, s); % Load (edited) list of regressors
    if strcmp(RegLists{s}(end), 'co');  RegLists{s}(end)=[]; end  % Get rid of constant
    
    % Set up contrasts: Weight every (non-nuisance) regressor
    matlabbatch{1}.spm.stats.con.spmmat = cellstr([ws.FLcon 'SPM.mat']);
    matlabbatch{1}.spm.stats.con.delete = 1;
    for r=1:size(RegLists{s},1)
        if isempty(RegLists{s}{r})==0 && strcmp(RegLists{s}{r}(1:2), 'n_')==0
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=RegLists{s}{r};
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=zeros(1,size(RegLists{s},1));
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec(r)=1;
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
            %
            Con.ConNames{s}{c,1}=matlabbatch{1}.spm.stats.con.consess{c}.tcon.name;
            Con.ConNames{s}{c,2}=r;
            c=c+1;
        end
    end
    if s==1;  disp('Contrasts for subjects #1: '); Con.matlabbatch=matlabbatch; disp(Con.ConNames{s}); input('Check for no nuisance regressors weighted. OK to continue?   '); end
    
    disp(['Running ' num2str(c-1) ' contrasts ---'])
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];ws=[];
end

% Output
Con.ConInstruc='Every non-nuisance regressor is weighted (1st derivative only)';
Con.SearchCriterias='1st derivatives of conditions that do not start with "n_"';
Con.RegLists=RegLists;

%% Second-level analyses

Con.SecondLevelAnalysis{1,1}='No 2nd level analyses planned for this FL model';
if  log.execute_contrasts~=1; error(Con.SecondLevelAnalysis{1,1}); end

end

