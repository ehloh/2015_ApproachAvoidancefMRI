function [Con log] = f3_1_Event(where, log)
% [Con log] = f3_1_Event(where, log)
% From flexible 1st level model, apply contrasts to every single trial
%
% Unlike other flexmod contrasts, this model just weights every
% non-nuisance regressor in the first-level model (1st derivative only).
% Nuisance regressors assumed prefixed n_;
%
%
%
% WARNING !! ########################################################
% 
% Because of the high number of contrasts + that every single contrast includes a representation of the entire matrix (x regressors * y scans)
% the SPM.mat variable becomes prohibitively large and cannot be saved as SPM runs through each contrast. To get around this issue, 
% one has to go directly into the function spm_contrasts, and either:
% 
%     (a) change the command used to save variables
%             save('SPM.mat', 'SPM', '-V6') --> save('SPM.mat', 'SPM', '-V7.3')
%             (but you may have to go through all the associated functions and change the save function repeatedly. 
%               failure to change back is non-fatal - v7.3 saving is merely slower. Also, SPM12 should not produce this error at all)
%               
% 	(b) remove data that you will not be using (chosen option here!)
%             these contrasts are used only to extract betas from, and are never fed into further processing in SPM.
%             therefore we can remove the variables SPM.xCon.X0 and SPM.xCon.X1o, since these are large data-holding
%             variables. To remove them, find the following lines of code in spm_contrasts, and use rmfield as follows (3 extra lines of code)
% 
%                     % place xCon back in SPM
%                     %--------------------------------------------------------------------------
%                     warning('MESSING UP WITH CONTRASTS');     % new lines of code
%                     xCon(end).X0=[];     % new lines of code
%                     xCon(end).X1o=[];      % new lines of code


%                     SPM.xCon = xCon;
%                     % ---------------------------------------------------------------------------------------
% 
%             PLEASE remember to remove these lines of extra code after you are done running these contrasts, as they will not work
%             normally in other SPM analyses otherwise. 
%
% ######################################################################                    
                    
if  log.execute_contrasts==1;
    edit spm_contrasts.m
    input('Running FL contrasts for this model requires that spm_contrasts.m is altered (see model dox for detail). Continue?   ');
end

% Folders
log.FLest=['2 First level' filesep log.onsetsmodel ' Estimated' filesep];
log.FLcon=['2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep];

RegLists=cell(log.n_subjs,1); Con.ConNames=cell(log.n_subjs,1);
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')  ---------'])
    ws.FLcon=[where.data_brain filesep log.subjects{s} filesep log.FLcon]; c=1;
    if isdir(ws.FLcon)==0 ; 
        ws.old=[where.data_brain filesep log.subjects{s} filesep log.FLest]; ws.new=ws.FLcon;
        eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
    end
    
    % Load edited list of regressors
    ws.s=load([where.data_brain filesep log.subjects{s} filesep log.FLcon  'SPM.mat']);
    ws.f_whatderiv=@(x)x(length(x)-2:length(x)); % Remove derivatives
    w.a=cellfun(ws.f_whatderiv, ws.s.SPM.xX.name','UniformOutput',0);
    ws.s.SPM.xX.name(strcmp(w.a,'(2)'))={'null'}; ws.s.SPM.xX.name(strcmp(w.a,'(3)'))={'null'};
    ws.f_readrealname=@(x)x(7:length(x)-6); % Simplify names
    RegLists{s}=cellfun(ws.f_readrealname,ws.s.SPM.xX.name,'UniformOutput',0)';
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



disp(' ##############################################')
disp('DONE with FL contrasts')
disp('If spm_contrasts has been adapted, REMOVE ALTERATIONS!!  ');

end

