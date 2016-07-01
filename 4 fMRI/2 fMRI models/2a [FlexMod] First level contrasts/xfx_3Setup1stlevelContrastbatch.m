function [matlabbatch] = xfx_3Setup1stlevelContrastbatch(where,log, ws, s, ConInstruc)
% [con] = f_setup1stlevelcontrastbatch(where,log, ws, s,con)
% Set up batch (all contrasts together), single subjecgt only
% 
% ConInstruc: Col 1= Contrast name, Col 3= Regressor Weights vector
%
% ------------------------------------------------------------------------------------

% Copy folder?
if isdir(ws.where_contrasts)==0; disp('Copying folder for contrasts . .'); copyfile([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Estimated'],ws.where_contrasts); disp('Done')
elseif s==1; input('Folder with the target name (onsets + contrasts type) already exists. Assume correct?   '); 
end

% Set up up batch
matlabbatch{1}.spm.stats.con.spmmat = cellstr([ws.where_contrasts 'SPM.mat' ]);
matlabbatch{1}.spm.stats.con.delete = 1;
for c=1:size(ConInstruc,1)
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =ConInstruc{c,1};
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=ConInstruc{c,3};
end

% Run batch
disp('Running contrasts -----------------')
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);

end

