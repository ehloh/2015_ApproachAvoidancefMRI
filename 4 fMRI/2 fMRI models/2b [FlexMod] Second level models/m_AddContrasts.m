function [ matlabbatch ] = m_AddContrasts(where, log, Des)
% [ matlabbatch ] = m_AddContrasts(where, log, Des)
%  Add contrasts to the (already specified & estimated) 2nd level model.
%  e.g. Can add sensible contrasts to an N-way factorial ANOVA
%
% Des.Directory:         Results folder (where's the 2nd level SPM.mat?)
%                                (Unlike other comparisons, this one doesn't create new folder)
% Des.ConInstruc:       Instructions for contrasts to be added 
%                                 (Col 1= Name, Col 2=Weights, Col 3=Contrast tyle)
%                                 - Col 3: T contrast=1, F contrast = 2
% Des.DeleteOldCon:   Delete existing contrasts from model? (1=Yes, 0=No)
%
% ------------------------------------------------------------------------------------

% Execute to debug: m=2; Des=Con.SecondLevelAnalysis{m,3};

% Compile instructions
% Des.ConInstruc={'Negative effect of Similarity' [-1 -1 1 1];
%             'Negative effect of Valence' [-1 1 -1 1];
%             'Negative Interaction: Similarity x Valence' [-1 1 1 -1];
%             'ME Sim_effect'          [1 1 0 0]; % Main effects (against baseline)
%             'ME Dis_effect'          [0 0 1 1];
%             'ME Val_effect'          [1 0 1 0];
%             'ME Neu_effect'         [0 1 0 1];
%             'SR_effect'                 [1 0 0 0];  % SR
%             'SR_vs_others'          [3 -1 -1 -1];
%             'SR_vs_SN'               [1 -1 0 0];
%             'SN_effect'                 [0 1 0 0];  % SN
%             'SN_vs_others'          [-1 3 -1 -1];
%             'SN_vs_SR'               [-1 1 0 0];
%             'DR_effect'                 [0 0 1 0]; % DR
%             'DR_vs_others'          [-1 -1 3 -1];
%             'DR_vs_DN'               [0 0 1 -1];
%             'DN_effect'                 [0 0 0 1]; % DN
%             'DN_vs_others'          [-1 -1 -1 3];
%             'DN_vs_DR'               [0 0 -1 1]};

% Compile batch
matlabbatch{1}.spm.stats.con.spmmat =  {[Des.Directory filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = Des.DeleteOldCon;
for c=1:size(Des.ConInstruc,1)
    switch Des.ConInstruc{c,3}
        case 1
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = Des.ConInstruc{c,1};
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = Des.ConInstruc{c,2};
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
        case 2
            matlabbatch{1}.spm.stats.con.consess{c}.fcon.name = Des.ConInstruc{c,1};
            matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec = {Des.ConInstruc{c,2}};
            matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none';
        otherwise
            error('Invalid contrast type (specifying new contrasts. Only t (1) or F (2) contrasts allowed.')
    end
end

end

