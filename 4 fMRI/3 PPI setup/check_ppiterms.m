clear all; close all hidden; clc

subjects={'p01_GV';'p02_YY'};
ppis={
%     'cellFam_roi_HPC_L_c3tc_psy_cF_Accept';
%     'cellFam_roi_HPC_L_c3tc_psy_cF_Explore';
%     'cellFam_roi_HPC_L_c3tc_psy_cF_Reject';
%     'cellFam_roi_HPC_L_c3tc_psy_ct_Bomb';
%     'cellFam_roi_HPC_L_c3tc_psy_ct_Explore';
%     'cellFam_roi_HPC_L_c3tc_psy_ct_NoBomb';
    'cellFam_roi_HPC_L_c7tct_psy_cF_Accept';
%     'cellFam_roi_HPC_L_c7tct_psy_cF_Explore';
    'cellFam_roi_HPC_L_c7tct_psy_cF_Reject';
%     'cellFam_roi_HPC_L_c7tct_psy_ct_Bomb';
%     'cellFam_roi_HPC_L_c7tct_psy_ct_Explore';
%     'cellFam_roi_HPC_L_c7tct_psy_ct_NoBomb';
};

%
where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\';

%%

dd=cell(length(subjects),1);

for s=1:length(subjects)
    disp(['Subject ' num2str(s) '  (' subjects{s} ') -----------------------------------------------------']);
    ws.wherePPIs=[where.data_brain subjects{s} filesep '2 First level' filesep 'm_c3_CompeteFull_XUVPEN PPIs' filesep];
    dd{s}=cell(length(ppis),1);
    
    for p=1:length(ppis)
        disp(['PPI #' num2str(p) ':  ' ppis{p}])
        dd{s}{p}=load([ws.wherePPIs ppis{p} filesep 'PPI_' ppis{p}]);
        disp(dd{s}{p}.PPI)
    end
end