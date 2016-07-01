% Count voxels of of overlap (and other things) 
clear all, close all hidden, clc



rois={
    
%     % [c13 ROIs ] ############################################
%     'G:\2 [Explore]\rois\c13 battery and amyg\HPC_aL_stc.nii';
%     'G:\2 [Explore]\rois\c13 battery and amyg\HPC_aR_stc.nii';
%     'G:\2 [Explore]\rois\c13 battery and amyg\HPC_aL_sc.nii';  % ME Choice 
%     'G:\2 [Explore]\rois\c13 battery and amyg\HPC_aR_sc.nii';
%     % Unmasked MEC HPC clusters 
%     'G:\2 [Explore]\2 Second level results s4Ants\m_c13_ChoiceFull_ULPEN_Basic\choice_2x3\ROI\HPC\HPC_aL_c.nii'
%     'G:\2 [Explore]\2 Second level results s4Ants\m_c13_ChoiceFull_ULPEN_Basic\choice_2x3\ROI\HPC\HPC_aR_c.nii'
    
    
%     % [v3g ROIs] ############################################
    'G:\2 [Explore]\rois\v3g HPC\HPC_aL_tv001.nii';
    'G:\2 [Explore]\rois\v3g HPC\HPC_aR_tv001.nii';

%     HPC vs AMYGDALA --- 
%     'G:\2 [Explore]\4a AntsAnatomical\2 A priori ROIs\HPC_bilat.nii'; 
%     'G:\2 [Explore]\4a AntsAnatomical\2 A priori ROIs\Amygdala_bilat.nii'; 
    
}; 


%% 

rvol= cell(size(rois,1),1); 
for r=1:size(rois,1)
    rvol{r} = spm_read_vols(spm_vol(rois{r} ));
end


sum( (rvol{1}(:)~=0))
sum( (rvol{2}(:)~=0))

disp(['No. overlapping voxels = '  num2str(sum( (rvol{1}(:)~=0).* (rvol{2}(:)~=0)))])
disp(['= ' num2str(100* sum( (rvol{1}(:)~=0).* (rvol{2}(:)~=0)) ./ sum( (rvol{1}(:)~=0)) ) '% of ROI #1'])




