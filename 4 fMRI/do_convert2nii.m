clear all; clc

roiname='Insula_L';
whereroi='C:\Users\eloh\Desktop\2 [Explore]\3 Anatomical\2 A priori ROIs\3 Unilateral Native Space';

%%

v=spm_vol([whereroi filesep roiname '.img']);
ima=spm_read_vols(v);
v.fname=[whereroi filesep roiname '.nii'];
spm_write_vol(v,ima);


cd(whereroi)