% Find the models with (a) fewer parameters that (b) do better than the
% existing fits
clear all, close all hidden, clc;

cd('/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/Simulation/1 Paramfit generative models/')
load('GenerativeFit cF (09-Dec-2014).mat');  bestbic=11088; best_npars=6; % cF
% load('GenerativeFit ct (09-Dec-2014).mat');  bestbic=11475; best_npars=6; % ct

%%

error


r_res=r_res(cell2mat(r_res(:, 3))<bestbic,:); % Better than 
for m=1:size(r_res,1)
    r_res{m,5}=size(r_res{m,2},2)-3;
end
r_res=r_res(cell2mat(r_res(:, 5))<best_npars,:);

%%

details.cF_top10={'bpjm10_fow';'bpjm16_feow';'bpjm08_fuw';'bpjm14_feuw';'bpmi14_feuw';'bpjm15_fevw';'bpmi08_fuw';'bjm16_feow';'bpm16_feow';'bpmi07_fe';};
details.cF_betterthantargetbic_fewerpars={'bpm08_fuw';'bjm08_fuw';'bpm10_fow';'bjm10_fow';'bpm07_fe';'bpmi02_f';'bm16_feow';'bmi08_fuw';'bm08_fuw';'bmi07_fe';'bm14_feuw';'bpj10_fow';'bpi08_fuw';'bmi10_fow';}; 
details.cF_hierar=unique([details.cF_top10; details.cF_betterthantargetbic_fewerpars]);
details.cF_hierar= {'bjm08_fuw';'bjm10_fow';'bjm16_feow';'bm08_fuw';'bm14_feuw';'bm16_feow';'bmi07_fe';'bmi08_fuw';'bmi10_fow';'bpi08_fuw';'bpj10_fow';'bpjm08_fuw';'bpjm10_fow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow';'bpm07_fe';'bpm08_fuw';'bpm10_fow';'bpm16_feow';'bpmi02_f';'bpmi07_fe';'bpmi08_fuw';'bpmi14_feuw';};
details.ct_top10={'bpmi11_euw';'bpjm11_euw';'bjm11_euw';'bmi11_euw';'bpjm04_uw';'bpm11_euw';'bm11_euw';'bpjm13_eow';'bpm13_eow';'bpjm12_evw';};
details.ct_betterthantargetbic_fewerpars={'bjm11_euw';'bmi11_euw';'bpjm04_uw';'bpm11_euw';'bm11_euw';'bpm13_eow';'bpmi04_uw';'bpm04_uw';'bpjm03_e';'bpm12_evw';}; 
details.ct_hierar=unique([details.ct_top10;  details.ct_betterthantargetbic_fewerpars]);
details.ct_hierar={'bjm11_euw';'bm11_euw';'bmi11_euw';'bpjm03_e';'bpjm04_uw';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpm04_uw';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpmi04_uw';'bpmi11_euw';};




