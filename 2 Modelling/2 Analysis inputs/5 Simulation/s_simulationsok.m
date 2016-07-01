% Check fitted parameters from simulated data: do we recover the parameters
% that were used to create the data in the first place?
clear all; close all; clc
% clear all; clc
for o1=1:1 % Model families
    details.modelfams.b={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow';'b17_yw';'b18_fyw';'b19_eyw';'b20_feyw'; 'b21_kw';'b22_fkw';'b23_ekw';'b24_fekw'};
    details.modelfams.bp={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow';'bp17_yw';'bp18_fyw';'bp19_eyw';'bp20_feyw'; 'bp21_kw';'bp22_fkw';'bp23_ekw';'bp24_fekw'};
    details.modelfams.bm={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow';'bm17_yw';'bm18_fyw';'bm19_eyw';'bm20_feyw'; 'bm21_kw';'bm22_fkw';'bm23_ekw';'bm24_fekw'};
    details.modelfams.bpm ={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow';'bpm17_yw';'bpm18_fyw';'bpm19_eyw';'bpm20_feyw'; 'bpm21_kw';'bpm22_fkw';'bpm23_ekw';'bpm24_fekw'};
    %
    details.modelfams.bi={'bi01';'bi02_f';'bi03_e';'bi04_uw';'bi05_vw';'bi06_ow';'bi07_fe';'bi08_fuw';'bi09_fvw';'bi10_fow';'bi11_euw';'bi12_evw';'bi13_eow';'bi14_feuw';'bi15_fevw';'bi16_feow';'bi17_yw';'bi18_fyw';'bi19_eyw';'bi20_feyw'; 'bi21_kw';'bi22_fkw';'bi23_ekw';'bi24_fekw'};
    details.modelfams.bpi ={'bpi01';'bpi02_f';'bpi03_e';'bpi04_uw';'bpi05_vw';'bpi06_ow';'bpi07_fe';'bpi08_fuw';'bpi09_fvw';'bpi10_fow';'bpi11_euw';'bpi12_evw';'bpi13_eow';'bpi14_feuw';'bpi15_fevw';'bpi16_feow';'bpi17_yw';'bpi18_fyw';'bpi19_eyw';'bpi20_feyw'; 'bpi21_kw';'bpi22_fkw';'bpi23_ekw';'bpi24_fekw'};
    details.modelfams.bmi ={'bmi01';'bmi02_f';'bmi03_e';'bmi04_uw';'bmi05_vw';'bmi06_ow';'bmi07_fe';'bmi08_fuw';'bmi09_fvw';'bmi10_fow';'bmi11_euw';'bmi12_evw';'bmi13_eow';'bmi14_feuw';'bmi15_fevw';'bmi16_feow';'bmi17_yw';'bmi18_fyw';'bmi19_eyw';'bmi20_feyw'; 'bmi21_kw';'bmi22_fkw';'bmi23_ekw';'bmi24_fekw'};
    details.modelfams.bpmi ={'bpmi01';'bpmi02_f';'bpmi03_e';'bpmi04_uw';'bpmi05_vw';'bpmi06_ow';'bpmi07_fe';'bpmi08_fuw';'bpmi09_fvw';'bpmi10_fow';'bpmi11_euw';'bpmi12_evw';'bpmi13_eow';'bpmi14_feuw';'bpmi15_fevw';'bpmi16_feow';'bpmi17_yw';'bpmi18_fyw';'bpmi19_eyw';'bpmi20_feyw'; 'bpmi21_kw';'bpmi22_fkw';'bpmi23_ekw';'bpmi24_fekw'};
    details.ct_modnums=[1 3 4 5 6 11 12 13 17 19]; % true models for ct
    %
    details.modelfams.bj={'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow';'bj17_yw';'bj18_fyw';'bj19_eyw';'bj20_feyw'; 'bj21_kw';'bj22_fkw';'bj23_ekw';'bj24_fekw'};
    details.modelfams.bpj={'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow';'bpj17_yw';'bpj18_fyw';'bpj19_eyw';'bpj20_feyw'; 'bpj21_kw';'bpj22_fkw';'bpj23_ekw';'bpj24_fekw'};
    details.modelfams.bjm={'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow';'bjm17_yw';'bjm18_fyw';'bjm19_eyw';'bjm20_feyw'; 'bjm21_kw';'bjm22_fkw';'bjm23_ekw';'bjm24_fekw'};
    details.modelfams.bpjm={'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm09_fvw';'bpjm10_fow';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow';'bpjm17_yw';'bpjm18_fyw';'bpjm19_eyw';'bpjm20_feyw'; 'bpjm21_kw';'bpjm22_fkw';'bpjm23_ekw';'bpjm24_fekw'};
end

req.tasktype=2;
%
req.generative_fittype={'Fit' '(18-Jul-2014) all20iter'};   % Spex for simulated data  ####         Basic
% req.generative_fittype={'Fit' '(02-Dec-2014) j'};   % Spex for simulated data  ####        j 
% req.generative_fittype={'Fit' '(05-Mar-2015) y k'};   % Spex for simulated data  ####               y k 

req.whichmod=details.modelfams.bpmi{11};

for o1=1:1 % Set up 
%     where.simfol='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/5 Simulation';
    where.where='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
    where.simfol=[where.where filesep '2 Analysis inputs' filesep '5 Simulation'];
    switch req.tasktype
        case 1; req.task='cF';
        case 2; req.task='ct';
    end
    details.generativeparfile=['Generative' req.generative_fittype{1} ' ' req.task ' ' req.generative_fittype{2}];
    
end

%% Load parameters for comparison 

% [Generative pars] Load details of requested model (gfit)
% gfit=load([where.simfol filesep '1 Paramfit generative models' filesep  details.generativeparfile '.mat']);
gfit=load([where.where filesep '2 Analysis inputs' filesep 'Det Jpower' filesep  'res_fitmodels_' req.task ' ' req.generative_fittype{2} '.mat']);
%
gfit.modnum=find(strcmp(gfit.r_res(:,1),req.whichmod)); if isempty(gfit.modnum); error('Could not find requested model in generative param file!'); end
gfit.par=gfit.r_res{gfit.modnum,2}(:, 4:end);
gfit.parnames=gfit.r_res{gfit.modnum,1}; gfit.parnames(   strfind(gfit.parnames, '_')-2  : strfind(gfit.parnames, '_') )=[]; gfit.parnames(strfind(gfit.parnames, 'u'))=[]; gfit.parnames(strfind(gfit.parnames, 'o'))=[]; gfit.parnames(strfind(gfit.parnames, 'v'))=[];

% [Refit pars] Load details of re-fit requested models (rfit)
rfit=load([where.simfol filesep '3 Paramfits refit' filesep   'Fitsim ' req.task ' - ' req.whichmod '.mat']);
% rfit=load([where.simfol filesep '3 Paramfits refit/Problematic' filesep   'Fitsim ' req.task ' - ' req.whichmod '.mat']);
% rfit=load([where.simfol filesep '3 Paramfits refit/Done' filesep   'Fitsim ' req.task ' - ' req.whichmod '.mat']);
rfit.modnum=find(strcmp(rfit.r_res(:,1),req.whichmod));
rfit.par=rfit.r_res{rfit.modnum,2}(:, 4:end);
rfit.parnames=rfit.r_res{rfit.modnum,1}; rfit.parnames(   strfind(rfit.parnames, '_')-2  : strfind(rfit.parnames, '_') )=[]; rfit.parnames(strfind(rfit.parnames, 'u'))=[]; rfit.parnames(strfind(rfit.parnames, 'o'))=[]; rfit.parnames(strfind(rfit.parnames, 'v'))=[];

%% Compare?
npar=size(gfit.par,2);


figure('color','w','Position', [500 200 npar*200 150],'Name',   ['[' req.task '] '  req.whichmod], 'NumberTitle', 'off');
subplot(1,npar+1,1)
scatter(gfit.r_res{gfit.modnum,2}(:,2), rfit.r_res{rfit.modnum,2}(:,2))
lsline,  title('nLL')
for p=1:npar
    subplot(1,npar+1,p+1)
    scatter(gfit.par(:,p),  rfit.par(:,p))
    lsline, title(rfit.parnames(p))
    
    axis square
end


% %%
% openvar gfit.par
% openvar rfit.par



