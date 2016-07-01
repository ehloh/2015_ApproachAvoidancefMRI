% Clean up runs (multiple options)
clear all; close all hidden; clc
er 

request.clearpartialfit=0;
request.whatsmissing=1;
request.removemods=0;
remove.choosemods=0;
request.checkparams_inspace =0;

for o=1:1 % General setup 
 % Folders
w=pwd;
if strcmp(w(1), '/')==0;  path(pathdef); addpath('D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models');
else addpath('/Users/EleanorL/Dropbox/sandisk/4 Explore experiment/3 Analysis/4 Fit computational models');
end

% Task names
tasknames={'cF';'ct'}; 
end
edit f_newname

%% % Take a partialfit or a workspace from a crashed fit and clean it up, turn  it into a proper results file

if request.clearpartialfit
    wherefile='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical';
    infilename= 'res_hierarfitmodels_cF (16-Apr-2015) almost';
    lastokrun=384;
    
    for o=1:1
        % ---------------------------------------------------------------------------------------------------
        ishierar=1;
        cd(wherefile);  load([infilename '.mat']); openvar r_res;
        input(['Check r_res:  Last ok run is :' num2str(lastokrun) '  ?   ' ]);
        
        r_res=r_res(1:lastokrun,:);
        for m=1:size(r_res,1); if isempty(r_res{m,2})==1;  r_res{m,3}=nan; end; end % unconverged
        if ishierar ==1;   [ r_res , r_iterations, r_iterd, missingres] = f_hierarmatch(r_res , r_iterations, r_iterd);  openvar missingres, disp('Models that need rerunning are listed in missingres'); end
        %
        rr={};  ri={}; rd={};  modnum=1; dm={}; 
        for m=1:size(r_res,1)
            ri(m,:)=r_iterations(strcmp(r_iterations(:,1) , r_res{m,1}),:);
            rd(m,:)=r_iterd(strcmp(r_iterd(:,1) , r_res{m,1}),:);
            dm(m,:)= details.models(strcmp(details.models(:, 1), r_res{m,1}), :);
        end
        details.models=dm; details.whichmodels=details.models(:,1); details.n_models=length(details.whichmodels);
        r_iterations=ri; if ishierar ==1;   r_iterd=rd; end
        for m=1:length(details.whichmodels) % CHECK
            if  strcmp(details.whichmodels{m,1}, r_res{m,1})+  strcmp(r_res{m,1}, r_iterd{m,1}) +  strcmp(r_res{m,1}, r_iterations{m,1})~=3
                error([r_res{m,1}   '     : Model matchup FAILS']);
            end
        end
        r_res=sortrows(r_res,3);
        
        % SAVE ---------------------------------------------------------------------------------------------------
        resfilename= f_newname([infilename ' CLEANpartial'], pwd); input(['Save as     ' resfilename  '  ?']);
        if ishierar ==1; save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');   % Save a hierarfit
        else  save(resfilename, 'details', 'r_iterations','r_res', 'rc','errorlog');  % Saved a fixed fit
        end
        
    end
end


%% Identify models that are still missing
% Checks for what other models need to be run. Assumes list of models here
% is the fulll list.
% Count only models where r_res, r_iterd and r_iterations match up.

if request.whatsmissing
    ishierar=0;
    
    wherefile='/Users/EleanorL/Dropbox/sandisk/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs';
    infilename= 'res_fitmodels_cF (09-Apr-2015) all';
%     wherefile='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical';
%     infilename= 'res_hierarfitmodels_ct (16-Apr-2015) all';
        
    % ---------------------------------------------------------------------------------------------------
    cd(wherefile);  load([infilename '.mat']); 
    
    for o=1:1 % All models (assume cF)
        details.ct_modnums=[1 3 4 5 6 11 12 13 17 19 21 23]; % true models for c
        if details.tasktype==1;  which=1:24; else which=details.ct_modnums; end
        
        modfams.b={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow';'b17_yw';'b18_fyw';'b19_eyw';'b20_feyw'; 'b21_kw';'b22_fkw';'b23_ekw';'b24_fekw'};
        modfams.bp={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow';'bp17_yw';'bp18_fyw';'bp19_eyw';'bp20_feyw'; 'bp21_kw';'bp22_fkw';'bp23_ekw';'bp24_fekw'};
        modfams.bm={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow';'bm17_yw';'bm18_fyw';'bm19_eyw';'bm20_feyw'; 'bm21_kw';'bm22_fkw';'bm23_ekw';'bm24_fekw'};
        modfams.bpm ={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow';'bpm17_yw';'bpm18_fyw';'bpm19_eyw';'bpm20_feyw'; 'bpm21_kw';'bpm22_fkw';'bpm23_ekw';'bpm24_fekw'};
        %
        modfams.bi={'bi01';'bi02_f';'bi03_e';'bi04_uw';'bi05_vw';'bi06_ow';'bi07_fe';'bi08_fuw';'bi09_fvw';'bi10_fow';'bi11_euw';'bi12_evw';'bi13_eow';'bi14_feuw';'bi15_fevw';'bi16_feow';'bi17_yw';'bi18_fyw';'bi19_eyw';'bi20_feyw'; 'bi21_kw';'bi22_fkw';'bi23_ekw';'bi24_fekw'};
        modfams.bpi ={'bpi01';'bpi02_f';'bpi03_e';'bpi04_uw';'bpi05_vw';'bpi06_ow';'bpi07_fe';'bpi08_fuw';'bpi09_fvw';'bpi10_fow';'bpi11_euw';'bpi12_evw';'bpi13_eow';'bpi14_feuw';'bpi15_fevw';'bpi16_feow';'bpi17_yw';'bpi18_fyw';'bpi19_eyw';'bpi20_feyw'; 'bpi21_kw';'bpi22_fkw';'bpi23_ekw';'bpi24_fekw'};
        modfams.bmi ={'bmi01';'bmi02_f';'bmi03_e';'bmi04_uw';'bmi05_vw';'bmi06_ow';'bmi07_fe';'bmi08_fuw';'bmi09_fvw';'bmi10_fow';'bmi11_euw';'bmi12_evw';'bmi13_eow';'bmi14_feuw';'bmi15_fevw';'bmi16_feow';'bmi17_yw';'bmi18_fyw';'bmi19_eyw';'bmi20_feyw'; 'bmi21_kw';'bmi22_fkw';'bmi23_ekw';'bmi24_fekw'};
        modfams.bpmi ={'bpmi01';'bpmi02_f';'bpmi03_e';'bpmi04_uw';'bpmi05_vw';'bpmi06_ow';'bpmi07_fe';'bpmi08_fuw';'bpmi09_fvw';'bpmi10_fow';'bpmi11_euw';'bpmi12_evw';'bpmi13_eow';'bpmi14_feuw';'bpmi15_fevw';'bpmi16_feow';'bpmi17_yw';'bpmi18_fyw';'bpmi19_eyw';'bpmi20_feyw'; 'bpmi21_kw';'bpmi22_fkw';'bpmi23_ekw';'bpmi24_fekw'};
        %
        modfams.bj={'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow';'bj17_yw';'bj18_fyw';'bj19_eyw';'bj20_feyw'; 'bj21_kw';'bj22_fkw';'bj23_ekw';'bj24_fekw'};
        modfams.bpj={'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow';'bpj17_yw';'bpj18_fyw';'bpj19_eyw';'bpj20_feyw'; 'bpj21_kw';'bpj22_fkw';'bpj23_ekw';'bpj24_fekw'};
        modfams.bjm={'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow';'bjm17_yw';'bjm18_fyw';'bjm19_eyw';'bjm20_feyw'; 'bjm21_kw';'bjm22_fkw';'bjm23_ekw';'bjm24_fekw'};
        modfams.bpjm={'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm09_fvw';'bpjm10_fow';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow';'bpjm17_yw';'bpjm18_fyw';'bpjm19_eyw';'bpjm20_feyw'; 'bpjm21_kw';'bpjm22_fkw';'bpjm23_ekw';'bpjm24_fekw'};
        %
        modfams.bji={'bji01';'bji02_f';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji07_fe';'bji08_fuw';'bji09_fvw';'bji10_fow';'bji11_euw';'bji12_evw';'bji13_eow';'bji14_feuw';'bji15_fevw';'bji16_feow';'bji17_yw';'bji18_fyw';'bji19_eyw';'bji20_feyw'; 'bji21_kw';'bji22_fkw';'bji23_ekw';'bji24_fekw'};
        modfams.bpji={'bpji01';'bpji02_f';'bpji03_e';'bpji04_uw';'bpji05_vw';'bpji06_ow';'bpji07_fe';'bpji08_fuw';'bpji09_fvw';'bpji10_fow';'bpji11_euw';'bpji12_evw';'bpji13_eow';'bpji14_feuw';'bpji15_fevw';'bpji16_feow';'bpji17_yw';'bpji18_fyw';'bpji19_eyw';'bpji20_feyw'; 'bpji21_kw';'bpji22_fkw';'bpji23_ekw';'bpji24_fekw'};
        modfams.bjmi={'bjmi01';'bjmi02_f';'bjmi03_e';'bjmi04_uw';'bjmi05_vw';'bjmi06_ow';'bjmi07_fe';'bjmi08_fuw';'bjmi09_fvw';'bjmi10_fow';'bjmi11_euw';'bjmi12_evw';'bjmi13_eow';'bjmi14_feuw';'bjmi15_fevw';'bjmi16_feow';'bjmi17_yw';'bjmi18_fyw';'bjmi19_eyw';'bjmi20_feyw'; 'bjmi21_kw';'bjmi22_fkw';'bjmi23_ekw';'bjmi24_fekw'};
        modfams.bpjmi={'bpjmi01';'bpjmi02_f';'bpjmi03_e';'bpjmi04_uw';'bpjmi05_vw';'bpjmi06_ow';'bpjmi07_fe';'bpjmi08_fuw';'bpjmi09_fvw';'bpjmi10_fow';'bpjmi11_euw';'bpjmi12_evw';'bpjmi13_eow';'bpjmi14_feuw';'bpjmi15_fevw';'bpjmi16_feow';'bpjmi17_yw';'bpjmi18_fyw';'bpjmi19_eyw';'bpjmi20_feyw'; 'bpjmi21_kw';'bpjmi22_fkw';'bpjmi23_ekw';'bpjmi24_fekw'};
        %
        allmods=[
            modfams.b(which)
            modfams.bp(which)
            modfams.bm(which)
            modfams.bpm(which)
            modfams.bi(which)   %
            modfams.bpi(which)
            modfams.bmi(which)
            modfams.bpmi(which)
            modfams.bj(which)   %
            modfams.bpj(which)
            modfams.bjm(which)
            modfams.bpjm(which)
            modfams.bji(which)  %
            modfams.bpji(which)
            modfams.bjmi(which)
            modfams.bpjmi(which)
            ];
    end
    % Check for models that are missing modelspaces
    if ishierar ==1;   [ r_res , r_iterations, r_iterd, missingres] = f_hierarmatch(r_res , r_iterations, r_iterd);  openvar missingres, disp('Models that need rerunning are listed in missingres'); end
    %
    r_modthere=zeros(length(allmods),1);
    for m=1:length(allmods)
        if sum(strcmp(r_res(:,1), allmods{m}))==1 && isempty(r_res(strcmp(r_res(:,1), allmods{m}),2))==0
            r_modthere(m)=1;
        else r_modthere(m)=0;
        end
    end
    if sum(r_modthere)==length(allmods); disp('All mods there! :)'); 
    else missingmods=allmods(r_modthere==0); openvar missingmods, disp('Missing mods:'), disp(missingmods)
    end
end

%% Compile a resuts file with ONLY requested parameters

if remove.choosemods
    ishierar=1;
    wherefile='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical\1 Manuscript fits\v2_2';
    infilename= 'res_hierarfitmodels_cF (21-Jul-2014) top10';
    
    
    whichmods={'bpm16_feow'};
    
    % ---------------------------------------------------------------------------------------------------
    cd(wherefile);  load([infilename '.mat']);
    
    if ishierar ==1;   [ r_res , r_iterations, r_iterd, resdiscrep] = f_hierarmatch(r_res , r_iterations, r_iterd);  openvar resdiscrep, disp('Models that need rerunning are listed in resdiscrep'); end
    rr={}; ri={}; rd={}; dm={}; notthere={}; modnum=1;
    for m=1:length(whichmods)
        rr(modnum, :)=r_res(strcmp(r_res(:,1), whichmods{m}), :);
        ri(modnum, :)=r_iterations(strcmp(r_iterations(:,1), whichmods{m}), :);
        rd(modnum, :)=r_iterd(strcmp(r_iterd(:,1), whichmods{m}), :);
        dm(modnum, :)=details.models(strcmp(details.models(:,1), whichmods{m}), :);        
    end 
    details.models=dm;  details.whichmodels=details.models(:,1);  details.n_models=length(details.whichmodels); 
    r_res =rr;  r_iterations=ri;  r_iterd=rd; 
    r_res =sortrows(r_res ,3);
    
    
     % SAVE ---------------------------------------------------------------------------------------------------
     if ishierar ==1  ; resfilename=['res_hierarfitmodels_' tasknames{details.tasktype} ' (' date ') NEWcompiled']; % Save a hierarfit
    else  resfilename=['res_fitmodels_'  tasknames{details.tasktype} ' (' date ') NEWcompiled'];    % Saved a fixed fit
    end
    resfilename= f_newname(resfilename, pwd); input(['Save as     ' resfilename  '  ?']);
    if ishierar ==1; save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');   % Save a hierarfit
    else  save(resfilename, 'details', 'r_iterations','r_res', 'rc','errorlog');  % Saved a fixed fit
    end
end

%% Check all models: Are parameters values converted to true parameter space? If not, convert

if request.checkparams_inspace  % hierarchical fits only!!
    request.wherefile='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical';
    request.filename='res_hierarfitmodels_ct (16-Apr-2015) all'; 
    request.datafile='All data (09-May-2014)';

%     % SLACK
%     request.wherefile='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling\1 Inputs\2b Hierar fits\New';
%     request.filename='res_hierarfitmodels_cF (16-Apr-2015) v1 all';
%     request.datafile='All data v1 (01-Apr-2015)';
    
    % --------------------------------------------------------------------------------
    for o=1:1  % Set up  
        cd(request.wherefile), load([request.filename '.mat']);
    
        
        % MARC 
        where.scripts='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
        path(pathdef); addpath(where.scripts); addpath(genpath([where.scripts filesep '1 Value functions']))
        sd=load([where.scripts filesep '2 Analysis inputs' filesep request.datafile '.mat']);
        d_design=nan(6*6, 10);  d_design(:,  [sd.details.col.EnvThreat sd.details.col.NTokens])=[sortrows(repmat((1:6)',6,1))/6 2*repmat((1:6)',6,1)];  d_design(:,sd.details.col.Task)=details.tasktype;
        
%         % SLACK 
%         where.scripts='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling';
%         where.origscripts='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
%         path(pathdef); addpath(genpath([where.scripts filesep '2 Val fxn']))
%         addpath(where.scripts);  addpath(where.origscripts); 
%         sd=load([where.scripts filesep '1 Inputs' filesep request.datafile '.mat']);
%         d_design=nan(4*4, 10);  d_design(:,  [sd.details.col.EnvThreat sd.details.col.NTokens])=[sortrows(repmat((1:4)',4,1))/4 2*repmat((1:4)',4,1)];  d_design(:,sd.details.col.Task)=details.tasktype;
        
        n_subs=size(sd.subjdata,1);
        
        % Generate design
        tasks={'conflict';'control'};  eval(['[d_design] = fpar_' tasks{details.tasktype} '(d_design, sd.details.col);'])

    end    
    
    
    % (a) Compile list of models that need correction ---------------------------------------
    request.mods2correct=1;   % empty to Find
    %
    if isempty(request.mods2correct)
        for o=1:1 
            input('Identify models that need correction?  '); 
            d_matchl=nan(size(r_res,1),3);
            r_res=sortrows(r_res,1); details.whichmodels=sortrows(details.whichmodels);
            for m=1:size(r_res,1)
                wm.parnames=details.models{strcmp(details.models(:,1), r_res{m,1}),3};
                %         wm.subs= datasample(1:details.n_subjs, 3, 'Replace', false);
                wm.subs=[3 6 12];
                wm.spars=r_res{m,2}(wm.subs, 4:end);
                wm.sL=r_res{m,2}(wm.subs, 2);
                
                for ss=1:length(wm.subs)
                    ws.subname=details.subjects{wm.subs(ss)};
                    ws.transpar=f_transpar(wm.parnames, wm.spars(ss, :) , 'from');
                    
                    %             ws.transpar=wm.spars(ss, :) ;  disp('NOT ARS')
                    
                    ws.subdata=sd.subjdata{strcmp(sd.subjdata(:,1), ws.subname),2};
                    
                    if isempty(strfind(r_res{m,1}, 'p'))==0
                        [nll,pch]=f_nllsoftmax_lapse(ws.transpar,  {r_res{m,1}  ws.subdata details.fixedpar sd.details.col},[], 0 );
                    else [nll,pch]=f_nllsoftmax(ws.transpar,  {r_res{m,1}  ws.subdata details.fixedpar sd.details.col},[], 0 );
                    end
                    
                    d_newl(m,ss)=nll;
                    d_recl(m,ss)=wm.sL(ss);
                    if wm.sL(ss)~=nll;  d_matchl(m, ss)=0;   % Compare Ls
                    else d_matchl(m, ss)=1;
                    end
                    ws=[];
                end
                wm=[];
            end
            d_matchl=[r_res(:,1) num2cell(d_matchl)];
            request.mods2correct=d_matchl(find(abs(3-sum(cell2mat(d_matchl(:,2:4)),2))),1);
            openvar d_newl, openvar d_recl, openvar d_matchl, disp('Mods to correct:'),  disp(request.mods2correct); input('Continue?');
            openvar request.mods2correct
        end
    end
    
    for o=1:1 % Settings for correction

        % v1 cF all (res_hierarfitmodels_cF (16-Apr-2015) all)
%         request.mods2correct={
%             'bji01';'bji02_f';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji07_fe';'bji08_fuw';'bji09_fvw';'bji10_fow';'bji11_euw';'bji12_evw';'bji13_eow';'bji14_feuw';'bji15_fevw';'bji16_feow';'bji17_yw';'bji18_fyw';'bji19_eyw';'bji20_feyw';'bji21_kw';'bji22_fkw';'bji23_ekw';'bji24_fekw';
%             'bpji01';'bpji02_f';'bpji03_e';'bpji04_uw';'bpji05_vw';'bpji06_ow';'bpji07_fe';'bpji08_fuw';'bpji09_fvw';'bpji11_euw';'bpji13_eow';'bpji14_feuw';'bpji15_fevw';'bpji16_feow';'bpji17_yw';'bpji18_fyw';'bpji19_eyw';'bpji20_feyw';'bpji21_kw';
%             };
        
        
        % v1 ct all (res_hierarfitmodels_ct (16-Apr-2015) all)
        request.mods2correct={'bji01';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji11_euw';'bji12_evw';'bji13_eow';'bji17_yw';'bji19_eyw';'bji21_kw';'bji23_ekw';'bpji01';'bpji03_e';'bpji04_uw';};
        
        for o1=1:1  % SLACK 
%             
%             % v1 cF all (res_hierarfitmodels_cF (16-Apr-2015) v1 all)
%             request.mods2correct={'bi15_fevw';'bi22_fkw';'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow';'bj17_yw';'bj18_fyw';'bj19_eyw';'bj20_feyw';'bj21_kw';'bj22_fkw';'bj23_ekw';'bj24_fekw';'bji01';'bji02_f';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji07_fe';'bji08_fuw';'bji09_fvw';'bji10_fow';'bji11_euw';'bji12_evw';'bji13_eow';'bji14_feuw';'bji15_fevw';'bji16_feow';'bji17_yw';'bji18_fyw';'bji19_eyw';'bji20_feyw';'bji21_kw';'bji22_fkw';'bji23_ekw';'bji24_fekw';'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow';'bjm17_yw';'bjm18_fyw';'bjm19_eyw';'bjm20_feyw';'bjm21_kw';'bjm22_fkw';'bjm23_ekw';'bjm24_fekw';'bjmi01';'bjmi02_f';'bjmi03_e';'bjmi04_uw';'bjmi05_vw';'bjmi07_fe';'bjmi08_fuw';'bjmi09_fvw';'bjmi11_euw';'bjmi12_evw';'bjmi13_eow';'bjmi14_feuw';'bjmi15_fevw';'bjmi16_feow';'bjmi17_yw';'bjmi18_fyw';'bjmi19_eyw';'bjmi20_feyw';'bjmi21_kw';'bjmi22_fkw';'bjmi23_ekw';'bjmi24_fekw';'bm09_fvw';'bp04_uw';'bp22_fkw';'bpi14_feuw';'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow';'bpj17_yw';'bpj18_fyw';'bpj19_eyw';'bpj20_feyw';'bpj21_kw';'bpj22_fkw';'bpj23_ekw';'bpj24_fekw';'bpji01';'bpji02_f';'bpji03_e';'bpji04_uw';'bpji05_vw';'bpji06_ow';'bpji07_fe';'bpji08_fuw';'bpji09_fvw';'bpji10_fow';'bpji11_euw';'bpji12_evw';'bpji13_eow';'bpji14_feuw';'bpji15_fevw';'bpji16_feow';'bpji17_yw';'bpji18_fyw';'bpji19_eyw';'bpji21_kw';'bpji22_fkw';'bpji23_ekw';'bpji24_fekw';'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm12_evw';'bpjm15_fevw';'bpjm16_feow';'bpjm18_fyw';'bpjm19_eyw';'bpjmi01';'bpjmi02_f';'bpjmi03_e';'bpjmi04_uw';'bpjmi05_vw';'bpjmi06_ow';'bpjmi07_fe';'bpjmi08_fuw';'bpjmi11_euw';'bpjmi12_evw';'bpjmi13_eow';'bpjmi14_feuw';'bpjmi16_feow';'bpjmi17_yw';'bpjmi18_fyw';'bpjmi19_eyw';'bpjmi23_ekw';'bpm02_f';'bpm05_vw';'bpm20_feyw';'bpmi02_f';};
%         
        
            % v2 cF all (res_hierarfitmodels_cF (16-Apr-2015) v2 all)
%             request.mods2correct={'bji01';'bji02_f';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji07_fe';'bji08_fuw';'bji09_fvw';'bji10_fow';'bji11_euw';'bji12_evw';'bji13_eow';'bji14_feuw';'bji15_fevw';'bji16_feow';'bji17_yw';'bji18_fyw';'bji19_eyw';'bji20_feyw';'bji21_kw';'bji22_fkw';'bji23_ekw';'bji24_fekw';'bjm06_ow';'bjm19_eyw';'bjmi01';'bjmi02_f';'bjmi03_e';'bjmi04_uw';'bjmi05_vw';'bjmi06_ow';'bjmi08_fuw';'bjmi09_fvw';'bjmi10_fow';'bjmi12_evw';'bjmi13_eow';'bjmi14_feuw';'bjmi15_fevw';'bjmi16_feow';'bjmi17_yw';'bjmi18_fyw';'bjmi19_eyw';'bjmi20_feyw';'bjmi21_kw';'bjmi22_fkw';'bjmi23_ekw';'bjmi24_fekw';'bp04_uw';'bp06_ow';'bp18_fyw';'bp20_feyw';'bpi01';'bpi08_fuw';'bpi11_euw';'bpi19_eyw';'bpj16_feow';'bpji01';'bpji02_f';'bpji03_e';'bpji04_uw';'bpji05_vw';'bpji06_ow';'bpji07_fe';'bpji08_fuw';'bpji09_fvw';'bpji10_fow';'bpji11_euw';'bpji12_evw';'bpji15_fevw';'bpji16_feow';'bpji20_feyw';'bpji23_ekw';'bpjm02_f';'bpjm03_e';'bpjm05_vw';'bpjm19_eyw';'bpjmi01';'bpjmi02_f';'bpjmi03_e';'bpjmi04_uw';'bpjmi05_vw';'bpjmi06_ow';'bpjmi13_eow';'bpm13_eow';'bpm20_feyw';'bpmi15_fevw';'bpmi17_yw';};
            
            
            
        end
        
        % Col 1=model name, col 2=Transform: +='to', -='from'
        request.mods_specifictransform=[request.mods2correct     num2cell( 1 * ones(length(request.mods2correct),1)       )];            
        input('Continue to implement corrections? Must be specified!!');
    end
    
    
    % (b) Implement correction ---------------------------------------
    d_corpar=[request.mods_specifictransform(:,1)  cell(size(request.mods_specifictransform,1),4)];  
    % Name, Original line, Re-computed line, simulated choice (group)
    for mm=1:length(request.mods2correct)
        disp(['Correcting  ' request.mods2correct{mm} ' -------------------']) 
        wm.res=r_res(  find(strcmp(r_res(:,1), request.mods2correct{mm}))  ,:);
        wm.ri=r_iterations(  find(strcmp(r_iterations(:,1), request.mods2correct{mm}))  ,:);
        wm.rd=r_iterd(  find(strcmp(r_iterd(:,1), request.mods2correct{mm}))  ,:);
        wm.dm=details.models(strcmp(details.models(:,1), request.mods2correct{mm}),:);
        d_corpar{mm, 2}= wm.res;   wm.newres([1 4])= wm.res([1 4]); 
        d_corpar{mm, 5}=1; % Is it possibe? Strike out if any betas are negative or epsilons are out of range 
        %
        wm.parnames= wm.dm{3};
        wm.transform=request.mods_specifictransform{strcmp(request.mods_specifictransform(:,1), request.mods2correct{mm}), 2};
        wm.predChoice=[cell(n_subs,3) ; repmat({zeros(6,6,1)},1,3)];
    % wm.predChoice=[cell(n_subs,3) ; repmat({zeros(4,4,1)},1,3)];  % SLACK 
        wm.design=d_design;

        for s=1:n_subs
            ws.subdata =sd.subjdata{s,2};
            ws.recordedpars=wm.res{2}(s, 4:end);
            
            % Recorded parameters are assumed to correspond to the correct
            % point in the parameter space (i.e. minimum). However, I'm not
            % sure what transformations have been applied/omitted. The only
            % way to know this is by looking at the simulations. 
            % Approach: (a) Apply probable transformations (b) Simulate
            % group-level choice for manual inspection (c) Record which
            % transforms are ok (d) Apply requested transformations and
            % save new parameters, archive old (original) parameters.
            
            switch wm.transform
                case 1, ws.newpars=f_transpar(wm.parnames, ws.recordedpars, 'to');
                case -1, ws.newpars=f_transpar(wm.parnames, ws.recordedpars, 'from');
                case 0, ws.newpars=ws.recordedpars; 
                case 2, 
                    ws.newpars=f_transpar(wm.parnames, ws.recordedpars, 'to');
                    ws.newpars=f_transpar(wm.parnames, ws.newpars, 'to');
                case -2, 
                    ws.newpars=f_transpar(wm.parnames, ws.recordedpars, 'from');
                    ws.newpars=f_transpar(wm.parnames, ws.newpars, 'from');
                otherwise error('Transform not yet defined!');
            end
                
            % Feed new params thru value fxn
            ws.newpars_trans=f_transpar(wm.parnames, ws.newpars, 'from');
            if isempty(strfind(wm.res{1}, 'p'))==0
                [wm.nll(s),pch]=f_nllsoftmax_lapse(ws.newpars_trans,  {wm.res{1}  ws.subdata details.fixedpar sd.details.col},[], 0 );
            else [wm.nll(s),pch]=f_nllsoftmax(ws.newpars_trans,  {wm.res{1}  ws.subdata details.fixedpar sd.details.col},[], 0 );
            end
            wm.newres{2}(s, 1)=s; 
            wm.newres{2}(s, 4: 3+length(wm.parnames))=ws.newpars;
            wm.newres{2}(s, 2)=wm.nll(s);   % Assumed r_res lists ACTUAL nLL
            
            % Generate simulated choice
            eval(['ws.d(:,1:3)= ' wm.res{1}   '(ws.newpars_trans, {[] wm.design  details.fixedpar sd.details.col});']) % v(A/R/E)
            ws.beta=ws.newpars(1); ws.softmaxbase=( exp(ws.beta.* ws.d(:,1))+exp(ws.beta.*ws.d(:,2))+exp(ws.beta.*ws.d(:,3)   )  );
            wm.predChoice{s,1}= exp(ws.beta.* ws.d(:,1)   ) ./ ws.softmaxbase;
            wm.predChoice{s,2}= exp(ws.beta.* ws.d(:,2)   ) ./ ws.softmaxbase;
            wm.predChoice{s,3}= exp(ws.beta.* ws.d(:,3)   ) ./ ws.softmaxbase;
            if isempty(strfind(wm.res{1}, 'p'))~=1;
                ws.epsilon=ws.newpars(2);
                wm.predChoice{s,1} =  ws.epsilon+ (1-3*ws.epsilon)*wm.predChoice{s,1};
                wm.predChoice{s,2} =  ws.epsilon+ (1-3*ws.epsilon)*wm.predChoice{s,2};
                wm.predChoice{s,3} =  ws.epsilon+ (1-3*ws.epsilon)*wm.predChoice{s,3};
            end
            wm.predChoice{s,1} = fliplr(reshape(wm.predChoice{s,1},6,6))';
            wm.predChoice{s,2} = fliplr(reshape(wm.predChoice{s,2},6,6))';
            wm.predChoice{s,3} = fliplr(reshape(wm.predChoice{s,3},6,6))';            
%             wm.predChoice{s,1} = fliplr(reshape(wm.predChoice{s,1},4,4))';  % SLACK 
%             wm.predChoice{s,2} = fliplr(reshape(wm.predChoice{s,2},4,4))';
%             wm.predChoice{s,3} = fliplr(reshape(wm.predChoice{s,3},4,4))';
            wm.predChoice{n_subs+1,1}= wm.predChoice{n_subs+1,1}+wm.predChoice{s,1} ;
            wm.predChoice{n_subs+1,2}= wm.predChoice{n_subs+1,2}+wm.predChoice{s,2} ;
            wm.predChoice{n_subs+1,3}= wm.predChoice{n_subs+1,3}+wm.predChoice{s,3} ;
            
            % Any insight looking at workspace variables?
            % ugh. 
         
            % Disqualifying?
            if ws.beta<0 , d_corpar{mm, 5}=0; end
            if isreal(wm.nll(s))==0, d_corpar{mm, 5}=0; end
            if isfield(ws, 'epsilon')==1; if ws.epsilon<0 || ws.epsilon>1/3; d_corpar{mm, 5}=0; end, end
            
            ws=[];
        end
        
        wm.predChoice{n_subs+1,1}=wm.predChoice{n_subs+1,1}./ n_subs;
        wm.predChoice{n_subs+1,2}=wm.predChoice{n_subs+1,2}./ n_subs;
        wm.predChoice{n_subs+1,3}=wm.predChoice{n_subs+1,3}./ n_subs;
        
        % Record new parameters 
        wm.newres(6:5+length(wm.parnames)) = num2cell(mean(wm.newres{2}(:, 4:end)));
        wm.newres{3}= sum( wm.newres{2}(:, 2));
        d_corpar{mm, 3} = wm.newres;        
        d_corpar{mm, 4}=wm.predChoice(end, 1:3);  % Group-level simulated saved only  
        
        wm=[];
    end
    
    % Plot
    close all hidden
    f.plotcols=4;   f.plotrows=10;    f.fontsize=15; f.fontsize_title=20; f.fontname='PT Sans Caption';
    f.subplot_VerHorz=[0.05 0.05]; f.fig_BotTop=[0.03 0.03]; f.fig_LeftRight=[0.02 0.02]; ff=1; f.size= [130 5 500 1100];
    figure('Name', ['Simulations ' num2str(ff)], 'NumberTitle', 'off', 'Position', f.size, 'Color', 'w');  k=1;
    for m=1:length(request.mods2correct)
        subtightplot(f.plotrows,  f.plotcols, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        text(0.5,0.5 , request.mods2correct{m}, 'FontSize', f.fontsize_title ); axis 'off'; k=k+1;
        
        for c=1:3
            subtightplot(f.plotrows,  f.plotcols, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            
            if isreal(d_corpar{m,4}{c})==1 && d_corpar{mm, 5}==1;
                imagesc(d_corpar{m,4}{c}, [0 1]); axis('square'); axis off
            else  text(0.5,0.5 , 'x'); axis 'off'; 
            end
        end
        if k == f.plotrows*f.plotcols+1, 
            ff=ff+1; figure('Name', ['Simulations ' num2str(ff)], 'NumberTitle', 'off', 'Position', f.size+[ff*70 0 0 0], 'Color', 'w');  k=1;
        end
    end    
    disp('FROM HERE, visually identify which models are correctly transformed now'); input('Continue to save? Only do so when all simulations look good.  '); 
    
    % Alter inputs. Originals (d_corpar, 2nd col) saved in details)
    for m=1:length(request.mods2correct)
        % REPLACE in results file
        
        wm.nr=find(strcmp( r_res(:,1), d_corpar{m,1}));
        %     wm.nri=find(strcmp( r_iterations(:,1), d_corpar{m,1}));
        %     wm.nrd=find(strcmp( r_iterd(:,1), d_corpar{m,1}));
        %
        %  Correct r_res ONLY. r_iterations and r_iterd are assumed to be ok
        r_res{wm.nr,2}= d_corpar{m,3}{2};
        r_res(wm.nr, 6:5+size(d_corpar{m,3}{2},2)-3) =  num2cell( mean(d_corpar{m,3}{2}(:, 4:end)) );
        
        wm=[];
    end
    details.d_corpar=d_corpar;
    
    % SAVE ---------------------------------------------------------------------------------------------------
    resfilename= f_newname([request.filename ' PARtransformed'], pwd); input(['Save as     ' resfilename  '  ?']);
    save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');
    
end


%% Remove models with certain parameters from the model space
% ALSO, remove duplicates


if request.removemods==1
    removepars={};
    
    wherefile='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical';
    filename='res_hierarfitmodels_ct (09-Apr-2015)1';
    ishierar=1;
    
    % ------------------------------------------------------
    cd(wherefile);    load([filename '.mat']);
    for p=1:length(removepars)
        
        for m=1:size(r_res,1)  % res
            if isempty(strfind(r_res{m,1}, removepars{p}))==0
                r_res{m,1}='zzzzzzzzzzzz';
            end
        end
        for m=1:size(r_iterations,1) % r_iterations
            if isempty(strfind(r_iterations{m,1}, removepars{p}))==0
                r_iterations{m,1}='zzzzzzzzzzzz';
            end
        end
        if ishierar ==1;
            for m=1:size(r_iterd,1) % r_iterd
                if isempty(strfind(r_iterd{m,1}, removepars{p}))==0
                    r_iterd{m,1}='zzzzzzzzzzzz';
                end
            end
        end
        
        % DETAILS #############
        for m=1:size(details.models,1)   % details.models
            if isempty(strfind(details.models{m,1}, removepars{p}))==0
                details.models{m,1}='zzzzzzzzzzzz';
            end
        end
        
        
    end
    for o=1:1  % Remove duplicates
        newfits={}; k=1;
        for m=1:size(r_res ,1)  % res
            if sum(strcmp(r_res(:,1), r_res{m,1}))>1 & strcmp( r_res{m,1},  'zzzzzzzzzzzz')==0
                wr.fits=r_res(find(strcmp(r_res(:,1), r_res{m,1})),:);
                %             wr.fits=
                r_res(find(strcmp(r_res(:,1), r_res{m,1})),1) = repmat(    {'zzzzzzzzzzzz'}, sum(strcmp(r_res(:,1), r_res{m,1})), 1);
                
                %
                wr.Ls=cell2mat(wr.fits(:, 3));
                if sum(wr.Ls==min(wr.Ls))~=1; error('MULTIPLE matches. Just delete one manually'); end
                wr.bestfit=wr.fits ( find(wr.Ls==min(wr.Ls)  ), :);
                newfits(k,  1:  size(wr.fits,2)) =  wr.bestfit;   k=k+1;
            end
        end
        r_res=[r_res; newfits];
        
        
        
        for m=1:size(details.models,1)  % res
            if sum(strcmp(details.models(:,1), details.models{m,1}))>1 & strcmp( details.models{m,1},  'zzzzzzzzzzzz')==0
                details.models{m,1}='zzzzzzzzzzzz';
            end
            
        end
        
    end
    
    % Get rid of
    r_res=sortrows(r_res,1); r_res(strcmp(r_res(:,1), 'zzzzzzzzzzzz'), :)=[];
    r_iterations=sortrows(r_iterations,1);  r_iterations(strcmp(r_iterations(:,1), 'zzzzzzzzzzzz'), :)=[];
    if ishierar ==1;   r_iterd= sortrows(r_iterd,1);  r_iterd(strcmp(r_iterd(:,1), 'zzzzzzzzzzzz'), :)=[]; end
    details.models=sortrows(details.models,1); details.models(strcmp(details.models(:,1), 'zzzzzzzzzzzz'), :)=[];
    details.whichmodels=sortrows(details.models(:,1),1);    details.n_models=length(details.whichmodels); r_res=sortrows(r_res,3);
    
    % SAVE ---------------------------------------------------------------------------------------------------
    if ishierar ==1  ; resfilename=['FIXED res_hierarfitmodels_' details.task ' (' date ')']; % Save a hierarfit
    else  resfilename=['FIXED res_fitmodels_' details.task ' (' date ')'];    % Saved a fixed fit
    end
    resfilename= f_newname(resfilename, pwd); input(['Save as     ' resfilename  '  ?']);
    if ishierar ==1; save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');   % Save a hierarfit
    else  save(resfilename, 'details', 'r_iterations','r_res', 'rc','errorlog');  % Saved a fixed fit
    end
end