% Correlate model parameters with  betas
clear all; close all hidden; clc

% Settings
request.fitcF_date='(10-Apr-2014)2'; % cF and ct have to be fit on same day (same data and procedures)
request.fitct_date=request.fitcF_date;
%
request.fitmethod='fit'; % fit/grid
request.dataset_date='(08-Apr-2014)';
%
request.roibeta_where='/Users/EleanorL/Dropbox/WorkPC/Explore beta/c3 choice betas/c3 ROI battery';
request.roibeta_where='D:\Dropbox\WorkPC\Explore beta\c3 choice betas\c3 ROI battery';
request.roibeta_file= '(03-Nov-2013) Extracted betas.txt';

for o1=1:1 % Request models
    log.whichmodel_fams{1}={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow'};
    log.whichmodel_fams{2}={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow'};
    log.whichmodel_fams{3}={'ba01';'ba02_f';'ba03_e';'ba04_uw';'ba05_vw';'ba06_ow';'ba07_fe';'ba08_fuw';'ba09_fvw';'ba10_fow';'ba11_euw';'ba12_evw';'ba13_eow';'ba14_feuw';'ba15_fevw';'ba16_feow'};
    %
    log.whichmodel_fams{4}={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow'};
    log.whichmodel_fams{5}={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow'};
    log.whichmodel_fams{6}={'bpa01';'bpa02_f';'bpa03_e';'bpa04_uw';'bpa05_vw';'bpa06_ow';'bpa07_fe';'bpa08_fuw';'bpa09_fvw';'bpa10_fow';'bpa11_euw';'bpa12_evw';'bpa13_eow';'bpa14_feuw';'bpa15_fevw';'bpa16_feow'};
end
log.whichmodel_cF='bpm16_feow';
log.whichmodel_ct='bm14_feuw';

for o1=1:1 %% Set up
    
    % Folders
    w.where=pwd;
    if strcmp(w.where(2), ':')==1;  where.scripts='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
    else  where.scripts='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
    end
    rand('state',sum(100*clock)); cd(where.scripts); 
    w.dd=load([where.scripts filesep '2 Analysis inputs' filesep 'All data ' request.dataset_date '.mat']); subjdata=w.dd.subjdata;
    
    % Load results
    w.dcf=load([where.scripts filesep '2 Analysis inputs' filesep 'res_'     request.fitmethod  'models_cF ' request.fitcF_date]); d_cf=w.dcf.r_res;
    w.dct=load([where.scripts filesep '2 Analysis inputs' filesep 'res_'     request.fitmethod  'models_ct ' request.fitct_date]); d_ct=w.dct.r_res;
    log.subjects=w.dcf.details.subjects; log.n_subjs=length(log.subjects);       if sum(strcmp(w.dcf.details.subjects, w.dct.details.subjects))~=log.n_subjs; error('List of subjects ~= for cF and ct modelfits. Check +implement workarounds'); end
    log.top5.cF= d_cf(1:5,1); log.top5.ct= d_ct(1:5,1); col=w.dcf.details.col;
    if sum(strcmp(d_cf(:,1), log.whichmodel_cF))~=1 || sum(strcmp(d_ct(:,1), log.whichmodel_ct))~=1;  w.w=cell(max([size(d_cf,1) size(d_ct,1)])+2,2); w.w(1:size(d_cf,1)+2,1)= [{'[ cF models ]'; '   '};  sortrows(d_ct(:,1))]; w.w(1:size(d_cf,1)+2,2)= [{'[ ct models ]'; '   '};   sortrows(d_ct(:,1))]; disp( [[{'# Available models'} {'##'}]; w.w]);    error(['Requested model is not in fit results ('  log.whichmodel_cF '  &  ' log.whichmodel_ct ')']);   end; % Check requested model
    [log.modspex  log.par_transformations log.cfmodspex] = f_modelsettings({log.whichmodel_cF});
    [log.modspex  log.par_transformations log.ctmodspex] = f_modelsettings({log.whichmodel_ct});
    
    disp('====================================')
    disp(['Requested: ' num2str(size(subjdata,1)) ' subjects']); disp(' ')
    disp('Running correlations for model:   '  ); disp(['                     cF     ' log.whichmodel_cF]); disp(['                    ct     ' log.whichmodel_ct]); disp(' ');
    disp(['Beta file :  ' request.roibeta_where])
    disp('====================================')
end

%% (1) Set up data

for o1=1:1 % Set up original data 
    
    % Data
    w.b=importdata([request.roibeta_where filesep request.roibeta_file]);
    d_beta= [strtrim(w.b.textdata(:,1))   [w.b.textdata(1,2:end);  num2cell(w.b.data)]];
    if sum(strcmp(d_cf(:,1), log.whichmodel_cF))+sum(strcmp(d_ct(:,1), log.whichmodel_ct))~=2; error(['Could not find requested model for beta correlation, for both cF and ct (' log.whichmodel_cF ')']); end
    log.cf_modpars=log.cfmodspex{find(strcmp(log.cfmodspex(:,1),log.whichmodel_cF)),3}';
    log.ct_modpars=log.ctmodspex{find(strcmp(log.ctmodspex(:,1),log.whichmodel_ct)),3}';
    d_behcf=[cellfun(@(x)['cF ' x], log.cf_modpars, 'UniformOutput',0);  num2cell(d_cf{find(strcmp(d_cf(:,1), log.whichmodel_cF)),2}(:,4:end))];
    d_behct=[cellfun(@(x)['ct ' x], log.ct_modpars, 'UniformOutput',0) ; num2cell(d_ct{find(strcmp(d_ct(:,1), log.whichmodel_ct)),2}(:,4:end))];
    d_beh=[[{'Subject'}; log.subjects] d_behcf d_behct];
    
    % Apply subject selection & ordering
    d_beh=[d_beh(1,:); sortrows(d_beh(2:end,:),1)]; d_beta=[d_beta(1,:); sortrows(d_beta(2:end,:),1)];
    if strcmp(char(d_beh(2:end,1)),char(d_beta(2:end,1)))~=1; error('Subjects from beh and beta tables do not immediately match. Implement selection procedure'); end
    log.subjects=d_beh(2:end,1);
    
    % Log of roicons, rois, cons, beh
    log.roicon=strtrim(d_beta(1,2:end)');d_beta(1,2:end)=log.roicon';
    log.beh=cellfun(@(x)[x ' - ' log.whichmodel_cF], strtrim(d_beh(1,2:end)'), 'UniformOutput',0);  d_beh(1,2:end)=log.beh';
    log.n_roicon=length(log.roicon); log.n_beh=length(log.beh);
    k=8; % characters for identifying rois
    log.firstcon_rois=find(1-cellfun(@isempty,strfind(log.roicon, log.roicon{1}(length(log.roicon{1})-k:end))));
    log.rois=cellfun(@(x)x(1:length(x)-length(log.roicon{log.firstcon_rois(1)}(find(1-(char(log.roicon{log.firstcon_rois(1)})==char(log.roicon{log.firstcon_rois(2)})), 1, 'last')+1:end))), log.roicon(log.firstcon_rois), 'UniformOutput',0);
    log.n_rois=length(log.rois);
    log.n_cons=sum(cell2mat(strfind(log.roicon, log.rois{1})));
    log.cons=cellfun(@(x)x(length(log.rois{1})+2:end), log.roicon(1:log.n_cons), 'UniformOutput',0); % change # here if the no of connecting characters between roi and contrast ~=1
    
    
    
    
    
    
    % Mean centre
    for r=1:log.n_rois
%         d_beta(2:end, 1+(r-1)*log.n_cons+1:  1+(r-1)*log.n_cons+log.n_cons)= num2cell(cell2mat(d_beta(2:end, 1+(r-1)*log.n_cons+1:  1+(r-1)*log.n_cons+log.n_cons))  -   repmat(mean(cell2mat(d_beta(2:end, 1+(r-1)*log.n_cons+1:  1+(r-1)*log.n_cons+log.n_cons)),2) , 1,log.n_cons));
    end
    
    
    
    
end
for o1=1:1 % Set up new beta values 
    orig.d_beta=d_beta; d_beta=orig.d_beta(:,1);
    log.firstcon_rois(log.n_rois+1)=length(log.roicon)+1;
    
    for r=1:log.n_rois
        wr.roicon_names=orig.d_beta(1,1+(log.firstcon_rois(r):log.firstcon_rois(r+1)-1));
        wr.roicon_betas=orig.d_beta(2:end,1+(log.firstcon_rois(r):log.firstcon_rois(r+1)-1));
        wr.roicon=[wr.roicon_names; wr.roicon_betas];
        wr.newcon=cell(log.n_subjs+1,0); wr.b=1;
        
        % Specify manual instructions for each new contrast (first line only) --------------
        %       xx= names of contrasts
        %       instruc={name of new contrast,  comparison expressed in terms of xx variables};
        
        
        % (1) cF Rej-Exp
        xx={'cF_Reject'; 'cF_Explore'}; instruc={'cF_Rej-Exp';  'x(:,1) - x(:,2)'};
        x=nan(log.n_subjs, length(xx));
        for i=1:length(xx)
            if sum(strcmp(log.cons, xx{i}))~=1; error(['[Creating new contrasts] could not find specified contrast (' xx{i} ')']); end
            x(:,i)= cell2mat(wr.roicon_betas(:, find(strcmp(log.cons, xx{i}))));
        end
        wr.newcon{1,wr.b}= [log.rois{r} '-' instruc{1}]; eval(['wr.newcon(2:end,wr.b)=num2cell('     instruc{2}   ');']); wr.b=wr.b+1;
        
        % (2) ct Rej-Exp
        xx={'ct_Bomb'; 'ct_Explore'}; instruc={'ct_Rej-Exp';  'x(:,1) - x(:,2)'};
        x=nan(log.n_subjs, length(xx));
        for i=1:length(xx)
            if sum(strcmp(log.cons, xx{i}))~=1; error(['[Creating new contrasts] could not find specified contrast (' xx{i} ')']); end
            x(:,i)= cell2mat(wr.roicon_betas(:, find(strcmp(log.cons, xx{i}))));
        end
        wr.newcon{1,wr.b}= [log.rois{r} '-' instruc{1}]; eval(['wr.newcon(2:end,wr.b)=num2cell('     instruc{2}   ');']); wr.b=wr.b+1;
        
        % (3) Rej_cF-ct
        xx={'cF_Reject'; 'ct_Bomb'}; instruc={'Rej_cF-ct';  'x(:,1) - x(:,2)'};
        x=nan(log.n_subjs, length(xx));
        for i=1:length(xx)
            if sum(strcmp(log.cons, xx{i}))~=1; error(['[Creating new contrasts] could not find specified contrast (' xx{i} ')']); end
            x(:,i)= cell2mat(wr.roicon_betas(:, find(strcmp(log.cons, xx{i}))));
        end
        wr.newcon{1,wr.b}= [log.rois{r} '-' instruc{1}]; eval(['wr.newcon(2:end,wr.b)=num2cell('     instruc{2}   ');']); wr.b=wr.b+1;
        
        % (4) Exp_cF-ct
        xx={'cF_Explore'; 'ct_Explore'}; instruc={'Exp_cF-ct';  'x(:,1) - x(:,2)'};
        x=nan(log.n_subjs, length(xx));
        for i=1:length(xx)
            if sum(strcmp(log.cons, xx{i}))~=1; error(['[Creating new contrasts] could not find specified contrast (' xx{i} ')']); end
            x(:,i)= cell2mat(wr.roicon_betas(:, find(strcmp(log.cons, xx{i}))));
        end
        wr.newcon{1,wr.b}= [log.rois{r} '-' instruc{1}]; eval(['wr.newcon(2:end,wr.b)=num2cell('     instruc{2}   ');']); wr.b=wr.b+1;
        
        %
        d_beta=[d_beta wr.roicon wr.newcon];
    end
    
    % Re-read beta log
    log.firstcon_rois=[];
    log.roicon=d_beta(1,2:end)';
    log.n_roicon=length(log.roicon);
    log.n_cons=log.n_roicon/log.n_rois;
    log.cons=cellfun(@(x)x(length(log.rois{1})+2:end),  log.roicon(1: log.n_cons) , 'UniformOutput',0); % change # here if the no of connecting characters between roi and contrast ~=1
end
for o1=1:1 % Set up new beh values 
    orig.d_beh=d_beh; 
    
    % For shared parameters, calculate difference
    log.sharedpars=intersect(log.cf_modpars, log.ct_modpars);
    for p=1:length(log.sharedpars)
        k=size(d_beh,2)+1;
        d_beh{1,k}= [log.sharedpars{p} '  cF-ct'];
        d_beh(2:end,k)=num2cell(cell2mat(d_behcf(2:end, find(strcmp(  log.sharedpars{p}, log.cf_modpars))))  - cell2mat(d_behct(2:end, find(strcmp(  log.sharedpars{p}, log.ct_modpars)))));
    end
    
    % Re-read details
    log.beh=d_beh(1,2:end)';
    log.n_beh=length(log.beh);
end

%%  Perform correlations (row=roicon, col=param)

logn.pthreshold=0.01;

rc_r=cell(log.n_roicon+1, log.n_beh+1); rc_r(2:end,1)=log.roicon; rc_r(1,2:end)=log.beh; rc_p=rc_r;
for rc=1:log.n_roicon
    for p=1:log.n_beh
        wp.beh=cell2mat(d_beh(2:end, find(strcmp(d_beh(1,:),log.beh{p}))));
        wp.beta=cell2mat(d_beta(2:end, find(strcmp(d_beta(1,:),log.roicon{rc}))));
        
        % Execute correlation
        [wp.r wp.p]=corr(wp.beh, wp.beta);
        if wp.p<=logn.pthreshold
            rc_p{rc+1,p+1}=wp.p;
            rc_r{rc+1,p+1}=wp.r;
        elseif wp.p<=0.1
            rc_p{rc+1,p+1}=wp.p;
        end
        wp=[];
    end
end
openvar rc_r
openvar rc_p