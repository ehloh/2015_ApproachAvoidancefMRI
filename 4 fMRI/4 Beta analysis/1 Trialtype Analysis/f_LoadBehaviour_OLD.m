function [ d_design ] = f_LoadBehaviour( request, col)
% d_design{subject, task}: trialstats for cF/ct (specified in col)
%
% EnvThreat-warped task space must be already computed (ETwarp conflict),
% stored where the modeling results are stored



keyboard

% d_betas:       d_betas{r,1}{task}: all data for single glm
%                     d_betas{r,2}{subject, task}:subject data for 2-level glms
%                     d_betas{r,3}{task}: betas from subject-specific first level glms (subject = row, col=iv)
% d_glmres{task}{roi+1, iv+1}: t stat if sig
% d_glmresp: p values



%% Compile canonical design

% Canonical design variables
desvars={'EnvThreat';'NTokens';'pLoss';'Entropy';'EV'; 'EnvThreatOrig'};
d_par=nan(36,3); d_par(:, [col.EnvThreat col.NTokens])=[(repmat(1:6, 1,6)'/6) sortrows(repmat(1:6, 1,6)')*2];
dcf=fpar_conflict(d_par, col);
dct=fpar_control(d_par, col);
dcf(:, col.EnvThreatOrig)= dcf(:, col.EnvThreat);  dct(:, col.EnvThreatOrig)= dct(:, col.EnvThreat);

    
    % Load task space as per subject's fitted parameters
    % cf=load([where.modres filesep 'res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968.mat']);
    % ct= load([where.modres filesep 'res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293.mat']);
    cfee=load([where.modres filesep 'ETwarp conflict.mat']); ctee=load([where.modres filesep 'ETwarp control.mat']); % Assumed to be the correct model
%     cfev=cfee.d_etwarpv; ctev=ctee.d_etwarpv;
    cfe= cfee.d_etwarp; cte= ctee.d_etwarp; 
    for s=1:log.n_subjs, 
        cfe{s,2}.NTokens= cfe{s,2}.NTok; cte{s,2}.NTokens= cte{s,2}.NTok; 
        cfe{s,2}.EnvThreatOrig= repmat((6:-1:1)',1,6)./6;  cte{s,2}.EnvThreatOrig= repmat((6:-1:1)',1,6)./6; 
    end
    
    
    
    % Compile subject-specific task pars according to requested
    %     d_design{s,1}=cf, d_design{s,2}=ct
    %     d_design{s+1,1}: entire concatenated string of pars
    d_design{log.n_subjs+1,1}=[]; d_design{log.n_subjs+1,2}=[];
    for s=1:log.n_subjs
        if request.DesignCanon
            for v=1:length(desvars)
                eval(['d_design{s,1}(:, col.' desvars{v} ')=dcf(:, col.' desvars{v} ');'])
                eval(['d_design{s,2}(:, col.' desvars{v} ')=dct(:, col.' desvars{v} ');'])
                d_design{s,1}(:, col.EnvThreatOrig)=dcf(:, col.EnvThreat);
                d_design{s,2}(:, col.EnvThreatOrig)=dct(:, col.EnvThreat);
            end
        else
            for v=1:length(desvars)
                eval(['ws.v=cfe{s,2}.' desvars{v} ';']), ws.v= flipud(ws.v);
                eval(['d_design{s,1}(:, col.' desvars{v} ')=ws.v(:);'])
                eval(['ws.v=cte{s,2}.' desvars{v} ';']), ws.v= flipud(ws.v);
                eval(['d_design{s,2}(:, col.' desvars{v} ')=ws.v(:);'])
                
                ws.v=cfe{s,2}.EnvThreatOrig; ws.v= flipud(ws.v);
                d_design{s,1}(:, col.EnvThreatOrig)=ws.v(:);
                ws.v= cte{s,2}.EnvThreatOrig; ws.v= flipud(ws.v);
                d_design{s,2}(:, col.EnvThreatOrig)=ws.v(:);
            end
        end
        
        % Compile the whole file
        d_design{log.n_subjs+1,1}=[d_design{log.n_subjs+1,1};   d_design{s,1}];
        d_design{log.n_subjs+1,2}=[d_design{log.n_subjs+1,2};   d_design{s,2}];
    end
    
    dcf=[];  dct=[]; cfe=[]; cte=[]; 
    



end

