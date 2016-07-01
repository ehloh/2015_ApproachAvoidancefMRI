function [logg d_roimatrix_cF d_roimatrix_ct] = f_LoadTrialtypeBetamatrix(request, logg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% OPTIONAL functionality
request.plotsubjecchecks=0;
request.check4nans=0;


%% (1) Load data and apply subject/ROI selections
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
%                       EnvThreat in reverse order (i.e. as per imagesc plot)

% Load data file
request.loadextractedbetas = spm_select('List', request.where_betas, 'Extracted beta matrix ');
if size(request.loadextractedbetas ,1)~=1; error(['Requested FIND extracted beta file, but none or >1 found (' request.where_betas ')']); end
load([request.where_betas filesep request.loadextractedbetas])

% Implement ROI requests/re-ordering
if isempty(request.ROIs)==0
    oldroin=cell2mat(cellfun(@(x)find(strcmp(log.ROInames, x)),  request.ROIs, 'UniformOutput',0));
    if length(oldroin)~=length(request.ROIs); error('Could not find all ROIs!!');  end 
    log.ROInames=request.ROIs;
    d_roimatrix_cF= d_roimatrix_cF(:,:, :, oldroin, :);
    d_roimatrix_ct= d_roimatrix_ct(:,:, :, oldroin, :);
    log.n_rois= length(log.ROInames);
    
%     % CHECK? 
%     newroimatrix_cF= d_roimatrix_cF(:,:, :, oldroin, :);
%     newroimatrix_ct= d_roimatrix_ct(:,:, :, oldroin, :);    
%     checksub=15;
%     for r=1:length(request.ROIs)
%         newroimatrix_cF(:,:,:, r, checksub) - d_roimatrix_cF(:,:,:, oldroin(r), checksub)
%         newroimatrix_ct(:,:,:, r, checksub) - d_roimatrix_ct(:,:,:, oldroin(r), checksub)
%     end
    
end

% Implement subject selections requested 
if isempty(request.specificsubjects)==0     
    subnums_inold=cellfun(@(x)find(strcmp(log.subjects, x)), logg.subjects);

    d_roimatrix_cF= d_roimatrix_cF(:,:, :, :, subnums_inold);
    d_roimatrix_ct= d_roimatrix_ct(:,:, :, :, subnums_inold);
    log.n_rois= length(log.ROInames);
    
%     % CHECK? 
%     newroimatrix_cF= d_roimatrix_cF(:,:, :, :, subnums_inold);
%     newroimatrix_ct= d_roimatrix_ct(:,:, :, :, subnums_inold);    
%     checkroi=2;
%     for r=1:length(subnums_inold)
%         newroimatrix_cF(:,:,:, checkroi, r) - d_roimatrix_cF(:,:,:, checkroi, subnums_inold(r))
%         newroimatrix_ct(:,:,:, checkroi, r) - d_roimatrix_ct(:,:,:, checkroi, subnums_inold(r))
%     end
    log.subjects=logg.subjects;  log.n_subjs =length(log.subjects);
end

%% Preprocessing if requested


% SETTINGS FOR PREPROCESSING ##################
request.flex_mintrials4inclusion=[]; % empty to omit trialwise filtering. Trialwise filtering applied first, then Subject-wise filtering
request.minsubs_cell=1;
request.filter_high=[];  % Empty to omit 
request.filter_low=[];

% High/Low-pass filter + Other optional transforms (e.g zscore)
if isempty(request.filter_high)==0
for s=1:log.n_subjs
    for r=1:log.n_rois
        
        % High/Low-pass filter (individual cells excluded only)
        if s==1 && r==1,  disp(['High/low-pass fltering between: ' num2str(request.filter_low) ' to ' num2str(request.filter_high) ]), end
        for e=1:log.nlevels
            for n=1:log.nlevels
                for c=1:request.n_choices
                    
                    % cF
                    if d_roimatrix_cF(e,n,c,r,s)>request.filter_high || d_roimatrix_cF(e,n,c,r,s)<request.filter_low
                        d_roimatrix_cF(e,n,c,r,s)=nan;
                    end
                    
                    % ct
                    if d_roimatrix_ct(e,n,c,r,s)>request.filter_high || d_roimatrix_ct(e,n,c,r,s)<request.filter_low
                        d_roimatrix_ct(e,n,c,r,s)=nan;
                    end
                    
                end
            end
        end
        
    end
end
else disp('NO high/low-pass filtering');
end

% Flexible models: no. of trials in each cell?
if isempty(request.flex_mintrials4inclusion)==0 && strcmp(log.FLmodel(1), 'f')==1
    for o=1:1 % Flex
        if s==1, input('FILTERING by no. trials etc. Continue?');  end
        
        % NTrials data to match s_roimatrix and s_roivect below
        %  s_ntrialmatrix {subject}  = Task x Choice x EnvThreat x NTokens
        %  s_ntrialvect= matrix of subject ntrials in horizontal vector  (row=subject, cols = task-choice-envthreat-ntokens)
        s_ntrialmatrix=cell(log.n_subjs+1, 1);  s_ntrialmatrix{log.n_subjs+1}=zeros(2, request.n_choices, log.nlevels, log.nlevels); s_ntrialvect=nan(log.n_subjs, 2*request.n_choices*(log.nlevels*log.nlevels+1)+ 4);
        w.beh=load([where.modellingdata filesep request.behfilename '.mat']);
        d_beh=w.beh.subjdata;  % Col 2=cF, col 3=ct, col 4=both, col 5=
        w.behcol=w.beh.details.col;
        w.behcol.Design_EnvThreatLevel=size(d_beh{1,2},2)+1;   w.behcol.Design_NTokensLevel=size(d_beh{1,2},2)+2;
        
        for s=1:log.n_subjs
            ws.cF=sortrows(d_beh{find(strcmp( d_beh(:,1),  log.subjects{s})),2}, [w.behcol.EnvThreat w.behcol.NTokens]);
            ws.ct=sortrows(d_beh{find(strcmp( d_beh(:,1),  log.subjects{s})),3}, [w.behcol.EnvThreat w.behcol.NTokens]);
            
            for o2=1:1
                % Mark design levels (EnvThreat, NTokens)
                if isempty(strfind(log.FLmodel, 'Chunk'))==1
                    ws.cF(:, w.behcol.Design_EnvThreatLevel) = ws.cF(:, w.behcol.EnvThreat)*6;
                    ws.cF(:, w.behcol.Design_NTokensLevel) = ws.cF(:, w.behcol.NTokens)/2;
                    ws.ct(:, w.behcol.Design_EnvThreatLevel) = ws.ct(:, w.behcol.EnvThreat)*6;
                    ws.ct(:, w.behcol.Design_NTokensLevel) = ws.ct(:, w.behcol.NTokens)/2;
                else
                    ws.cF(:, w.behcol.Design_EnvThreatLevel) = 10+ws.cF(:, w.behcol.EnvThreat)*6;
                    ws.cF(:, w.behcol.Design_NTokensLevel) = 10+ws.cF(:, w.behcol.NTokens)/2;
                    ws.ct(:, w.behcol.Design_EnvThreatLevel) = 10+ws.ct(:, w.behcol.EnvThreat)*6;
                    ws.ct(:, w.behcol.Design_NTokensLevel) = 10+ws.ct(:, w.behcol.NTokens)/2;
                    ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==11  | ws.cF(:, w.behcol.Design_EnvThreatLevel)==12, w.behcol.Design_EnvThreatLevel)=1;
                    ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==13 | ws.cF(:, w.behcol.Design_EnvThreatLevel)==14, w.behcol.Design_EnvThreatLevel)=2;
                    ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==15 | ws.cF(:, w.behcol.Design_EnvThreatLevel)==16, w.behcol.Design_EnvThreatLevel)=3;
                    ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==11 |  ws.cF(:, w.behcol.Design_NTokensLevel)==12, w.behcol.Design_NTokensLevel)=1;
                    ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==13 |  ws.cF(:, w.behcol.Design_NTokensLevel)==14, w.behcol.Design_NTokensLevel)=2;
                    ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==15 |  ws.cF(:, w.behcol.Design_NTokensLevel)==16, w.behcol.Design_NTokensLevel)=3;
                    ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==11  | ws.ct(:, w.behcol.Design_EnvThreatLevel)==12, w.behcol.Design_EnvThreatLevel)=1;
                    ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==13 | ws.ct(:, w.behcol.Design_EnvThreatLevel)==14, w.behcol.Design_EnvThreatLevel)=2;
                    ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==15 | ws.ct(:, w.behcol.Design_EnvThreatLevel)==16, w.behcol.Design_EnvThreatLevel)=3;
                    ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==11 |  ws.ct(:, w.behcol.Design_NTokensLevel)==12, w.behcol.Design_NTokensLevel)=1;
                    ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==13 |  ws.ct(:, w.behcol.Design_NTokensLevel)==14, w.behcol.Design_NTokensLevel)=2;
                    ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==15 |  ws.ct(:, w.behcol.Design_NTokensLevel)==16, w.behcol.Design_NTokensLevel)=3;
                end
            end
            
            % Scroll through choice x envthreat x ntokens, nan-out trials with insufficient n trials
            for ch=1:request.n_choices
                for ee=1:log.nlevels
                    e=log.nlevels+1-ee;
                    for n=1:log.nlevels
                        
                        % No. of trials
                        s_ntrialmatrix{s}(1, ch, e, n)= sum(ws.cF(ws.cF(:,w.behcol.Design_EnvThreatLevel)==ee &  ws.cF(:,w.behcol.Design_NTokensLevel)==n , w.behcol.Choice)==ch);
                        s_ntrialmatrix{s}(2, ch, e, n)= sum(ws.ct(ws.ct(:,w.behcol.Design_EnvThreatLevel)==ee &  ws.ct(:,w.behcol.Design_NTokensLevel)==n , w.behcol.Choice)==ch);
                        
                        % Nan-out cF and cts cells with insufficient trials
                        if s_ntrialmatrix{s}(1, ch, e, n)<request.flex_mintrials4inclusion;
                            d_roimatrix_cF(e,n, ch, :, s)=nan;
                        end
                        if s_ntrialmatrix{s}(2, ch, e, n)<request.flex_mintrials4inclusion;
                            d_roimatrix_cF(e,n, ch, :, s)=nan;
                        end
                        
                    end
                end
            end
            s_ntrialmatrix{log.n_subjs+1}=s_ntrialmatrix{log.n_subjs+1}+s_ntrialmatrix{s};
        end
        s_ntrialmatrix{log.n_subjs+1}=s_ntrialmatrix{log.n_subjs+1}/log.n_subjs;
        
        % Assemble vector.  Conversion of EnvThreat x NTokens matrix to a vector
        %   Vector order: EnvThreat (descending threat) followed by NTokens (ascending N)
        for s=1:log.n_subjs+1
            k=1;
            for t=1:2
                for ch=1:request.n_choices
                    wv=(squeeze(s_ntrialmatrix{s}(t,ch,:,:)))';
                    s_ntrialvect(s, k:k+length(wv(:))-1)=(wv(:))';
                    k=k+length(wv(:))+1;
                end
                k=k+2;
            end
        end
    end
end

% Mean Centre
if request.meancentre
    for s=1:log.n_subjs,
        for r=1:log.n_rois
            
            % Mean centre across both tasks 
            wsr.mean=mean2(nanmean(nanmean(([d_roimatrix_cF(:,:,:,r,s); d_roimatrix_ct(:,:,:,r,s)]))));
            d_roimatrix_cF(:,:,:,r,s)=d_roimatrix_cF(:,:,:,r,s)-wsr.mean;
            d_roimatrix_ct(:,:,:,r,s)=d_roimatrix_ct(:,:,:,r,s)-wsr.mean;
        end
    end
end

% Check betas for nans
if request.check4nans
    for s=1:log.n_subjs, k=1;
        for c=1:request.n_choices
            for r=1:log.n_rois
                ws.r(k)=sum(sum(isnan(d_roimatrix_cF(:,:,1,r,s)))); k=k+1;
            end
        end
        disp([log.subjects{s} '  : '  num2str(sum(ws.r)) ' bad cells'])
        ws=[];
    end
end

%% Subject sanity checks (usually not needed for Trialtype models)
%  s_roimatrix {subject, roi}  = Task x Choice x EnvThreat x NTokens
%  s_roivect{roi}= matrix of subject betas in horizontal vector  (row=subject, cols = task-choice-envthreat-ntokens)

if request.subchecks_forflex
    
    % Assemble subject beta vectors
    s_roimatrix=cell(log.n_subjs, log.n_rois); s_roivect=repmat({nan(log.n_subjs, 2*request.n_choices*(log.nlevels*log.nlevels+1)+ 4)}, log.n_rois,1);
    for s=1:log.n_subjs
        for r=1:log.n_rois
            s_roimatrix{s, r}=nan(2, request.n_choices, log.nlevels, log.nlevels); s_roimatrix{s, r}=nan(2, request.n_choices, log.nlevels, log.nlevels);
            for ch=1:request.n_choices
                s_roimatrix{s, r}(1, ch,:,:)=d_roimatrix_cF(:,:, ch, r, s);
                s_roimatrix{s, r}(2, ch,:,:)=d_roimatrix_ct(:,:, ch, r, s);
            end
            
            % Assemble vector.  Conversion of EnvThreat x NTokens matrix to a vector
            %   Vector order: EnvThreat (descending threat) followed by NTokens (ascending N)
            k=1;
            for t=1:2
                for ch=1:request.n_choices
                    wv=(squeeze(s_roimatrix{s,r}(t,ch,:,:)))';
                    s_roivect{r}(s, k:k+length(wv(:))-1)=(wv(:))';
                    k=k+length(wv(:))+1;
                end
                k=k+2;
            end
        end
    end
    
    % PLOT
    log.plotsubjbetarange=[]; % -10 10]; % Plotting range only, not actually data filtering
    f.subplotcols=3; k=1; % How many columns?
    f.figwidth= 1400; f.figheight=1000;  f.subplot_VerHorz=[0.05 0.1];
    f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.01];  f.nancol=[0 0 0];
    fs=figure('Name', [log.FLmodel ': Subject betas in each cell'], 'NumberTitle', 'off', 'Position',[200,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
    if isempty(log.plotsubjbetarange)==0; disp(['Beta range for subject plot: ' num2str(log.plotsubjbetarange(1)) ' to '  num2str(log.plotsubjbetarange(2)) ]); end
    for r=1:log.n_rois
        subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        if isempty(log.plotsubjbetarange)==1;  imagescnan(s_roivect{r}, 'NaNColor', f.nancol);   else imagescnan(s_roivect{r}, log.plotsubjbetarange, 'NaNColor', f.nancol);  end
        colorbar; set(gca, 'XTick', []); title(log.rois{r})
        k=k+1;
    end
    
    if request.n_choices ==3
        % Total no. of trials in each cell for each subject
        subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagescnan(s_ntrialvect, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
        title('Total no. of trials in each cell'); k=k+1;
        
        % Mean no. of trials in each cell for each subject
        subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagescnan(s_ntrialvect(log.n_subjs+1,:), 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
        title('Mean no. of trials in each cell'); k=k+1;
        
        % Cell inclusion
        subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagescnan(s_ntrialvect>=request.flex_mintrials4inclusion, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
        title(['Cell included or not for each subject (min ntrials=' num2str(request.flex_mintrials4inclusion) ')']); k=k+1;
    end
end

%% OUTPUT 

logg.n_rois = log.n_rois; 
logg.rois=log.ROInames;



%% CHECK 
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
% %                       EnvThreat in reverse order (i.e. as per imagesc plot)
% s=4;
% for r=1:log.n_rois
%     a{r,1}=log.rois{r};
%     a(r,2:7)=num2cell( squeeze(mean(abs(d_roimatrix_cF(:,:,1,r,s)))) );
% end


end

