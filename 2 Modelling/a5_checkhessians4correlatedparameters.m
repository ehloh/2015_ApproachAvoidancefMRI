% Use hessian matrices to check for correlated parameters
%           Covariance matrix =matrix inverse (hessian matrix); 
%           Correlation between 2 parameters= covariance/( SD(1st parameter) *  SD(2nd parameter) )
clear all; close all hidden; clc

% where='I:\4 Explore experiment\3 Analysis\4 Fit computational models';
where='/Volumes/PENNYDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
load([where filesep '(16-Apr-2013) res_fitmodels_control.mat'])

plot=1;

%% Compute correlation matrices
%       r_corr: Col1=model, col2=subject-specific correlation matrices, Col 3=Mean correlation matrix

% Obtain correlation matrices
r_corr=cell(size(r_res,1),3);
for m=1:size(r_res,1);
    disp(['Model ' num2str(m) '  - ' r_res{m,1} ' ############## '])
    r_corr{m,1}=r_res{m,1};
    if size(r_res{m,4}{1},1)==1 % FOR all models ####
        disp('Skipped (single parameter)')
    else
        r_corr{m,2}=cell(size(r_res{m,2},1),1);
        for s=1:size(r_res{m,2},1) % FOR all subjects ####
            disp(['Subject ' num2str(s)])
            
            
            % Begin computation --------------------------------
            ws.hess=r_res{m,4}{s}; 
            ws.hess(ws.hess(:,:)==0)=1*1e-10;   % Correct for Hessian values of 0?
            ws.cov=inv(ws.hess);
            ws.sds=nan*zeros(1,size(ws.hess,1)); ws.cor=nan*zeros(size(ws.hess,1),size(ws.hess,1));  % Correlation between 2 parameters= covariance/( SD(1st parameter) *  SD(2nd parameter) )
            for i=1:size(ws.hess,1)  % Calculate SD: SD=
                ws.sds(i)=sqrt(ws.cov(i,i));
            end
            for i=1:size(ws.hess,1)  % Compute correlation matrix 
                for j=1:size(ws.hess,1) % 
                    ws.cor(i,j)=ws.cov(i,j)/(ws.sds(i)*ws.sds(j));
                end
            end
            
            % Checks ----------------------------------------
            ws.diagonal=[];for i=1:size(ws.hess,1); ws.diagonal=vertcat(ws.diagonal, ws.cor(i,i)); end;
            r_corr{m,2}{s,1}=[];
             if isreal(ws.cor)~=1 
                 disp('Subject skipped - Computed Correlation matrix has non-real numbers')
             elseif isreal(ws.cov)~=1 
                 disp('Subject skipped - Covariance matrix has non-real numbers')
%              elseif mean(abs(ws.diagonal))~=1
%                  disp('Subject skipped - wtf ?')
             else
                 r_corr{m,2}{s,1}=ws.cor;
             end
            diag{m}{s}=ws.diagonal;
            ws=[];
            
            
            % Using function
%             [ws.sd ws.fcor]=cov2corr(inv(r_res{m,4}{s}));
%             r_corr{m,2}{s,1}=ws.fcor;

        end
    end
end

% input('Done computing correlations. Plot correlation matrices?    ')

%% Plot correlation matrices

plotTopN=5;
if plot==1
    
    % Plot mean correlation coefficients -----------------
    figure('Position', [200 200 500 800], 'Name', 'Mean correlation plots for each model'); set(gcf,'color','w')
    for m=1:plotTopN 
        if isempty(r_corr{m,2})==0
            
            % Calculate average (not absoluted)
            wm.s=zeros(size(r_corr{m,2}{1})); k=0;
            for s=1:size(r_corr{m,2},1)
                if isreal(r_corr{m,2}{s})==1 && isempty(r_corr{m,2}{s})==0
                    wm.s=wm.s+r_corr{m,2}{s}; k=k+1;
                end
            end
            wm.s=wm.s/size(r_corr{m,2},1);
            r_corr{m,3}=wm.s;
            
            % Plot 
            subplot(plotTopN,2, (m-1)*2+1)
            axis 'off' ; text(0,0, r_corr{m,1})
            subplot(plotTopN,2,(m-1)*2+2)
            imagesc(r_corr{m,3}); axis square; axis 'off'; colorbar
        end
    end
    
    
    % Plot Individual-subject's correlation coefficients -----------------
    figure('Position', [205 100 1500 1200], 'Name', 'Correlation plots for each subject (possibly limited)'); set(gcf,'color','w')
    w.nsubs=size(r_corr{m,2},1);
    for m=1:plotTopN 
        % Labels
        subplot(plotTopN+1, w.nsubs+1 , (m-1)*(w.nsubs+1)+1)
        text(0,0, r_res{m,1}); axis 'off';
        
        % Plot correlation coefficients
        if isempty(r_corr{m,2})==0
            for s=1:size(r_corr{m,2},1)
                if isreal(r_corr{m,2}{s})==1 && isempty(r_corr{m,2}{s})==0
                    subplot(plotTopN+1, w.nsubs+1 , (m-1)*(w.nsubs+1)+s+1)
                    imagesc(r_corr{m,2}{s}, [-0.5 0.5]);axis square; axis 'off'; colorbar
                end
            end
        end
    end
    
end

