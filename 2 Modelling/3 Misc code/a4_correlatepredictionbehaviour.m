% analyse4_correlateplots
%     Perform 2D correlations on plots of behaviour vs model predictions (individual subjects & overall)
clear all; close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Conflict Exploration/3 Analysis';
where.where='H:\5 Conflict Exploration\3 Analysis';
where.plots=[where.where filesep 'Fit computational models'];

% Select data
details.behaviour='res_21-Feb-2013_plotbehaviour.mat';
details.predictions='res_21-Feb-2013_plotpredictions.mat';


for o1=1:1 %% Set up


% Load
beh=load([where.plots filesep details.behaviour]);
pred=load([where.plots filesep details.predictions]); 

% Models to be run are determined by settings of prediction-plots
details.whichmodels=pred.details.whichmodels;
details.n_subs=pred.details.n_subjs;
details.n_models=size(details.whichmodels,2);

disp('====================================')
disp('Run correlations on choice contingency plots, comparing')
disp('predictions & behaviour (group & individual level)')
disp(['No. of subjects: ' num2str(details.n_subs)])
disp('Models:')
disp(details.whichmodels')
disp('====================================')

end

%%

r_corr=cell(details.n_models,1); r_pval=cell(details.n_models,1);
for m=1:details.n_models
    % r_corr{m} shows 'r' for each model, r_corr{m} shows associated p value
    % for each model
    %       (row=subject, col 1=accept, col 2= reject, col 3=explore, col 4=mean over choice)
    r_corr{m}=zeros(details.n_subs+1,4);
    for s=1:details.n_subs+1
        % Group means are in last row
        for c=1:3
            [ws.r ws.p]=corrcoef(beh.plots{c}{s}, pred.plots{c}{s,m});
            r_corr{m}(s,c)=ws.r(1,2);
            r_pval{m}(s,c)=ws.p(1,2);
        end
        r_corr{m}(s,4)=mean(r_corr{m}(s,1:3));
        r_pval{m}(s,4)=mean(r_pval{m}(s,1:3));
    end
end

%% Plot R values

% Format for plotting
plots=cell(4,1); mplots=cell(4,1);
for m=1:details.n_models
    % plots{c}: row=subject, col=model
    for c=1:3
        plots{c}(:,m)=r_corr{m}(1:end-1, c);
        mplots{c}(:,m)=r_corr{m}(end, c);  % Organize data for mean plots
    end
    % Overall choice
    plots{4}(:,m)=r_corr{m}(1:end-1, 4);
    mplots{4}(:,m)=r_corr{m}(end, 4);
end

% Plot
figure('Name', 'Mean r values (correlating predictions & behaviour on the subject level)','NumberTitle', 'off', 'Position', [600,400,1200,400]); 
c=1;
subplot(1,4,c); boxplot(plots{c}, 'labels', details.whichmodels, 'plotstyle','compact'); title('Accept');
c=2;
subplot(1,4,c); boxplot(plots{c}, 'labels', details.whichmodels, 'plotstyle','compact'); title('Reject')
c=3;
subplot(1,4,c); boxplot(plots{c}, 'labels', details.whichmodels, 'plotstyle','compact'); title('Explore')
c=4;
subplot(1,4,c); boxplot(plots{c}, 'labels', details.whichmodels, 'plotstyle','compact'); title('Mean r values')



