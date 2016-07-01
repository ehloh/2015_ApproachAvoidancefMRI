clear all; clc; close all hidden; path(pathdef)

tasktype=2;    % cF or ct
WhereRes=[filesep '2 High iter']; 
% WhereRes=[];
% WhereRes=[filesep 'CompareFits'];
% WhereRes=[filesep 'WorkingFits'];

% ##################################################
% winmod={'bpm16_feow'; 'bpm11_euw';}; % Which models?
winmod={'bpmi16_feow'; 'bpmi11_euw';}; % Which models?
% winmod={'bpma01';};
fmod=winmod{tasktype}; gmod=fmod; switch tasktype; case 1; tn='cF';  case 2; tn='ct'; end; cd(['2 Analysis inputs' WhereRes])

% % cF 
% f=load(['res_fitmodels_' tn ' (22-Jul-2014) bpmi16 1000iter.mat']);
% g=load(['res_fitmodels_' tn ' (22-Jul-2014) bpmi16 bpm16 1000iter.mat']);

% % ct 
% f=load(['res_fitmodels_' tn ' (22-Jul-2014) bpmi11 bpm11 1000iter.mat']);
% g=load(['res_fitmodels_' tn ' (22-Jul-2014) top10 200iter.mat']);

% cF vs ct #########################
f=load(['2 Analysis inputs' WhereRes filesep 'res_fitmodels_cF (22-Jul-2014) bpmi16 2000iter recomb.mat']);
g=load(['2 Analysis inputs' WhereRes filesep 'res_fitmodels_ct (22-Jul-2014) bpmi11 bpm11 1000iter.mat']);
fmod=winmod{1}; gmod=winmod{2};

for o1=1:1 % Set up data (f=r_res, fd=details.models, fm=model num, fsp=subject par vals)
f.r_res=sortrows(f.r_res,1);
fm=find(strcmp(f.r_res(:,1), fmod));    if isempty(fm)==1; error(['Could not find requested model (' fmod ') in fit res! ']); end
fd=sortrows(f.details.models,1);
g.r_res=sortrows(g.r_res,1);
gm=find(strcmp(g.r_res(:,1), gmod));   if isempty(gm)==1; error(['Could not find requested model (' gmod ') in grid res! ']); end
gd=sortrows(g.details.models,1);
%
f=f.r_res;
g=g.r_res;
fsp=f{fm,2}(:,4:end);
gsp=g{gm,2}(:,4:end);

end

%%

% for cF and ct, adjust ct params to line up
gsp(:, 6:7)=gsp(:, 5:6);
gsp(:,5)=0;
gd{2, 3}(6:7)=gd{2, 3}(5:6);


% Range(fit)
for p=1:fd{fm,2}
    disp([fd{fm,3}{p} '   : ' num2str(min(fsp(:,p)),2) ' to '  num2str(max(fsp(:,p)),2)])
end

% Plot distribution of parameters 
figure; set(gcf, 'Color', 'w');
% labels={'fit';'grid'};
labels={'cF';'ct'};
for p=1:fd{fm,2}
    subplot(4,  fd{fm,2},  p)  % Hist f
    hist(fsp(:,p)); 
    title(['[' labels{1} '] ' fmod ' -  ' fd{fm,3}{p}])
    
    subplot(4, fd{fm,2},   p  + fd{fm,2}*1)  % Hist g
    hist(gsp(:,p)); 
    title(['[' labels{2} '] ' gmod ' -  ' gd{gm,3}{p}])
    
    subplot(4, fd{fm,2},   p  + fd{fm,2}*2)  % Plot f vs g
    scatter(fsp(:,p), gsp(:,p))
    xlabel(labels{1});  ylabel(labels{2});
    
    subplot(4, fd{fm,2},   p  + fd{fm,2}*3)  % Range
    text(0,1.2, [fd{fm,3}{p} '   : ' num2str(min(fsp(:,p)),2) ' to '  num2str(max(fsp(:,p)),2)])
    axis off;
end


% nLLs
% nlls=[f{fm,2}(:,2)   g{gm,2}(:,2)  f{fm,2}(:,2)-   g{gm,2}(:,2)];
% disp(nlls(1,1)-nlls(1,2))



figure,scatter(fsp(:, 3).*fsp(:, 4),gsp(:, 3).*gsp(:, 4))
