function fpar_plot( data, col)
% function fpar_plot( data, col)
%   Plot requested quantities (detailed in col)
%   To alter titles etc other variables, see script itself (fpar_plot, in
%   comp modelling folder)
% ------------------------------------



figure('Name', 'Task space variables', 'NumberTitle', 'off', 'Position', [100 100 1000 600], 'Color', 'w'); k=1;
f.FontSize=20; f.FontName='PT Sans Caption';   
f.FontSize_Title=35; 

    
% Translate data from trialstats to 6x6 cells  
dmat=cell(6,6); f.plotrows=2; f.plotcol=5;
for e=1:6
    for n=1:6
        dmat{7-e, n}=data(data(:,col.EnvThreat)== e/6 & data(:,col.NTokens)== n*2,:);
        dmat{7-e, n}=dmat{7-e, n}(1,:);
        
        % Individual variables
        et(7-e,n)=dmat{7-e, n}(:, col.EnvThreat);
        nt(7-e,n)=dmat{7-e, n}(:, col.NTokens);
        pl(7-e,n)=dmat{7-e, n}(:, col.pLoss);
        en(7-e,n)=dmat{7-e, n}(:, col.Entropy);
        ev(7-e,n)=dmat{7-e, n}(:, col.EV);
        if isfield(col, 'VExplore');  vex(7-e,n)=dmat{7-e, n}(1, col.VExplore); end
        if isfield(col, 'EntropyNTok'); entok(7-e,n)=dmat{7-e, n}(1, col.EntropyNTok); end
        if isfield(col, 'EntropyEV'); enev(7-e,n)=dmat{7-e, n}(1, col.EntropyEV); end
        if isfield(col, 'Conflict'); cf(7-e,n)=dmat{7-e, n}(1, col.Conflict); end
        %
        if isfield(col, 'StanDev'); sd(7-e,n)=dmat{7-e, n}(1, col.StanDev); end
        if isfield(col, 'vMeanVar'); vmv(7-e,n)=dmat{7-e, n}(1, col.vMeanVar); end
        if isfield(col, 'BinomVar'); bnv(7-e,n)=dmat{7-e, n}(1, col.BinomVar); end
        %
        if isfield(col, 'EVGain'); evg(7-e,n)=dmat{7-e, n}(1, col.EVGain); end
        if isfield(col, 'EVLoss'); evl(7-e,n)=dmat{7-e, n}(1, col.EVLoss); end
        if isfield(col, 'EVConflict'); evc(7-e,n)=dmat{7-e, n}(1, col.EVConflict); end
        if isfield(col, 'PavConflict'); pc(7-e,n)=dmat{7-e, n}(1, col.PavConflict); end
        if isfield(col, 'Misc'); misc(7-e,n)=dmat{7-e, n}(1, col.Misc); end
    end
end


%% Plot

% EnvThreat
subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(et,[0 1]); title('Env Threat', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)

% NTokens
subplot(f.plotrows, f.plotcol,k), k=k+1;  imagesc(nt,[2 12]);  title('N Tokens', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)

% pLoss
subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(pl ,[0 1]); title('pLoss', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)

% Entropy
subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(en); title('Entropy', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)

% #################################################################

% VExplore
if isfield(col, 'VExplore');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(vex); title('VExplore (pprob)', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
    % subplot(f.plotrows, f.plotcol,7); imagesc(vex, [min(min(vex)) max(max(vex))]); title('VExplore (min=0)')
    %         set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
end

% EntropyNTok
if isfield(col, 'EntropyNTok');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(entok,[min(min(entok)) max(max(entok))]); title('Entropy * NTok', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)

end

% #################################################################

% EV
subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(ev,[min(min(ev)) max(max(ev))]); title('EV', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)

if isfield(col, 'EVGain');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(evg); title('EV Gain', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end
if isfield(col, 'EVLoss');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(evl); title('EV Loss', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end
if isfield(col, 'EVConflict');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(evc); title('EV Conflict', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end

% Pavlovian conflict
if isfield(col, 'PavConflict');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(pc); title('Pav Conflict', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end

% #################################################################

% [Mean Variance theory quantities]

% Binomial variance
if isfield(col, 'BinomVar');
    subplot(f.plotrows, f.plotcol,8); imagesc(bnv); title('Binomial variance', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end

% Theoretical standev + mean-variance value
if isfield(col, 'StanDev');
    subplot(f.plotrows, f.plotcol,9); imagesc(sd); title('Standard Dev', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end
if isfield(col, 'vMeanVar');
    subplot(f.plotrows, f.plotcol,10); imagesc(vmv); title('Mean-Var value', 'FontSize', f.FontSize_Title, 'FontName', f.FontName)
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end

% #################################################################

if isfield(col, 'Misc');
    subplot(f.plotrows, f.plotcol,k), k=k+1; imagesc(misc); title('Misc')
    set(gca,'YTick', 1:6, 'YTickLabel', {'6/6' '5/6' '4/6' '3/6' '2/6' '1/6'},'XTick', 1:6, 'XTickLabel', 2:2:12); axis square; colorbar;
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
end

input('Plot OK. Continue?'); 

end

