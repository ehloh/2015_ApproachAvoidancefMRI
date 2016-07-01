% Format betas for analysis (mean centre, append behaviour scores)
clear all; close all hidden; clc


% Requested
log.subjects={};
request.meancentre=1;
request.addbeh=1;

% Settings
ncols_perroi=4;



something is wrong with this script for processing c7 models. c3 model appears to be okay. 


% Data details
behfile_where='I:\5 Explore fMRI\1 Behavioural data\Group behaviour files';
behfile_name='Behavioural profile for correlation'; % excel
%
betafile_where='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\m_c7_Cluster6CompeteFull_XUVPEN\choice_cluster2x2\ROI\2 Selected c7 rois';
betafile_name='(03-Nov-2013) Extracted betas'; % txt

%% Mean centre


% Load betafile --------------------
%   Assumed format: col=beta row, row=subject, 1st col and 1st row subject and beta name
w.betafile=importdata([betafile_where filesep betafile_name '.txt']);
betafile=[w.betafile.textdata(1,:); [w.betafile.textdata(2:end,1) num2cell(w.betafile.data)]];
[log.subjects log.n_subjs betafile]=f_selectsubjects(betafile, log.subjects, [betafile [{'ok'}; num2cell(ones(size(betafile,1)-1,1))]], 'ok');

ncols=size(betafile,2); nrows=size(betafile,1); % Count out parameters
nrois=(ncols-1)/ncols_perroi; 
betalist=strtrim(betafile(1,2:end)'); openvar betalist

if request.meancentre
    
    disp(['Assumed format for mean centreing:   ' num2str(nrois) ' rois,  ' num2str(ncols_perroi) ' beta columns per roi']);
    for r=1:nrois
        
        % Load details
        wr.startcol=1+(r-1)*ncols_perroi+1;
        wr.endcol=1+r*ncols_perroi;
        wr.betatitles=betafile(1,wr.startcol:wr.endcol);
        disp(['ROI #' num2str(r) ' includes betas:']);
        disp(char(wr.betatitles')); disp(' ')
        
        wr.dat=betafile(2:end,wr.startcol:wr.endcol);
        wr.dat=num2cell(    cell2mat(wr.dat)-repmat(mean(cell2mat(wr.dat),2), 1, ncols_perroi));
        betafile(2:end,wr.startcol:wr.endcol)=wr.dat;
        
        
        wr=[];
    end

end

%% Add behavioural scores

if request.addbeh==1;
    % Load profiles
    [w.n w.t w.behinhibit]=xlsread([behfile_where filesep behfile_name],'BehInhibition');
    [w.n w.t w.behexplore]=xlsread([behfile_where filesep behfile_name],'Explore');
    
    % Apply subject selection according to beta text betafile
    [log.subjects log.n_subjs w.behinhibit]=f_selectsubjects(w.behinhibit, strtrim(betafile(2:end,1)), [w.behinhibit [{'ok'}; num2cell(ones(size(w.behinhibit,1)-1,1))]], 'ok');
    [log.subjects log.n_subjs w.behexplore]=f_selectsubjects(w.behexplore, strtrim(betafile(2:end,1)), [w.behexplore [{'ok'}; num2cell(ones(size(w.behexplore,1)-1,1))]], 'ok');
    behfile=[w.behinhibit w.behexplore(:,2:end)];
    
    % Add beh to beta (check subject order)
    file=[behfile betafile(:,2:end)];
end

%% Export

if request.meancentre==1 && request.addbeh==1
    printok=print2txt(betafile_where, [betafile_name ' mean centred w beh'], file);
elseif request.meancentre==1
    printok=print2txt(betafile_where, [betafile_name ' mean centred'], file);
elseif request.addbeh==1
    printok=print2txt(betafile_where, [betafile_name ' w beh'], file);
end
disp(printok)



