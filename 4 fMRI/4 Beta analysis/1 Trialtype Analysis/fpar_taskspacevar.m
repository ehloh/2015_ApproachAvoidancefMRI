function [ d_design] = f_LoadBehaviour( request, col, subpars)
% [ d_design] = f_LoadBehaviour(request, col, subpars)
% Load task space per subject, including EnvThreat warped if specified
%
% Inputs:
%         request	 n_subjs
%         col           EnvThreat, NTokens
%         subpars    Results from behavioural modelling 
%                           subpars=[] to use canonical design
%                           subpars{1} = {parnames, subpars for cF}
%                           subpars{2} = {parnames, subpars for ct}
% Output:
%           d_design{subject, task}: trialstats for cF/ct (specified in col)
%                             nsubjs+1: average of all sub variables
%                             nsubjs+2: mean param vals.
%
% ##################################

request.whichvar={'EnvThreat';'NTokens'; 'pLoss';'Entropy';'EntropyNTok';'EV';'EVGain'; 'EVLoss'; 'EVConflict'; 'EnvThreatOrig'; 'PavConflict';};

%% Canonical task space

canon.d_par=nan(36,3); canon.d_par(:, [col.EnvThreat col.NTokens])=[(repmat(1:6, 1,6)'/6) sortrows(repmat(1:6, 1,6)')*2];
canon.dcf=fpar_conflict(canon.d_par, col);
canon.dct=fpar_control(canon.d_par, col);
canon.dcf(:, col.EnvThreatOrig)= canon.dcf(:, col.EnvThreat);  canon.dct(:, col.EnvThreatOrig)= canon.dct(:, col.EnvThreat);


%% Load task space + warp for each subject

if isempty(subpars)==0
    mod.parvals={subpars{1}{3} subpars{2}{3}};
    mod.parname={subpars{1}{2} subpars{2}{2}};
    mod.modname={subpars{1}{1} subpars{2}{1}};
    %
    for t=1:2 % Calculate warped task space
        if isempty(strcmp(mod.parname{t}, 'j'))==1;
            disp('No j param requested. Continuing without warp!');
            mod.parvals{t} = [mod.parvals{t}  ones(size(mod.parvals{t},1),1)];
            mod.parname{t}=[mod.parname{t}; {'j'}];
        end
        
        % Set up data holders
        d_etwarp{t}=num2cell( mod.parvals{t}(:,  find(strcmp(mod.parname{t}, 'j'))) ); % j param, all subjects
        d_etwarp{t}{request.n_subjs +2}= mean(cell2mat(d_etwarp{t})); % nsubjs+1: average of all sub variables, nsubjs+2: mean param vals.
        for p=1:length(request.whichvar)  % for mean
            eval(['d_etwarp{t}{request.n_subjs +1,2}.' request.whichvar{p} '=zeros(6,6);'])
        end
        
        % Initial task space
        whichsubs=[1:request.n_subjs request.n_subjs+2]';
        wt.des=nan(6*6, 10); wt.des(:,  [col.EnvThreat col.NTokens])=[sortrows(repmat((1:6)',6,1))/6 2*repmat((1:6)',6,1)];
        wt.des(:,  col.EnvThreatOrig) = wt.des(:,  col.EnvThreat);
        switch t
            case 1,   [wt.des] = fpar_conflict(wt.des, col);
            case 2,   [wt.des] = fpar_control(wt.des, col);
        end
        
        for ss=1: length(whichsubs)
            s=whichsubs(ss);
            
            % Fill in fixed params
            ws.mv=[];
            if sum(strcmp(mod.parname{t}, 'f'))==1,
                if s==request.n_subjs  +2,  ws.mv.FixedLoss=     mean(  mod.parvals{t}(:,  find(strcmp(mod.parname{t}, 'f')))  );
                else ws.mv.FixedLoss=     mod.parvals{t}(s,  find(strcmp(mod.parname{t}, 'f'))) ;
                end
            end
            
            % Apply distortion
            switch t
                case 1, [ ws.ov] = fcf_changeEnvThreat(power(wt.des(:, col.EnvThreat), d_etwarp{t}{s,1}), wt.des(:, col.NTokens), ws.mv);
                case 2, [ ws.ov] = fct_changeEnvThreat(power(wt.des(:, col.EnvThreat), d_etwarp{t}{s,1}), wt.des(:, col.NTokens), ws.mv);
            end
            
            % Record
            for p=1:length(request.whichvar)
                if strcmp(request.whichvar{p}, 'EnvThreatOrig')==0
                    eval(['ws.ov.d=ws.ov.' request.whichvar{p} ';'])
                    ws.ov.d= fliplr(reshape(ws.ov.d, 6,6))'   ;
                    
                    eval(['d_etwarp{t}{s,2}.' request.whichvar{p} '=ws.ov.d;']);
                    if s<request.n_subjs+1, eval(['d_etwarp{t}{request.n_subjs+1,2}.' request.whichvar{p} '=d_etwarp{t}{request.n_subjs+1,2}.' request.whichvar{p} ' + ws.ov.d;']); end
                end
            end
            d_etwarp{t}{s,2}.EnvThreatOrig = fliplr(reshape(wt.des(:,  col.EnvThreatOrig),6,6))';
            ws=[];
        end
        for p=1:length(request.whichvar) % nsubs+2: Average of all subjects' taskspace
            eval(['d_etwarp{t}{request.n_subjs+1,2}.' request.whichvar{p} '=d_etwarp{t}{request.n_subjs+1,2}.' request.whichvar{p} './request.n_subjs;']);
        end
    end
    
    % Output
    d_design=cell(request.n_subjs ,2);
    for s=1:request.n_subjs +2
        for t=1:2
            for v=1:length(request.whichvar)
                eval(['wv.v= d_etwarp{t}{s,2}.' request.whichvar{v} ';'])
                eval(['d_design{s, t}(:, col.' request.whichvar{v} ')=wv.v(:);']);
            end
        end
    end
    
else % Construct canonical task space per subject
    d_design=  [repmat({canon.dcf}, request.n_subjs+2 , 1) repmat({canon.dct}, request.n_subjs+2 , 1)] ;
end

%% Calculate additional values

% % WHICH quantities?
% %     vChoice, predChoice, pBestChoice, vChosenAnd, vGamble, vGamblePosNeg
% %     vBestUnchosenPosNeg, RejOrvBU, vSeeNosee
% %
% 
% % if isempty(subpars)==0, 
%     whichsubs= size(subpars{1}{3},1);
% % else whichsubs=1;
% % end
% 
% % Which if the most often chosen choice
% for t=1:2
%     for s=1:size(d_beh,1)
%         for e=1:6
%             for n=1:6 
%                 % Not set up to deal with multiple modes
%                 d_beh{s, 2+t}(e, n)=  mode(d_beh{s,t}(d_beh{s,t}(:, col.EnvThreat).*6==e & d_beh{s,t}(:, col.NTokens)/2==n, col.Choice));
%                 d_design{s,t}(d_design{s,t}(:, col.EnvThreatOrig).*6==e & d_design{s,t}(:, col.NTokens)./2==n,  col.Modalchoice) = d_beh{s, 2+t}(e, n);
%             end
%         end
%     end
% end
% col.vChoice=[col.vAccept  col.vReject  col.vExplore];
% 
% for t=1:2
%     wt.fixedpar.cF_FL=-12; wt.fixedpar.cF_EC=-2; wt.fixedpar.ct_FL=0; wt.fixedpar.ct_EC=-2;
%     
%     for s=1:whichsubs
% %         ws.par=subpars{t}{3}(s,:)
%         [ ws.transpar] = f_transpar(subpars{t}{2}, subpars{t}{3}(s,:), 'from');  % Apply inverse transformation to params
%         
%         % (1) Choice values
%         ws.des=d_design{s,t}; ws.des(:, col.Task)=t;  ws.des(:, col.Trialnum)=1:size(ws.des,1); 
%         eval(['[ ws.vchoice] = ' subpars{t}{1} '(ws.transpar, {[] ws.des wt.fixedpar, col});'])
%         ws.des(:, [col.vAccept  col.vReject  col.vExplore]) = squeeze(ws.vchoice);
%         
%         % (2) V(Modal choice) & V(Best Nonmodal)
%         for tn=1:size(ws.des,1)
%             ws.des(tn, col.vModalchoice) = ws.des(tn,  col.vChoice( ws.des(tn, col.Modalchoice) )); 
%             ws.des(tn, col.vBestNonmodal) = max( ws.des(tn,  col.vChoice(find(1- (ws.des(tn, col.Modalchoice) == [1 2 3])))));
%         end
%         
%         
%         
%         % vGamble
% % WHICH quantities?
% %     vChoice, predChoice, pBestChoice, vChosenAnd, vGamble, vGamblePosNeg
% %     vBestUnchosenPosNeg, RejOrvBU, vSeeNosee



end

