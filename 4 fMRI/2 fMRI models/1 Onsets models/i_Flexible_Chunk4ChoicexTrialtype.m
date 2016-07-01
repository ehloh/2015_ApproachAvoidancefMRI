function [variables c] = i_Flexible_Chunk4ChoicexTrialtype(data, prefix, choicenames,variables, col, c )
% [variables c] = i_Flexible_Chunk4ChoicexTrialtype(data, prefix, choicenames,variables, col, c )
% 'Trial Type' first-level model - allows for flexible analysis at 2nd level
%  4 cells chunked, i.e. 3x3 rather than 6x6 design
%
% ------------------------------------------------------------------------------------------

% Evaluate to debug: data=datastruc{1}, prefix='cF_'; choicenames={'Accept';'Reject';'Explore'}; 


%% Set up design fow new cells     
%       3x3 cell; In each cell, 4 further cells each = [EnvThreatLevel NTokensLevel] of original 6x6-cell to be included

chunk=cell(3,3); 
chunk{1,1}= {[1 1]   [1 2]   [2 1]   [2 2]};  % i.e. Cell 1-1: t1-1, t1-2, t2-1, t2-2
chunk{1,2}= {[1 3]   [1 4]   [2 3]   [2 4]}; 
chunk{1,3}= {[1 5]   [1 6]   [2 5]   [2 6]}; 
chunk{2,1}= {[3 1]   [3 2]   [4 1]   [4 2]};
chunk{2,2}= {[3 3]   [3 4]   [4 3]   [4 4]}; 
chunk{2,3}= {[3 5]   [3 6]   [4 5]   [4 6]}; 
chunk{3,1}= {[5 1]   [5 2]   [6 1]   [6 2]};
chunk{3,2}= {[5 3]   [5 4]   [6 3]   [6 4]}; 
chunk{3,3}= {[5 5]   [5 6]   [6 5]   [6 6]}; 

%% Create regressor variables

for choice=1:3 % Choices
%     checkvar=repmat({nan(3,3)}, 10,1);
    
    for ee=1:3
        for nn=1:3
            
            % Compile sample
            datatrial=[];
            for cc=1:length(chunk{ee,nn})
                datatrial=[datatrial;  data(data(:,col.EnvThreat).*6==chunk{ee,nn}{cc}(1) & data(:,col.NTokens)./2==chunk{ee,nn}{cc}(2),:)];
            end
            datatrial=datatrial(datatrial(:,col.Resp1)==choice,:);
            
            if isempty(datatrial)==0
                variables.names{c}=[prefix choicenames{choice} '_e' num2str(ee) '-n' num2str(nn)];
                variables.onsets{c}=datatrial(:, col.Onset_Offer);
                variables.durations{c}=datatrial(:, col.Duration_Offer);
                variables.durations{c}=zeros(size(variables.onsets{c})); % Reset durations to 0
                c=c+1;
                
%                 % Check
%                 checkvar{1}(4-ee, nn)=mean(datatrial(:, col.Resp1)==1);
%                 checkvar{2}(4-ee, nn)=mean(datatrial(:, col.Resp1)==2);
%                 checkvar{3}(4-ee, nn)=mean(datatrial(:, col.Resp1)==3);
%                 %
%                 checkvar{4}(4-ee, nn)=mean( datatrial(:, col.EnvThreat));
%                 checkvar{5}(4-ee, nn)=mean( datatrial(:, col.NTokens));
%                 checkvar{6}(4-ee, nn)=mean( datatrial(:, col.pLoss));
%                 checkvar{7}(4-ee, nn)=mean( datatrial(:, col.Entropy));
%                 checkvar{8}(4-ee, nn)=mean( datatrial(:, col.VExplore));
%                 checkvar{9}(4-ee, nn)=mean( datatrial(:, col.EV));
%                 checkvar{10}(4-ee, nn)=mean( datatrial(:, col.OutcomeMean));
%                 checkvar{11}(4-ee, nn)=mean( datatrial(:, col.RT1));
            end
        end
    end
    
%     cchoices{choice}=checkvar;
end


% % Part of checks
% for cc=1:3
%     figure
%     for p=1:11
%         subplot(4,3, p);
%         imagescnan(cchoices{cc}{p}, 'NaNColor', [0 0 0]);axis square; colorbar
%     end
% end



end

