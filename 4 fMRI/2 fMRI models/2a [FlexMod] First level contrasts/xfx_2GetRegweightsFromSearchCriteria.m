function [RegWeights] = xfx_2GetRegweightsFromSearchCriteria(SearchCriteria,RegList)
%  [RegWeights] = f_2GetRegweightsFromSearchCriteria(SearchCriteria,RegList)
% Given search critera + names of regressors in a 1st-level model, returns
%   vector of regressor weights to be applied (in SPM'sc format) for each contrast
%
% Search criteria is subject-specific, depending on what choices subjects
% show per cell/trial type
%
%   SearchCriteria/RegWeights:
%         Col 1: Name of contrast
%         Col 2: Details for weights
%                   Col 1: Weight value (numeric)
%                   Col 2: Search criteria (inclusive) for this weight value [cell]
%                   Col 3: Names of matching regressors, that receive this weight value  [cell]
%                   Col 4: Weights vector, for regressors that receive this weight value  [double]
%         Col 3: Overall weights vector 
%
%
%  Note: 
%   (1) Search criteria is inclusive, returns a hit for all items in the
%   list that return non-empty strfind. 
%       - criteria specifying 'Bomb' and 'ct' will produce overinclusions
%         (e.g. criteria 'Bomb' returns positive hit for 'NoBomb')    
%
% 
% ----------------------------------------------------------------------------------------------

% Execute to debug: RegList=ws.RegList; 
% OR execute: RegList=RegLists{sc}; SearchCriteria=SearchCriterias{sc};
% SearchCriteria={'test'      {-1 {'cF_Accept'; 'cF_Reject'}; 2 {'cF_Explore'}};
%                           'test2'     {1 {'cF_Accept'; 'cF_Reject'}; -2 {'cF_Explore'}}     };  

RegWeights=SearchCriteria;
for c=1:size(SearchCriteria,1) % For all contrasts
    disp(['Contrast ' num2str(c) '  ' SearchCriteria{c,1} ' -----------------------------------------------' ])
    disp(['No. of different weight values to apply: ' num2str(size(SearchCriteria{c,2},1))])
    %
    RegWeights{c,1}=SearchCriteria{c,1};
    RegWeights{c,2}=SearchCriteria{c,2};
    RegWeights{c,3}=zeros(length(RegList),1);
    
    for v=1:size(SearchCriteria{c,2},1) % (1) Apply all distinct weight-values separately
        disp(['Weight value #' num2str(v) ': ' num2str(SearchCriteria{c,2}{v,1}) '  --------'])
        disp('Search criteria for this weight value:'); disp(SearchCriteria{c,2}{v,2}); 
        wv.weights=zeros(length(RegList),1);
        if SearchCriteria{c,2}{v,1}==0;   SearchCriteria{c,2}{v,1}=0.001;  RegWeights{c,2}{v,1}=0.001; end % 0-weighted contrasts changed to 0.001
        
        for sc=1: length(SearchCriteria{c,2}{v,2}) % (2) Apply all search criteria (find target regressors)
            ws.w=strfind(RegList, char(SearchCriteria{c,2}{v,2}(sc)));
            ws.w=double(~cellfun(@isempty,ws.w));  % Temporary Non-zero weights applied to all regressors that find a string match
            %   
            wv.weights=wv.weights+ws.w;
            disp(['Search criteria #' num2str(sc) ':' char(SearchCriteria{c,2}{v,2}(sc)) ': ' num2str(sum(ws.w)) ' regressors weighted'])
            ws=[];
        end
        
        % Convert temporay Non-zero weight to specified real weight value
        wv.weights=wv.weights.*SearchCriteria{c,2}{v,1}; % If weight value is 0, this will flag an error. One way of resolving this is to mean-centre all parameters
        
        % Check for errors
        if sum(wv.weights.*RegWeights{c,3})~=0  
            % Is there an overlap between these weighted regressors and those
            %   that are already weighted? There really shouldn't be.
            CulpritReg=RegList(find(sum(wv.weights.*RegWeights{c,3})~=0));
            disp('Culprit Reg:'); disp(CulpritReg); error('Search criteria tries to weight a certain regressor more than once (CulpritReg)')
        elseif sum(wv.weights)==0 && SearchCriteria{c,2}{v,1}~=0
            % No matches for this search criteria?
            CulpritReg=RegList(find(sum(wv.weights.*RegWeights{c,3})~=0));
            disp('Bad search criteria:'); disp(CulpritReg); error('Search criteria finds no matching regressors (CulpritReg)')
        end
        
        % Record
        RegWeights{c,2}{v,3}=RegList(find(wv.weights));
        RegWeights{c,2}{v,4}=wv.weights; 
        disp([num2str(sum(wv.weights~=0)) ' regs with this weight-value ' num2str(SearchCriteria{c,2}{v,1}) ' applied:']); disp(RegWeights{c,2}{v,3});disp(' ')
        
        % Integrate regressor-weights of this value with regressor weights
        % for this contrast
        RegWeights{c,3}=RegWeights{c,3}+RegWeights{c,2}{v,4};
    
        wv=[]; 
    end % end of this weight-value

end

end

