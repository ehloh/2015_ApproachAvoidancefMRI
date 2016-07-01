function [variables c] =  i_choiceRT_conditions_cluster(data, variables, prefix,  choicenames, col, c , InnerCluster)
% Create regressors for Choice, each modelled as a condition
% But design split into the N cells where Explore peaks, and the rest of
% the design
%
%   choicenames: cell with names, should correspond to code for choice
%                           (1=Accept/NoBomb, 2=Reject/Bomb, 3=Explore)
% ---------------------------------------------------------------

% Evaluate to debug: datastruc=ws.d; data=datastruc{1}; prefix='cF_'; choicenames={'Accept';'Reject';'Explore'}; % InnerCluster=[2 5; 2 6; 3 5; 3 6];

% Define Clusters (Col=EnvThreat, Row=NTokens)
OuterCluster=vertcat([1*ones(6,1) (1:6)'],[2*ones(6,1) (1:6)'],[3*ones(6,1) (1:6)'],[4*ones(6,1) (1:6)'],[5*ones(6,1) (1:6)'],[6*ones(6,1) (1:6)']);
for i=1:size(InnerCluster,1) % Define outer cluster - excluding inner cluster
    OuterCluster(find(OuterCluster(:,1)==InnerCluster(i,1) & OuterCluster(:,2)==InnerCluster(i,2)),:)=[];
end

% Note: Clusters are specified with NTokens as 1:6, but data's NTokens column in in NTokens (2:2:12) 
% rather than pairs. In compiling the samples below (data_inner=, data_outer= ), 
% Cluster instructions' NToken value must be multiplied by 2.
% 

%% InnerCluster

% Compile data for this cluster
data_inner=[];
for i=1:size(InnerCluster,1)
    data_inner=vertcat(data_inner, data(data(:,col.EnvThreat)*6==InnerCluster(i,1) & data(:,col.NTokens)==InnerCluster(i,2)*2,:));
end
data_inner=sortrows(data_inner,col.Onset_Offer);
if isempty(data_inner)==1; error('No trials in the in zone!'); end 


% keyboard

% Construct regerssors
for p=1:length(choicenames)
    choiced=data_inner(data_inner(:,col.Resp1)==p,:);
    
    variables.names{c}=['in_' prefix choicenames{p}];
    variables.onsets{c}=choiced(:,col.Onset_Motor1); % Onset= Locked to time of decision
    variables.durations{c}=choiced(:,col.Duration_Choice); 
    variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0
    
    % Parametric modulators - Choice RT
    variables.pmod(c).name{1}=['in_' prefix choicenames{p} '_RT'];
    variables.pmod(c).param{1}=choiced(:,col.RT1);
    variables.pmod(c).poly{1}=1;
    
    c=c+1;
end


%% OuterCluster

% Compile data for this cluster
data_outer=[];
for i=1:size(OuterCluster,1)
    data_outer=vertcat(data_outer, data(data(:,col.EnvThreat)==OuterCluster(i,1) & data(:,col.NTokens)==OuterCluster(i,2)*2,:));
end
data_outer=sortrows(data_outer,col.Onset_Offer);

% Construct regerssors
for p=1:length(choicenames)
    choiced=data_outer(data_outer(:,col.Resp1)==p,:);
    
    variables.names{c}=['out_' prefix choicenames{p}];
    variables.onsets{c}=choiced(:,col.Onset_Motor1); % Onset= Locked to time of decision
    variables.durations{c}=choiced(:,col.Duration_Choice); 
    variables.durations{c}=zeros(size(variables.durations{c})); % Reset  durations to 0
    
    % Parametric modulators - Choice RT
    variables.pmod(c).name{1}=['out_' prefix choicenames{p} '_RT'];
    variables.pmod(c).param{1}=choiced(:,col.RT1);
    variables.pmod(c).poly{1}=1;
    
    c=c+1;
end

end

