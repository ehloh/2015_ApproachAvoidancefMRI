function [variables c] = i_nointerest(data, NoInterestVar, variables, col, c )
% [variables c] = i_nointerest(data, NoInterestVar, variables, c )
%   Create variables for no-interest regressors
%
%   All variables: {'Error'; 'Motor'; 'OutcomePresented'; 'Info'}
% ------------------------------------------------------------------------------------

% Execute to debug: data=datastruc{3}; NoInterestVar={'Error'; 'Motor'; 'OutcomePresented'; 'Info'}; 

%% Format data

% Errors
d.ErrorOnset=data(data(:,col.TrialValid)==0,col.Onset_Offer); % ERROR: Error marked from start of offer to end of trial
d.ErrorDuration=data(data(:,col.TrialValid)==0,col.Onset_TrialEnd)-data(data(:,col.TrialValid)==0,col.Onset_Offer);

% Motor
d.Motor1Onset=data(isnan(data(:,col.Onset_Motor1))==0, col.Onset_Motor1);
d.Motor2Onset=data(isnan(data(:,col.Onset_Motor2))==0, col.Onset_Motor2);
d.Motor1Choice=data(isnan(data(:,col.Onset_Motor1))==0, col.Resp1);
d.Motor2Choice=data(isnan(data(:,col.Onset_Motor2))==0, col.Resp2);
d.MotorOnset=vertcat(d.Motor1Onset,d.Motor2Onset);
d.MotorDuration=zeros(size(d.MotorOnset,1),1);
d.MotorChoice=vertcat(d.Motor1Choice,d.Motor2Choice);

% Information from Exploring
d.InfoOnset=data(data(:,col.OutcomePresented)==1 & data(:,col.Resp1)==3, col.Onset_Info);
d.InfoDuration=data(data(:,col.OutcomePresented)==1 & data(:,col.Resp1)==3, col.Duration_Info);
d.InfoxEntropy=data(data(:,col.OutcomePresented)==1 & data(:,col.Resp1)==3, col.Entropy);

% Outcome (Presented + Pmod Amount)
d.OutcomePresentedOnset=data(data(:,col.OutcomePresented)==1, col.Onset_Outcome);
d.OutcomePresentedDuration=data(data(:,col.OutcomePresented)==1, col.Duration_Outcome);
d.OutcomePresentedMagnitude=data(data(:,col.OutcomePresented)==1, col.OutcomeMagnitude);
d.OutcomeOnlyOnset=d.OutcomePresentedOnset;
d.OutcomeOnlyDuration=d.OutcomePresentedDuration;

% Adjust ALL durations = 0
d.ErrorDuration=zeros(length(d.ErrorDuration),1);
d.MotorDuration=zeros(length(d.MotorDuration),1);
d.OutcomePresentedDuration=zeros(length(d.OutcomePresentedDuration),1);
d.InfoDuration=zeros(length(d.InfoDuration),1);

%% Construct no interest regressors

for p=1:length(NoInterestVar)
    variables.names{c}=['n_' NoInterestVar{p}];
    eval(['variables.onsets{c}=d.' NoInterestVar{p} 'Onset;'])
    eval(['variables.durations{c}=d.' NoInterestVar{p} 'Duration;'])
    
    % Parametric modulators (for particular no-interest variables)
    if strcmp(NoInterestVar{p}, 'OutcomePresented')==1
        variables.pmod(c).name{1}='n_OutcomePresentedMagnitude';
        variables.pmod(c).param{1}=d.OutcomePresentedMagnitude;
        variables.pmod(c).poly{1}=1;
    elseif strcmp(NoInterestVar{p}, 'Motor')==1
        %             variables.pmod(c).name{1}='n_MotorChoice';   % Choices are now put in as discrete regressors, to allow for direct comparison
        %             variables.pmod(c).param{1}=d.MotorChoice;
        %             variables.pmod(c).poly{1}=1;
    elseif  strcmp(NoInterestVar{p}, 'Info')==1
        variables.pmod(c).name{1}='InfoxEntropy';
        variables.pmod(c).param{1}=d.InfoxEntropy;
        variables.pmod(c).poly{1}=1;
    end
    %
    c=c+1;
    
end

end

