function [variables c] = m_f3_Event(datastruc, col, c)
% 'Event' first-level model - Each individual trial is modelled as an
% individual event. 
%
% 	datastruc: {1}=Conflict data, {2}=Control data, {3}=All data
%
%  -----------------------------------------------------------------------------
% Evaluate to debug: datastruc=ws.d;

%% 

% Each individual trial is represented by one regressor
data=datastruc{3}; c=1; 
for t=1:size(data,1)
    variables.names{c}=['tr' num2str( t  )];
    variables.onsets{c}=data(t, col.Onset_Offer);
    variables.durations{c}=data(t, col.Duration_Offer);
    variables.durations{c}=0; % Reset durations to 0
    %
    c=c+1;
end

% NO-INTEREST regressors
%       Info from exploration is NOT modelled
[variables c] = i_nointerest(datastruc{3}, {'Motor'}, variables, col, c );


end

