function [ V ] = bpmi04_uw(x, modelinputs)
%   Inoptimal value function (post-prob) + distortion of V(Accept|~See) 
%     + Variable exploration bonus (w) x Uncertainty (u)
%
%  - Inputs
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[] data   fPar   col}
%                                   data:   Data from all trials     (row=trialnum, col=fixed params in col)
%                                   col:      Column specifications for data
%                                   fPar:    Fixed parameters
%
%  -  Output 
%        V:        Value associated with Accept (1), Reject (2), Explore (3)
%                        (row=parampoint, col=trialnum, 3rd dimension=Accept/Reject/Explore)
%
% ---------------------------------------------------------------------------------------------------------------------------------
% Free parameters 'x':        x(1)= softmax beta        (Applied outside of value fxn)
%                                        x(2)= softmax epsilon  (p)
%                                        x(3)= Unoptimal posterior probability multipler  (m)
%                                        x(4)= V(Accept|~See) distortion (a)
%                                        x(5)= Variable exploration bonus (w)
%-------------------------------------------------------------------------------------------------------------------------------

data=modelinputs{2}; fPar=modelinputs{3}; col=modelinputs{4};

ExploreBonus=x(:,5);

%% EV calculation using basic value function

% Get parameters in model
modname=mfilename; modname(strfind(modname, '_')-2:strfind(modname, '_'))=[];

[ V ] = bpmi01(x, {modname data, fPar,col});

% Add variable Exploration bonus x Uncertainty/Entropy
nTrials=size(data,1); nPP=size(x,1); 
V(:,:,3)=V(:,:,3)  +   repmat(ExploreBonus, [1 nTrials]).*repmat(data(:,col.Entropy)', [nPP 1]);

end

