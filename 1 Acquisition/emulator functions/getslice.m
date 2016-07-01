function [slices, times] = getslice( port )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WARNING - This version of getslice is MODIFIED %%
%%%         DO NOT USE this version in MRI         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [slices, times] = getslice( port )
%
% Returns the most recent slice numbers with time stamps,
% NB This WILL WAIT for complete slice information to arrive in the serial port buffer!
%
% See also:
%   config_serial, waitslice, logslice.
%
% Cogent2000 function
% 
% Version 2.0 09-10-2009

% Version History
% 2.0, 09-10-2009, E.F. - Added TRIO-STIM to list of "banned" computers
% 1.2, 22-07-2008, E.F. : emulscan mode altered for emulscan v?? and later
% 1.1, 23-08-2005? E.F. : emulscan mode added
% 1.0, 24-09-2002? E.F. : Creation

% or returns [-1, -1] if no data available.
% NB This will NOT WAIT for slice information if the serial port buffer is empty!
% It WILL WAIT if one of the two slice information bytes has been received.

drawnow
global cogent

switch port
    case 0 % Emulscan...
    if any( strcmpi( {'ALLEGRA-STIMULU' 'SONATA-STIM' 'TRIO-STIM'}, getenv( 'COMPUTERNAME' )))
        stop_cogent
        clear global cogent
        error( 'DO NOT use this version of getslice on MRI stim PCs')
    end % if any
    global emulscan
    slices = emulscan.slice;
    times = emulscan.time;
    
otherwise  % Standard "serial port" scanner pulses...
    temp = cogent.serial{port};
    ix = []; values = []; times = [];
    
    while isempty( ix ) % MUST HAVE some closely spaced bytes 
        values0 = []; times0 = [];
        while length( values0 ) < 2 % MUST HAVE 2 bytes or more, to calculate their spacing
            values1 = []; times1 = [];
            while isempty( values1 ) % MUST HAVE some bytes to work with (please?)
                [ values1, times1 ] = CogSerial( 'GetEvents', temp.hPort );
                drawnow; % hopefuly giving Ctrl-C more chance!
            end % while isempty(values1)
            values0 = [values0; values1];
            times0 = [times0; times1];
        end % while length(values0)<2
        values = [values; values0];
        times = [times; times0];
        ix = find( diff( times ) < 0.010 ); % i.e. < 10ms
    end % while isempty(ix)
    % return values...
    slices = ( values( ix )*256 + values( ix+1 ) ); % convert byte pairs into slice numbers
    times = times( ix+1 ); % pick out one time stamp per byte pair
    times = floor( times * 1000 ); % convert from sec to ms
end % switch port

% and store final values for logging...
scans.slices = slices;
scans.times = times;
scans.number_of_slices = length( scans.times );
cogent.scanner = scans; % copy data into the cogent structure