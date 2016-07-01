function fs_config_emulscan( scanner, slices, volumes )
% Configure Emulscan - Emulscan function
%
% Syntax:  emulscan = config_emulscan( scanner, slices, volumes );
%           scanner - 'sonata', 'allegra', 'trio' or a custom TR in ms
%            slices - number of slices per volume
%           volumes - number of volumes
%
% e.g.
% config_emulscan( 'sonata', 32, 100 )
%    will start and set-up emulscan to emulate the sonata running 100 volumes of 32 slices each
% config_emulscan( 80, 16, 50 )
%    will start and set-up emulscan to emulate a scanner running 50 volumes
%    of 16 slices each and 80ms per slice
% 
% Version 2.0 09-10-2009

% Version History
% 2.0, 09-10-2009, E.F. - Version numbers made consistant
% 1.1, 22-07-2008, E.F. Changed from using seperate Matlab session via pnet.dll to a Matlab timer
% 1.0, 09-12-2005, E.F.

global emulscan

switch class( scanner )
    case { 'char' }
        switch lower( scanner )
            case { 'trio' }
                TR = 0.068; % 68ms
            case { 'allegra' }
                TR = 0.065; % 65ms
            case { 'sonata' }
                TR = 0.090; % 90ms
            otherwise
                error('"scanner" must be "trio", "allegra", "sonata" or a numeric TR')
        end % switch tolower( scanner )
    case { 'double' }
        TR = scanner/1000; % convert from s to ms
    otherwise
        error('"scanner" must be "allegra" or "sonata" or a numeric TR')
end % switch class( scanner )

nSlices = slices*volumes;
emulscan.timer = timer( ...
    'BusyMode',      'queue', ...
    'ExecutionMode', 'fixedRate', ...
    'Name',          'emulscan', ...
    'Period',        TR, ...
    'StartDelay',    0, ...
    'StartFcn',      'emulscan_start', ...
    'StopFcn',       'stop_emulscan', ...
    'TasksToExecute',nSlices, ...
    'TimerFcn',     'emulscan_next_slice' ...
    );