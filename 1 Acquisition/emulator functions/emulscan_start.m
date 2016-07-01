function emulscan_start
% Scanner emulator StartFcn - Emulscan helper function
% 
% Version 2.0 09-10-2009

% Version History
% 2.0, 09-10-2009, E.F. - Version numbers made consistant
% 1.0, 22-07-2008, E.F.

global emulscan

emulscan.slice = 1;
emulscan.time = time;
disp('Scanner emulator started...')
drawnow