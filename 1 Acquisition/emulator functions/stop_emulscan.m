function stop_emulscan
% scanner emulator StopFcn - Emulscan helper function
% 
% Version 2.0 09-10-2009

% Version History
% 2.0, 09-10-2009, E.F. - Version numbers made consistant
% 1.0, 22-07-2008, E.F.

global emulscan

t = emulscan.timer;
disp( sprintf( 'Scanner emulator stopped : Average period %dms, Total pulses %d.', 1000*t.AveragePeriod, t.TasksExecuted ) )

if abs( (t.AveragePeriod-t.Period)/t.Period ) > 0.1
    warning('Emulscan timing error more than 10%!')
end
if t.TasksExecuted ~= t.TasksToExecute
    warning('Emulscan interupted or incomplete')
end