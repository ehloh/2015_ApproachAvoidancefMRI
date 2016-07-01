% Emulscan Demo - An example of how to use Emulscan
%
% config_emulscan( 'sonata',16,2 ) % sonata (i.e. 90ms per slice), 16 slices (per volume), 2 volumes 
% start_emulscan
% slicelist = 5:5:30;
% for n = 1:length(slicelist)
%     disp('waiting...')
%     [s,t] = waitslice( 0,slicelist(n) );
%     disp( ['   ...found slice ' num2str( s(end) ) ' @ ' num2str(t(end))] )
%     drawnow
% end
%
% Version 2.0 09-10-2009

% Version History
% 2.0, 09-10-2009, E.F. - Corrected bug in recording of timestamps. Added more comments.

% Configure emulator
config_emulscan( 'allegra',16,2 ) % sonata (i.e. 90ms per slice), 16 slices (per volume), 2 volumes 
% Start emulator
start_emulscan

% Wait for slices 5 10 15 20 25 & 30
slicelist = 5:5:30; clear t
for n = 1:length(slicelist)
    disp('waiting...')
    % Wait for a slice.
    % NB look on port ZERO !
    [s,t] = waitslice( 0,slicelist(n) );
    disp( ['   ...found slice ' num2str( s(end) ) ' @ ' num2str(t(end))] )
    drawnow
end