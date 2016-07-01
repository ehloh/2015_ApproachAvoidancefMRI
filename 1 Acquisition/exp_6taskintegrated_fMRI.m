% Interleaved Conflict & Control blocks (MRI)
clear all; close all hidden; clc;

for o1=1:1 % Documentation 
% Separate par & data files for Conflict & Control task    
%
%     Col 1:      Trial type (read # from grid, starting from bottom-left corner)
%     Col 2:      # Token-pairs offered (Trial type X, 1-6)
%     Col 3:      Set pBomb level  (Trial type Y, in ascending order of probability)
%     Col 4:      Explored half (if applicable)
%     Col 5:      Task type (1=Conflict task, 2=Control task)
%     Col 6:      Bomb present in Activated tokens?
%     Col 7:      Trial number
% [BEHAVIOUR]    
%     Col 8:      [1st Response] 
%                             Conflict task (1=Accept, 2=Reject, 3=Explore) 
%                             Control task (1=No bomb, 2=Bomb, 3=Explore) 
%     Col 9:      [1st Response] RT
%     Col 10:     [2nd Response] 
%                             Conflict task (1=Accept, 2=Reject) 
%                             Control task (1=No bomb, 2=Bomb) 
%     Col 11:     [2nd Response] RT
%     Col 12:     Outcome/Accuracy 
%                                 Conflict task (1=Gain, 0=Null, -1=Loss)
%                                 Control task (1=Correct, 0=Wrong) 
%     Col 13:     Trial valid? (0=Aborted)
%     Col 14:     Outcome presented? (1=Yes)
%     Col 15:     Outcome magnitude
%     Col 16:     
%     Col 17:     Block
%     Col 18:     [Explore info] Did exploration reveal a bomb?
%     Col 19:     [Explore info] Position of revealed bomb (if applicable)
% [TIMINGS] 
%     Col 20:     START TRIAL
%     Col 21:     Fixation onset
%     Col 22:     [1st Offer] Onset ##
%     Col 23:     [Response] Onset
%     Col 24:     [Response] RT (not duration)
%     Col 25:     [Explored-offer] Onset # #
%     Col 26:     [2nd Response] Onset
%     Col 27:     [2nd Response] RT (not duration)
%     Col 28:     Outcome onset
%     Col 29:     END TRIAL
%
end 
for o1=1:1 % Execution settings: Testing/Coding? 
    
    testing=0;
    w.scanning=0;
    if testing==1
        w.subjname=input('Subject ID: ','s');
        w.screenmode=1;
        block=input('Block number:   ');
    else
        w.subjname='t1';
        w.screenmode=0;
        block=1;
    end
    %
    exec.plotfigs=0;
    exec.checkstruc=0; % Check statistical structure of task (derived via probabilistic sampling)
    where=pwd;
    if w.scanning==1 && testing==1
        w.screenmode=2; % IF SCANNING: 2
        w.resolution=3;
    else
        w.screenmode=0;
        w.resolution=3;
    end
end
for o1=1:1 % Paramters for this block
    load(['Data\' w.subjname '_file_4fMRIparameters.mat'])
    rdata=blockpar{block};
    if p.counterbal_startask==1
        rdata=sortrows(rdata,5);
    elseif p.counterbal_startask==2
        rdata=sortrows(rdata,-5);
    else
        input('ERROR: Invalid start-task counterbalancing specified')
    end
    
    p.block=block;
    if w.scanning==1
        switch p.counterbal_hands
            case 1
                p.buttons.key1=28; % Accept/Reject keys
                p.buttons.key2=29;
                p.buttons.key3=30;
                p.buttons.key4=35; % No bomb/bomb keys
                p.buttons.key5=34;
                p.buttons.key6=33;
            case 2
                p.buttons.key1=35; % Accept/Reject keys
                p.buttons.key2=34;
                p.buttons.key3=33;
                p.buttons.key4=28; % No bomb/bomb keys
                p.buttons.key5=29;
                p.buttons.key6=30;
        end
    else
        switch p.counterbal_hands
            case 1
                p.buttons.key1=97; % Accept/Reject keys
                p.buttons.key2=100;
                p.buttons.key3=98;
                p.buttons.key4=26; % No bomb/bomb keys
                p.buttons.key5=24;
                p.buttons.key6=3;
            case 2
                p.buttons.key1=26; % Accept/Reject keys
                p.buttons.key2=24;
                p.buttons.key3=3;
                p.buttons.key4=97; % No bomb/bomb keys
                p.buttons.key5=100;
                p.buttons.key6=98;
        end
    end
    
    %
    diary([where filesep 'Data' filesep w.subjname '_diary_6intergratedfMRI_b' num2str(block)]); diary on
    w.line=50; w.back=p.disp.settokencolor;
    input(['Number of blocks= ' num2str(p.nblocks) '. OK?   ']);
end
for o1=1:1 % Parameters for scanning
scan.NumDummies=4; % N dummy volumes 
scan.SlicesVol=48; % n slices/vol (N slices + oversampling)
scan.DummySlicesWait = scan.NumDummies*scan.SlicesVol;
scan.STWF = scan.DummySlicesWait+1;  % Slice To Wait For = No. 

switch w.scanning
    case 1 % Scanning
        config_serial(1);   
        scan.port=1; 
    case 0 % Simulation
        config_serial(1);   
        scan.port=0;         
        %
        scan.sim.TR=70.0; % in ms
        scan.sim.slices=48;
        scan.sim.volumes=100;
end
end

% COGENT
config_display(w.screenmode,w.resolution, [0 0 0], [1 1 1], 'Helvetica', 40, 9,0, 0); % marker 5g start
config_keyboard;
cd([where filesep 'Data'])
config_log([w.subjname '_cogentlog_6integratedfMRI_b' num2str(block)]); 
start_cogent % if using testing laptops, w.resolution in config dis must be 3! also, w.lineoptions =-330;
cgtext('EXPERIMENTER PRESS SPACEBAR', 0, w.line*2)
t.startcogent=cgflip(0.5,0.5,0.5);
waitkeydown(inf,71);

%% SET UP STOCK DISPLAYS

for o1=1:1 % Create stock displays 
% Sprite details
%         Sprite 11-16: Empty tokens for all pBomb levels (1-6)
%         Sprite 6: [Offer] White circle (Offered token)
%         Sprite 7: [Explore] - Bomb
%         Sprite 8: [Explore] - Not a bomb
%         Sprite 9: [Outcome] Green circle (Won token)
%         Sprite 10: [Outcome] Red circle (Loss token - Bomb)
%         Sprite 17: Empty tokens for tutorial
%         Sprite 18: [Outcome] Yellow circle (No bomb token)
for i=1:p.task.envthreat_nlevels % Sprites 11-16 (Empty tokens for all levels)
    w.backcol=p.disp.contingencycolours{i,1};
    cgmakesprite(10+i,1200,800, w.backcol)
    cgsetsprite(10+i) 
    for j=1:p.task.envthreat_nlevels 
        cgpencol(p.disp.settokencolor)
        cgpenwid(8)
        cgellipse(p.disp.tokenpos(j),p.disp.tokenY_row1,p.disp.tokensize,p.disp.tokensize)
        cgellipse(p.disp.tokenpos(j),p.disp.tokenY_row2,p.disp.tokensize,p.disp.tokensize)
    end
    cgsetsprite(0)
end
for i=1:5 % Create sprites 6-10
cgmakesprite(5+i,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
cgtrncol(5+i,'n') 
end
cgmakesprite(18,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
cgtrncol(18,'n') 
cgsetsprite(6) % Sprite 6: [Offer] Offered token
cgpencol([p.disp.settokencolor])
cgellipse(0,0,p.disp.tokensize,p.disp.tokensize,'f')
cgsetsprite(7) % Sprite 7: [Explore] Bomb token
cgpenwid(10)
cgpencol(p.disp.color_bomb)
cgellipse(0,0,p.disp.tokensize,p.disp.tokensize)
cgsetsprite(8) % Sprite 8: [Explore] Not a bomb
cgpenwid(10)
cgpencol(p.disp.color_nobomb)
cgellipse(0,0,p.disp.tokensize,p.disp.tokensize)
cgsetsprite(9) % Sprite 9: [Outcome] Win token
cgpencol([p.disp.color_nobomb])
cgellipse(0,0,p.disp.tokensize+8,p.disp.tokensize+8,'f')
cgsetsprite(10) % Sprite 10: [Outcome] Loss token
cgpencol([p.disp.color_bomb])
cgellipse(0,0,p.disp.tokensize+8,p.disp.tokensize+8,'f')
cgsetsprite(18) % Sprite 18: [Outcome] Loss token
cgpencol([p.disp.color_nowinnolose])
cgellipse(0,0,p.disp.tokensize+8,p.disp.tokensize+8,'f')
cgmakesprite(4,700,100, [1 0 0]) % Sprite 4: Response options (Conflict task)
cgtrncol(4,'r') 
cgsetsprite(4)
cgrect(-230,0,170, 60, [0.3 0.3 0.3])
cgrect(0,0,170, 60, [0.3 0.3 0.3])
cgrect(230,0,170, 60, [0.3 0.3 0.3])
cgpencol([0 0 0])
cgfont('Helvetica',30)
cgtext('ACCEPT',-230,0)
cgtext('REJECT',0,0)
cgtext('EXPLORE',230,0)
cgfont('Helvetica',40)
cgmakesprite(5,700,100, [1 0 0]) % Sprite 5: Response options (Control task)
cgtrncol(5,'r') 
cgsetsprite(5)
cgrect(-230,0,170, 60, [0.3 0.3 0.3])
cgrect(0,0,170, 60, [0.3 0.3 0.3])
cgrect(230,0,170, 60, [0.3 0.3 0.3])
cgpencol([0 0 0])
cgfont('Helvetica',30)
cgtext('NO BOMB',-230,0)
cgtext('BOMB',0,0)
cgtext('EXPLORE',230,0)
cgfont('Helvetica',40)
cgsetsprite(0)
for o2=1:1 % TEST DISPLAYS --------------------
w.testdisplay=0;
if w.testdisplay==1
    w.width=100; % Test colours
    for i=1:p.task.envthreat_nlevels 
        w.col=p.disp.contingencycolours{i,1};
        cgmakesprite(20+i,w.width,600,w.col)
        cgdrawsprite(20+i,-(w.width*3)-40+i*w.width,0)
    end
    cgflip(p.disp.settokencolor)
    for i=1:p.task.envthreat_nlevels % Test displays
        cgdrawsprite(10+i,0,0)
        cgflip(0,0.1,0)
        waitkeydown(inf)
    end
    cgtext('START TOKENS',0,0)
    cgflip(0,0.1,0)
    waitkeydown(inf) 
    for i=1:p.task.envthreat_nlevels 
    cgdrawsprite(10+i,0,0)
    cgdrawsprite(5+i,p.disp.tokenpos(i),p.disp.tokenY_row1)
    cgflip(0,0,0)
    waitkeydown(inf)
    end 
end
end
end
for o1=1:1 % Buttons 
    cgtext('Please position your fingers on the keys', 0, w.line*5)
    cgtext('NO BOMB/BOMB TASK',0,w.line*3)
    switch p.counterbal_hands
        case 1
            cgtext('(Left hand)',0,w.line*2)
            cgtext('(Right hand)',0,w.line*-3-20)
        case 2
            cgtext('(Right hand)',0,w.line*2)
            cgtext('(Left hand)',0,w.line*-3-20)
    end
    cgdrawsprite(5,0,w.line*0+40)
    cgtext('ACCEPT/REJECT TASK',0,w.line*-2-20)
    cgdrawsprite(4,0,w.line*-4-40)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf);
    
    % Read buttons
    cgtext('Press the button for NO BOMB', 0, w.line*4)
    cgtext('(in the No bomb/Bomb task)',0,w.line*2)
    cgdrawsprite(5,0,w.line*-4)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf,p.buttons.key4);
%     p.buttons.key4=waitkeydown(inf);
    cgtext('Press the button for BOMB', 0, w.line*4)
    cgtext('(in the No bomb/Bomb task)',0,w.line*2)
    cgdrawsprite(5,0,w.line*-4)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf,p.buttons.key5);
%     p.buttons.key5=waitkeydown(inf);
    cgtext('Press the button for EXPLORE', 0, w.line*4)
    cgtext('(in the No bomb/Bomb task)',0,w.line*2)
    cgdrawsprite(5,0,w.line*-4)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf, p.buttons.key6);
    cgtext('Press the button for ACCEPT', 0, w.line*4)
    cgtext('(in the Accept/Reject task)',0,w.line*2)
    cgdrawsprite(4,0,w.line*-4)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf,p.buttons.key1);
    cgtext('Press the button for REJECT', 0, w.line*4)
    cgtext('(in the Accept/Reject task)',0,w.line*2)
    cgdrawsprite(4,0,w.line*-4)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf,p.buttons.key2);
    cgtext('Press the button for EXPLORE', 0, w.line*4)
    cgtext('(in the Accept/Reject task)',0,w.line*2)
    cgdrawsprite(4,0,w.line*-4)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf,p.buttons.key3);
    
    % Done
    cgtext('Experimenter, set up SPIKE and', 0, w.line*3)
    cgtext('press spacebar when ready', 0, w.line*2)
    cgflip(0.5,0.5,0.5);
    waitkeydown(inf,71);
end

%% ########### SCANNER STARTUP ############################################

if rdata(1,5)==1
   wt.taskname='ACCEPT/REJECT';
elseif rdata(1,5)==2
   wt.taskname='NO BOMB/BOMB';
end
cd(where)
cgtext(['First you do a block of the ' wt.taskname ' trials'],0,w.line*2)
cgtext('Press any key to start',0,w.line*-3)
times.starttask=cgflip(w.back);
waitkeydown(inf)
clock
cgtext('+',0,0)           
logstring('##########################################################################')
for o1=1:1;
    if w.scanning==1 % Scanner setup
        cd(where)
        cgtext('GET READY TO START', 0, w.line*4)
        cgtext('+', 0, w.line*0)
        times.present_startscreen=cgflip(w.back);
        
        % -----> SCANNER STARTS (outside this script)
        
        % Wait for pulses & trigger startpoint
        if w.scanning==1 % Real scanning
            [scan.StartSlice1, scan.SliceTime1] = waitslice(scan.port,scan.STWF);
            times.start=scan.SliceTime1; % Crucial start point
            logstring('Scanner status: Trio scanning')
            logstring(['START POINT - ' num2str(times.start)])
            disp(['Start point: ' num2str(times.start)])
        elseif w.scanning==0 % Simulator
            fs_config_emulscan(scan.sim.TR,scan.sim.slices,scan.sim.volumes) % sonata (i.e. 90ms per slice), 16 slices (per volume), 2 volumes
            fs_start_emulscan % Start emulator
            [scan.StartSlice1, scan.SliceTime1] = waitslice(scan.port,scan.STWF);
            times.start_sim=scan.SliceTime1; % Crucial start point
            logstring('Scanner status: Simulation')
            logstring(['START POINT - ' num2str(times.start_sim)])
        elseif w.scanning==-1 % Manual trigger
            wait(2000)
            cgtext('Press any key to start trials', 0, w.line*2)
            cgflip(0,0,0);
            times.start=waitkeydown(inf);
            logstring('Scanner status: Manual (no scanning or simulation)')
            logstring(['START POINT(Manually triggered) - ' num2str(times.start)])
        else
            input('ERROR: w.scanning is not set!')
        end
    end; 
end

%%  EXECUTION LOOP
w.alltrials=0;
w.omissions=0;
w.explorations=0;
w.block=1;

for trialnum=1: size(rdata,1)
    rdata(trialnum,7)=trialnum;

        for o1=1:1 % Set up displays & details for this trial
            wt=[];   
            wt.pos=[p.disp.tokenpos  zeros(p.task.envthreat_nlevels,1)];  % Assemble (shuffled) positions
            wt.startpos=randi(p.task.envthreat_nlevels+1-rdata(trialnum,2),1);
            for i=wt.startpos:(wt.startpos+rdata(trialnum,2)-1)
                wt.pos(i,2)=1;
            end
            wt.posX(1:rdata(trialnum,2),1)=wt.pos((wt.pos(:,2)==1),1);
            if rdata(trialnum,5)==1
                wt.keys=[p.buttons.key1;p.buttons.key2;p.buttons.key3];
            else
                wt.keys=[p.buttons.key4;p.buttons.key5;p.buttons.key6];
            end
            % Resets
            wk=[];
            if rdata(trialnum,5)==1
                wt.task='Accept, Reject, or Explore?';
            else
                wt.task='Is there an activated bomb?';
            end
            rdata(trialnum,13)=999;
        end
        cgtext(wt.task,0, 220)
        cgdrawsprite(3+rdata(trialnum,5),0,-250)
        cgtext('+',0,0)           
        t.fixate=cgflip(w.back);
        rdata(trialnum,17)=w.block;
        rdata(trialnum,20)=t.fixate;
        rdata(trialnum,21)=t.fixate;
        rdata(trialnum,30)=nan;
        
       % Verify parameters-display match-up
       logstring(['############# Trial #: ' num2str(trialnum) ' ########################################### '])
       logstring(['Task: ' num2str(rdata(trialnum,5)) ' (1=Conflict, 2=Control)'])
       logstring(['Offer: Colour: ' p.disp.contingencycolours{rdata(trialnum,3),4} ' (' num2str(rdata(trialnum,3)) '), ' num2str(2*rdata(trialnum,2)) ' activated'])
       logstring(['Activated bomb: ' num2str(rdata(trialnum,6))])
       logstring(['Outcome shown? '  num2str(rdata(trialnum,14))])
       %             wtt.goon=waitkeydown(inf,60);
        
        % [OFFER] -----------
        cgdrawsprite(rdata(trialnum,3)+10,0,0)
        cgtext('+',0,0)           
        cgtext(wt.task,0, 220)
        for j=1:rdata(trialnum,2)
            cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
        end
        cgdrawsprite(3+rdata(trialnum,5),0,-250)
        waituntil(t.fixate*1000+800+rand*500)
        clearkeys
        wt.offer1=cgflip(w.back);
        [wk.key1press wk.key1time wk.a]=waitkeydown(p.disp.firstoffer,[wt.keys]);
        if wk.a>0
            rdata(trialnum,9)= wk.key1time(1)-wt.offer1*1000;
            switch wk.key1press(1)
                case wt.keys(1)
                    rdata(trialnum,8)=1;
                case wt.keys(2)
                    rdata(trialnum,8)=2;
                case wt.keys(3)
                    rdata(trialnum,8)=3;
            end
            rdata(trialnum,13)=1;
            rdata(trialnum,23)=wk.key1time(1);
            rdata(trialnum,24)=rdata(trialnum,9);
            wt.responseok=1;
        else 
            rdata(trialnum,8)=0;
            rdata(trialnum,13)=0;
            rdata(trialnum,23)=nan;
            rdata(trialnum,24)=0;
            wt.responseok=0;
        end
        rdata(trialnum,22)=wt.offer1;
        
        % [If Explored + Outcome shown] 2ND CHOICE --------------
        if rdata(trialnum,14)==1 && rdata(trialnum,8)==3
            logstring('EXPLORATION PHASE')
            wt.responseok=0;
            cgdrawsprite(rdata(trialnum,3)+10,0,0)
            cgtext('+',0,0)
            cgtext(wt.task,0,220)
            for j=1:rdata(trialnum,2) 
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
            end
            cgdrawsprite(3+rdata(trialnum,5),0,-250)
            cgrect(230,-250,120, 30, [0.1 0.1 0.1])
            eval(['wt.revealed_Y=p.disp.tokenY_row' num2str(rdata(trialnum, 4)) ';'])  % Reveal status of half
            for j=1:rdata(trialnum,2) 
                cgdrawsprite(8,wt.posX(j),wt.revealed_Y);
            end
            wt.x=randi(2);
            if rdata(trialnum, 6)==1
                wt.posxB=wt.posX(randi(size(wt.posX,1)));
                if wt.x<2
                    cgdrawsprite(7,wt.posxB,wt.revealed_Y);
                    rdata(trialnum,18)=1;
                    rdata(trialnum,19)=wt.posxB;
                    logstring('Explore Info: Bomb activated, and shown')
                else
                    rdata(trialnum,18)=0;
                    rdata(trialnum,19)=999;
                    logstring('Explore Info: Bomb activated, not shown')
                end
            else
                logstring('Explore Info: No bomb, no bomb shown')
            end
            cgtext(['Information cost: - ' num2str(p.task.explorationcost*10) ' p'],0,130)
            waituntil(wt.offer1*1000+p.disp.firstoffer)
            clearkeys
            wt.offer2=cgflip(w.back);
            [wk.key2press wk.key2time wk.a2]=waitkeydown(p.disp.secondoffer,[wt.keys(1) wt.keys(2)]);       
            if wk.a2>0
                rdata(trialnum,11)= wk.key2time(1)-wt.offer2*1000;
                switch wk.key2press(1)
                    case wt.keys(1)
                        rdata(trialnum,10)=1;
                    case wt.keys(2)
                        rdata(trialnum,10)=2;
                end
                rdata(trialnum,13)=1;
                rdata(trialnum,26)=wk.key2time(1);
                rdata(trialnum,27)=rdata(trialnum,11);
                wt.responseok=1;
            else 
                rdata(trialnum,10)=0;
                rdata(trialnum,13)=0;
                rdata(trialnum,26)=nan;
                rdata(trialnum,27)=nan;
                wt.responseok=0;
            end
            rdata(trialnum,25)=wt.offer2;
            
            % Workspace for continuation
            wt.respcol=10;
            wt.lastofferonset=wt.offer2;
        else
            % Workspace for continuation
            wt.respcol=8;
            wt.lastofferonset=wt.offer1;
            rdata(trialnum,25)=nan;
            rdata(trialnum,26)=nan;
            rdata(trialnum,27)=nan;
        end
        
        
        % [EVALUATE OUTCOME] -----------
        if  wt.responseok==1 % Mark outcome magnitude
             if rdata(trialnum, 8)==3 % Explored
                    wt.explcost=p.task.explorationcost;
                    wt.infocost=['Information cost: - ' num2str(10*p.task.explorationcost) ' p'];
                else
                    wt.explcost=0;
                    wt.infocost=' ';
             end
            if rdata(trialnum,5)==1 % CONFLICT TASK
                if rdata(trialnum,wt.respcol)==2 % reject
                    rdata(trialnum,12)=0;
                    rdata(trialnum,15)=0-wt.explcost;
                    wt.won= '?';
                    wt.tokenspritenum=6;
                elseif rdata(trialnum,6)==0 && rdata(trialnum,wt.respcol)==1 % accept + no bomb
                    rdata(trialnum,12)=1;
                    rdata(trialnum,15)=rdata(trialnum, 2)*2-wt.explcost;
                    wt.tokenspritenum=9;
                    wt.won=['+ ' num2str(10*rdata(trialnum, 15)) ' p'];
                    if rdata(trialnum, 8)==3
                        wt.infocost=['(Information cost: ' num2str(10*p.task.explorationcost) ' p)'];
                    end
                elseif rdata(trialnum,6)==1 && rdata(trialnum,wt.respcol)==1 % accept + bomb
                    rdata(trialnum,12)=-1;
                    rdata(trialnum,15)=p.task.fixedloss-wt.explcost;
                    wt.tokenspritenum=10;
                    wt.won=[num2str(rdata(trialnum,15)*10) ' p'];
                    if rdata(trialnum, 8)==3
                        wt.infocost=['(Information cost: ' num2str(10*p.task.explorationcost) ' p)'];
                    end
                else
                    rdata(trialnum,15)=0-wt.explcost;
                end
            elseif rdata(trialnum,5)==2 % CONTROL TASK
                if rdata(trialnum,6)+1==rdata(trialnum,wt.respcol)
                    rdata(trialnum,12)=1;
                    rdata(trialnum,15)=rdata(trialnum, 2)*2-wt.explcost;
                    wt.tokenspritenum=9;
                    wt.won=['+ ' num2str(rdata(trialnum,15)*10) ' p'];
                else
                    rdata(trialnum,12)=0;
                    rdata(trialnum,15)=0-wt.explcost;
                    wt.tokenspritenum=18;
                    wt.won= '0 p';
                end
            end
        else
            rdata(trialnum,28)=nan;
            rdata(trialnum,12)=0;
            rdata(trialnum,15)=0;
            wt.won= 'No response';
            wt.tokenspritenum=6;
            wt.infocost=' ';
        end

        % [DISPLAY OUTCOME IF NECESSARY] ---------------
        if rdata(trialnum,14)==1
            cgdrawsprite(rdata(trialnum,3)+10,0,0)
            cgtext('+',0,0)
            cgdrawsprite(3+rdata(trialnum,5),0,-250)
            for j=1:rdata(trialnum,2)
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
            end
            for j=1:rdata(trialnum,2) % win/loss tokens 
                cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row2);
            end
            cgfont('Helvetica',70)
%             cgtext(wt.won,0,160)
            cgtext(wt.won,0, 220)
            cgfont('Helvetica',25)
%             cgtext(wt.infocost,0,-130)
            cgtext(wt.infocost,0,130)
            cgfont('Helvetica',40)
            waituntil(wt.lastofferonset*1000+p.disp.firstoffer)
            wt.outcome=cgflip(w.back);
            % Workspace for continuation
           wt.lastoftrialonset=wt.outcome;
           wt.timebeforenexttrial=p.disp.firstoffer/2;
            rdata(trialnum,28)=wt.outcome;
        else
            % Workspace for continuation
            wt.lastoftrialonset=wt.lastofferonset;
           wt.timebeforenexttrial=p.disp.firstoffer;
            rdata(trialnum,28)=nan;
        end
        
       % Verify choices recorded
       logstring(' ---------------------------------------------------')
       logstring(['Fixation onset: ' num2str(t.fixate)])
       logstring(['1st Offer:  ' num2str(rdata(trialnum,22))])
       logstring(['2nd Offer:  ' num2str(rdata(trialnum,25))])
       logstring(['Outcome: ' num2str(rdata(trialnum,15)) '        Onset '  num2str(rdata(trialnum,28))])
       logstring(' ==============================')
       logstring(['Response valid: ' num2str(rdata(trialnum,13))])
       logstring(['Response 1: ' num2str(rdata(trialnum,8)) '  (1=Accept/No bomb, 2=Reject/Bomb, 3=Explore)'])
       logstring(['                ' num2str(rdata(trialnum,23)) '    '   num2str(rdata(trialnum,24))])
       logstring(['Response 2: ' num2str(rdata(trialnum,10)) ' (1=Accept/No bomb, 2=Reject/Bomb)'])
       logstring(['                ' num2str(rdata(trialnum,26)) '    '   num2str(rdata(trialnum,27))])
       logstring('Trial summary:')
       logstring([rdata(trialnum,1) rdata(trialnum,2) rdata(trialnum,3) rdata(trialnum,4) rdata(trialnum,5) rdata(trialnum,6) rdata(trialnum,7) rdata(trialnum,8) rdata(trialnum,9)  rdata(trialnum,10) rdata(trialnum,11) rdata(trialnum,12) rdata(trialnum,13) rdata(trialnum,14)  rdata(trialnum,15) rdata(trialnum,16) rdata(trialnum,17) rdata(trialnum,18) rdata(trialnum,19) ]) 
       logstring([rdata(trialnum,20) rdata(trialnum,21) rdata(trialnum,22)  rdata(trialnum,23) rdata(trialnum,24) rdata(trialnum,25) rdata(trialnum,26) rdata(trialnum,27) rdata(trialnum,28) rdata(trialnum,29)])
       disp(num2str([rdata(trialnum,1) rdata(trialnum,2) rdata(trialnum,3) rdata(trialnum,4) rdata(trialnum,5) rdata(trialnum,6) rdata(trialnum,7) rdata(trialnum,8) rdata(trialnum,9) rdata(trialnum,10) rdata(trialnum,11) rdata(trialnum,12) rdata(trialnum,13) rdata(trialnum,14)  rdata(trialnum,15) rdata(trialnum,16) rdata(trialnum,17) rdata(trialnum,18) rdata(trialnum,19) rdata(trialnum,20) rdata(trialnum,21) rdata(trialnum,22) rdata(trialnum,23) rdata(trialnum,24) rdata(trialnum,25) rdata(trialnum,26) rdata(trialnum,27) rdata(trialnum,28) rdata(trialnum,29)]));
%           wtt.goon=waitkeydown(inf,60);        
        
       % [BREAKS] --------------------------------
       if trialnum~=size(rdata,1) && rdata(trialnum,5)~=rdata(trialnum+1,5) % Task transition
           logstring('########################################')
           w.thisblock=rdata(rdata(:,17)==w.block,:); 
           w.score=sum(w.thisblock(w.thisblock(:,13)==1, 15))*10;
           disp(['Score for block: ' num2str(w.score)]);
           clock
           %
           wait(p.disp.firstoffer)
           clock
           if rdata(trialnum+1,5)==1  
               wt.taskname='ACCEPT/REJECT';
           else
               wt.taskname='NO BOMB/BOMB';
           end
           times.breakstart(w.block)=cgflip(w.back);
           logstring(['Break start: ' num2str(times.breakstart(w.block))])
           cgtext(['You have just completed ' num2str(w.block) '/' num2str(round(12/p.nblocks)) ' of this block'],0,w.line*4)
           cgtext('Now you do a some trials of the ',0,w.line*2)
           cgtext([wt.taskname ' task'],0,w.line*1)
           cgtext(['Score for section just completed: ' num2str(w.score) ' p'], 0, w.line*-2)
           cgtext('The task will continue in 30 seconds',0,w.line*-4)
           times.breakend(w.block)=cgflip(w.back);
           wait(30*1000)
           cgtext('+',0,0)
           logstring(['Break end: ' num2str(times.breakend(w.block))])
           w.block=w.block+1;
           logstring('########################################')
           
       end

    % Save workspace 
    save(['Data' filesep w.subjname '_workspace_6integratedfMRI_b' num2str(block)]);
    
    cgtext('+',0,0)
	rdata(trialnum,29)=wt.lastoftrialonset*1000+ wt.timebeforenexttrial;
    waituntil(wt.lastoftrialonset*1000+ wt.timebeforenexttrial)   
       
       
end

clock

%%  FINISH UP 

try
    [times.aftertrials_slicenum times.aftertrials_slicetime]=getslice(scan.port);
end
try
    w.time=clock;
    f_sendemail('kurzlich',['[SCANNING] (' num2str(w.time(4)) ':', num2str(w.time(5)) 'hrs) ' w.subjname ' has finished scanning (block ' num2str(block) ')'], ['Finished block ' num2str(block)])
catch
end

cgtext('+',0,0)           
times.end_lastfixationcross=cgflip(w.back);
disp('##########################################################################')
disp('END OF BLOCK')
disp('WAIT 30 seconds')
wait(30*1000)
disp('DONE - Experimenter press spacebar')
waitkeydown(inf,71);

stop_cogent

% Post-hoc calcs
outcomes.n_omissions=w.omissions;
outcomes.n_explorations=w.explorations;
outcomes.net=sum(rdata((rdata(:,13)==1),15))-w.explorations*p.task.explorationcost;

% Save
w.dat.rdata=rdata;
w.dat.data=rdata(rdata(:,13)==1,:);
w.dat.conflictdata=w.dat.data(w.dat.data(:,5)==1,:);
w.dat.controldata=w.dat.data(w.dat.data(:,5)==2,:);
w.dat.par=par;
p.log.subject=w.subjname; p.log.date=date; p.log.clock=clock;
w.dat.settings=p;
w.dat.outcomes=outcomes;
w.dat.times=times;

try
    w.dat.abortedtrials=abortedtrials;
catch
    w.dat.abortedtrials='No aborted trials';
end
integratedfMRI=w.dat;
filename=[w.subjname '_file_6integratedfMRI_b' num2str(block) '.mat'];
save(filename, 'integratedfMRI')

try 
    movefile(filename, [where filesep 'Data'])
catch
    disp('EXPERIMENTER MESSAGE: Data saved in curerent directory, not in rdata folder')
end
try % Transfer file to DeletedDaily
movetowhere_drive='\\Asia\DeletedDaily';
movetowhere_folder='EL Goal-Conflict Data';
cd(movetowhere_drive)
w.dirs=dir;
w.dir_there=0;
for i=3:length(w.dirs) % Data folder in Drive already there?
    if strcmp(w.dirs(i,1).name,movetowhere_folder)==1
        w.dir_there=1;
    end
end
if w.dir_there==0
    mkdir(movetowhere_folder)
else
end
cd(movetowhere_folder)
copyfile([where filesep 'Data' filesep filename])
copyfile([where filesep 'Data' filesep w.subjname '_workspace_6integratedfMRI_b' num2str(block)])
copyfile([where filesep 'Data' filesep w.subjname '_diary_6intergratedfMRI_b' num2str(block)])
copyfile([where filesep 'Data' filesep w.subjname '_cogentlog_6integratedfMRI_b' num2str(block)])
w.warning=' ';
catch
    w.warning='WARNING: data not transfered to Deleted Daily. Transfer files manually';
end

%
disp('#################################################')
try
    disp(['Duration: ' num2str((t.end-t.startcogent)/60) ' minutes'])
catch
end

disp(['Score: ' num2str(sum(rdata(:,15)))])
disp(['# of Omissions: ' num2str(size(rdata,1)-sum(rdata(:,13)))])
disp('(See variable ''w.dat.behaviour'' for behavioural performance)')
disp(w.warning)
disp('#################################################')

cd(where)
[w.dat.beh_conflict] = f_plotindividualbehaviour(p,w.dat.conflictdata, [w.subjname ' (Conflict)']); % Plot subject's behaviour
[w.dat.beh_control] = f_plotindividualbehaviour(p,w.dat.controldata, [w.subjname ' (Control)']); 
diary off

