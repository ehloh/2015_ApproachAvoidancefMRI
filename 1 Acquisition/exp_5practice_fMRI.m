% Interleaved Conflict & Control blocks (MRI)
clear all;close all hidden;clc;

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
%     Col 13:     Trial aborted? (1=Aborted)
%     Col 14:     Outcome presented? (1=Yes)
%     Col 15:     Outcome magnitude
%     Col 16:     Position of revealed bomb (if applicable)
%
end
for o1=1:1 % Execution settings: Testing/Coding? 
    
    testing=1;
    exec.plotfigs=0;
    exec.checkstruc=0; % Check statistical structure of task (derived via probabilistic sampling)
    choicetask=1;
    %
    if testing==1
        w.subjname=input('Subject ID: ','s');
        w.screenmode=0;
    else
        w.subjname='test3';
        w.taskrep=1;
%         w.taskrep=input('# times this stage repeated (1=1st time, 2=2nd, etc): ');
        w.screenmode=0;
    end
    w.taskstage='taskconflictcontrol';
    where=pwd;
    w.res=2;
end
for o1=1:1 % Parameters for practice
     load(['Data\' w.subjname '_file_4fMRIparameters.mat'])
     p.npracticetrials=20;
     %
     par=sortrows(par,7);
     rdata=par(1:p.npracticetrials,:);
     rdata(1:round(p.npracticetrials/2),5)=1;
     rdata(round(p.npracticetrials/2)+1:p.npracticetrials,5)=2;
     rdata(:,7)=rand( p.npracticetrials,1); rdata=sortrows(rdata,7);
     rdata(:,14)=0;
     rdata(1:round(p.npracticetrials/3),14)=1;
     rdata(:,7)=rand( p.npracticetrials,1); rdata=sortrows(rdata,7);
     rdata=sortrows(rdata,5);
    %
    w.line=50;
    rand('state',sum(100*clock));
    w.back=p.disp.settokencolor;
end
for o1=1:1 % Checks for running 
disp('------------------------------')
disp('[SESSION 4: INSTRUCTIONS & PRACTICE FOR INTEGRATED FMRI TASK')
% disp(['No. levels pBombEnv: ' num2str(p.task.envthreat_nlevels)])
% disp(['No. reps per trial type: ' num2str(p.n.reps_pertrialtype)])
disp(['First block task type: ' num2str(rdata(1,5))])
disp(['Exploration cost: ' num2str(p.task.explorationcost*10) ' p'])
disp('------------------------------')
% input('Setup OK?  ');
end

% COGENT
config_display(w.screenmode,w.res, [0 0 0], [1 1 1], 'Helvetica', 40, 9,0, 0); % marker 5g start
config_keyboard;
start_cogent % if using testing laptops, w.res in config dis must be 3! also, w.lineoptions =-330;
cgtext('Now you will do a practice task', 0, w.line*4)
cgtext('to prepare you for the fMRI session', 0, w.line*3)
cgtext('The computer will now give you instructions', 0, w.line*0)
cgtext('Press any key to continue with the instructions', 0, w.line*-1)
t.startcogent=cgflip(0.8,0.8,0.8);
waitkeydown(inf);

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
    w.backcol=p.disp.contingencycolours{i};
    cgmakesprite(10+i,1000,700, w.backcol)
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
        w.col=p.disp.contingencycolours{i};
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

%% INSTRUCTIONS FOR INTEGRATED TASK

for o1=1:1
    if testing==1
        cgpencol([0 0 0])
        cgtext('Please read the instructions carefully, as this',0,w.line*3)
        cgtext('stage of the experiment is very different from', 0, w.line*2)
        cgtext('the previous ones',0,w.line*1)
        cgflip(w.back)
        waitkeydown(inf)

        cgtext('In the previous session,  you won & lost money',0,w.line*5+25)
        cgtext('playing TWO different decision-making games', 0, w.line*4+25)
        cgtext('In both tasks, you saw ''offers'' on each trial,',0,w.line*3+10)
        cgtext('consisting of a background colour & a number of',0,w.line*2+10)
        cgtext(['activated tokens (out of ' num2str(p.task.envthreat_nlevels*2) '). Each offer tells you how'],0,w.line*1+10)
        cgtext('likely there is to be an ''activated bomb'', on that trial',0,w.line*0+10)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgtext('In the previous session,  you won & lost money',0,w.line*5+25)
        cgtext('playing TWO different decision-making games', 0, w.line*4+25)
        cgtext('In both tasks, you saw ''offers'' on each trial,',0,w.line*3+10)
        cgtext('consisting of a background colour & a number of',0,w.line*2+10)
        cgtext(['activated tokens (out of ' num2str(p.task.envthreat_nlevels*2) '). Each offer tells you how'],0,w.line*1+10)
        cgtext('likely there is to be an ''activated bomb'', on that trial',0,w.line*0+10)
        cgtext('In the ACCEPT/REJECT task, you won money',0,w.line*-1-25)
        cgtext('by Accepting offers without an activated bomb',0,w.line*-2-25)
        cgtext('In the BOMB/NO BOMB task, you won money by',0,w.line*-4)
        cgtext('predicting if there was an activated bomb or not',0,w.line*-5)
        cgflip(w.back)
        waitkeydown(inf)

        % Create displays for instructions
        w.egback=[0.3 0.1 0.0];
        if p.task.envthreat_nlevels==6
            w.egpos=[175 105 35 -35 -105 -175];
        elseif p.task.envthreat_nlevels==4
            w.egpos=[105 35 -35 -105];
        else
            input('Error: no positions for tutorial!')
        end
        cgmakesprite(17,500,300, w.egback)
        cgsetsprite(17) 
        for j=1:p.task.envthreat_nlevels
            cgpencol(w.back)
            cgpenwid(6)
            cgellipse(w.egpos(j),40,60,60)
            cgellipse(w.egpos(j),-40,60,60)
            cgpencol(0,0,0)
        end
        cgpencol(w.back)
        cgellipse(w.egpos(3),40,60,60,'f')
        cgellipse(w.egpos(3),-40,60,60,'f')
        cgellipse(w.egpos(4),40,60,60,'f')
        cgellipse(w.egpos(4),-40,60,60,'f')
        cgpencol(0,0,0)
        cgsetsprite(0)

        % RECAP: PROBABILISTIC STRUCTURE
        cgdrawsprite(17,0,140)
        cgtext('Remember: the likelihood of a there being',0,w.line*-1)
        cgtext('an ''activated'' bomb or not on each trial is',0,w.line*-2)
        cgtext('signalled by the combination of the background',0,w.line*-3)
        cgtext('colour & the number of activated tokens',0,w.line*-4)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgdrawsprite(17,0,140)
        cgfont('Helvetica',35)
        cgtext('The background colour signals how likely there is to be a',0,w.line*-1+15)
        cgtext('bomb (activated or not) in the offer ',0,w.line*-2+17)
        cgtext(['- i.e. under any one of the ' num2str(p.task.envthreat_nlevels*2) ' tokens'],0,w.line*-3+20)
        cgfont('Helvetica',40)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgdrawsprite(17,0,140)
        cgfont('Helvetica',35)
        cgtext('The background colour signals how likely there is to be a',0,w.line*-1+15)
        cgtext('bomb (activated or not) in the offer ',0,w.line*-2+17)
        cgtext(['- i.e. under any one of the ' num2str(p.task.envthreat_nlevels*2) ' tokens'],0,w.line*-3+20)
        cgtext('If there IS a bomb in the offer, it is only considered',0,w.line*-4+15)
        cgtext('''activated'' if it falls under one of the activated tokens',0,w.line*-5+17)
        cgtext('(i.e. one of the tokens coloured in white)',0,w.line*-6+20)
        cgfont('Helvetica',40)
        cgflip(w.back)
        waitkeydown(inf)

        cgdrawsprite(17,0,140)
        cgfont('Helvetica',35)
        cgtext('Remember -  if the bomb falls under one of the',0,w.line*-1+15)
        cgtext(' INACTIVATED tokens, it will NOT cause you to lose money', 0, w.line*-2+15)
        cgtext('in the ACCEPT/REJECT task',0,w.line*-3+15)
        cgtext('And in the BOMB/NO BOMB task, you will either win money',0,w.line*-4)
        cgtext('or win nothing - you WON''T lose money for being wrong',0,w.line*-5)
        cgfont('Helvetica',40)
        cgflip(w.back)
        waitkeydown(inf)

        cgdrawsprite(17,0,140)
        cgtext('For both tasks, the money you can win on each trial',0,w.line*-1)
        cgtext('is proportionate to the number of activated tokens',0,w.line*-2)
        cgtext('So an offer like the one shown here means that you',0,w.line*-4)
        cgtext('you could potentially win 40 p',0,w.line*-5)
        cgflip(w.back)
        waitkeydown(inf)

        cgdrawsprite(17,0,140)
        cgpencol(0.8, 0, 0) % Explore
        cgellipse(w.egpos(4),40+140,60,60)
        cgpencol(0, 0.7, 0)
        cgellipse(w.egpos(3),40+140,60,60)
        cgpencol(0,0,0)
        cgfont('Helvetica',35)
        cgtext(['Information cost: - ' num2str(p.task.explorationcost*10) ' p'],0,w.line*1-10)
        cgtext('In BOTH tasks, you will be able to EXPLORE',0,w.line*-1+10)
        cgtext('the offer before deciding how you want to respond', 0, w.line*-2+10)
        cgtext(['Each time you ''Explore'', it will cost you ' num2str(p.task.explorationcost*10) ' p,'],0,w.line*-3-15)
        cgtext('- but it can help you decide what response to make,',0,w.line*-4-15)
        cgtext('to help you win more money overall',0,w.line*-5-15)
        cgfont('Helvetica',40)
        cgflip(w.back)
        waitkeydown(inf)

        % BOTH TASKS TOGETHER
        if rdata(1,5)==1
           w.task1='ACCEPT/REJECT';
           w.task2='BOMB/NO BOMB';
        else
           w.task1='BOMB/NO BOMB';
           w.task2='ACCEPT/REJECT';
        end
        cgtext('This stage of the experiment will have',0,w.line*5)
        cgtext('trials from BOTH of these tasks', 0, w.line*4)
        cgtext('You will start with a block of trials from the', 0, w.line*2)
        cgtext([w.task1 ' task, followed by a block of trials'], 0, w.line*1)
        cgtext(['from the ' w.task2 ' task'], 0, w.line*0)
        cgtext(['You will do several blocks in the scanner,'], 0, w.line*-2)
        cgtext('alternating between the two different tasks', 0, w.line*-3)
        cgflip(w.back)
        waitkeydown(inf)
        
        % Sign-posted
        cgdrawsprite(17,0,140)
        cgtext('Accept, Reject, Explore?',0,w.line*5)
        cgtext('You will need to concentrate hard to not get',0,w.line*-2)
        cgtext('confused about what type of task you are doing,', 0, w.line*-3)
        cgtext('at each point in time!', 0, w.line*-4)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgdrawsprite(17,0,140)
        cgtext('Accept, Reject, Explore?',0,w.line*5)
        cgfont('Helvetica',35)
        cgtext('To help you keep track of which task you are currently doing,',0,w.line*-2)
        cgtext('the task question will be onscreen, for each offer',0,w.line*-3)
        cgtext('- just like in the display above',0,w.line*-4)
        cgfont('Helvetica',40)
        cgflip(w.back)        
        waitkeydown(inf)
       
        % Hands -------------------------
        switch p.counterbal_hands
            case 1
                cgdrawsprite(17,0,140)
                cgtext('Accept, Reject, Explore?',0,w.line*5)
                cgfont('Helvetica',35)
                cgtext('=Left arrow          =Down arrow       =Right arrow',0,w.line*-2)
                cgdrawsprite(4,0,w.line*-1)
                cgtext('Now you will do a quick practice of the task',0,w.line*-3-15)
                cgtext('For this practice session, please use your',0,w.line*-4-15)
                cgtext('Use your RIGHT hand for the ACCEPT/REJECT task',0,w.line*-5-15)
                cgfont('Helvetica',40)
                cgflip(w.back)
                waitkeydown(inf)
                
                cgdrawsprite(17,0,140)
                cgtext('Is there an activated bomb?',0,w.line*5)
                cgfont('Helvetica',35)
                cgtext('= Z key                  = X key                    = C key',0,w.line*-2)
                cgdrawsprite(5,0,w.line*-1)
                cgtext('For the BOMB/NO BOMB task, use your left hand',0,w.line*-3-15)
                cgtext('on the Z, X and C keys',0,w.line*-4-15)
                cgfont('Helvetica',20)
                cgtext('The experimenter will tell you which buttons to use in the scanner later on!',0,w.line*-5-20)
                cgfont('Helvetica',40)
                cgflip(w.back)
                waitkeydown(inf)
            case 2
                
                cgdrawsprite(17,0,140)
                cgtext('Accept, Reject, Explore?',0,w.line*5)
                cgfont('Helvetica',35)
                cgtext('= Z key                  = X key                    = C key',0,w.line*-2)
                cgdrawsprite(4,0,w.line*-1)
                cgtext('Now you will do a quick practice of the task',0,w.line*-3-15)
                cgtext('For this practice session, please use your',0,w.line*-4-15)
                cgtext('Use your LEFT hand for the ACCEPT/REJECT task',0,w.line*-5-15)
                cgfont('Helvetica',40)
                cgflip(w.back)
                waitkeydown(inf)
                
                cgdrawsprite(17,0,140)
                cgtext('Is there an activated bomb?',0,w.line*5)
                cgfont('Helvetica',35)
                cgtext('=Left arrow          =Down arrow       =Right arrow',0,w.line*-2)
                cgdrawsprite(5,0,w.line*-1)
                cgtext('For the BOMB/NO BOMB task, use your RIGHT hand',0,w.line*-3-15)
                cgtext('and the above keys',0,w.line*-4-15)
                cgfont('Helvetica',20)
                cgtext('The experimenter will tell you which buttons to use in the scanner later on!',0,w.line*-5-20)
                cgfont('Helvetica',40)
                cgflip(w.back)
                waitkeydown(inf)
        end
        
        
        
        % Outcomes omitted + do later
        cgtext('This stage of the experiment will be different',0,w.line*5)
        cgtext('from the previous stages in one more important way', 0, w.line*4)
        cgtext('In the previous stages of the experiment, the',0,w.line*3-15)
        cgtext('computer would tell you how much money you', 0, w.line*2-15)
        cgtext('would have won or lost on every single trial', 0, w.line*1-15)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgtext('This stage of the experiment will be different',0,w.line*5)
        cgtext('from the previous stages in one more important way', 0, w.line*4)
        cgtext('In the previous stages of the experiment, the',0,w.line*3-15)
        cgtext('computer would tell you how much money you', 0, w.line*2-15)
        cgtext('would have won or lost on every single trial', 0, w.line*1-15)
        cgtext('In this stage, the computer will only show you the ', 0, w.line*-1)
        cgtext('outcome on 50% of the trials. On the other 50%,', 0, w.line*-2)
        cgtext('the computer will take note of your outcome,', 0, w.line*-3)
        cgtext('but it won''t show you how much you got -', 0, w.line*-4)
        cgtext('it will just go straight to the next trial ', 0, w.line*-5)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgtext('Even though you are not seeing the outcome, the',0,w.line*5)
        cgtext('computer is still keeping track of your winnings',0,w.line*4)
        cgtext('We will tally your winnings at the end',0,w.line*3)
        cgtext('You also should have a good guess of how much',0,w.line*1)
        cgtext('money you are winning anyway, since you know',0,w.line*0)
        cgtext('how likely there is to be an activated bomb or not',0,w.line*-1)
        cgtext('on that trial',0,w.line*-2)
        cgflip(w.back)
        waitkeydown(inf)
        
        
        cgtext('On those 50% of trials where the computer is',0,w.line*5)
        cgtext('not showing you the outcome, the Exploring will', 0, w.line*4)
        cgtext('also work slightly differently', 0,w.line*3)
        cgflip(w.back)
        waitkeydown(inf)        
        
        cgtext('On those 50% of trials where the computer is',0,w.line*5)
        cgtext('not showing you the outcome, the Exploring will', 0, w.line*4)
        cgtext('also work slightly differently', 0,w.line*3)
        cgtext('If you choose to Explore on those trials,',0,w.line*1)
        cgtext('you will NOT immediately go through the procedure', 0, w.line*0)
        cgtext('of getting more  information (and making your,', 0,w.line*-1)
        cgtext(' response) during this session', 0,w.line*-2)
        cgtext('You will do the Exploring for those trials LATER',0,w.line*-4)
        cgflip(w.back)
        waitkeydown(inf)        
        
        cgtext('After you are done with this session, you will do',0,w.line*5+20)
        cgtext('another session where you see the trials in which',0,w.line*4+20)
        cgtext('you chose to explore. THEN, you will go through',0,w.line*3+20)
        cgtext('the process of getting more information, and',0,w.line*2+20)
        cgtext(' then making your choice',0,w.line*1+20)
        cgdrawsprite(17,0,-120)
        cgpencol(p.disp.color_nobomb)
        cgellipse(w.egpos(4),-120+40,60,60)
        cgellipse(w.egpos(3),-120+40,60,60)
        cgpencol(0,0,0)
        cgfont('Helvetica', 30)
        cgtext('You chose EXPLORE',0,w.line*-1+30)
        cgtext('Accept or Reject?',0,w.line*-4-25)
        cgfont('Helvetica', 40)
        cgflip(w.back)
        waitkeydown(inf)        
        
        cgtext('This should not affect your strategy in THIS task',0,w.line*5)
        cgtext('You should play every trial however you think is best', 0, w.line*4)
        cgtext('for you to win money or avoid losing money',0,w.line*3)
        cgtext('You should choose to Explore whenever you think',0,w.line*2-20)
        cgtext('that it would be worthwhile, given the type of offer', 0, w.line*1-20)
        cgtext('that you are facing', 0, w.line*0-20)
        cgtext('Don''t hold back on Exploring, just because you are', 0, w.line*-2)
        cgtext('not immediately getting the information! All that will', 0, w.line*-3)
        cgtext('be done later on', 0, w.line*-4)
%         cgfont('Helvetica',30)
%         cgtext('Note: You will NOT lose money for the Explore decisions on trials aside', 0, w.line*-2+10)
%         cgtext('from the 10 trials that we choose at the end of the experiment,', 0, w.line*-2-25)
%         cgtext('to calculate your total winnings. ', 0,w.line*-2-60)
%         cgtext(' You should make the decision about whether or not to Explore based', 0, w.line*-2-105)
%         cgtext('on whether or not you think it would be useful - don''t hold back on', 0, w.line*-2-140)
%         cgtext('exploring just because you are not immediately getting the information!', 0,w.line*-2-175)
%         cgfont('Helvetica',40)
        cgflip(w.back)
        waitkeydown(inf)
        
        cgtext('You shouldn''t hold back from exploring just so',0,w.line*5)
        cgtext('that you will be able to skip the next session', 0, w.line*4)
        cgtext('- we will keep you here for a set amount of time,',0,w.line*3)
        cgtext('no matter how much exploring you do in this session',0,w.line*2)
        cgtext('Focus on exploring when you think the type of offer', 0, w.line*0)
        cgtext('you are facing makes it worth-while!', 0, w.line*-1)
        cgfont('Helvetica', 30)
        cgtext('Note: This is effectively the same as in the previous stages! ', 0, w.line*-3-25)
        cgtext('We just have you do the ''Explore'' bit later, rather than immediately,', 0, w.line*-3-60)
        cgtext('just to save time for this particular session (for analysis reasons)', 0, w.line*-3-95)
        cgfont('Helvetica', 40)
        cgflip(w.back)
        waitkeydown(inf)
        
%         cgtext('This should not affect your strategy, however',0,w.line*5)
%         cgtext('You should play every trial however you think is best', 0, w.line*4)
%         cgtext('for you to win money or avoid losing money',0,w.line*3)
%         cgtext('You should choose to Explore whenever you think',0,w.line*2-20)
%         cgtext('that it would be worthwhile, given the type of offer', 0, w.line*1-20)
%         cgtext('that you are facing', 0, w.line*0-20)
%         cgfont('Helvetica',30)
%         cgtext('Note: You will NOT lose money for the Explore decisions on trials aside', 0, w.line*-2+10)
%         cgtext('from the 10 trials that we choose at the end of the experiment,', 0, w.line*-2-25)
%         cgtext('to calculate your total winnings. ', 0,w.line*-2-60)
% %         cgtext(' You should make the decision about whether or not to Explore based', 0, w.line*-2-105)
% %         cgtext('on whether or not you think it would be useful - don''t hold back on', 0, w.line*-2-140)
% %         cgtext('exploring just because you are not immediately getting the information!', 0,w.line*-2-175)
%         cgfont('Helvetica',40)
%         cgflip(w.back)
%         waitkeydown(inf)
%         
%         cgtext('This should not affect your strategy, however',0,w.line*5)
%         cgtext('You should play every trial however you think is best', 0, w.line*4)
%         cgtext('for you to win money or avoid losing money',0,w.line*3)
%         cgtext('You should choose to Explore whenever you think',0,w.line*2-20)
%         cgtext('that it would be worthwhile, given the type of offer', 0, w.line*1-20)
%         cgtext('that you are facing', 0, w.line*0-20)
%         cgfont('Helvetica',30)
%         cgtext('Note: You will NOT lose money for the Explore decisions on trials aside', 0, w.line*-2+10)
%         cgtext('from the 10 trials that we choose at the end of the experiment,', 0, w.line*-2-25)
%         cgtext('to calculate your total winnings. ', 0,w.line*-2-60)
%         cgtext(' You should make the decision about whether or not to Explore based', 0, w.line*-2-105)
%         cgtext('on whether or not you think it would be useful - don''t hold back on', 0, w.line*-2-140)
%         cgtext('exploring just because you are not immediately getting the information!', 0,w.line*-2-175)
%         cgfont('Helvetica',40)
%         cgflip(w.back)
%         waitkeydown(inf)

        
        % Win money: subset of trials
        cgtext('Just like in previous stages, a proportion of',0,w.line*4)
        cgtext('the money you win in this game will go towards ', 0, w.line*3)
        cgtext('the extra money we pay you at the end of this ', 0, w.line*2)
        cgtext('experiment', 0, w.line*1)
%         cgtext('To determine how much extra money you', 0, w.line*3)
%         cgtext('will get at the end of the experiment, we will choose', 0, w.line*2)
%         cgtext('10 trials at random, and you will get however much', 0, w.line*1)
%         cgtext('money you won on those trials', 0, w.line*0)
%         cgtext('This means that every trial counts for real money,',0, w.line*-2)
%         cgtext('and you should try your best to do whatever ', 0, w.line*-3)
        cgtext('You should try your best to do whatever ', 0, w.line*-1)
        cgtext('would be the most likely to win you money', 0, w.line*-2)
        cgtext(' (or stop you losing money) on every single trial', 0, w.line*-3)
        cgflip(w.back)        
        waitkeydown(inf)        

        % This practice
        cgtext(['In this practice session, you will do ' num2str(p.npracticetrials/2) ' trials'],0,w.line*5)
        cgtext('of each task',0,w.line*4)
        cgtext('Try to get used to the speed of the task, as well as',0,w.line*2+20)
        cgtext('to the idea of not seeing the outcome on each trial',0,w.line*1+20)
%         cgtext('You should also try to Explore a little bit during this',0,w.line*-1+20)
%         cgtext('practice just to see how it will work later on',0,w.line*-2+20)
%         cgtext('(This practice session will not count for real money)',0,w.line*-3+20)
%         cgfont('Helvetica',20)
%         cgtext('NOTE: The probability of winning and losing in this practice may not match what you have learned',0,w.line*-4)
%         cgtext('about the colours - but this is just because it is a low number of trials. ',0,w.line*-5+25)
%         cgtext('When you are doing the real task later, the probabilities will be the same as what you learned',0,w.line*-5)
%         cgtext('in the previous session of the experiment',0,w.line*-5-25)
%         cgfont('Helvetica',40)
        cgflip(w.back)    
        waitkeydown(inf)
        
        
        % This practice
        cgtext(['In this practice session, you will do ' num2str(p.npracticetrials/2) ' trials'],0,w.line*5)
        cgtext('of each task',0,w.line*4)
        cgtext('Try to get used to the speed of the task, as well as',0,w.line*2+20)
        cgtext('to the idea of not seeing the outcome on each trial',0,w.line*1+20)
        cgtext('You should also try to Explore a little bit during this',0,w.line*-1+20)
        cgtext('practice just to see how it will work later on',0,w.line*-2+20)
        cgtext('(This practice session will not count for real money)',0,w.line*-3+20)
        cgfont('Helvetica',20)
        cgtext('NOTE: The probability of winning and losing in this practice may not match what you have learned',0,w.line*-4)
        cgtext('about the colours - but this is just because it is a low number of trials. ',0,w.line*-5+25)
        cgtext('When you are doing the real task later, the probabilities will be the same as what you learned',0,w.line*-5)
        cgtext('in the previous session of the experiment',0,w.line*-5-25)
        cgfont('Helvetica',40)
        cgflip(w.back)    
        waitkeydown(inf)


        % Call experimenter
        cgtext('Please call the experimenter now',0,w.line*2)
        cgfont('Helvetica',30)
        cgtext('Note: It is VERY important that you understand exactly how to win',0,w.line*0)
        cgtext('money on this task. Please think now about how the two tasks work,',0,w.line*-1)
        cgtext('and ask the experimenter if there is anything about the task',0,w.line*-2)
        cgtext('that you are unsure about',0,w.line*-3)
        cgfont('Helvetica',40)
        cgflip(w.back)    
        waitkeydown(inf)
        
        cgfont('Helvetica',30)
        cgtext('Experimenter message:',0,w.line*1)
        cgtext('(1) Remember both tasks?',0,w.line*0)
        cgtext('(2) Outcomes omitted makes sense? Exploring the same?',0,w.line*-1)
        cgfont('Helvetica',40)
        cgflip(w.back)    
        waitkeydown(inf)
    end
    
    % Last instruction before task starts
    cgtext('Please position your fingers on the response keys',0,w.line*4)
    cgtext('Response keys:',0,w.line*2)
    cgfont('Helvetica',30)
    switch p.counterbal_hands
        case 1
            cgtext('NO BOMB/BOMB TASK', -230, w.line*0)
            cgtext('ACCEPT/REJECT TASK', 230, w.line*0)
            cgfont('Helvetica',20)
            cgtext('(Left hand)', -230, w.line*0+35)
            cgtext('(Right hand)', 230, w.line*0+35)
            cgtext('[No bomb]            [Bomb]             [Explore]', -230, w.line*-1)
            cgtext('  = Z key                = X key              = C key', -230, w.line*-1-25)
            cgtext('  [Accept]            [Reject]                [Explore]', 230, w.line*-1)
            cgtext('= Left arrow      = Down arrow    = Right arrow', 230, w.line*-1-25)
        case 2
            cgtext('NO BOMB/BOMB TASK', 230, w.line*0)
            cgtext('ACCEPT/REJECT TASK', -230, w.line*0)
            cgfont('Helvetica',20)
            cgtext('(Left hand)', -230, w.line*0+35)
            cgtext('(Right hand)', 230, w.line*0+35)
            cgtext('[No bomb]            [Bomb]             [Explore]', 230, w.line*-1)
            cgtext('= Left arrow      = Down arrow    = Right arrow', 230, w.line*-1-25)
            cgtext('  [Accept]            [Reject]                [Explore]', -230, w.line*-1)
            cgtext('  = Z key                = X key              = C key', -230, w.line*-1-25)
    end
    
    cgfont('Helvetica',40)
    cgtext('Press the down arrow to start', 0, w.line*-4)
    cgflip(w.back)    
    waitkeydown(inf,100)    
   
end

%% PRACTICE: INTEGRATED FMRI TASK

if rdata(1,5)==1
   wt.taskname='ACCEPT/REJECT';
elseif rdata(1,5)==2
   wt.taskname='NO BOMB/BOMB';
end
cgtext(['Now you do a block of the ' wt.taskname ' trials'],0,w.line*2)
cgtext('Press any key to start',0,w.line*-3)
t.taskstart=cgflip(w.back);
cgtext('+',0,0)           

for trialnum=1:size(rdata,1)
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
        
       % Verify parameters-display match-up
       disp('--------------------------')
       disp(['Trial #: ' num2str(trialnum)])
       disp(['Task type: ' num2str(rdata(trialnum,5)) ' (1=Conflict, 2=Control)'])
       disp(['Offer: Colour ' num2str(rdata(trialnum,3)) ', ' num2str(2*rdata(trialnum,2)) ' activated'])
       disp(['Activated bomb: ' num2str(rdata(trialnum,6))])
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
            wt.responseok=1;
        else 
            rdata(trialnum,8)=0;
            rdata(trialnum,13)=0;
            wt.responseok=0;
        end
        
        % [If Explored + Outcome shown] 2ND CHOICE --------------
        if rdata(trialnum,14)==1 && rdata(trialnum,8)==3
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
                wt.posx=wt.posX(randi(size(wt.posX,1)));
                rdata(trialnum,16)=wt.posx;
                if wt.x<2
                    cgdrawsprite(7,wt.posx,wt.revealed_Y);
                    rdata(trialnum,17)=1;
                    rdata(trialnum,16)=wt.posx;
                else
                    rdata(trialnum,17)=0;
                    rdata(trialnum,16)=999;
                end
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
                wt.responseok=1;
            else 
                rdata(trialnum,10)=0;
                rdata(trialnum,13)=0;
                wt.responseok=0;
            end
            % Workspace for continuation
            wt.respcol=10;
            wt.lastofferonset=wt.offer2;
        else
            % Workspace for continuation
            wt.respcol=8;
            wt.lastofferonset=wt.offer1;
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
        else
            % Workspace for continuation
            wt.lastoftrialonset=wt.lastofferonset;
           wt.timebeforenexttrial=p.disp.firstoffer;
        end
        
       % [BREAKS] --------------------------------
       if trialnum~=size(rdata,1) && rdata(trialnum,5)~=rdata(trialnum+1,5) % Task transition
           wait(p.disp.firstoffer)
           if rdata(trialnum+1,5)==1  
               wt.taskname='ACCEPT/REJECT';
           else
               wt.taskname='NO BOMB/BOMB';
           end
           cgtext('Now you do a some practice trials of the ',0,w.line*2)
           cgtext([wt.taskname ' task'],0,w.line*1)
           cgtext('The trials will start in 10 seconds',0,w.line*-4)
           cgflip(w.back);
           wait(10*1000)
       end
       
       cgtext('+',0,0)
       waituntil(wt.lastoftrialonset*1000+ wt.timebeforenexttrial)
       
       % Verify choices recorded
       disp(['Response valid: ' num2str(rdata(trialnum,13))])
       disp(['Response 1: ' num2str(rdata(trialnum,8)) '  (1=Accept/No bomb, 2=Reject/Bomb, 3=Explore)'])
       disp(['Response 2: ' num2str(rdata(trialnum,10)) ' (1=Accept/No bomb, 2=Reject/Bomb)'])
       disp(['Outcome: ' num2str(rdata(trialnum,15))])
       disp('------------------------------------------------------')
%         wtt.goon=waitkeydown(inf,60);        
             
end

%% INTRUCTIONS FOR OMITTED OUTCOMES

cgflip(w.back)
cgtext('Thank you - you have finished the first practice task', 0, w.line*5)
cgtext('After you are done with this practice session, you ', 0, w.line*4)
cgtext('will do this task in the scanner', 0, w.line*3)
cgtext('After you are done with the task in the scanner,', 0, w.line*1)
cgtext('you will then do another short task', 0, w.line*0)
cgtext('In this last after-scanning task, you will do' , 0, w.line*-2)
cgtext('the EXPLORING for those trials that you had ', 0, w.line*-3)
cgtext('previously seen in the scanner, but weren''t', 0, w.line*-4)
cgtext(' immediately given a hint for', 0, w.line*-5)
cgflip(w.back)
waitkeydown(inf);

cgtext('Now we will show you how this last session goes', 0, w.line*2)
cgtext('Press any key to continue with the instructions', 0, w.line*1)
cgflip(w.back)
waitkeydown(inf);
        
        
for o1=1:1
    if testing==1
        
        % Create displays for instructions
        w.egback=[0.3 0.1 0.0];
        if p.task.envthreat_nlevels==6
            w.egpos=[175 105 35 -35 -105 -175];
        elseif p.task.envthreat_nlevels==4
            w.egpos=[105 35 -35 -105];
        else
            input('Error: no positions for tutorial!')
        end
        cgmakesprite(17,500,300, w.egback)
        cgsetsprite(17) 
        for j=1:p.task.envthreat_nlevels
            cgpencol(w.back)
            cgpenwid(6)
            cgellipse(w.egpos(j),40,60,60)
            cgellipse(w.egpos(j),-40,60,60)
            cgpencol(0,0,0)
        end
        cgpencol(w.back)
        cgellipse(w.egpos(3),40,60,60,'f')
        cgellipse(w.egpos(3),-40,60,60,'f')
        cgellipse(w.egpos(4),40,60,60,'f')
        cgellipse(w.egpos(4),-40,60,60,'f')
        cgpencol(0,0,0)
        cgsetsprite(0)

        % Trial structure
        cgdrawsprite(17,0,140)
        cgfont('Helvetica',30)
        cgtext('Accept, Reject, Explore?',0,265)
        cgtext('You chose EXPLORE',0,230)
        cgfont('Helvetica',40)
        cgtext('For each trial, the computer will show you a', 0, w.line*-1)
        cgtext('trial from earlier on which you wanted to EXPLORE', 0, w.line*-2)
        cgtext('Just like in the above display', 0, w.line*-3)
        cgflip(w.back)
        waitkeydown(inf);
        
%         cgdrawsprite(17,0,140)
%         cgfont('Helvetica',30)
%         cgtext('Accept, Reject, Explore?',0,265)
%         cgtext('You chose EXPLORE',0,230)
%         cgfont('Helvetica',40)
%         cgtext('If you chose to ACCEPT/REJECT, the computer', 0, w.line*-1)
%         cgtext('will then immediately show you what outcome', 0, w.line*-2)
%         cgtext('you got on that trial', 0, w.line*-3)
%         cgflip(w.back)
%         waitkeydown(inf);
%         cgdrawsprite(17,0,140)
%         cgtext('+ 40 p',0,268)
%         cgpencol(0, 0.8, 0)
%         cgellipse(w.egpos(3),40+140,66,66,'f')
%         cgellipse(w.egpos(3),-40+140,66,66,'f')
%         cgellipse(w.egpos(4),40+140,66,66,'f')
%         cgellipse(w.egpos(4),-40+140,66,66,'f')
%         cgpencol(0,0,0)
%         cgfont('Helvetica',30)
%         cgtext('You chose ACCEPT',0,230)
%         cgfont('Helvetica',40)
%         cgtext('If you chose to ACCEPT/REJECT, the computer', 0, w.line*-1)
%         cgtext('will then immediately show you what outcome', 0, w.line*-2)
%         cgtext('you got on that trial', 0, w.line*-3)
%         cgflip(w.back)
%         waitkeydown(inf);
        
        cgdrawsprite(17,0,140)
        cgfont('Helvetica',30)
        cgtext('Accept, Reject, Explore?',0,265)
        cgtext('You chose EXPLORE',0,230)
        cgpencol(0, 0.8, 0)
        cgellipse(w.egpos(3),40+140,64,64)
        cgellipse(w.egpos(4),40+140,64,64)
        cgpencol(0, 0, 0)
        cgfont('Helvetica',40)
        cgtext('Then the computer will give you more information', 0, w.line*-1+10)
        cgtext('about the trial, and you can make your decision', 0, w.line*-2+10)
        cgtext('Press the same keys as before, to make', 0, w.line*-3-10)
        cgtext('your choice. Unlike in previous stages,', 0, w.line*-4-10)
        cgtext('you can take however long you''d like, to decide', 0, w.line*-5-10)
        cgflip(w.back)
        waitkeydown(inf);
        
        cgdrawsprite(17,0,140)
        cgtext('+ 40 p',0,268)
        cgpencol(0, 0.8, 0)
        cgellipse(w.egpos(3),40+140,66,66,'f')
        cgellipse(w.egpos(3),-40+140,66,66,'f')
        cgellipse(w.egpos(4),40+140,66,66,'f')
        cgellipse(w.egpos(4),-40+140,66,66,'f')
        cgpencol(0,0,0)
        cgfont('Helvetica',30)
        cgtext('You chose EXPLORE',0,230)
        cgtext('Your second choice was ACCEPT',0,30)
        cgfont('Helvetica',40)
        cgtext('After you make your decision, the computer', 0, w.line*-1)
        cgtext('will then tell you what outcome', 0, w.line*-2)
        cgtext('you got on this trial', 0, w.line*-3)
        cgflip(w.back)
        waitkeydown(inf);
        
        % Instruction for BOMB/NO BOMB task
        cgdrawsprite(17,0,140)
        cgfont('Helvetica',30)
        cgtext('Is there an activated bomb?',0,265)
        cgpencol(0, 0.8, 0)
        cgellipse(w.egpos(3),40+140,64,64)
        cgellipse(w.egpos(4),40+140,64,64)
        cgpencol(0, 0, 0)
        cgfont('Helvetica',40)
        cgtext('Same with the NO BOMB/BOMB task:', 0, w.line*-1)
        cgtext('the computer will show you a trial on which you', 0, w.line*-2)
        cgtext('chose to EXPLORE. You will then get more ', 0, w.line*-3)
        cgtext('information before you make your choice', 0, w.line*-4)
        cgflip(w.back)
        waitkeydown(inf);
        
        % Distinguish the two tasks
%         cgdrawsprite(17,0,140)
%         cgtext('Is there an activated bomb?',0,260)
%         cgtext('Like before, you need to pay attention', 0, w.line*-1+15)
%         cgtext('to which task you are currently doing,', 0, w.line*-2+15)
%         cgtext('and use the correct keys to respond!', 0, w.line*-3+15)
%         cgtext('The computer will show the task question', 0, w.line*-4)
%         cgtext('on the screen, to help you out', 0, w.line*-5)
%         cgflip(w.back)
%         waitkeydown(inf);
    end
    
    
    % Last instruction before task starts
    cgtext('Please position your fingers on the response keys',0,w.line*4)
    cgtext('Response keys:',0,w.line*2)
    cgfont('Helvetica',30)
    switch p.counterbal_hands
        case 1
            cgtext('NO BOMB/BOMB TASK', -230, w.line*0)
            cgtext('ACCEPT/REJECT TASK', 230, w.line*0)
            cgfont('Helvetica',20)
            cgtext('(Left hand)', -230, w.line*0+35)
            cgtext('(Right hand)', 230, w.line*0+35)
            cgtext('[No bomb]            [Bomb]             [Explore]', -230, w.line*-1)
            cgtext('  = Z key                = X key              = C key', -230, w.line*-1-25)
            cgtext('  [Accept]            [Reject]                [Explore]', 230, w.line*-1)
            cgtext('= Left arrow      = Down arrow    = Right arrow', 230, w.line*-1-25)
        case 2
            cgtext('NO BOMB/BOMB TASK', 230, w.line*0)
            cgtext('ACCEPT/REJECT TASK', -230, w.line*0)
            cgfont('Helvetica',20)
            cgtext('(Left hand)', -230, w.line*0+35)
            cgtext('(Right hand)', 230, w.line*0+35)
            cgtext('[No bomb]            [Bomb]             [Explore]', 230, w.line*-1)
            cgtext('= Left arrow      = Down arrow    = Right arrow', 230, w.line*-1-25)
            cgtext('  [Accept]            [Reject]                [Explore]', -230, w.line*-1)
            cgtext('  = Z key                = X key              = C key', -230, w.line*-1-25)
    end
    cgfont('Helvetica',40)
    cgtext('Press the down arrow to start', 0, w.line*-4)
    cgflip(w.back)    
    waitkeydown(inf,100)     
    
end

%% PRACTICE: OMITTED OUTCOMES

pracdata=rdata;
rdata=pracdata(pracdata(:,8)==3 & pracdata(:,14)==0,:);
rdata(:,14)=1;

if rdata(1,5)==1
    w.task1='ACCEPT/REJECT';
    w.task2='BOMB/NO BOMB';
else
    w.task2='ACCEPT/REJECT';
    w.task1='BOMB/NO BOMB';
end

% Any Exploring - Omitted trials to practice with? 
if size(rdata,1)>0
    cgtext('Now you will do the omitted EXPLORING for',0,w.line*3)
    cgtext(['the '  w.task1 ' trials'],0,w.line*1)
    cgtext('Press any key to start',0,w.line*-3)
    t.taskstart=cgflip(w.back);
    waitkeydown(inf)
    cgtext('+',0,0)
else
    cgtext('Please call the experimenter',0,w.line*3)
    cgtext('Experimenter message: No omitted explorations',0,w.line*1)
    cgflip(w.back);
    waitkeydown(inf)
end

%

for trialnum=1:size(rdata,1)
    rdata(trialnum,15)=nan;
    rdata(trialnum,7)=trialnum;
    
    % Verify parameters-display match-up
    disp('--------------------------')
    disp(['Trial #: ' num2str(trialnum)])
    disp(['# of tokens: ' num2str(2*rdata(trialnum,2))])
    disp(['pBomb category/COLOUR: ' num2str(rdata(trialnum,3))])
    disp(['Activated bomb?: ' num2str(rdata(trialnum,6))])
    %             wtt.goon=waitkeydown(inf,60);
    
        for o1=1:1 % set up displays & details for this trial
            wt=[];
            wt.pos=[p.disp.tokenpos  zeros(p.task.envthreat_nlevels,1)]; % assemble (shuffled) positions
            wt.startpos=randi(p.task.envthreat_nlevels+1-rdata(trialnum,2),1);
            for i=wt.startpos:(wt.startpos+rdata(trialnum,2)-1)
                wt.pos(i,2)=1;
            end
            wt.posx(1:rdata(trialnum,2),1)=wt.pos((wt.pos(:,2)==1),1);
            % resets
            wk=[];
            % displays
            if rdata(trialnum,5)==1
                wt.task='Accept, reject or explore?';
                wt.keys=[p.buttons.key1;p.buttons.key2];
                wt.resps=4;
                switch rdata(trialnum,8)
                    case 1
                        wt.resp='You chose accept';
                    case 2
                        wt.resp='You chose reject';
                    case 3
                        wt.resp='You chose explore';
                end
            else
                wt.task='Is there an activated bomb?';
                wt.keys=[p.buttons.key4;p.buttons.key5];
                wt.resps=5;
                switch rdata(trialnum,8)
                    case 1
                        wt.resp='You chose bomb';
                    case 2
                        wt.resp='You chose no bomb';
                    case 3
                        wt.resp='You chose explore';
                end
            end
        end
        cgtext(wt.task,0,260)
        wt.fixate=cgflip(0.3,0.3,0.3);
        
        % [offer + first choice] -----------
        cgdrawsprite(rdata(trialnum,3)+10,0,0)
        cgtext(wt.task,0,260)
        cgtext('+',0,0)
        for j=1:rdata(trialnum,2)
            cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row1);
            cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row2);
        end
        waituntil(wt.fixate*1000+500+rand*500)
        wt.offer1=cgflip(w.back);
        cgdrawsprite(rdata(trialnum,3)+10,0,0)
        cgtext(wt.task,0,260)
        cgtext(wt.resp,0,160)
        cgtext('+',0,0)
        for j=1:rdata(trialnum,2)
            cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row1);
            cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row2);
        end
        waituntil(wt.offer1*1000+p.disp.firstoffer)
        wt.whichchoice=cgflip(w.back);
        
        % [explore if necessary]----------
        if rdata(trialnum,8)==3
            cgdrawsprite(rdata(trialnum,3)+10,0,0) % offer + explore
            cgtext(wt.task,0,260)
            cgtext(wt.resp,0,160)
            cgtext('+',0,0)
            for j=1:rdata(trialnum,2)
                cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row1);
                cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row2);
            end
            eval(['wt.revealed_y=p.disp.tokenY_row' num2str(rdata(trialnum, 4)) ';'])  % reveal status of half
            for j=1:rdata(trialnum,2)
                cgdrawsprite(8,wt.posx(j),wt.revealed_y);
            end
            wt.x=randi(2);
            if rdata(trialnum, 6)==1
                wt.posxB=wt.posx(randi(size(wt.posx,1)));
                if wt.x<2
                    cgdrawsprite(7,wt.posxB,wt.revealed_y);
                    rdata(trialnum,14)=1;
                else
                    rdata(trialnum,14)=0;
                end
            end
            cgtext(['(Information cost: - ' num2str(p.task.explorationcost*10) ' p)'],0,-130)
            cgdrawsprite(wt.resps,0,-200)
            cgrect(230,-200,120, 30, [0.1 0.1 0.1])
            waituntil(wt.whichchoice*1000+p.disp.firstoffer)
            clearkeys
            wt.choice=cgflip(w.back);
            wt.keyok=0;
            while wt.keyok==0
                clearkeys
                [wk.key1press wk.key1time wk.a]=waitkeydown(inf,[wt.keys]);
                if wk.a==1
                    rdata(trialnum,11)= wk.key1time(1)-wt.choice*1000;
                    switch wk.key1press
                        case wt.keys(1)
                            rdata(trialnum,10)=1;
                        case wt.keys(2)
                            rdata(trialnum,10)=2;
                    end
                    wt.keyok=1;
                end
            end
        end
        
        % [calculate outcome]-------
        if rdata(trialnum,8)==3
            wt.explorecost=p.task.explorationcost;
            wt.col=10;
            wt.infocost=['(Information cost: - ' num2str(p.task.explorationcost*10) ' p)'];
            wt.preoutcome=wt.keyok;
        else
            wt.explorecost=0;
            wt.col=8;
            wt.infocost=' ';
            wt.preoutcome=wt.whichchoice;
        end
        if rdata(trialnum,5)==1 % accept/reject task
            if rdata(trialnum,wt.col)==1 && rdata(trialnum,6)==0 % accept + no bomb
                rdata(trialnum,12)=1;
                rdata(trialnum,15)=rdata(trialnum,12)*rdata(trialnum,2)*2-wt.explorecost;
                wt.spritenum=9;
                wt.won=['+ ' num2str((rdata(trialnum,15))*10) ' p'];
            elseif rdata(trialnum,wt.col)==1 && rdata(trialnum,6)==1 % accept + bomb
                rdata(trialnum,12)=-1;
                rdata(trialnum,15)=p.task.fixedloss-wt.explorecost;
                wt.spritenum=10;
                wt.won=[num2str((rdata(trialnum,15))*10) ' p'];
            elseif rdata(trialnum,wt.col)==2 % reject
                rdata(trialnum,12)=0;
                rdata(trialnum,15)=-wt.explorecost;
                wt.spritenum=6;
                wt.won='?';
            else
                input('error')
            end
        elseif rdata(trialnum,5)==2 % bomb/no bomb task
            if rdata(trialnum,wt.col)-1== rdata(trialnum,6)
                rdata(trialnum,12)=1;
                rdata(trialnum,15)=2*rdata(trialnum,2)-wt.explorecost;
                wt.spritenum=9;
                wt.won=['+ ' num2str((rdata(trialnum,15))*10) ' p'];
            else
                rdata(trialnum,12)=0;
                rdata(trialnum,15)=0-wt.explorecost;
                wt.spritenum=18;
                wt.won=' 0 p';
            end
        end
        
        wt.spritenum
        
        % display outcome
        cgdrawsprite(rdata(trialnum,3)+10,0,0)
        cgtext('+',0,0)
        for j=1:rdata(trialnum,2)
            cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row1);
            cgdrawsprite(6,wt.posx(j),p.disp.tokenY_row2);
        end
        for j=1:rdata(trialnum,2) % win/loss tokens
            cgdrawsprite(wt.spritenum,wt.posx(j),p.disp.tokenY_row1);
            cgdrawsprite(wt.spritenum,wt.posx(j),p.disp.tokenY_row2);
        end
        cgfont('helvetica',70)
        cgtext(wt.won,0,160)
        cgfont('helvetica',40)
        cgtext(wt.infocost,0,-130)
        waituntil(wt.preoutcome*1000+p.disp.firstoffer)
        wt.outcome=cgflip(w.back);
        
        cgtext('+',0,0)
        waituntil(wt.outcome*1000+p.disp.outcome)
        
        % verify choices recorded
        disp(['1st response:' num2str(rdata(trialnum,8)) '  (1=bomb, 2=no bomb, 3=explore)'])
        if rdata(trialnum,8)==3
            disp('----')
            disp(['2nd response:' num2str(rdata(trialnum,10)) '  (1=bomb, 2=no bomb)'])
            disp(['did exploration reveal a bomb?: ' num2str(rdata(trialnum,14))])
            disp('----')
        end
        disp(['outcome: ' num2str(rdata(trialnum,12))])
        disp(['trial aborted?: ' num2str(rdata(trialnum,13))])
        if rdata(trialnum,13)==1
            disp(['how many trials:' num2str(size(rdata,1))])
        end
        %         wtt.goon=waitkeydown(inf,60);
        
        if trialnum~=size(rdata,1) && rdata(trialnum,5)~=rdata(trialnum+1,5)
            cgtext('Now you will do the omitted EXPLORING for',0,w.line*3)
            cgtext(['the '  w.task2 ' trials'],0,w.line*1)
            cgtext('Press any key to start',0,w.line*-3)
            cgflip(w.back);
            waitkeydown(inf)
            cgtext('+',0,0)
        end

end

%%  FINISH UP 

cgfont('Helvetica',40)
cgflip(0.5,0.5,0.5)
cgtext('Thank you - you have now completed ',0,w.line*5)
cgtext('the practice session',0,w.line*4)
cgtext('You will now go on to do the ACCEPT/REJECT', 0, w.line*2)
cgtext('and NO BOMB/BOMB tasks inside the scanner', 0, w.line*1)
cgtext('- just like you did in the FIRST practice session', 0, w.line*0)
cgflip(0.5,0.5,0.5);
waitkeydown(inf)

cgtext('After you come out of the scanner, you will',0,w.line*2)
cgtext('then do the omitted EXPLORING from the ',0,w.line*1)
cgtext('scanner session', 0, w.line*0)
cgtext('- just like you did in the SECOND practice session', 0, w.line*-1)
cgflip(0.5,0.5,0.5);
waitkeydown(inf)

cgtext('If the practice counted for money, you would',0,w.line*2)
cgtext([' have won + ' num2str(sum(rdata(:,15))*10+sum(pracdata(:,15))*10) '  p'],0,w.line*1)
cgtext('Please call the experimenter',0, w.line*-3)
cgflip(0.5,0.5,0.5);
waitkeydown(inf)

stop_cogent

% Save
practice.rdata=rdata;
practice.alldata=pracdata;
practice.settings=p;
filename=[w.subjname '_file_5practice.mat'];
save(filename, 'practice')

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
w.warning=' ';
catch
    w.warning='WARNING: data not transfered to Deleted Daily. Transfer files manually';
end
%
disp('-----------------------------------')
try
    disp(['Duration: ' num2str((t.end-t.startcogent)/60) ' minutes'])
catch
end
disp(['Winnings: ' num2str(sum(rdata(:,15))*10+sum(pracdata(:,15))*10)])
disp(['No. of omitted explorations: ' num2str(sum(rdata(:,5)==1 & rdata(:,8)==3 & rdata(:,14)==0)) '  Conflict, '  num2str(sum(rdata(:,5)==2 & rdata(:,8)==3 & rdata(:,14)==0)) '  Control'])
disp(w.warning)
disp('-----------------------------------')
   