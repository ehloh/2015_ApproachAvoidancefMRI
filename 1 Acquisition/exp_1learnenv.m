% Learn p(BombEnv), via forced acceptace + reveal inactivated bombs
clear all; close all hidden; clc

for o1=1:1 % Documentation for par 
%     Col 1:      Trial type (read # from grid, starting from bottom-left corner)
%     Col 2:      # Token-pairs offered (Trial type X, 1-6)
%     Col 3:      Set pBomb level  (Trial type Y, 1-6 in ascending order of probability)
%     Col 4:      Monitoring task trial?
%     Col 5:      Bomb present in set?
%     Col 6:      Bomb present in Activated tokens?
%     Col 7:      Trial number
%     Col 8:      [Monitor] Responded correctly? (1=Hit, 0=Miss, 999=n/a)
%     Col 9:      [Monitor] RT
%     Col 10:     [BLANK]
%     Col 11:     [BLANK]
%     Col 12:     Outcome (1=Gain, -1=Loss, 0=No effect) 
%     Col 13:     [BLANK]
%     Col 14:     [BLANK]
%     Col 15:     Position of first (leftmost) token (1-6)
%
%
% CHOICE TASK ================
%
%     Col 1: Trial number
%     Col 2: Left Colour #
%     Col 3: Right Colour #
%     Col 4: Better colour (correct choice - 1=Left, 2=Right)
%     Col 5: Choice (1=Left, 2=Right)
%     Col 6: Accuracy
end
for o1=1:1 % Execution settings: Testing/Coding?  

    testing=0;
    exec.plotfigs=0;
    exec.checkstruc=1; % Check statistical structure of task (derived via probabilistic sampling)
    choicetask=1;
    %    
    if testing==1
        w.subjname=input('Subject ID: ','s');
        p.counterbal_hands=input('Counterbalanced hands (1 or 2) :  ');
        w.screenmode=1;
    else
        w.subjname='t2';
        w.screenmode=0;
        p.counterbal_hands=1;
    end
    w.res=2;
    w.taskstage='learnenv';
    where=pwd;
end
for o1=1:1 % Generate general parameters
    
    % TASK DESIGN ------------------------
    p.task.envthreat_nlevels=6; % How many levels of p(Bomb present)
    p.n.reps_pertrialtype=12; % If you're going to change this, CHECK the prob struc generated!
    p.task.pBombEnv=nan*zeros(p.task.envthreat_nlevels,1);
    for i=1:p.task.envthreat_nlevels
        p.task.pBombEnv(p.task.envthreat_nlevels+1-i,1)=i/p.task.envthreat_nlevels; % Contingency, in descending probability
    end
    p.task.explorationcost=2; % # of tokens
    p.task.fixedloss=-2*p.task.envthreat_nlevels;
    p.n.reps_choicetask=6;
    p.applysubjectivecolourorders=1;
    
    % Pragmatics -----------------------------
    where=pwd;
    
    % ----------------------------------------
    rand('state',sum(100*clock));
    p.monitoringtask_percent=0.1;
    p.n.ntrials=p.n.reps_pertrialtype*p.task.envthreat_nlevels*p.task.envthreat_nlevels;
    %
    if p.task.envthreat_nlevels==6
        p.disp.tokenpos=[-250; -150; -50;50;150;250];
    elseif p.task.envthreat_nlevels==4
        p.disp.tokenpos=[ -150; -50;50;150];
    else
        input('Error: Token positions incorrect');
    end
    p.disp.tokenY_row1=50; % Display
    p.disp.tokenY_row2=-50;
    p.disp.tokensize=80;
    p.disp.tokengap=20;
    p.disp.settokencolor=[0.8 0.8 0.8];    
    p.disp.set_contingencycolours={[0.1 0.3 0.7] 'blue'; [0.1 0.4 0.4] 'green'; [0.7 0.4 0] 'orange';[0.4 0 0.6] 'purple';[0.7 0.7 0.0] 'yellow';[0.8 0.5 0.5] 'pink'}; 
    p.disp.color_bomb=[0.7 0 0];
    p.disp.color_nobomb=[0 0.7 0];
    p.disp.color_nowinnolose=[1 1 0]; 
    p.disp.monitorlength=1300; % Timings
    p.disp.itifixation=750;
    p.disp.firstoffer=2000;
    p.disp.firstresponse=000;
    p.disp.secondoffer=p.disp.firstoffer;
    p.disp.secondresponse=000;
    p.disp.outcome=1000;
    switch p.counterbal_hands
        case 1 % Right hand =Accept/Reject/Explore
            p.buttons.key1=97; % Buttons
            p.buttons.key2=100;
            p.buttons.key3=98;
            p.buttons.key4=26; % z
            p.buttons.key5=24; % x
            p.buttons.key6=3; % c
        case 2 % Left hand =No Bomb/Bomb/Explore
            p.buttons.key4=97; % Buttons
            p.buttons.key5=100;
            p.buttons.key6=98;
            p.buttons.key1=26; % z
            p.buttons.key2=24; % x
            p.buttons.key3=3; % c
    end
    p.buttons.monitor=100;
    % Misc
    w.line=50;
    w.back=p.disp.settokencolor;
end
for o1=1:1 % Generate subject-specific parameters
[p par Norm check]=f_generate_taskstruc(p,exec);
% Randomize background pBomb colours
for i=1:size(p.disp.set_contingencycolours,1)
    p.disp.set_contingencycolours{i,3}=rand;
end
p.disp.set_contingencycolours=sortrows(p.disp.set_contingencycolours,3);
for i=1:p.task.envthreat_nlevels
   p.disp.contingencycolours{i,1}=p.disp.set_contingencycolours{i,1};
   p.disp.contingencycolours{i,4}=p.disp.set_contingencycolours{i,2};
   p.disp.colournames{i,1}=p.disp.set_contingencycolours{i,2};
end
% Mark monitoring task
par(:,4)=0;
par(:,7)=rand(size(par,1),1); 
par=sortrows(par,7);
par(1:ceil(size(par,1)*p.monitoringtask_percent),4)=1; % Somthing is not an integer here
% Last randomization & set up for execution
par(:,7)=rand(size(par,1),1); 
par=sortrows(par,7);
par(:,7)=1:size(par,1);
data=par;
end
for o1=1:1 % Checks for running 
diary([where filesep 'Data' filesep w.subjname '_diary_1learnenv']); diary on;
disp('------------------------------')
disp('[SESSION 1: LEARNING STAGE (pBombEnv)]')
disp(['No. levels pBombEnv: ' num2str(p.task.envthreat_nlevels)])
disp(['No. reps per trial type: ' num2str(p.n.reps_pertrialtype)])
disp(['Choice task on/off: ' num2str(choicetask)])
disp(['Estimated time: ' num2str(5.75*p.task.envthreat_nlevels*p.task.envthreat_nlevels*p.n.reps_pertrialtype/60) ' min'])
disp('------------------------------')
input('EXPERIMENTER: Subject turn away  ');
disp(['Reshuffle colours after learning, according to subjective order?   : ' num2str(p.applysubjectivecolourorders)])
disp('Experimenter colours')
disp(p.disp.colournames)
input('Hit enter to start           ');
end

% COGENT
config_display(w.screenmode,w.res, [0 0 0], [1 1 1], 'Helvetica', 40, 9,0, 0); % marker 5g start
config_keyboard;
start_cogent % if using testing laptops, w.res in config dis must be 3! also, w.lineoptions =-330;
cgtext('The experimenter should have given you', 0, w.line*4)
cgtext('instructions for this task. ', 0, w.line*3)
cgtext(' If the trials start without you getting', 0, w.line*1)
cgtext('instructions, please tell the experimenter', 0, w.line*0)
cgtext(' immediately', 0, w.line*-1)
cgtext('Press any key to start the experiment', 0, w.line*-3)
t.startcogent=cgflip(0.8,0.8,0.8);
waitkeydown(inf);

%% Set up stock displays

for o1=1:1 % Create stock displays 
% Sprite details
%         Sprite 11-16: Empty tokens for all pBomb levels (1-6)
%         Sprite 6: [Offer] White circle (Offered token)
%         Sprite 7: [Explore] - Bomb
%         Sprite 8: [Explore] - Not a bomb
%         Sprite 9: [Outcome] Green circle (Won token)
%         Sprite 10: [Outcome] Red circle (Loss token - Bomb)
%         Sprite 17: Empty tokens for tutorial
for i=1:p.task.envthreat_nlevels % Sprites 11-16 (Empty tokens for all levels)
    w.backcol=p.disp.contingencycolours{i,1};
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
cgmakesprite(5,700,100, [1 0 0]) % Sprite 5: Response options
cgtrncol(5,'r') 
cgsetsprite(5)
cgrect(-230,0,170, 60, [0.3 0.3 0.3])
cgrect(0,0,170, 60, [0.3 0.3 0.3])
cgrect(230,0,170, 60, [0.3 0.3 0.3])
cgpencol([0 0 0])
cgfont('Helvetica',30)
cgtext('ACCEPT',-230,0)
cgtext('REJECT',0,0)
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
    cgpencol(1,1,1)
    cgfont('Helvetica', 80)
    cgtext('HOW DOES THIS LOOK DFL:KDFFASDF?',0,0)
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
    cgfont('Helvetica', 40)
end
end
w.back=p.disp.settokencolor;
end

%% Instructions
for o1=1:1 
if testing==1111 % Instructions off
    cgpencol(0,0,0)
    w.egback=[0.8 0.5 0.8];
    
    cgtext('STAGE 1: LEARNING PHASE',0,w.line*3)
    cgtext('Press any key to scroll through the instructions',0,w.line*-2)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Create displays for the instructions
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
    cgtext('+',0,0)
    cgsetsprite(0)
    
    % Introductory
    cgtext('In this experiment, you will have the chance to',0,w.line*5)
    cgtext('win money for yourself, by playing a',0,w.line*4)
    cgtext('decision-making game',0,w.line*3)    
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('In this experiment, you will have the chance to',0,w.line*5)
    cgtext('win money for yourself, by playing a',0,w.line*4)
    cgtext('decision-making game',0,w.line*3)    
    cgtext('There are several stages in this experiment',0,w.line*1+20)
    cgtext('In this 1st stage, you will learn how the game works,',0,w.line*0)
    cgtext('and you will learn to predict the likelihood of winning',0,w.line*-1)
    cgtext('on each trial, based on different variables in the task',0,w.line*-2)
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('In this experiment, you will have the chance to',0,w.line*5)
    cgtext('win money for yourself, by playing a',0,w.line*4)
    cgtext('decision-making game',0,w.line*3)    
    cgtext('There are several stages in this experiment',0,w.line*1+20)
    cgtext('In this 1st stage, you will learn how the game works,',0,w.line*0)
    cgtext('and you will learn to predict the likelihood of winning',0,w.line*-1)
    cgtext('on each trial, based on different variables in the task',0,w.line*-2)
    cgtext('All this will help you win more money later on in',0,w.line*-4)
    cgtext('the task, so pay close attention!',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Aim for 1st stage
    cgtext('On each trial, you will either win or lose money from',0,w.line*5+20)
    cgtext('the computer. Whether you win or lose is decided',0,w.line*4+20)
    cgtext('entirely by the computer, for now',0,w.line*3+20)
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('On each trial, you will either win or lose money from',0,w.line*5+20)
    cgtext('the computer. Whether you win or lose is decided',0,w.line*4+20)
    cgtext('entirely by the computer, for now',0,w.line*3+20)
    cgtext('Although you cannot control whether you win or lose',0,w.line*2)
    cgtext('money in this stage of the experiment, you''ll be able',0,w.line*1)
    cgtext('to influence your winnings later on in the experiment',0,w.line*0)
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('On each trial, you will either win or lose money from',0,w.line*5+20)
    cgtext('the computer. Whether you win or lose is decided',0,w.line*4+20)
    cgtext('entirely by the computer, for now',0,w.line*3+20)
    cgtext('Although you cannot control whether you win or lose',0,w.line*2)
    cgtext('money in this stage of the experiment, you''ll be able',0,w.line*1)
    cgtext('to influence your winnings later on in the experiment',0,w.line*0)
    cgtext('For now, you should try to figure out when you are',0,w.line*-1-20)
    cgtext('likely to win money and when you are likely to lose',0,w.line*-2-20)
    cgtext('This will help you figure out what you should do',0,w.line*-4)
    cgtext('to maximise your winnings later on',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Trial sequence for 1st stage
    cgdrawsprite(17,0,148)
    cgtext('For every trial, the likelihood of winning or losing is',0,w.line*-1+20)
    cgtext('predicted by the picture that is shown on the screen',0,w.line*-2+20)
    cgtext(['This picture tells you how many tokens (out of ' num2str(2*p.task.envthreat_nlevels) ')'],0,w.line*-3)
    cgtext('are activated on each trial',0,w.line*-4)
    cgflip(w.back);    
    waitkeydown(inf);
    cgdrawsprite(17,0,148)
    cgtext('For every trial, the likelihood of winning or losing is',0,w.line*-1+20)
    cgtext('predicted by the picture that is shown on the screen',0,w.line*-2+20)
    cgtext(['This picture tells you how many tokens (out of ' num2str(2*p.task.envthreat_nlevels) ')'],0,w.line*-3)
    cgtext('are activated on each trial',0,w.line*-4)
    cgpencol(0,0.8,0)
    cgfont('Helvetica',30)
    cgtext('4 tokens',330,260)
    cgtext('activated',330,230)
    cgtext(['(out of ' num2str(2*p.task.envthreat_nlevels) ')'],330,200)
    cgpenwid(4)
    cgdraw(270,260,0,260)
    cgdraw(0,260,0,235)
    cgpencol(0,0,0)
    cgfont('Helvetica',40) 
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148)
    cgpencol(0,0.8,0)
    cgfont('Helvetica',30)
    cgtext('4 tokens',330,260)
    cgtext('activated',330,230)
    cgtext(['(out of ' num2str(2*p.task.envthreat_nlevels) ')'],330,200)
    cgpenwid(4)
    cgdraw(270,260,0,260)
    cgdraw(0,260,0,235)
    cgpencol(0,0,0)
    cgfont('Helvetica',40) 
    cgtext('The tokens tell you how much money you can win,',0,w.line*-1+20)
    cgtext('with more tokens indicating more money at stake',0,w.line*-2+20)
    cgtext(['e.g.  4 tokens activated (out of ' num2str(2*p.task.envthreat_nlevels) ') means that '],0,w.line*-3)
    cgtext('you will win 4 x 10 p, if the computer decides',0,w.line*-4)
    cgtext('to give you money on this trial',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) % Introduce loss
    cgpencol(0,0.8,0)
    cgfont('Helvetica',30)
    cgtext('4 tokens',330,260)
    cgtext('activated',330,230)
    cgtext(['(out of ' num2str(2*p.task.envthreat_nlevels) ')'],330,200)
    cgpenwid(4)
    cgdraw(270,260,0,260)
    cgdraw(0,260,0,235)
    cgpencol(0,0,0)
    cgfont('Helvetica',40) 
    cgtext('BUT, the number of tokens also changes how likely',0,w.line*-1+20)
    cgtext('you are to LOSE money on a given trial',0,w.line*-2+20)
    cgtext('On some trials, there will be a bomb under one of the',0,w.line*-3)
    cgtext('tokens. If the bomb is underneath one of the activated',0,w.line*-4)    
    cgtext('tokens, you will LOSE money instead of winning',0,w.line*-5)    
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgtext('Only bombs that are under the activated tokens WILL',0,w.line*-1+20)
    cgtext('cause you to lose money',0,w.line*-2+20)
    cgflip(w.back);
    waitkeydown(inf);
    cgdrawsprite(17,0,148) 
    cgpencol(0,0.8,0)
    cgfont('Helvetica',23)
    cgtext('If there is a',325,275)
    cgtext('bomb under one',325,253)
    cgtext('one of these',325,231)
    cgtext('tokens, you',325,209)
    cgtext('will lose money',325,187)
    cgpenwid(3)
    cgdraw(260,275,-30,275)
    cgdraw(-30,275,-30,235)
    cgpencol(0,0,0)
    cgfont('Helvetica',40) 
    cgtext('Only bombs that are under the activated tokens will',0,w.line*-1+20)
    cgtext('cause you to lose money',0,w.line*-2+20)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgtext('Only bombs that are under the activated tokens will',0,w.line*-1+20)
    cgtext('cause you to lose money',0,w.line*-2+20)
    cgtext('If there IS a bomb in the offer, but it falls underneath',0,w.line*-3)
    cgtext('one of the inactivated tokens, you will NOT',0,w.line*-4)
    cgtext('lose any money',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    cgdrawsprite(17,0,148) 
    cgpencol(0,0.8,0)
    cgfont('Helvetica',23)
    cgtext('If there is a',325,120)
    cgtext('bomb under one',325,98)
    cgtext('one of these',325,76)
    cgtext('tokens, you will',325,54)
    cgtext('NOT lose money',325,32)
    cgpenwid(3)
    cgdraw(253,52,100,52)
    cgdraw(100,52,100,67)
    cgpencol(0,0,0)
    cgfont('Helvetica',40) 
    cgtext('Only bombs that are under the activated tokens will',0,w.line*-1+20)
    cgtext('cause you to lose money',0,w.line*-2+20)
    cgtext('If there IS a bomb in the offer, but it falls underneath',0,w.line*-3)
    cgtext('one of the inactivated tokens, you will NOT',0,w.line*-4)
    cgtext('lose any money',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgtext('IF there is a bomb in the offer, it will be under',0,w.line*-1+25)
    cgtext('any one of the tokens, active or inactive',0,w.line*-2+25)
    cgtext('There is nothing in the task that will help you',0,w.line*-3)
    cgtext('figure out where a bomb might be',0,w.line*-4)
    cgtext('- that is randomly decided by the computer',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgtext('Thus, more activated tokens indicate greater',0,w.line*-1)
    cgtext('potential winnings, but also a greater chance',0,w.line*-2)
    cgtext('of losing money - since the odds of  encountering',0,w.line*-3)
    cgtext('an ''activated bomb'' are more likely with more',0,w.line*-4)
    cgtext('activated tokens',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Outcomes + run trial sequence
    cgdrawsprite(17,0,148)    
    cgpencol(0.8,0,0)
    cgellipse(w.egpos(3),40+148,66,66,'f')
    cgellipse(w.egpos(3),-40+148,66,66,'f')
    cgellipse(w.egpos(4),40+148,66,66,'f')
    cgellipse(w.egpos(4),-40+148,66,66,'f')
    cgpencol(0,0,0)
    cgtext(['- ' num2str(-10*p.task.fixedloss) ' p'],0,255)
    cgtext('On every trial, the computer will tell you whether',0,w.line*-1+25)
    cgtext('you won or lost money',0,w.line*-2+25)
    cgtext(['For each trial where you lose, you will lose ' num2str(-10*p.task.fixedloss) ' p,'],0,w.line*-3)
    cgtext('and you will know that there was a bomb',0,w.line*-4)
    cgtext('under one of the activated tokens',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgpencol(0,0,0)
    cgtext('But if you win money, there may or may not be',0,w.line*-1+25)
    cgtext('a bomb under one of the inactivated tokens',0,w.line*-2+25)
    cgtext('For win trials, the computer will first tell you how',0,w.line*-3)    
    cgtext('much money you won. Then, it will tell you whether',0,w.line*-4)    
    cgtext('or not there was an inactivated bomb present',0,w.line*-5)
    cgflip(w.back)
    waitkeydown(inf); % Trigger
    cgdrawsprite(17,0,148) 
    cgtext('But if you win money, there may or may not be',0,w.line*-1+25)
    cgtext('a bomb under one of the inactivated tokens',0,w.line*-2+25)
    cgtext('For win trials, the computer will first tell you how',0,w.line*-3)    
    cgtext('much money you won. Then, it will tell you whether',0,w.line*-4)    
    cgtext('or not there was an inactivated bomb present',0,w.line*-5)
    wt.trig1=cgflip(w.back);
    cgtext('But if you win money, there may or may not be',0,w.line*-1+25)
    cgtext('a bomb under one of the inactivated tokens',0,w.line*-2+25)
    cgtext('For win trials, the computer will first tell you how',0,w.line*-3)    
    cgtext('much money you won. Then, it will tell you whether',0,w.line*-4)    
    cgtext('or not there was an inactivated bomb present',0,w.line*-5)
    waituntil(wt.trig1*1000+100)
    wt.trig11=cgflip(w.back);
    cgdrawsprite(17,0,148) 
    cgtext('But if you win money, there may or may not be',0,w.line*-1+25)
    cgtext('a bomb under one of the inactivated tokens',0,w.line*-2+25)
    cgtext('For win trials, the computer will first tell you how',0,w.line*-3)    
    cgtext('much money you won. Then, it will tell you whether',0,w.line*-4)    
    cgtext('or not there was an inactivated bomb present',0,w.line*-5)
    waituntil(wt.trig11*1000+700)
    wt.trig2=cgflip(w.back);    
    cgdrawsprite(17,0,148) 
    cgpencol(0,0.8,0)
    cgpenwid(6)
    cgellipse(w.egpos(3),40+148,66,66,'f')
    cgellipse(w.egpos(3),-40+148,66,66,'f')
    cgellipse(w.egpos(4),40+148,66,66,'f')
    cgellipse(w.egpos(4),-40+148,66,66,'f')
    cgpencol(0.8,0,0)
    cgpencol(0,0,0)
    cgtext('+ 40 p ',0,255)
    cgtext('But if you win money, there may or may not be',0,w.line*-1+25)
    cgtext('a bomb under one of the inactivated tokens',0,w.line*-2+25)
    cgtext('For win trials, the computer will first tell you how',0,w.line*-3)    
    cgtext('much money you won. Then, it will tell you whether',0,w.line*-4)    
    cgtext('or not there was an inactivated bomb present',0,w.line*-5)
    waituntil(wt.trig2*1000+1500)
    wt.trig3=cgflip(w.back);    
    cgdrawsprite(17,0,148) 
    cgpencol(0,0.8,0)
    cgpenwid(6)
    cgellipse(w.egpos(3),40+148,66,66,'f')
    cgellipse(w.egpos(3),-40+148,66,66,'f')
    cgellipse(w.egpos(4),40+148,66,66,'f')
    cgellipse(w.egpos(4),-40+148,66,66,'f')
    cgpencol(0.8,0,0)
    cgellipse(w.egpos(2),-40+148,60,60) % inactive bomb
    cgpencol(0,0,0)
    cgtext('+ 40 p ',0,255)
    cgtext('Inactivated bomb present ',0,40)
    cgtext('But if you win money, there may or may not be',0,w.line*-1+25)
    cgtext('a bomb under one of the inactivated tokens',0,w.line*-2+25)
    cgtext('For win trials, the computer will first tell you how',0,w.line*-3)    
    cgtext('much money you won. Then, it will tell you whether',0,w.line*-4)    
    cgtext('or not there was an inactivated bomb present',0,w.line*-5)
    waituntil(wt.trig3*1000+1500)
    cgflip(w.back);    
    waitkeydown(inf);
    
    % pBomb in set (background colours)
    cgdrawsprite(17,0,148) 
    cgtext('There won''t always be a bomb in the set, however',0,w.line*-1)
    cgtext('On each trial, the probability that there will be a',0,w.line*-2)
    cgtext('bomb SOMEWHERE in the set is indicated by ',0,w.line*-3)
    cgtext('the background colour',0,w.line*-4)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgtext('You will have to figure out how likely a bomb',0,w.line*-1+25)
    cgtext('is to be in the set (activated or not), given each',0,w.line*-2+25)
    cgtext('of the different background colours',0,w.line*-3+25)
    cgtext('There''s nothing about the colours that intrinsically',0,w.line*-4+25)
    cgtext('predicts the likelihood of a bomb being around -',0,w.line*-5+25)
    cgtext('you just have to figure it out, based on experience',0,w.line*-6+25)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Incidental task
    cgdrawsprite(17,0,148) 
    cgpencol(1,1,1)
    cgtext('+',0,150)
    cgpencol(0,0,0)
    cgtext('While you try to figure out which colours are more',0,w.line*-1+25)
    cgtext('likely to indicate a bomb than others, you will',0,w.line*-2+25)
    cgtext('also have to watch the black cross in the middle',0,w.line*-3+25)
    cgtext('of the screen. Every once in a while, this black',0,w.line*-4+25)
    cgtext('cross will turn white in colour',0,w.line*-5+25)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,148) 
    cgpencol(1,1,1)
    cgtext('+',0,150)
    cgpencol(0,0,0)
    cgtext('You should press the DOWN ARROW',0,w.line*-1)
    cgtext('whenever you see the black cross turn white',0,w.line*-2)
    cgtext('This will only happen on some trials',0,w.line*-4+25)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Summary
    cgtext('To summarize, your task for this stage of the',0,w.line*5+20)
    cgtext('experiment is to:',0,w.line*4+20)
    cgtext('(1) Press the down arrow every time you see the',0,w.line*3)
    cgtext('black cross turn white',0,w.line*2)
    cgtext('(2) Figure out how likely a bomb is to be present,',0,w.line*0+20)
    cgtext('given each of the background colours',0,w.line*-1+20)
    cgfont('Helvetica',30)
    cgtext(['(Note: There are only ' num2str(p.task.envthreat_nlevels) ' different background colours, in total)'],0,w.line*-2+30)
    cgfont('Helvetica',40)
    cgtext('(3) Figure out how good your chances of winning are,',0,w.line*-3)
    cgtext('given different combinations of background colours',0,w.line*-4-20)
    cgtext('+ number of activated tokens',0,w.line*-5-20)
    cgflip(w.back);
    waitkeydown(inf)
    
    % Performance/attention warning
    cgtext('A proportion of the money that you win and lose',0,w.line*5+20)
    cgtext('in this task WILL count for real money that we pay',0,w.line*4+20)
    cgtext('you at the end of the study',0,w.line*3+20)
    cgtext('BUT, if you miss  too many of the white crosses,',0,w.line*2)
    cgtext('OR, if, by the end of this session, you do not have',0,w.line*1)
    cgtext('some idea of which background colours indicate a',0,w.line*0)
    cgtext('high likelihood of a bomb and which ones don''t,',0,w.line*-1)
    cgtext('you will NOT receive all of the extra money that',0,w.line*-2)
    cgtext('you win on the task',0,w.line*-3)
    cgflip(w.back);
    waitkeydown(inf)
    
    cgtext('It''s very important that you are 100% clear on how',0,w.line*2)
    cgtext('this task works and what you have to do',0,w.line*1)
    cgtext('If you have any questions, please ask the',0,w.line*-1)
    cgtext('researcher now!',0,w.line*-2)
    cgflip(w.back);
    waitkeydown(inf)
    
    cgtext('Please call the experimenter now',0,w.line*1)
    cgflip(w.back);
    waitkeydown(inf);
    
end

    cgtext('The task is about to start',0,w.line*2)
    cgtext('Please position your index finger on the down arrow',0,w.line*0)
    cgtext('Press the down arrow to start the task',0,w.line*-2)
    cgflip(w.back);
    waitkeydown(inf,100)
end

%%  EXECUTION LOOP
w.omissions=0;
w.block=1;

cgfont('Helvetica',70)
cgtext('+',0,0)

for trialnum=1:size(data,1);
        data(trialnum,14)=nan;
        data(trialnum,7)=trialnum;
        t.fixate=cgflip(w.back);
        
        for o1=1:1 % Set up displays & details for this trial
            wt=[];   
            wt.pos=[p.disp.tokenpos  zeros(p.task.envthreat_nlevels,1)]; % Assemble (shuffled) positions
            wt.startpos=randi(p.task.envthreat_nlevels+1-data(trialnum,2),1);
            for i=wt.startpos:(wt.startpos+data(trialnum,2)-1)
                wt.pos(i,2)=1;
            end
            wt.posX(1:data(trialnum,2),1)=wt.pos((wt.pos(:,2)==1),1);
            if data(trialnum,5)==1 && data(trialnum,6)==0 % Inactivated bomb displays (for training stage)
                wtb=[];
                wtb.pos=wt.pos((wt.pos(:,2)==0),1);
                wtb.choose=randi(size(wtb.pos,1));
                wtb.chooserow=randi(2);
                eval(['wtb.whereinactivebomb=[wtb.pos(wtb.choose,1); p.disp.tokenY_row' num2str(randi(2)) '];' ]);
            end
            % Resets
            wk=[];
        end
        
       % Verify parameters-display match-up ----------
       disp('######################################################################')
       disp(['Trial #: ' num2str(trialnum)])
       disp(['# of token-pairs: ' num2str(2*data(trialnum,2))])
       disp(['pBomb category/COLOUR: ' num2str(data(trialnum,3)) ' (' p.disp.colournames{data(trialnum,3)} ')' ])
       disp(['Activated bomb?: ' num2str(data(trialnum,6))])
       disp(['Inactivated bomb?: ' num2str(data(trialnum,5))])
%         wtt.goon=waitkeydown(inf,90);

        % [FIRST OFFER] -----------
        cgdrawsprite(data(trialnum,3)+10,0,0)
        cgtext('+',0,0)
        for j=1:data(trialnum,2)
            cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
        end
        waituntil(t.fixate*1000+700+rand*500)
        clearkeys
        t.offer1=cgflip(w.back);
        
        % [OUTCOME 1 (Real outcome)] -------------------
            cgdrawsprite(data(trialnum,3)+10,0,0)
            cgtext('+',0,0)
            if data(trialnum,6)==1 % Lose
                for j=1:data(trialnum,2)
                    cgdrawsprite(10,wt.posX(j),p.disp.tokenY_row1);
                    cgdrawsprite(10,wt.posX(j),p.disp.tokenY_row2);
                end
                cgfont('Helvetica',70)
                cgtext([num2str(p.task.fixedloss*10) ' p'],0,160)
                data(trialnum,12)=-1;
            else % Win
                for j=1:data(trialnum,2)
                cgdrawsprite(9,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(9,wt.posX(j),p.disp.tokenY_row2);
                end
                cgfont('Helvetica',70)
                cgtext(['+ ' num2str(2*data(trialnum, 2)*10) ' p'],0,160)
                data(trialnum,12)=1;
            end
           waituntil(t.offer1*1000+p.disp.firstoffer)
           t.outcome1=cgflip(w.back);
           
        % [OUTCOME 2 (Inactivated bombs)] -------------------
            cgdrawsprite(data(trialnum,3)+10,0,0)
            cgtext('+',0,0)
            if data(trialnum,6)==1 % Lose
                for j=1:data(trialnum,2)
                    cgdrawsprite(10,wt.posX(j),p.disp.tokenY_row1);
                    cgdrawsprite(10,wt.posX(j),p.disp.tokenY_row2);
                end
                cgfont('Helvetica',70)
                cgtext([num2str(p.task.fixedloss*10) ' p'],0,160)
                data(trialnum,12)=-1;
                w.secondoutcome=p.disp.firstoffer/2;
            else % Win
                for j=1:data(trialnum,2)
                cgdrawsprite(9,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(9,wt.posX(j),p.disp.tokenY_row2);
                end
                cgfont('Helvetica',70)
                cgtext(['+ ' num2str(2*data(trialnum, 2)*10) ' p'],0,160)
                data(trialnum,12)=1;
                if data(trialnum, 5)==1 % Inactive bomb
                    cgfont('Helvetica',40)
                    cgdrawsprite(7, wtb.whereinactivebomb(1), wtb.whereinactivebomb(2));
                    cgtext('Inactive bomb present',0,-130)
                    cgfont('Helvetica',70)
                else
                    cgfont('Helvetica',40)
                    cgtext('No inactive bomb',0,-130)
                    cgfont('Helvetica',70)
                end
                w.secondoutcome=p.disp.firstoffer/2; 
            end
           waituntil(t.outcome1*1000+w.secondoutcome)
           t.outcome2=cgflip(w.back);

        if data(trialnum, 4)==1
            % Display from previous screen
            cgdrawsprite(data(trialnum,3)+10,0,0)
            if data(trialnum,6)==1 % Lose
                for j=1:data(trialnum,2)
                    cgdrawsprite(10,wt.posX(j),p.disp.tokenY_row1);
                    cgdrawsprite(10,wt.posX(j),p.disp.tokenY_row2);
                end
                cgfont('Helvetica',70)
                cgtext([num2str(p.task.fixedloss*10) ' p'],0,160)
                data(trialnum,12)=-1;
            else % Win
                for j=1:data(trialnum,2)
                cgdrawsprite(9,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(9,wt.posX(j),p.disp.tokenY_row2);
                end
                cgfont('Helvetica',70)
                cgtext(['+ ' num2str(2*data(trialnum, 2)*10) ' p'],0,160)
                data(trialnum,12)=1;
                if data(trialnum, 5)==1 % Inactive bomb
                    cgfont('Helvetica',40)
                    cgdrawsprite(7, wtb.whereinactivebomb(1), wtb.whereinactivebomb(2));
                    cgtext('Inactive bomb present',0,-130)
                    cgfont('Helvetica',70)
                else
                    cgfont('Helvetica',40)
                    cgtext('No inactive bomb',0,-130)
                    cgfont('Helvetica',70)
                end
            end
            cgpencol(1,1,1)
            cgtext('+',0,0)
            cgpencol(0,0,0)
            wm.jit=randi(1000);
            waituntil(t.outcome2*1000+1000+wm.jit)
            clearkeys
            t.monitoron=cgflip(w.back);
            [wm.key wm.keytime wm.a]=waitkeydown(p.disp.monitorlength,p.buttons.monitor);
            if isempty(wm.key)==0
                data(trialnum,8)=1;
                data(trialnum,9)=wm.keytime-t.monitoron*1000;
            else
                data(trialnum,8)=0;
                data(trialnum,9)=999;
            end
        else
            data(trialnum,8)=999;
            data(trialnum,9)=999;
        end
           
       % Breaks
       if trialnum==floor(size(data,1)*w.block/5)       
           disp('######################################################################')
           disp('BREAK SCREEN: START --------------')
           cgfont('Helvetica',40)
           wait(1000)
           cgflip(0.5,0.5,0.5)
           cgtext('Please take a break if you like',0,200)
           cgtext('Please ask the experimenter if you are',0,100)
           cgtext('still confused about the task',0,50)
           cgtext(['You have completed ' num2str(w.block) '/5 of this stage'],0,-100)
           cgtext('Press any key to continue',0,-200)
           w.block=w.block+1;
           t.break=cgflip(0.5,0.5,0.5);
           waitkeydown(inf)           
           disp('BREAK SCREEN: END --------------')
           disp('######################################################################')
       end
       
       % Display performance
       if data(trialnum,4)==1
            disp(['Monitoring task: ' num2str(data(trialnum,4)) ' ------ Performance: ' num2str(data(trialnum,8)) '  (' num2str(data(trialnum,9)) ')'])
       end
       disp(['Outcome: ' num2str(data(trialnum, 12))])     
       
       % Save workspace
       try
           save([where filesep 'Data' filesep w.subjname '_workspace_1learnenv'])
       catch
           save([w.subjname '_workspace_1learnenv'])
       end
       
       % Email researcher if almost done
       if trialnum==size(data,1)-15
            try
               w.time=clock;
               f_sendemail('learnreward',['[GoalExplore] (' num2str(w.time(4)) ':', num2str(w.time(5)) 'hrs) ' w.subjname ' has 1.5 min left (' w.taskstage ')'], ['Session almost complete: ' w.taskstage])        
            catch
            end
       end
      
       % Wait for next trial
       cgtext('+',0,0)
       waituntil(t.outcome2*1000+p.disp.outcome)
       clock 
end
   
cgfont('Helvetica',40)

%%  FINISH UP 

% CHOICE TASK
if choicetask==1
    disp('######################################################################')
    disp('START CHOICE TASK')
    [w.dat.choice]=f_choicetask_envcolours(p, p.n.reps_choicetask);
    outcomes.learningaccuracy=mean(w.dat.choice.data(:,6));
    if p.applysubjectivecolourorders==1
        % Reshuffle background colours according to subjective value
        w.colourscores=zeros(6,1);  % Score colours
        for i=1:size(w.dat.choice.data,1)
            switch w.dat.choice.data(i,5)
                case 1
                    w.colourscores(w.dat.choice.data(i,2))=w.colourscores(w.dat.choice.data(i,2))+1;
                case 2
                    w.colourscores(w.dat.choice.data(i,3))=w.colourscores(w.dat.choice.data(i,3))+1;
            end
        end
        p.disp.contingencycolours_original=p.disp.contingencycolours;
        w.dat.choice.colourscores=w.colourscores;
        for i=1:6;  % p.disp.contingencycolours:
            p.disp.contingencycolours{i,2}=i; % Col 2: Original order
            p.disp.contingencycolours{i,3}=w.colourscores(i);
            p.disp.contingencycolours{i,4}=p.disp.colournames{i};
        end
        % Implement reshuffling
        p.disp.contingencycolours=sortrows(p.disp.contingencycolours,3);
        p.disp.colournames_original=p.disp.colournames; p.disp.colournames=[];
        for i=1:6; p.disp.colournames{i,1}=p.disp.contingencycolours{i,4}; end
        p.disp.contingencycolours=sortrows(p.disp.contingencycolours,-3);
    end
else
    outcomes.learningaccuracy=999;
end
disp('######################################################################')

cgflip(0.5,0.5,0.5)
cgtext('You have now completed',0,100)
cgtext('this part of the task',0,50)
cgtext('Please call the experimenter',0,-50)
t.end=cgflip(0.5,0.5,0.5);
waitkeydown(inf)
stop_cogent

% Post-hoc calcs
outcomes.n_monitoromissions=mean(data((data(:,8)<5),8));
outcomes.win=sum(data((data(:,12)>0),12).*data((data(:,12)>0),2))*2;
outcomes.loss=size(find(data(:,12)<0),1)*p.task.fixedloss;
outcomes.net=outcomes.win+outcomes.loss;

% Save
w.dat.data=data;
w.dat.par=par;
p.log.subject=w.subjname;
p.log.date=date;
p.log.clock=clock;
w.dat.settings=p;
w.dat.outcomes=outcomes;
w.dat.check=check;
try     
    w.dat.abortedtrials=abortedtrials;
catch
    w.dat.abortedtrials='No aborted trials';
end
eval([w.taskstage '=w.dat;'])
try
    w.savecomm=['save([where filesep ''Data'' filesep w.subjname ''_file_1learnenv.mat''], ''' w.taskstage ''');'];
    eval(w.savecomm)
catch
    disp('EXPERIMENTER MESSAGE: Data saved in curerent directory, not in Data folder')
    w.savecomm=['save([where filesep w.subjname ''_file_1learnenv.mat''], ''' w.taskstage ''');'];
    eval(w.savecomm)
end
try % File transfer to Deleted Daily
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
copyfile([where filesep 'Data' filesep w.subjname '_file_1learnenv.mat'])
w.warning=' ';
catch
    w.warning='WARNING: Data not transfered to Deleted Daily. Transfer files manually';
end

%
disp('-----------------------------------')
try
    disp(['Duration: ' num2str((t.end-t.startcogent)/60) ' minutes'])
catch
end
disp(['Score: ' num2str(outcomes.net)])
disp(['Monitoring score: ' num2str(outcomes.n_monitoromissions) ' (<80% = poor)'])
disp(['Learning accuracy: ' num2str(outcomes.learningaccuracy)])
disp(w.warning)
disp('Experimenter, check colour orders :   p.disp.contingencycolours')
disp('(Col 2=original, Col 3=Score, Col 4=Name)')
disp('-----------------------------------'); diary off