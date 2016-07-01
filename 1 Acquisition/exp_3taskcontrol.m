% Behavioural test/training for control task
clear all;close all hidden;clc;

for o1=1:1 % Documentation 
%     Col 1:      Trial type (read # from grid, starting from bottom-left corner)
%     Col 2:      # Token-pairs offered (Trial type X, 1-6)
%     Col 3:      Set pBomb level  (Trial type Y, 1-6 in ascending order of probability)
%     Col 4:      Explored half (if applicable)
%     Col 5:      Explored position (if applicable)
%     Col 6:      Bomb present in Activated tokens?
%     Col 7:      Trial number
%     Col 8:      [1st Response] Yes bomb, No bomb, Explore?  (1=Yes, 2=No, 3=Explore)
%     Col 9:      [1st Response] RT
%     Col 10:     [2nd Response] Yes bomb or No bomb? 
%     Col 11:     [2nd Response] RT
%     Col 12:     Outcome/Accuracy (1=Gain/Correct, -0=Wrong) 
%     Col 13:     Trial aborted (to be repeated later) (1=Aborted)  
%     Col 14:     Did exploration reveal a bomb? (1=Yes, 0=No, 999=Not explored)
%     Col 15:     Outcome in amount (payment)
end
for o1=1:1 % Execution settings: Testing/Coding? 
    
    testing=1;
    exec.plotfigs=0;
    exec.checkstruc=1; % Check statistical structure of task (derived via probabilistic sampling)
    choicetask=1;
    %
    if testing==1
        w.subjname=input('Subject ID: ','s');
        w.screenmode=1;
    else
        w.subjname='t1';
        w.screenmode=0;
    end
    w.taskstage='taskcontrol';
    where=pwd;
    w.res=2;
end
for o1=1:1 % General parameters (Load from learning) 
    load ([where filesep 'Data' filesep w.subjname '_file_1learnenv.mat'])
    p=learnenv.settings;
    % New settings for Test phase
    p.n.reps_pertrialtype=12; % If you change this, CHECK the probabilistic structure!
    p.subsetpercent_control=1;
    %
    w.line=50;
    rand('state',sum(100*clock));
end
for o1=1:1 % Subject-specific parameters 
[p par Norm check]=f_generate_taskstruc(p,exec);
% Mark which side to reveal if explored
par=sortrows(par,1);
for i=1:size(par,1)
    par(i,4)=3-randi(2);
end
% Choose subset of trials
par(:,7)=rand(size(par,1),1); 
par=sortrows(par,7);
w.ntrials=ceil(p.subsetpercent_control*size(par,1));
par_original=par;
par=[];
par_original(:,7)=rand(size(par_original,1),1);
par_original(:,7)=0;
par=par_original(1:w.ntrials,:);
% Last randomization & set up for execution
par(:,7)=rand(size(par,1),1); 
par=sortrows(par,7);
par(:,7)=1:size(par,1);
rdata=par;
% Misc
w.back=p.disp.settokencolor;
w.task='Is there an activated bomb?';
end
for o1=1:1 % Checks for running 
diary([where filesep 'Data' filesep w.subjname '_diary_3taskcontrol']); diary on;
disp('------------------------------')
disp('[SESSION 3: TASK CONTROL]')
% % disp(['No. levels pBombEnv: ' num2str(p.task.envthreat_nlevels)])
% disp(['No. reps per trial type: ' num2str(p.n.reps_pertrialtype)])
% disp(['Exploration cost: ' num2str(p.task.explorationcost*10) ' p'])
disp(['Choice task on/off: ' num2str(choicetask)])
disp(['No. of trials: ' num2str(size(rdata,1)) ' (Subset=' num2str(p.subsetpercent_control) ')'])
disp(['Estimated time: ' num2str(4.25*size(rdata,1)/60) ' min'])
disp('------------------------------')
% input('Setup OK?  ');
end

% COGENT
config_display(w.screenmode,w.res, [0 0 0], [1 1 1], 'Helvetica', 40, 9,0, 0); % marker 5g start
config_keyboard;
start_cogent % if using testing laptops, w.res in config dis must be 3! also, w.lineoptions =-330;
cgtext('The computer will now give you', 0, w.line*4)
cgtext('instructions for this task. ', 0, w.line*3)
cgtext(' If the trials start without you getting', 0, w.line*1)
cgtext('instructions, please tell the experimenter', 0, w.line*0)
cgtext(' immediately', 0, w.line*-1)
cgtext('Press any key to continue with the instructions', 0, w.line*-3)
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
%         Sprite 10: [Outcome] Yellow circle (Loss token - Bomb)
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
cgpencol([p.disp.color_nowinnolose])
cgellipse(0,0,p.disp.tokensize+8,p.disp.tokensize+8,'f')
cgmakesprite(5,700,100, [1 0 0]) % Sprite 5: Response options
cgtrncol(5,'r') 
cgsetsprite(5)
% cgrect(-200,0,170, 60, [0.3 0.3 0.3])
% cgrect(200,0,170, 60, [0.3 0.3 0.3])
cgrect(-230,0,170, 60, [0.3 0.3 0.3])
cgrect(0,0,170, 60, [0.3 0.3 0.3])
cgrect(230,0,170, 60, [0.3 0.3 0.3])
cgpencol([0 0 0])
cgfont('Helvetica',30)
% cgtext('NO BOMB',-230,0)
% cgtext('BOMB',0,0)xc
% cgtext('EXPLORE',230,0)
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

%% Instructions

for o1=1:1 
if testing==1
    cgpencol(0,0,0)
    w.egback=[0.8 0.5 0.8];
    
    cgtext('Please read these instructions carefully, as they will',0,w.line*3)
    cgtext('help you figure out how to win money in the task',0,w.line*2)
    cgtext('Press any key to scroll through the instructions',0,w.line*0)
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
    cgsetsprite(0)
    
    % Probabilistic structure is the same as in previous stage
    cgtext('By now, you should have a good idea of how',0,w.line*5)
    cgtext('likely you are to encounter an activated bomb, given',0,w.line*4)
    cgtext('each type of background, as well as each ',0,w.line*3)
    cgtext('combination of background colour',0,w.line*2)
    cgtext('and number of activated tokens',0,w.line*1)
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('By now, you should have a good idea of how',0,w.line*5)
    cgtext('likely you are to encounter an activated bomb, given',0,w.line*4)
    cgtext('each type of background, as well as each ',0,w.line*3)
    cgtext('combination of background colour',0,w.line*2)
    cgtext('and number of activated tokens',0,w.line*1)
    cgtext('In this stage of the experiment, these',0,w.line*-2)
    cgtext('same likelihoods apply',0,w.line*-3)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgtext('i.e. If a certain combination of background colour &',0,w.line*3)
    cgtext('number of activated tokens meant a good chance',0,w.line*2)
    cgtext('encountering an activated bomb in the previous',0,w.line*1)
    cgtext('stage of the experiment, it will indicate the same', 0, w.line*0)
    cgtext('thing in this stage of the experiment',0,w.line*-1)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgtext('In this stage of the experiment, an activated',0,w.line*5)
    cgtext('bomb will NOT always lose you money',0,w.line*4)
    cgtext('Your job in this stage of the experiment is to',0,w.line*2)
    cgtext('decide whether or not you think there is an',0,w.line*1)
    cgtext('activated bomb or not, on each trial',0,w.line*0)
    cgtext('You will win money for correctly predicting',0,w.line*-2)
    cgtext('whether there is an activated bomb on each trial',0,w.line*-3)
    cgflip(w.back);
    waitkeydown(inf);
    
    
    % TRIAL SEQUENCE
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext('On each trial, the computer will show you the',0,w.line*-2)
    cgtext('offer, consisting of a number of activated tokens',0,w.line*-3)
    cgtext('and a coloured background',0,w.line*-4)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    switch p.counterbal_hands
        case 1
            cgtext('Press the Z key if you think there is NO activated',0,w.line*-2)
            cgtext('bomb in the offer, or the X key if you think there IS',0,w.line*-3)
            cgtext('an activated bomb in the offer',0,w.line*-4)
            cgtext('Use your LEFT HAND for this task',0,w.line*-5-15)
        case 2
            cgtext('Press the Left arrow key if you think there is NO',0,w.line*-2)
            cgtext('activated bomb in the offer, or the DOWN arrow if',0,w.line*-3)
            cgtext('you think there IS an activated bomb in the offer',0,w.line*-4)
            cgtext('Use your RIGHT HAND for this task',0,w.line*-5-15)
    end
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext('+ 40 p', 0, 240)
    cgpencol(0, 0.8, 0)
    cgellipse(w.egpos(3),40+140,66,66,'f')
    cgellipse(w.egpos(3),-40+140,66,66,'f')
    cgellipse(w.egpos(4),40+140,66,66,'f')
    cgellipse(w.egpos(4),-40+140,66,66,'f')
    cgpencol(0,0,0)
    cgtext('For every trial where you correctly predict whether',0,w.line*-2)
    cgtext('or not there is an activated bomb, you will win money',0,w.line*-3)
    cgtext('The amount of money you win is proportionate',0,w.line*-4)
    cgtext('to the number of activated tokens (10p per token)',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext(' 0 p', 0, 240)
    cgpencol(p.disp.color_nowinnolose)
    cgellipse(w.egpos(3),40+140,66,66,'f')
    cgellipse(w.egpos(3),-40+140,66,66,'f')
    cgellipse(w.egpos(4),40+140,66,66,'f')
    cgellipse(w.egpos(4),-40+140,66,66,'f')
    cgpencol(0,0,0)
    cgtext('For every trial where you are INCORRECT in your',0,w.line*-2)
    cgtext('guess, you will NEITHER win NOR lose money',0,w.line*-3)
    cgtext('Unlike in the last stage of the experiment, you will',0,w.line*-4-20)
    cgtext('NOT lose any money in this task',0,w.line*-5-20)
    cgflip(w.back);
    waitkeydown(inf);
    
    % EXPLORE
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext('You can also choose to Explore, if you''d like more',0,w.line*-2)
    cgtext('more information before making your guess',0,w.line*-3)
    switch p.counterbal_hands
        case 1
            cgtext('Press the C key with your left hand, if',0,w.line*-4-15)
        case 2
            cgtext('Press the RIGHT arrow key with your right hand, if',0,w.line*-4-15)
    end
    cgtext('you''d like to Explore',0,w.line*-5-15)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgfont('Helvetica',20)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgtext('RED OUTLINE =',-330,285)
    cgtext('There is a bomb',-330,260)
    cgtext('underneath this',-330,230)
    cgtext('token',-330,200)
    cgdraw(-255,260,-140,235)
    cgdraw(-150,252,-140,235)
    cgdraw(-155,225,-140,235)
%     cgdraw(-255,260,-40,260)
%     cgdraw(-40,260,-40,235)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgtext('GREEN OUTLINE =',325,285)
    cgtext('There is NO bomb',330,260)
    cgtext('underneath this',330,230)
    cgtext('token',330,200)
    cgdraw(255,260,140,235)
    cgdraw(150,252,140,235)
    cgdraw(155,225,140,235)
%     cgdraw(255,260,40,260)
%     cgdraw(40,260,40,235)
    cgpencol(0,0,0)
    cgfont('Helvetica',40)
    cgtext('If you choose to Explore, the computer will reveal',0,w.line*-1+20)
    cgtext('the ''status'' of 50% of the activated tokens, showing',0,w.line*-2+20)
    cgtext('whether there are bombs or not, underneath',0,w.line*-3+20)
    cgtext('A green outline indicates NO bomb, while a red',0,w.line*-4)
    cgtext('outline indicates a bomb, underneath each token',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgdrawsprite(5,0,-45)
    cgrect(230,-45,120, 30, [0.1 0.1 0.1])
    cgtext('After you learn the status of half of the activated',0,w.line*-2)
    cgtext('tokens, you will then have to decide whether you ',0,w.line*-3)
    cgtext('think there is likely to be an activated bomb ',0,w.line*-4)
    cgtext('or not, on that trial',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgtext('The information gained by Exploring can help you',0,w.line*-1+17)
    cgtext('decide which response to make',0,w.line*-2+17)
    cgflip(w.back);
    waitkeydown(inf);
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgtext('The information gained by Exploring can help you',0,w.line*-1+17)
    cgtext('decide which response to make',0,w.line*-2+17)
    cgtext('But this information is not free - each ''Explore''',0,w.line*-3)
    cgtext(['choice will cost you ' num2str(p.task.explorationcost*10) ' p, regardless of what'  ],0,w.line*-4)
    cgtext('choice you make after that',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext('You will have TWO seconds, for each response',0,w.line*-2-15)
    cgtext('you have to make.  If you do not respond in time,',0,w.line*-3-15)
    cgtext('the trial will be stopped, and you will have to',0,w.line*-4-15)
    cgtext('redo that trial again later in the game',0,w.line*-5-15)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Money: random sample
    cgtext('In this stage, every single trial has equal influence',0,w.line*5+20)
    cgtext('on how much extra money you win on the task', 0, w.line*4+20)
    cgtext('The amount of extra money we pay you at the end', 0, w.line*2)
    cgtext('of the experiment is proportionate to the total amount', 0, w.line*1)
    cgtext('of money you win on these tasks', 0, w.line*0)
    cgtext('This means that every trial counts for real money,',0, w.line*-2)
    cgtext('and you should try your best to do whatever ', 0, w.line*-3)
    cgtext('would be the most likely to win you money', 0, w.line*-4)
    cgtext(' on every single trial', 0, w.line*-5)
    cgflip(w.back)        
    waitkeydown(inf)
    
    cgtext('Please call the experimenter now',0,w.line*1)
    cgflip(w.back);
    waitkeydown(inf);
end

    cgtext('Please position your fingers on the response keys',0,w.line*3)
    cgdrawsprite(5,0,-5)
    switch p.counterbal_hands
        case 1
            cgtext('Use your LEFT HAND for this task',0,w.line*2)
            cgfont('Helvetica',30)
            cgtext('= Z key                  = X key                  = C key', 0, w.line*-1)
        case 2
            cgtext('Use your RIGHT HAND for this task',0,w.line*2)
            cgfont('Helvetica',30)    
            cgtext('= Left arrow        = Down arrow         = Right arrow', 0, w.line*-1)
    end
    cgfont('Helvetica',40)
    cgtext('Press the spacebar to start the task',0,w.line*-3-20)
    cgflip(w.back);
    waitkeydown(inf,71)

end

%%  EXECUTION LOOP
w.alltrials=0;
w.omissions=0;
abortedtrials=zeros(1,15);
w.explorations=0;
trialnum=1;
w.block=1;
cgtext('+',0,0)

while w.alltrials~=1
    rdata(trialnum,15)=nan;
    rdata(trialnum,7)=trialnum;
    
    % Verify parameters-display match-up
       disp('--------------------------')
       disp(['Trial #: ' num2str(trialnum)])
       disp(['Offer: Colour ' num2str(rdata(trialnum,3)) '(' p.disp.contingencycolours{rdata(trialnum,3),4} '), ' num2str(2*rdata(trialnum,2)) ' activated'])
       disp(['Bomb present=' num2str(rdata(trialnum,5)) ', Bomb active=' num2str(rdata(trialnum,6))])
%             wtt.goon=waitkeydown(inf,60);

        for o1=1:1 % Set up displays & details for this trial
            wt=[];   
            wt.pos=[p.disp.tokenpos  zeros(p.task.envthreat_nlevels,1)]; % Assemble (shuffled) positions
            wt.startpos=randi(p.task.envthreat_nlevels+1-rdata(trialnum,2),1);
            for i=wt.startpos:(wt.startpos+rdata(trialnum,2)-1)
                wt.pos(i,2)=1;
            end
            wt.posX(1:rdata(trialnum,2),1)=wt.pos((wt.pos(:,2)==1),1);
            % Resets
            wk=[];
        end
        cgtext(w.task,0,160)
        wt.fixate=cgflip(0.3,0.3,0.3);
        
        % [FIRST OFFER] -----------
        cgdrawsprite(rdata(trialnum,3)+10,0,0)
        cgtext(w.task,0,160)
        cgtext('+',0,0)
        for j=1:rdata(trialnum,2)
            cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
        end
        cgdrawsprite(5,0,-250)
        waituntil(wt.fixate*1000+500+rand*500)
        clearkeys
        wt.offer1=cgflip(w.back);
        [wk.key1press wk.key1time wk.a1]=waitkeydown(p.disp.firstoffer,[p.buttons.key4 p.buttons.key5 p.buttons.key6]);
        if wk.a1==0 % Aborted
            cgdrawsprite(rdata(trialnum,3)+10,0,0)
            for j=1:rdata(trialnum,2)
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
            end
            cgdrawsprite(5,0,-250)
            cgtext('NO RESPONSE',0,160)
            waituntil(wt.offer1*1000+p.disp.firstoffer)
            wt.outcome=cgflip(w.back);
            w.omissions=w.omissions+1;
            abortedtrials=vertcat(abortedtrials, rdata(trialnum,:));
            rdata(trialnum,13)=0;
            rdata(trialnum,8:12)=nan;
        else
            rdata(trialnum,13)=1;
            switch wk.key1press(1) % log responses
                case p.buttons.key4
                    rdata(trialnum,8)=1;
                case p.buttons.key5
                    rdata(trialnum,8)=2;
                case p.buttons.key6
                    rdata(trialnum,8)=3;
            end
            rdata(trialnum,9)= wk.key1time(1)-wt.offer1*1000;
        end
        
        % [EXPLORE/2ND OFFER] -------------
        if rdata(trialnum,8)==3 
            cgdrawsprite(rdata(trialnum,3)+10,0,0)
            cgtext('+',0,0)
            cgtext(w.task,0,160)
            for j=1:rdata(trialnum,2) 
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
            end
            cgdrawsprite(5,0,-250)
            cgrect(230,-250,120, 30, [0.1 0.1 0.1])
            eval(['wt.revealed_Y=p.disp.tokenY_row' num2str(rdata(trialnum, 4)) ';'])  % Reveal status of half
            for j=1:rdata(trialnum,2) 
                cgdrawsprite(8,wt.posX(j),wt.revealed_Y);
            end
            wt.x=randi(2);
            if rdata(trialnum, 6)==1
                wt.posx=wt.posX(randi(size(wt.posX,1)));
                rdata(trialnum,5)=wt.posx;
                if wt.x<2
                    cgdrawsprite(7,wt.posx,wt.revealed_Y);
                    rdata(trialnum,14)=1;
                else
                    rdata(trialnum,14)=0;
                end
            end
            cgtext(['Information cost: - ' num2str(p.task.explorationcost*10) ' p'],0,-130)
            waituntil(wt.offer1*1000+p.disp.firstoffer)
            clearkeys
            wt.offer2=cgflip(w.back);
            [wk.key2press wk.key2time wk.a2]=waitkeydown(p.disp.secondoffer,[p.buttons.key4 p.buttons.key5]);
            w.explorations=w.explorations+1;
            if wk.a2==0 % Aborted
                cgdrawsprite(rdata(trialnum,3)+10,0,0)
                for j=1:rdata(trialnum,2) 
                    cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
                    cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
                end
                cgdrawsprite(5,0,-250)
                cgrect(230,-200,120, 30, [0.1 0.1 0.1])
                for j=1:rdata(trialnum,2) 
                    cgdrawsprite(8,wt.posX(j),wt.revealed_Y);
                end
                if rdata(trialnum,14)==1
                    cgdrawsprite(7,wt.posx,wt.revealed_Y);
                else
                end
                cgtext(['Information cost: - ' num2str(p.task.explorationcost*10) ' p'],0,-130)
                cgtext('NO RESPONSE',0,160)
                waituntil(wt.offer2*1000+p.disp.secondoffer)
                wt.outcome=cgflip(w.back);
                w.omissions=w.omissions+1;
                abortedtrials=vertcat(abortedtrials, rdata(trialnum,:));
                rdata(trialnum,13)=0;
                rdata(trialnum,10:12)=nan;
            else
                rdata(trialnum,13)=1;
                switch wk.key2press(1) % log responses
                    case p.buttons.key4
                        rdata(trialnum,10)=1;
                    case p.buttons.key5
                        rdata(trialnum,10)=2;
                end
                rdata(trialnum,11)= wk.key2time(1)-wt.offer2*1000;
            end
        else
            wt.offer2=wt.offer1;
            rdata(trialnum,14)=999;
        end
        
        % [Calculate outcome]-------
        if rdata(trialnum,13)==1
            if rdata(trialnum,8)==3 % Explore?
                wt.explorecost=p.task.explorationcost;
%                 wt.infocost=' ';
                wt.infocost=['Information cost: - ' num2str(p.task.explorationcost*10) ' p'];
            else
                wt.explorecost=0;
                wt.infocost=' ';
            end
            if rdata(trialnum,8)-1==rdata(trialnum,6) % Correct?
                rdata(trialnum,12)=1;
            elseif rdata(trialnum,8)==3 && rdata(trialnum,10)-1==rdata(trialnum,6) 
                rdata(trialnum,12)=1;
            else
                rdata(trialnum,12)=0;
            end
            rdata(trialnum,15)=2*rdata(trialnum,2)*rdata(trialnum,12)-wt.explorecost;
            if rdata(trialnum,15)>0 % Display colours
                wt.won=['+ ' num2str((rdata(trialnum,15))*10) ' p'];
                wt.tokenspritenum=9;
            else
                wt.won='0 p';
                wt.tokenspritenum=10;
            end
        end
        
        % Display outcome
        if rdata(trialnum,13)==1
            cgdrawsprite(rdata(trialnum,3)+10,0,0)
            cgtext('+',0,0)
            cgdrawsprite(5,0,-250)
            for j=1:rdata(trialnum,2)
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(6,wt.posX(j),p.disp.tokenY_row2);
            end
            for j=1:rdata(trialnum,2) % win/loss tokens 
                cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row2);
            end
            cgfont('Helvetica',70)
            cgtext(wt.won,0,160)
            cgfont('Helvetica',40)
            cgtext(wt.infocost,0,-130)
           waituntil(wt.offer2*1000+p.disp.firstoffer)
           wt.outcome=cgflip(w.back);
        end
        
       % Breaks
       if trialnum==floor(size(par,1)*w.block/5)
           wait(2000)
           cgflip(0.5,0.5,0.5)
           cgtext('Please take a break if you like',0,100)
           cgtext(['You have completed ' num2str(w.block) '/5 of this task'],0,0)
           cgtext('Press any key to continue',0,-100)
           w.block=w.block+1;
           t.break=cgflip(0.5,0.5,0.5);
           waitkeydown(inf)
       end
       
       cgtext('+',0,0)
       waituntil(wt.outcome*1000+p.disp.outcome)
       
      % Email researcher if almost done
       if trialnum==size(rdata,1)-15
           try
                w.time=clock;
                f_sendemail('learnreward',['[GoalExplore] (' num2str(w.time(4)) ':', num2str(w.time(5)) 'hrs) ' w.subjname ' has 1.5 min left (' w.taskstage ')'], ['Session almost complete: ' w.taskstage])        
           catch
            end
       end
       
       % Verify choices recorded
       disp(['Response valid: ' num2str(rdata(trialnum,13))])
       disp(['Response 1:' num2str(rdata(trialnum,8)) '  (1=Accept, 2=Reject, 3=Explore)'])
       disp(['Response 2:' num2str(rdata(trialnum,10)) '  (1=Accept, 2=Reject, 3=Explore)'])
       disp(['Outcome: ' num2str(rdata(trialnum,12))])
%         wtt.goon=waitkeydown(inf,60);        
       
       % Done with trials yet?
       if rdata(trialnum,13)==0
            rdata=vertcat(rdata, rdata(trialnum,:));
       end
       if trialnum==size(rdata,1)
           w.alltrials=1;
       else
           trialnum=trialnum+1;
       end
       
       
       % Save workspace
       try
           save([where filesep 'Data' filesep w.subjname '_workspace_3taskcontrol'])
       catch
           save([w.subjname '_workspace_3taskcontrol'])
       end
end
   
%%  FINISH UP 

cgfont('Helvetica',40)
if choicetask==1
    [w.dat.choice]=f_choicetask_envcolours(p, p.n.reps_choicetask);
    outcomes.learningaccuracy=mean(w.dat.choice.data(:,6));
else
    outcomes.learningaccuracy=999;
end

cgflip(0.5,0.5,0.5)
cgtext('You have now completed',0,100)
cgtext('this part of the task',0,50)
cgtext('Please call the experimenter',0,-50)
t.end=cgflip(0.5,0.5,0.5);
waitkeydown(inf)
stop_cogent

% Post-hoc calcs
outcomes.n_omissions=w.omissions;
outcomes.n_explorations=w.explorations;
outcomes.net=sum(rdata((rdata(:,13)==1),15))-w.explorations*p.task.explorationcost;

% Save
w.dat.rdata=rdata;
w.dat.data=rdata(rdata(:,13)==1,:);
w.dat.par=par;
p.log.subject=w.subjname;
p.log.date=date;
p.log.clock=clock;
w.dat.settings=p;
w.dat.outcomes=outcomes;
w.dat.check=check;
[w.dat.behaviour] = f_plotindividualbehaviour(p,w.dat.data, w.subjname); % Plot subject's behaviour
try
    w.dat.abortedtrials=abortedtrials;
catch
    w.dat.abortedtrials='No aborted trials';
end
eval([w.taskstage '=w.dat;'])

filename=[w.subjname '_file_3' w.taskstage '.mat'];
eval(['save(filename, ''' w.taskstage ''')']) % Save file
try 
    movefile(filename, [where filesep 'Data'])
catch
    disp('EXPERIMENTER MESSAGE: Data saved in curerent directory, not in rdata folder')
end

try % Transfer file to DeletedDaily
    cd([where filesep 'Data'])
    files2transfer=dir([w.subjname '*']);
    movetowhere_drive='\\Asia\DeletedDaily'; %
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
    for i=1:size(files2transfer,1)
        copyfile([where filesep 'Data' filesep files2transfer(i).name])
    end
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
disp(['Score: ' num2str(outcomes.net) ' (' num2str(outcomes.net/(2520/2)) '%  = ' num2str(7*outcomes.net/(2520/2)) ')'])
disp(['# of Omissions: ' num2str(w.omissions)])
disp(['Learning accuracy: ' num2str(outcomes.learningaccuracy)])
disp('(See variable ''w.dat.behaviour'' for behavioural performance)')
disp(w.warning)
disp('-----------------------------------'); diary off