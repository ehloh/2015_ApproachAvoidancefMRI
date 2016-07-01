% Behavioural test/training for conflict task
clear all;close all hidden;clc;
taskversion='fMRI';

for o1=1:1 % Documentation 
%     Col 1:      Trial type (read # from grid, starting from bottom-left corner)
%     Col 2:      # Token-pairs offered (Trial type X, 1-6)
%     Col 3:      Set pBomb level  (Trial type Y, 1-6 in ascending order of probability)
%     Col 4:      Explored half (if applicable)
%     Col 5:      Explored position (if applicable)
%     Col 6:      Bomb present in Activated tokens?
%     Col 7:      Trial number
%     Col 8:      [1st Response] Accept, Reject, or Explore?  (1=Accept, 2=Reject, 3=Explore)
%     Col 9:      [1st Response] RT
%     Col 10:     [2nd Response] Accept or Reject? 
%     Col 11:     [2nd Response] RT
%     Col 12:     Outcome (1=Gain, -1=Loss, 0=No effect) 
%     Col 13:     Trial aborted (to be repeated later) (1=Aborted)  
%     Col 14:     Did exploration reveal a bomb? (1=Yes, 0=No, 999=Not explored)
%     Col 15:     Outcome magnitude
end
for o1=1:1 % Execution settings: Testing/Coding? 
    
    testing=0;
    exec.plotfigs=0;
    exec.checkstruc=1; % Check statistical structure of task (derived via probabilistic sampling)
    choicetask=0;
    %
    if testing==1
        w.subjname=input('Subject ID: ','s');
        w.taskrep=input('# times this stage repeated (1=1st time, 2=2nd, etc): ');
        w.screenmode=1;
    else
        w.subjname='t1';
        w.taskrep=1;
%         w.taskrep=input('# times this stage repeated (1=1st time, 2=2nd, etc): ');
        w.screenmode=0;
    end
    if w.taskrep==1
        w.taskstage='taskconflict';
%         w.subset=1;
    else
        w.taskstage=['taskconflict' num2str(w.taskrep)];
%         w.subset=1;
    end 
%     input(['Subset of trials: ' num2str(w.subset) ' OK?']);
    where=pwd;
    w.res=2;
end
for o1=1:1 % General parameters (Load from learning) 
    load ([where filesep 'Data' filesep w.subjname '_file_1learnenv.mat'])
    p=learnenv.settings;
    p.taskversion=taskversion;
    % New settings for Test phase
    p.n.reps_pertrialtype=12; % If you change this, CHECK the probabilistic structure!
    w.subset=1;
    p.subsetpercent_conflict=w.subset;
    input(['Subset of trials: ' num2str(w.subset) ' OK?']);
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
w.ntrials=ceil(p.subsetpercent_conflict*size(par,1));
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
w.task='Accept, Reject, or Explore?';
end
for o1=1:1 % Checks for running 
diary([where filesep 'Data' filesep w.subjname '_diary_2taskconflict']); diary on;
disp('------------------------------')
disp(['[SESSION 2: TASK CONFLICT (repetition # ' num2str(w.taskrep) ']'])
% disp(['No. levels pBombEnv: ' num2str(p.task.envthreat_nlevels)])
% disp(['No. reps per trial type: ' num2str(p.n.reps_pertrialtype)])
disp(['Exploration cost: ' num2str(p.task.explorationcost*10) ' p'])
disp(['Choice task on/off: ' num2str(choicetask)])
disp(['No. of trials: ' num2str(size(rdata,1)) ' (Subset=' num2str(p.subsetpercent_conflict) ')'])
disp(['Estimated time: ' num2str(4.25*size(rdata,1)/60) ' min'])
disp(['Task version: ' p.taskversion])
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
p.log.start=clock;
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
cgpencol(0,0,0)
for o1=1:1 
if testing==1
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
    cgtext('likely there is to be an ''activated'' bomb, given',0,w.line*4)
    cgtext('combination of background colour and number',0,w.line*3)
    cgtext(' of tokens',0,w.line*2)
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('By now, you should have a good idea of how',0,w.line*5)
    cgtext('likely there is to be an ''activated'' bomb, given',0,w.line*4)
    cgtext('combination of background colour and number',0,w.line*3)
    cgtext(' of tokens',0,w.line*2)
    cgtext('In this stage of the experiment, these',0,w.line*0+15)
    cgtext('same likelihoods apply',0,w.line*-1+15)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgtext('i.e. If a certain combination of background colour &',0,w.line*5)
    cgtext('number of activated tokens meant a good chance',0,w.line*4)
    cgtext('of there being an activated bomb in the previous',0,w.line*3)
    cgtext('stage of the experiment, it will indicate the same', 0, w.line*2)
    cgtext('thing in this stage of the experiment',0,w.line*1)
    cgflip(w.back);
    waitkeydown(inf);
    cgtext('i.e. If a certain combination of background colour &',0,w.line*5)
    cgtext('number of activated tokens meant a good chance',0,w.line*4)
    cgtext('of there being an activated bomb in the previous',0,w.line*3)
    cgtext('stage of the experiment, it will indicate the same', 0, w.line*2)
    cgtext('thing in this stage of the experiment',0,w.line*1)
    cgtext('In this stage of the experiment, however,',0,w.line*-1)
    cgtext('you will perform a slightly different task',0,w.line*-2)
    cgtext('Now you will learn how this task works',0,w.line*-4)
    cgflip(w.back);
    waitkeydown(inf);
   
    % Trial sequence
    cgdrawsprite(17,0,140)
    cgtext('On each trial, the computer will make you an',0,w.line*-1+15)
    cgtext('offer, consisting of a number of activated tokens',0,w.line*-2+15)
    cgtext('and a coloured background',0,w.line*-3+15)
    cgtext('Like before, the background colour indicates the',0,w.line*-4)
    cgtext('probability of there being a bomb in the offer',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgtext('Just like in the previous stage, more activated',0,w.line*-1+15)
    cgtext('tokens means greater potential winnings, but also a',0,w.line*-2+15)
    cgtext('greater likelihood of encountering an activated bomb',0,w.line*-3+15)
    cgfont('Helvetica',30)
    cgtext('Note: The offers indicate both (a) the amount of money you can win',0,w.line*-4)
    cgtext('as well as (b) how likely you are to win or lose. Please ask the ',0,w.line*-4-30)
    cgtext('experimenter if you''d like a reminder of how the game works!',0,w.line*-4-60)
    cgfont('Helvetica',40)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgtext('When you see the offer, you should have a rough',0,w.line*-1)
    cgtext('idea about whether there is likely to be an activated',0,w.line*-2)
    cgtext('bomb or not - since you''ve learned this in',0,w.line*-3)
    cgtext('the previous stage of the experiment',0,w.line*-4)
    cgflip(w.back);
    waitkeydown(inf);
   
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext('For each offer, you have the option of Accepting,',0,w.line*-2)
    cgtext('Rejecting, or Exploring the offer',0,w.line*-3)
    switch p.counterbal_hands
        case 1 
            cgtext('Make your choice by pressing the Left, Down or',0,w.line*-4-15)
            cgtext('Right arrow respectively (with your RIGHT hand)',0,w.line*-5-15)
        case 2            
            cgtext('Make your choice by pressing the Z key, X key',0,w.line*-4-15)
            cgtext('or C key respectively (with your LEFT hand)',0,w.line*-5-15)
    end
    cgflip(w.back);
    waitkeydown(inf);
    
    % You can lose money
    cgdrawsprite(17,0,140)
    cgdrawsprite(5,0,-45)
    cgtext('+ 40 p', 0, 240)
    cgpencol(0, 0.8, 0)
    cgellipse(w.egpos(3),40+140,66,66,'f')
    cgellipse(w.egpos(3),-40+140,66,66,'f')
    cgellipse(w.egpos(4),40+140,66,66,'f')
    cgellipse(w.egpos(4),-40+140,66,66,'f')
    cgpencol(0,0,0)
    cgtext('If you ACCEPT an offer with no activated bomb',0,w.line*-2)
    cgtext('you will win money - 10p per activated token',0,w.line*-3)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgtext([num2str(p.task.fixedloss*10) ' p'], 0, 240)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(3),40+140,66,66,'f')
    cgellipse(w.egpos(3),-40+140,66,66,'f')
    cgellipse(w.egpos(4),40+140,66,66,'f')
    cgellipse(w.egpos(4),-40+140,66,66,'f')
    cgpencol(0,0,0)
    cgtext('But if you Accept an offer WITH an activated',0,w.line*-1+20)
    cgtext('bomb, you will LOSE money',0,w.line*-2+20)
    cgtext(['Each time this happens, you will lose ' num2str(p.task.fixedloss*10) ' p'],0,w.line*-3)
    cgtext('The computer will tell you whether you have won or',0,w.line*-4-20)
    cgtext('lost money, on each trial',0,w.line*-5-20)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgtext('?', 0, 240)
    cgtext('You can also Reject offers, to avoid losing money',0,w.line*-1)
    cgtext('If you Reject the offer, you will not win or lose any',0,w.line*-2)
    cgtext('money, no matter what outcome goes with the offer',0,w.line*-3)
    cgtext('Also, the computer will NOT tell you what outcome',0,w.line*-4-20)
    cgtext('you would have gotten, if you had accepted the offer',0,w.line*-5-20)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgtext('If you always Reject offers, however, you will never',0,w.line*-1+10)
    cgtext('win any money. Sometimes it''s better to Explore',0,w.line*-2+10)
    cgtext('If you choose to Explore, the computer',0,w.line*-3-10)
    cgtext('will give you more information, before you choose',0,w.line*-4-10)
    cgtext('to Reject or Accept the offer',0,w.line*-5-10)
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
    cgtext('The computer will reveal the ''status'' of 50% of the',0,w.line*-1+20)
    cgtext('activated tokens offered, showing whether there are',0,w.line*-2+20)
    cgtext('bombs or not, underneath these tokens',0,w.line*-3+20)
    cgtext('For each token, a green outline indicates NO bomb, ',0,w.line*-4)
    cgtext('while a red outline indicates a bomb',0,w.line*-5)
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
    cgtext('tokens, you will then have to decide whether to ',0,w.line*-3)
    cgtext('Accept or Reject the offer',0,w.line*-4)
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
    cgtext('Once again, you will receive the associated outcome',0,w.line*-2)
    cgtext('if you Accept the offer, but you will not win or lose',0,w.line*-3)
    cgtext('money (nor learn what the outcome would have',0,w.line*-4)
    cgtext('been) if you choose to Reject the offer',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
 
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgtext('The computer won''t tell you about inactivated bombs',0,w.line*-1+15)
    cgtext('either, in all the rest of the experiment',0,w.line*-2+15)
    cgtext('They won''t have any effect at all on your winnings,',0,w.line*-3)
    cgtext('throughout the experiment',0,w.line*-4)
    cgflip(w.back);
    waitkeydown(inf);   
    
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgtext('The information gained by Exploring can help you',0,w.line*-1+17)
    cgtext('decide if it is worth risking loss, to gain more money',0,w.line*-2+17)
    cgflip(w.back);
    waitkeydown(inf);
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgtext('The information gained by Exploring can help you',0,w.line*-1+17)
    cgtext('decide if it is worth risking loss, to gain more money',0,w.line*-2+17)
    cgtext('This information is not free, however! Each ''Explore''',0,w.line*-3)
    cgtext(['choice will cost you ' num2str(p.task.explorationcost*10) ' p, regardless of whether'  ],0,w.line*-4)
    cgtext('you decide to Accept or Reject the offer later on',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgdrawsprite(17,0,140)
    cgpencol(0.8, 0, 0)
    cgellipse(w.egpos(4),40+140,60,60)
    cgpencol(0, 0.7, 0)
    cgellipse(w.egpos(3),40+140,60,60)
    cgpencol(0,0,0)
    cgtext('Don''t be afraid to explore! If you use the Explorations',0,w.line*-1)
    cgtext('well, you will end up winning more money',0,w.line*-2)
    cgtext('overall, so it''s worth your while to try it out if you',0,w.line*-3)
    cgtext('feel like you''d like more information before',0,w.line*-4)
    cgtext('making your decision',0,w.line*-5)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Timing
    cgtext('You will have TWO seconds, for each response',0,w.line*5)
    cgtext('you have to make.  If you do not respond in time,',0,w.line*4)
    cgtext('the trial will be stopped, and you will have to',0,w.line*3)
    cgtext('redo that trial again later in the game',0,w.line*2)
    cgflip(w.back);
    waitkeydown(inf);
    
    cgtext('You will have TWO seconds, for each response',0,w.line*5)
    cgtext('you have to make.  If you do not respond in time,',0,w.line*4)
    cgtext('the trial will be stopped, and you will have to',0,w.line*3)
    cgtext('redo that trial again later in the game',0,w.line*2)
    cgtext('Don''t worry if it seems like it''s going very',0,w.line*0)
    cgtext('quickly at first! You will get used to the pace,',0,w.line*-1)
    cgtext('and if you are too slow on a trial, you won''t lose,',0,w.line*-2)
    cgtext('any money because of it. All that will happen is that',0,w.line*-3)
    cgtext('you will have to do it later',0,w.line*-4)
    cgflip(w.back);
    waitkeydown(inf);
    
    % Money
    cgtext('Remember, you are playing for real money', 0, w.line*4)
    cgtext('in this game - the amount of extra money we', 0, w.line*3)
    cgtext('pay you at the end of this experiment will be', 0, w.line*2)
    cgtext('proportional to your winnings in this game', 0, w.line*1)
    cgtext('Therefore, you should whatever you think ', 0, w.line*-1)
    cgtext('would be the most likely to win you money, ', 0, w.line*-2)
    cgtext('on every single trial!', 0, w.line*-3)
    cgflip(w.back)        
    waitkeydown(inf)
    
    cgtext('Please call the experimenter now',0,w.line*1)
    cgflip(w.back);
    waitkeydown(inf);
end

    cgtext('Please position your fingers on the response keys',0,w.line*3)
    switch p.counterbal_hands
        case 1
            cgtext('Use your RIGHT HAND for this task',0,w.line*2)
        case 2
            cgtext('Use your LEFT HAND for this task',0,w.line*2)
    end
    cgfont('Helvetica',30)
    cgdrawsprite(5,0,-5)
    switch p.counterbal_hands
        case 1 
            cgtext('= Left arrow              = Down arrow             = Right arrow', 0, w.line*-1)    
        case 2            
            cgtext('=  Z Key                   = X Key                 = C Key', 0, w.line*-1)
    end
    cgfont('Helvetica',40)
    cgtext('Press the down arrow to start the task',0,w.line*-3-20)
    cgflip(w.back);
    waitkeydown(inf,100)

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
    disp('######################################################################')
   disp(['Trial #: ' num2str(trialnum)])
   disp(['Offer: Colour ' num2str(rdata(trialnum,3)) '(' p.disp.contingencycolours{rdata(trialnum,3),4} '), ' num2str(2*rdata(trialnum,2)) ' activated'])
   disp(['Active Bomb=' num2str(rdata(trialnum,6))])
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
        [wk.key1press wk.key1time wk.a1]=waitkeydown(p.disp.firstoffer,[p.buttons.key1 p.buttons.key2 p.buttons.key3]);
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
                case p.buttons.key1(1)
                    rdata(trialnum,8)=1;
                case p.buttons.key2
                    rdata(trialnum,8)=2;
                case p.buttons.key3
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
            [wk.key2press wk.key2time wk.a2]=waitkeydown(p.disp.secondoffer,[p.buttons.key1 p.buttons.key2]);
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
                % Print tokens as before
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
                    case p.buttons.key1
                        rdata(trialnum,10)=1;
                    case p.buttons.key2
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
                wt.respcol=10;
            else
                wt.explorecost=0;
                wt.infocost=' ';
                wt.respcol=8;
            end
            if rdata(trialnum,wt.respcol)==2 % Reject: nothing
                rdata(trialnum,12)=0;
                rdata(trialnum,15)=0-wt.explorecost;
            elseif rdata(trialnum,wt.respcol)==1 && rdata(trialnum,6)==1 % Accept +Bomb
                rdata(trialnum,12)=-1;
                rdata(trialnum,15)=p.task.fixedloss-wt.explorecost;
            elseif rdata(trialnum,wt.respcol)==1 && rdata(trialnum,6)==0 % Accept +No Bomb
                rdata(trialnum,12)=1;
                rdata(trialnum,15)=2*rdata(trialnum,2)-wt.explorecost;
            end
            if rdata(trialnum,15)>0 % Display colours
                wt.won=['+ ' num2str((rdata(trialnum,15))*10) ' p'];
                wt.tokenspritenum=9;
            elseif rdata(trialnum,15)==0
                if rdata(trialnum,wt.respcol)==1
                    wt.won='+ 0 p';
                    wt.tokenspritenum=6;
                else
                    wt.won='?';
                    wt.tokenspritenum=6;
                end
            else
                wt.won=[num2str(rdata(trialnum,15)*10) ' p'];
                wt.tokenspritenum=10;
                if rdata(trialnum,wt.respcol)==2
                    wt.won='?';
                    wt.tokenspritenum=6;
                end
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
       
      % Email researcher if almost done
       if trialnum==size(rdata,1)-15
           try
                w.time=clock;
                f_sendemail('learnreward',['[GoalExplore] (' num2str(w.time(4)) ':', num2str(w.time(5)) 'hrs) ' w.subjname ' has 1.5 min left (' w.taskstage ')'], ['Session almost complete: ' w.taskstage])        
           catch
            end
       end
       
       % Verify choices recorded
       disp(['Response OK? ' num2str(rdata(trialnum,13))])
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
           save([where filesep 'Data' filesep w.subjname '_workspace_2taskconflict'])
       catch
           save([w.subjname '_workspace_2taskconflict'])
       end
       
       % Wait for next trial
       cgtext('+',0,0)
       waituntil(wt.outcome*1000+p.disp.outcome)
       
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
p.log.end=clock;
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
filename=[w.subjname '_file_2' w.taskstage '.mat'];
eval(['save(filename, ''' w.taskstage ''')']) % Save file
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
disp(['Score: ' num2str(outcomes.net) ' (' num2str(outcomes.net/(2520/2)) '%  = ' num2str(7*outcomes.net/(2520/2)) ')'])
disp(['# of Omissions: ' num2str(w.omissions)])
disp(['Learning accuracy: ' num2str(outcomes.learningaccuracy)])
disp('(See variable ''w.dat.behaviour'' for behavioural performance)')
disp(w.warning)
disp('-----------------------------------'); diary off