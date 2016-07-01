% Interleaved Conflict & Control blocks (MRI)
clear all;close all hidden;clc;
nblocks=6;

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
        w.screenmode=1;
    else
        w.subjname='t1';
        w.taskrep=1;
%         w.taskrep=input('# times this stage repeated (1=1st time, 2=2nd, etc): ');
        w.screenmode=0;
    end
    w.taskstage='taskconflictcontrol';
    where=pwd;
    w.res=2;
end
for o1=1:1 % General parameters (Load from learning) 
    load ([where filesep 'Data' filesep w.subjname '_file_1learnenv.mat'])
    p=learnenv.settings;
    % New settings for Test phase
    p.n.reps_pertrialtype=18; % If you change this, CHECK the probabilistic structure!
    p.n.nblocks=6; % per task
    if round(p.n.reps_pertrialtype*2*p.task.envthreat_nlevels*p.task.envthreat_nlevels/p.n.nblocks)~=(p.n.reps_pertrialtype*2*p.task.envthreat_nlevels*p.task.envthreat_nlevels/p.n.nblocks)
        input('Error: Invalid number of blocks. Must be a denominator of total number of trials');
    end
    %
    w.line=50;
    rand('state',sum(100*clock));
end
for o1=1:1 % Subject-specific parameters 
[p w.par Norm check]=f_generate_taskstruc(p,exec);
% Combine conflict/pavlovian & control trials (alternating blocks)
w.ntrials=size(w.par,1);
w.ntrialsblock=w.ntrials/p.n.nblocks;
for i=1:size(w.par,1)  % Mark whether outcome is presented
    if i/2==round(i/2)
        w.par(i,14)=1;
    else
        w.par(i,14)=0;
    end
end
parP=w.par;
parP(:,5)=1;
parP(:,7)=rand(w.ntrials,1);
parP=sortrows(parP,7);
parC=w.par;
parC(:,5)=2;
parC(:,7)=rand(w.ntrials,1);
parC=sortrows(parC,7);
for i=1:p.n.nblocks:w.ntrials-1
    parP(i:i+p.n.nblocks-1, 7)=1:p.n.nblocks;
    parC(i:i+p.n.nblocks-1, 7)=1:p.n.nblocks;
end
par=vertcat(parP,parC);
p.whichtaskfirst=randi(2);
switch p.whichtaskfirst
    case 1
        par=sortrows(par,[7 5]);
    case 2
        par=sortrows(par,[7 -5]);
end
% Mark which side to reveal if explored
par(:,4)=randi(2,size(par,1),1);
rdata=par;
% Misc
w.back=p.disp.settokencolor;
end
for o1=1:1 % Checks for running 
disp('------------------------------------------------------------------')
disp('CREATING PARAMETERS FOR FMRI SESSION (Interleaved Conflict & Control)')
% disp(['No. levels pBombEnv: ' num2str(p.task.envthreat_nlevels)])
% disp(['No. reps per trial type: ' num2str(p.n.reps_pertrialtype)])
disp(['First block task type: ' num2str(rdata(1,5))])
disp(['Exploration cost: ' num2str(p.task.explorationcost*10) ' p'])
disp(['Number of blocks :' num2str(nblocks)])
disp('------------------------------------------------------------------')
% input('Setup OK?  ');
end


%% Format for fMRI
filename=[w.subjname '_file_4fMRIparameters.mat'];
p.counterbal_startask=input('Counterbalancing which task to start?   (1=Conflict first, 2=Control first):   ');
    
switch nblocks
    case 3
        input('ERROR: No. of blocks should be 4!')
%         % 3 Blocks
%         blockpar{1}=par(1:round(size(par,1)/3),:);
%         blockpar{2}=par(round(size(par,1)/3)+1:round(size(par,1)*2/3),:);
%         blockpar{3}=par(round(size(par,1)*2/3)+1:size(par,1),:);
%         if size(blockpar{1},1) + size(blockpar{2},1)+size(blockpar{3},1) ~=size(rdata,1)
%             disp('ERROR: No. of trials in each block does not equal no. of trials overall!')
%         end
%         p.nblocks=3;
%     case 4
        % 4 Blocks
%         blockpar{1}=par(1:round(size(par,1)/4),:);
%         blockpar{2}=par(round(size(par,1)/4)+1:round(size(par,1)*2/4),:);
%         blockpar{3}=par(round(size(par,1)*2/4)+1:round(size(par,1)*3/4),:);
%         blockpar{4}=par(round(size(par,1)*3/4)+1:size(par,1),:);
%         if size(blockpar{1},1) + size(blockpar{2},1)+size(blockpar{3},1)+size(blockpar{4},1) ~=size(rdata,1)
%             disp('ERROR: No. of trials in each block does not equal no. of trials overall!')
%         end
%         p.nblocks=4;
    case 6
        % 6 blocks
        blockpar{1}=par(1:round(size(par,1)/6),:);
        blockpar{2}=par(round(size(par,1)/6)+1:round(size(par,1)*2/6),:);
        blockpar{3}=par(round(size(par,1)*2/6)+1:round(size(par,1)*3/6),:);
        blockpar{4}=par(round(size(par,1)*3/6)+1:round(size(par,1)*4/6),:);
        blockpar{5}=par(round(size(par,1)*4/6)+1:round(size(par,1)*5/6),:);
        blockpar{6}=par(round(size(par,1)*5/6)+1:size(par,1),:);
        if size(blockpar{1},1) + size(blockpar{2},1)+size(blockpar{3},1)+size(blockpar{4},1)+size(blockpar{5},1)+size(blockpar{6},1) ~=size(rdata,1)
            disp('ERROR: No. of trials in each block does not equal no. of trials overall!')
        end
        p.nblocks=6;
end



%% Save & put to Deleted Daily
p.datablock=rdata;
save(filename, 'par','parC','parP','p', 'blockpar');


try % Save file 
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
