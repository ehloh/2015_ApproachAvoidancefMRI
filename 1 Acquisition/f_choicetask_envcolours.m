function [c]=f_choicetask_envcolours(p, numreps_perpair)
% [c]=f_choicetask_envcolours(p, numreps_perpair)
%
%------------------------------------------------------

for o1=1:1 % Documentation for choice task
%
% CHOICE TASK ================
%     Col 1: Trial number
%     Col 2: Left Colour #
%     Col 3: Right Colour #
%     Col 4: Better colour (correct choice - 1=Left, 2=Right)
%     Col 5: Choice (1=Left, 2=Right)
%     Col 6: Accuracy
%    
end

w.line=50;
w.back=p.disp.settokencolor;
c.numreps_perpair=numreps_perpair;

% Instructions
cgflip(w.back)
cgtext('Now, you have a quick choice task to complete',0,w.line*4)
cgtext('The computer will show you two different colours',0,w.line*2)
cgtext('at a time. You will have seen these colours before',0,w.line*1)
cgtext('in the task',0,w.line*-0)
for i=1:p.task.envthreat_nlevels
p.disp.contingencycolours{i,2}=i;
p.disp.contingencycolours{i,3}=rand;
end
p.disp.contingencycolours=sortrows(p.disp.contingencycolours,3);
for i=1:p.task.envthreat_nlevels
    w.col=p.disp.contingencycolours{i};
    cgrect(-305+(i-1)*125,w.line*-2,120, 80, w.col)
end
cgflip(w.back);
waitkeydown(inf);
cgtext('For each pair of colours, you have to decide which',0,w.line*5)
cgtext('of the two colours you prefer - which is less likely',0,w.line*4)
cgtext('to indicate a bomb',0,w.line*3)
switch p.counterbal_hands
    case 1
        cgtext('Press the LEFT arrow if you prefer the colour on',0,w.line*-3)
        cgtext('the left, or the RIGHT arrow ',0,w.line*-4)
        cgtext('if you prefer the colour on the right',0,w.line*-5)
    case 2
        cgtext('Press the Z key if you prefer the colour on',0,w.line*-3)
        cgtext('the left, or the C key ',0,w.line*-4)
        cgtext('if you prefer the colour on the right',0,w.line*-5)
end
for i=1:2
    w.col=p.disp.contingencycolours{i};
    cgrect(-185+(i-1)*355,w.line*0,240, 170, w.col)
end

% Set up
p.disp.contingencycolours=sortrows(p.disp.contingencycolours,2);
c.colours=sortrows(p.disp.contingencycolours,2);
wc.choicecount=1;
for j=1:size(p.disp.contingencycolours,1)
    for i=1:size(p.disp.contingencycolours,1)
        if i>j
            if randi(2)==1
                c.data(wc.choicecount,2)=c.colours{i,2};
                c.data(wc.choicecount,3)=c.colours{j,2};
                c.data(wc.choicecount,4)=2;
            else
                c.data(wc.choicecount,2)=c.colours{j,2};
                c.data(wc.choicecount,3)=c.colours{i,2};
                c.data(wc.choicecount,4)=1;
            end
            wc.choicecount=wc.choicecount+1;
        end
    end
end
c.dat=c.data;
for i=1:c.numreps_perpair-1 % As many reps as requested
    c.data=vertcat(c.dat,c.data);
end
c.data(:,1)=rand(size(c.data,1),1);
c.data=sortrows(c.data,1);


cgflip(w.back);
waitkeydown(inf);
cgtext(['You will do ' num2str(size(c.data,1)) ' trials of this task'],0,w.line*150)
cgtext('Press any key to start',0,w.line*-100)
cgflip(w.back);
waitkeydown(inf);

% Execute
for i=1:size(c.data,1)
    c.data(i,1)=i;
    cgtext(['Trial number: ' num2str(i)],0,200)
    cgtext('Which of these two colours indicates',0,-100)
    cgtext('a LOWER probability of a bomb being present?',0,-150)
    cgrect(-200,50,240,170,c.colours{c.data(i,2),1})
    cgrect(200,50,240,170,c.colours{c.data(i,3),1})
    cgflip(w.back);
    wm.keyright=0;
    while wm.keyright==0
        clearkeys
        [wm.key wm.keytime wm.a]=waitkeydown(inf,[p.buttons.key1 p.buttons.key3]);
        if size(wm.key)==1
            wm.keyright=1;
        end
    end
        switch wm.key(1)
            case p.buttons.key1
                c.data(i,5)=1;
            case p.buttons.key3
                c.data(i,5)=2;
        end
        if c.data(i,5)==c.data(i,4)
            c.data(i,6)=1;
        else
            c.data(i,6)=0;
        end
        waituntil(500)
        disp([num2str(i) ' --- Correct? ' num2str(c.data(i,6))])
        wm=[];
end

end