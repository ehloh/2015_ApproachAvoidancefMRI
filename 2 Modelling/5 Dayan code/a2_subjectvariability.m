% Subject variability? Whats the top model for all subjects?
clear all; clc

task=1;

for o=1:1
    where.res='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\Simulation\1 Paramfit generative models';
%     where.res='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/Simulation/1 Paramfit generative models';
    if task==1; rr=load([where.res filesep 'GenerativeFit cF (09-Dec-2014).mat']);
    else rr=load([where.res filesep 'GenerativeFit ct (09-Dec-2014).mat']);end
    rr.r_res=sortrows(rr.r_res,1); rr.mods=rr.r_res(:,1); rr.subjects=rr.details.subjects;
    
    
%    where.where='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
% %     where.where='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
%     where.inputs=[where.where filesep '2 Analysis inputs'];    
%     path(pathdef); addpath(where.where), addpath([where.where filesep '1 Value functions Det Jpower']);
%     addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bjm']);
%     addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bpm']);
%     addpath([where.where filesep '1 Value functions Det Jpower' filesep 'bm']);
    
end

%%

% Assemble master matrix: model #, subject, BIC, L
r_mat=[];
for m=1:size(rr.r_res,1)
    wm.row=(m-1)*length(rr.subjects)+1:m*length(rr.subjects);
    r_mat(wm.row,1)=m;
    r_mat(wm.row,2)=1:length(rr.subjects);
    r_mat(wm.row,3)=rr.r_res{m,2}(:,1);
    r_mat(wm.row,4)=rr.r_res{m,2}(:,2);
end

% Best model on a subjectwise level
r_modbictally=[rr.mods num2cell(zeros(length(rr.mods),1))];
for m=1:size(r_modbictally,1)
    r_modbictally{m,3}=length(r_modbictally{m,1})-3;  if isempty(strfind(r_modbictally{m,1},'w'))==0; r_modbictally{m,3}=r_modbictally{m,3}-1; end
end
r_modLtally=r_modbictally;
r_partally=[{'b';'p';'j';'m'; 'i';'f';'e';'ow';'uw';'vw'} num2cell(zeros(10,1))];
for s=1:length(rr.subjects);
    ws.res=r_mat(r_mat(:,2)==s, :);
    r_modbictally{ws.res(ws.res(:,3)==min(ws.res(:,3)), 1),2}=r_modbictally{ws.res(ws.res(:,3)==min(ws.res(:,3)), 1),2}+1;
    r_modLtally{ws.res(ws.res(:,4)==min(ws.res(:,4)), 1),2}=r_modLtally{ws.res(ws.res(:,4)==min(ws.res(:,4)), 1),2}+1;
    
    ws.mod=r_modbictally{ws.res(ws.res(:,3)==min(ws.res(:,3)), 1),1};  % Parameters in best-BIC model?
    for p=1:size(r_partally,1)
        if isempty(strfind(ws.mod,r_partally{p,1}))==0, r_partally{p,2}=r_partally{p,2}+1; end
    end
end

% r_modLtally
% r_modbictally
% r_partally
openvar r_partally
openvar r_modbictally
% mean([cell2mat(r_modbictally(1:20,3))])
r_modbictally=sortrows(r_modbictally,-2);



