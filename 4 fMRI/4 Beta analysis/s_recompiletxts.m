% Compile several text documents into one. Assume 1st column of each text
% document is subject.
clear all; clc

% where.inputs='I:\Conflict Beta analysis\PPI Beta files\Rej-Exp';
where.inputs='/Volumes/PENNYDISK/Conflict Beta analysis/PPI Beta files/cF_Rej-Exp';

log.text_inputs=cellstr(spm_select('List', where.inputs, '.*.txt'));
disp('Order of subjects not checked - assumed the same')

%%

nchar4roi=[16 19 2 21 22 21 22];
%     
% k=20;
% j=26;
% a{j}(length(a{j})-k:end)

dat=[];
k=2;
for t=1:length(log.text_inputs)
    wt.file=importdata([where.inputs filesep log.text_inputs{t}]);
    wt.d=[wt.file.textdata(1,:);    [strtrim(wt.file.textdata(2:end,1)) num2cell(wt.file.data)]];
    disp(['Doc ' num2str(t)])
    
    % Flip titles LR order of roi & contrast
    for c=2:9
        disp(num2str(c))
        wc.name=wt.d{1,c}(1:length(wt.d{1,c})-nchar4roi(t)-1);
        disp(wc.name)
        
        
        
    end
    
    
    %
    
    if t==1; dat=wt.d; 
    else
        dat=[dat wt.d(:,2:end)];
    end
    wt=[];
    
end





%%


input('Continue to print?')
printok=print2txt(where.inputs, ['(' date ') Combined inputs'], dat);
disp(printok)

