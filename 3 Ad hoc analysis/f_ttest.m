function [ h, p ] = f_ttest( d)
% T-test, but return a 0.5 if there is a trend  

[h, p]= ttest(d); 
if p<0.1 & p>0.05, h=0.5; end

end

