% Load modelling pars and look at 
clear all, close all hidden, clc

cf = load('3 Hierarchical/res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968');  % [ Current, Manuscript v2-6 ]
ct = load('3 Hierarchical/res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293');

%%
cf_sp = cf.r_res{1,2}(:, 4:end);  % Subject params 
ct_sp = ct.r_res{1,2}(:, 4:end);
cf_set = cf.details.models( find(strcmp(cf.details.models(:,1), cf.r_res{1,1})), :); 
ct_set = ct.details.models( find(strcmp(ct.details.models(:,1), ct.r_res{1,1})), :); 
  
req.sharedpar={'b','p','j','i','w'};

close all hidden
figure('color','w'), f.FontSize=25; f.FontName='PT Sans Caption';    
for p=1:length(req.sharedpar)
    wp.cf= cf_sp(:, find(strcmp( cf_set{3},  req.sharedpar{p}))); 
    wp.ct= ct_sp(:, find(strcmp( ct_set{3},  req.sharedpar{p})));
    wp.d = [wp.cf wp.ct];
    
       
    [h pp ci st] = ttest(wp.cf-wp.ct);
    disp([req.sharedpar{p} ': t=' num2str(st.tstat) ' p='  num2str(pp)])
    
    subplot(1, 5,p)
    barwitherr(std(wp.d)/sqrt(size(cf_sp,1)),  mean(wp.d ))
    title([req.sharedpar{p} ' parameter'])
    set(gca,'FontSize', f.FontSize, 'FontName', f.FontName )
    ylabel('Parameter value')
    
end 

  