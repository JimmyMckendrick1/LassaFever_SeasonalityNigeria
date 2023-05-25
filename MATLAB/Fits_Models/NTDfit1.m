mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1); 
load(strcat(newdir,'/TestFits/RatModel/RatEpi4.mat'))
Folderdir = strcat(newdir,'/TestFits/RatModel/');
FileName = 'RatEpi4_res_%d';

sv_file = strcat(Folderdir,FileName);

[All_gen_results_cell,Total_iterations,tolVec,seed] = Algorithm_Modelselection_AQ_SMC_P...
    (Models,Prior_models,Priors_given_m,Inputs,data,K,Quantile,T,error_method,NoPart,MaxIter,PKernel,sv_file);
