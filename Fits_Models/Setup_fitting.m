mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1); 
load(strcat(newdir,'/Inputs/Manuscript_input.mat')) % Load an input file
Folderdir = strcat(newdir,'/Results/'); % Save location
FileName = 'Manuscript_input_res_%d'; % Save name

sv_file = strcat(Folderdir,FileName);

[All_gen_results_cell,Total_iterations,tolVec,seed] = Algorithm_AQ_SMC_P_epi...
    (K,Quantile,T,MaxIter,NoPart,model,Input,Priors_given_m,data,error_method,PKernel,sv_file);
