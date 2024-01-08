function [] = Setup_fitting_sensitivity(i,j)
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1); 
load(strcat(newdir,'/Inputs/Manuscript_input.mat')) % Load an input file
Folderdir = strcat(newdir,'/TestFits/Review/');
FileName = 'Sensitivity_';
sv_file = strcat(Folderdir,FileName);
Recovered_arr = [0.3,0.2,0.1];
data_arr = [1,2,3];


if isequal(j,1) 
    load(strcat(newdir,'/TestFits/Review/Manuscript_results.mat'),'tolVec')
    Input.S_h_0 = (1-Recovered_arr(i))*Input.N_h_0 - Input.E_h_0-Input.A_h_0-Input.I_h_0;
    tolVec = tolVec(1:end-1);
    sv_file = strcat(Folderdir,FileName,string(Recovered_arr(i)*100),"Recov_x",string(j),"data_%d");
    [All_gen_results_cell,Total_iterations,seed] = Algorithm_ABCSMC_P_epi...
    (tolVec,MaxIter,NoPart,model,Input,Priors_given_m,data,error_method,PKernel,sv_file);
else
    load(strcat(newdir,'/TestFits/Review/Manuscript_results.mat'),'tolVec')
    tolVec = tolVec(1:end-1);
    data(1:16,2) = data_arr(j)*data(1:16,2); %2018 data scaled for possible under-reporting
    Input.E_h_0 = data_arr(j)*Input.E_h_0; 
    Input.A_h_0 = data_arr(j)*Input.A_h_0;
    Input.I_h_0 = data_arr(j)*Input.I_h_0;
    Input.C_h_0 = data_arr(j)*Input.C_h_0;

    Input.S_h_0 = (1-Recovered_arr(i))*Input.N_h_0 - Input.E_h_0-Input.A_h_0-Input.I_h_0;
    sv_file = strcat(Folderdir,FileName,string(Recovered_arr(i)*100),"Recov_x",string(j),"data_%d");
    [All_gen_results_cell,Total_iterations,seed] = Algorithm_ABCSMC_P_epi...
    (tolVec,MaxIter,NoPart,model,Input,Priors_given_m,data,error_method,PKernel,sv_file);


end