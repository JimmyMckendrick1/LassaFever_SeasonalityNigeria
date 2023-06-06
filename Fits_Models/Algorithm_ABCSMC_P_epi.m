function [All_gen_results_cell,Total_iterations,seed] = Algorithm_ABCSMC_P_epi...
    (tolVec,MaxIter,M,model,Input,Priors_given_m,data,error_method,PKernel,save_name)
%%%%%%%%%%%%%%%%%%%%%%%%McKendrick et al%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximate Bayesian Computation SMC parallelised scheme for epidemiological models.
% This algorithm applies my adapted SMC for a single model.
% This is the SMC version from Toni et al 2009 where the posterior converges interatively with
% a sequence of tolerances. It produces posterior distributions.
% This function takes the inputs of number of particles, maximum number of
% iterations, the inputs for the model, and of them which parameters are to be fitted,
% along with the data and error calculaton method.
% This will essentially use the rejection scheme for the first generation
% and then subsequent generations will be created by perturbing the
% previous with which ever method you choose. This will then tell you the 
% rng seed used, the total number of parameter samples and the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Priors_given_m: The priors of the tested parameters for each model (Table)
% Input: The inputs for each of the functions, needs to be in the same order for 
% the functions in Models (Cell, table entries which have their first entry as 
% the inputs of the distribution, second entry is the distribution name)
% Data: The data, with time in days (Array)
% tolVec: The decreasing sequence of tolerances used for each generation (Array)
% M: Total number of particles generated (int64)
% model: The function handle for the model being fitted (Function handle)
% MaxIter: Maximum number of parameters sampled before running too long (int64)
% Pkernel: How the perturbation kernel works. A cell of the method and whatever
% is needed, such as a matrix variance for a multivariate Gauss dist. (Cell,string) 
% save_name: String path of where the results will be saved. This includes
% folder path and the actual file name. The algorithm will save iteratively
% after each generation and don't include .mat as this done automatically (string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All_gen_results_cell: A cell of all generations results. In order of gen
% seed: The rng used in this run
% Total_iterations: The total number of iterations
% tolVec: tolerance vector that was computed during the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For convenience, running in the file this will have these dummy parameters.
if nargin == 0
    M = 200; MaxIter = inf;
    model = @Model_SEAIR_const;  sigma=1;
    beta = [1]; nu = [0.1]; eps = [1]; p = 0.8; gamma = [0.2]; B = [1/(53.5*365.0)]; mu  = [0.016]; muI = [0.02];
    N = [1000]; S0 = [950]; E0 = [2]; I0 = [1]; A0 = [1]; C0  = [1]; D0 = [0]; MaxTime = [80]; 
    Input = table(beta,sigma,nu,eps,p,gamma,B,mu,muI,N,S0,E0,A0,I0,C0,D0,MaxTime);
    beta = {{0;4};'Uniform'};
    eps = {{0;10};'Uniform'};
    gamma = {{0.3;0.1};'Normal'};
    Priors_given_m = table(beta,eps,gamma);
    PKernel = {'covarianceMVN'};
    data = [0:1:15;1:2:31]';
    error_method = {'SquareDiff',{'C'}};
    tolVec = [1e4,1e3,1e2];
    save_name = 'Test_%g';
end
tic;
seed = rng; % Save the rng used so that the results can be repeated if necessary
T = length(tolVec);
All_gen_results_cell = cell(1,T); % Preset the results cell
Dates = data(:,1); % The dates of each data point
Bool_data = data(:,2:end) ~= 0;
Total_iterations = zeros(1,T); % Number of parameter samples drawn in a generation
Num_iterations = zeros(1,M); % Number of parameter samples drawn in a generation
PKernel_matrix = cell(1,T); % In this current formulation, the covariance of 
% the previous generation of particles is used to as the perturbation to sample and
% produce the next generation of parameter particles. No other methods are
% supported currently

Max_params = width(Priors_given_m); % Find the maximum number of parameters tested for a model
p_set_dim = zeros(1,Max_params+1); 
p_set_dim(1) = width(Priors_given_m);
% The number of parameters per named parameter that will be tested for each model

% Set up the weightings which will be used to sample particles from the
% previous generation for that model
weights = cell(1,T);
for i = 1:T
    weights{1,i} = NaN(1,M);
end 

% Get the name of the parameters to be tested to label the results
Template_gen_results_T = table();
Names_model_i = cell(1,p_set_dim(1)+1);
for Name_check_vars = 1:p_set_dim(1)
    var_name = Priors_given_m.Properties.VariableNames{Name_check_vars};
    Names_model_i{Name_check_vars} = var_name;
    
    J = length(Input.(Priors_given_m.Properties.VariableNames{Name_check_vars}));
    Template_gen_results_T{:,Name_check_vars} = NaN(M,J);
    p_set_dim(Name_check_vars+1) = J;
end
Names_model_i{end} = 'error';

% Set up the All_gen_results_cell entries to be tables of appropriately
% labelled columns for the results of each sample, repeated for each of the T generations.
% Also, record the number of subparameters there are under that name in
% the inputs of the model

Template_gen_results_T{:,p_set_dim(1)+1} = NaN(M,1);
Template_gen_results_T.Properties.VariableNames = Names_model_i;
for i=1:T
    % Every results table for that model is identicle
    All_gen_results_cell{1,i} = Template_gen_results_T;
end 
True_max_pdim = sum(p_set_dim(1,2:end),2); % the number of params once we don't segregate by name
saved_results = NaN(M,True_max_pdim+1); % Array for sampled parameters of all models, 
% these are not grouped by name as they are in the results cell
y_vec = NaN(1,M); % indexes for reordering parameters by error

clear Curr_gen_results Name_check_vars Names_model_i True_max_pdim J
% Main iteration of the ABC algorithm
g = 1;
%%% First generation
parfor i_particle=1:M
    ParamInd = 0; % Indicator as to whether we've found a reasonable parameter set
    % While we haven't gone over the maximum number of iterations 
    % or found an acceptable parameter set do:
    while Num_iterations(i_particle) < MaxIter && ParamInd == 0 
        % Sample from the prior distributions of each parameter being tested and
        % add it to the input vector.
        Sampling_cell_args = {Priors_given_m, p_set_dim, Input};
        Input_and_param = AlgorithmModule('SampleFromPrior',Sampling_cell_args);
        Prop_prm = Input_and_param{2}; 

        % Calculate the error between the data and the proposed model, at
        % the times that we have data for.
        Output_T = model(table2cell(Input_and_param{1}));
        Err_cell_args = {Dates, data(:,2:end), Output_T, error_method{2},Bool_data};
        err =  AlgorithmModule(error_method{1},Err_cell_args);    
        if err < tolVec(g)
            saved_results(i_particle,:) = [Prop_prm.Variables err];  
            % Save parameters, leaving gaps if we the model needs fewer params the
            % the most needed
            ParamInd = 1;
        end
        Num_iterations(i_particle) = Num_iterations(i_particle)+1;
    end
end
% Organise in order of error 
[~,idy] = sort(saved_results(:,end));
saved_results = saved_results(idy,:);
saved_results = saved_results(1:M,:);
Total_iterations(1) = sum(Num_iterations);
Num_iterations = zeros(1,M);
% Unpack the information, due to parallelisation cannot slice arrays
% and cells and so we need to do it in a slightly clunky way
cell_args = {All_gen_results_cell{g}, saved_results,p_set_dim};
All_gen_results_cell{g} = AlgorithmModule('Unpack',cell_args);
% Calculate the weight of the particle for future generation sampling
weights{g}(:) = 1;
% Normalise the weight
weights{g}=weights{g}/norm(weights{g},1);
Kernel_cell_args = {PKernel,All_gen_results_cell{g}(:,1:end-1)};

PKernel_matrix{1} = AlgorithmModule("PerturbationKernMatrix",Kernel_cell_args);

save_file_g = strcat(sprintf(save_name, g),'.mat');
save(save_file_g,'All_gen_results_cell','Total_iterations','tolVec','seed','M','Priors_given_m','Input','PKernel','data','error_method','model','MaxIter','weights')

disp(g)
disp(toc)
% Set the next tolerance entry by the quantile of errors in the generation
% current generation
g = 2;
while g <=T
    % While we haven't run through T generations nor has the tolerance
    % between each generation become too small, we run this main block
    Num_iterations = zeros(1, M);
    parfor i_particle = 1:M
        ParamInd = 0; % Indicator as to whether we've found a reasonable parameter set
        % While we haven't gone over the maximum number of iterations 
        % or found an acceptable parameter set do:
        while Num_iterations(i_particle) < MaxIter && ParamInd == 0 
            % Indicator for where to save the parameters for this particular model and 
            % generation in our results cell, amongst other things
            Prev_Ind = g-1; % Find where you were previously
            % For subsequent generations sample from our previous 
            % generation of particles with the weightings calculated 

            y_vec(i_particle) = datasample(1:length(weights{Prev_Ind}),...
                1,'Weights', weights{1,Prev_Ind}); 

            Sampling_cell_args = cell(1,6);
            %sampled particle
            Sampling_cell_args{1} = All_gen_results_cell{Prev_Ind}(y_vec(i_particle),1:end-1);
            Sampling_cell_args{2} = Priors_given_m; %Priors of this model
            Sampling_cell_args{3} = p_set_dim; %Total number of parameters tested
            Sampling_cell_args{4} = Input; %The input for the model
            Sampling_cell_args{5} = PKernel_matrix{Prev_Ind};
            Sampling_cell_args{6} = PKernel{1}; %Information on the peturbation kernel
            Input_and_param = AlgorithmModule('SampleFromLastGen',Sampling_cell_args);
            Prop_prm = Input_and_param{2};

            % Once an acceptable parameter set is found simulate the data

            Output_T = model(table2cell(Input_and_param{1}));
            Err_cell_args = {Dates,data(:,2:end),Output_T,error_method{1,2},Bool_data};
            err =  AlgorithmModule(error_method{1},Err_cell_args); 

            if err<=tolVec(g) && isreal(err) %If error is acceptable, record it
                saved_results(i_particle,:) = [Prop_prm.Variables err];
                ParamInd = 1;
            else %Otherwise resample a parameter set
                ParamInd = 0;
            end
         
            Num_iterations(i_particle) = Num_iterations(i_particle) + 1;
        end
    end
    [~,idy] = sort(saved_results(:,end));
    saved_results = saved_results(idy,:);
    Total_iterations(g) = sum(Num_iterations); 

        % Unpack the information
    cell_args = {All_gen_results_cell{g}, saved_results,p_set_dim};
    All_gen_results_cell{g} = AlgorithmModule('Unpack',cell_args);    
    
    % Calculate the weight of the particle for future generation
    if g < T 
        Prev_Ind = g-1;
        wght_cell_args = cell(1,6);
        wght_cell_args{1} = Priors_given_m;
        wght_cell_args{2} = transpose(weights{Prev_Ind});
        wght_cell_args{4} = All_gen_results_cell{Prev_Ind}{:,1:end-1};
        wght_cell_args{5} = PKernel{1};
        wght_cell_args{6} = PKernel_matrix{Prev_Ind};
        for i_particle = 1:height(All_gen_results_cell{g})
            wght_cell_args{3} = All_gen_results_cell{g}(i_particle,1:end-1);
            weights{g}(i_particle) = AlgorithmModule('WeightCalc',wght_cell_args);  
        end
    end

    weights{g} = rmmissing(weights{g});        
    weights{g} = weights{g}/norm(weights{g},1);

    % Remove NaN values from the tables in the results cell for this
    % generation, which exist because I have no way of knowing how many
    % particles will be accepted for each model and I do not want to
    % increase the size of the tables if we overfill.

    Kernel_cell_args = {PKernel,All_gen_results_cell{g}(:,1:end-1)};

    PKernel_matrix{g} = AlgorithmModule("PerturbationKernMatrix",Kernel_cell_args);

    save_file_g = strcat(sprintf(save_name, g),'.mat');
    save(save_file_g,'All_gen_results_cell','Total_iterations','tolVec','PKernel','weights')
    g = g+1;
    disp(toc)
end

if isequal(i_particle,MaxIter)
    disp('Generated the maximum number of parameter sets, may not have generated sufficient acceptable sets')
end
end