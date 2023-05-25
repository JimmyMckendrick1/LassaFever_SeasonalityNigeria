function resp = AlgorithmModule(str_function, cell_args)
% This function is a python-like module, containing sub functions of the
% Approximate Bayesian Computation algorithms that I have for fitting my
% epidemiological models.
allowed_functions = {'ErrorCalculation_Index', 'CustomDiscreteDistribution',...
    'SampleFromPrior','SampleFromPrior_Adaptive','SampleFromLastGen','Error','WeightCalc',...
    'PerturbationKernMatrix','MCMC_alpha','Unpack','SquareDiff','SquareDiffMeta','qt_adaptive','W_hat';
                     @ErrorCalculation_Index, @CustomDiscreteDistribution,...
    @SampleFromPrior,@SampleFromPrior_Adaptive,@SampleFromLastGen,@Error,@WeightCalc,...
    @PerturbationKernMatrix,@MCMC_alpha,@Unpack,@SquareDiff,@SquareDiffMeta,@qt_adaptive,@W_hat};

[found, where] = ismember(str_function, allowed_functions(1, :));  
if ~found
    error('Function %s is not a valid function');
end
resp = allowed_functions{2, where}(cell_args);  
end
%%

function [Model_Index] = CustomDiscreteDistribution(cell_args)
% This function selects an index from a custom discrete uniform in a binary
% tree search. Also just does the 

% Input: cell_args is in this case just an array, with entries increasing 
% to 1 as a custom Discrete Uniform distribution
if strcmp(cell_args{1},"Discrete uniform")
    Model_Index = randi([1,cell_args{2}]);
else
    Model_Index = datasample(1:length(cell_args),1,'Weights',[cell_args{:}]);
end
end

function [Ans] = SampleFromPrior(cell_args)
% This function looks through the inputs of a model and the priors of the
% parameters for that model and then finds the index of the parameters to
% be fitted so that they can be easily switched during the main body of an
% Approximate Bayesian Computation scheme

[Priors_given_m,p_set_dim,input] = cell_args{:};
Prop_prm = table(NaN(1,p_set_dim(1)));
Prop_input = input;
Off_set_dists = {'Gamma','Logistic'};
for i = 1:p_set_dim
    J = p_set_dim(i+1);
    Pp_entry = NaN(1,J);
    if isequal(width(Priors_given_m.(Priors_given_m.Properties.VariableNames{i})),J) && J>1
        for j = 1:J        
            if ismember(Priors_given_m{2,i}{1},Off_set_dists)
                Pp_entry(j) = Priors_given_m{1,i}{1}{end,j}+random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{1:end-1,j});
            else
                Pp_entry(j) = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:});
                while strcmp(Priors_given_m{2,i}{1},'Normal') && Pp_entry(j)< 0
                    Pp_entry(j) = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:,j});
                end        
            end        
        end        
    else
        for j = 1:J        
            if ismember(Priors_given_m{2,i}{1},Off_set_dists)
                Pp_entry(j) = Priors_given_m{1,i}{1}{end}+random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{1:end-1});
            else
                Pp_entry(j) = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:});
                while strcmp(Priors_given_m{2,i}{1},'Normal') && Pp_entry(j)< 0
                    Pp_entry(j) = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:});
                end        
            end        
        end
    end
    Prop_prm.(i)= Pp_entry;
    Prop_input.(Priors_given_m.Properties.VariableNames{i}) = Prop_prm.(i);

end
Ans = {Prop_input,Prop_prm};
end

function [Prop_prm] = SampleFromPrior_Adaptive(cell_args)
% This function looks through the inputs of a model and the priors of the
% parameters for that model. For the adaptive scheme this alternative
% function is needed to only provide the proposed parameters and not a new
% input for the model.

[Priors_given_m,p_set_dim] = cell_args{:};
Prop_prm = NaN(1,p_set_dim);

Off_set_dists = {'Gamma','Logistic'};

for i = 1:p_set_dim
    for j = 1:length(Prop_input.Properties.VariableNames{i})
        if ismember(Priors_given_m{2,i}{1},Off_set_dists)
            Prop_prm(i) = Priors_given_m{1,i}{1}{end}+random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{1:end-1});
        else
            Prop_prm(i) = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:});
            while strcmp(Priors_given_m{2,i}{1},'Normal') && Prop_prm(i) < 0
                Prop_prm(i) = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:});
            end        
        end
    end
end

end

function [Ans] = SampleFromLastGen(cell_args)
% This function will perturb a particle that was chosen by a weighting given
% from the previous generation of particles. The function will
% then change the input args of the model
% Currently has work in progress if different kernel types are added
% other than multivariate Gauss. The check for this is commented out at the moment. 

[Sampled_Part,Priors_given_m,p_set_dim,...
    Input,Perturbation,PKernel_type] = cell_args{:};

Prop_prm = table(NaN(1,p_set_dim(1)));
Prop_input = Input;

Off_set_dists = {'Gamma','Logistic'};

Param_Ind = 0;
count = 0;
if strcmp(PKernel_type,"covarianceMVN") || strcmp(PKernel_type,"MVN")
    while Param_Ind < sum(p_set_dim(2:end))
        % Test whether each sampled parameter from the
        % perturbation of the previous generation is still
        % within the prior distribution, and other acceptable
        % criterion. Probably functionify it too
        Param_Ind = 0;
        Prop_prm_arr = mvnrnd(Sampled_Part.Variables, Perturbation);

        for i = 1:p_set_dim(1)
            J = p_set_dim(i+1);   
            if ismember(Priors_given_m{2,i}{1},Off_set_dists)
                y = pdf(Priors_given_m{2,i}{1},Prop_prm_arr(Param_Ind+1:Param_Ind+J)-Priors_given_m{1,i}{1}{end},Priors_given_m{1,i}{1}{1:end-1});
            else
                y = pdf(Priors_given_m{2,i}{1},Prop_prm_arr(Param_Ind+1:Param_Ind+J),Priors_given_m{1,i}{1}{:});
            end
            Pp_entry = Prop_prm_arr(Param_Ind+1:Param_Ind+J);
            if strcmp(Priors_given_m{2,i}{1},'Normal') && any(Pp_entry<=0)
                break
            elseif any(y<=0)
                break
            end
            Param_Ind= Param_Ind+J;
            Prop_prm.(i)= Pp_entry;
            Prop_input.(Priors_given_m.Properties.VariableNames{i}) = Prop_prm.(i);
        end
        count = count+1;
    end
elseif strcmp(PKernel_type,"Piecewise")
    % Test whether each sampled parameter from the
    % perturbation of the previous generation is still
    % within the prior distribution, and other acceptable
    % criterion. Probably functionify it too
    p_i = 1;
    while p_i < p_set_dim(1)+1
        if strcmp(Perturbation{2,p_i},'Discrete Uniform')
            Prop_prm.(p_i) = Sampled_Part.(p_i) + random(Perturbation{2,p_i}, Perturbation{1,p_i})-floor(Perturbation{1,p_i}/2);
        else
            Prop_prm.(p_i) = random(Perturbation{2,p_i}, Sampled_Part.(p_i), Perturbation{1,p_i});
        end
        if ismember(Priors_given_m{2,p_i}{1},Off_set_dists)
            y = pdf(Priors_given_m{2,p_i}{1},Prop_prm.(p_i)-Priors_given_m{1,p_i}{1}{end},Priors_given_m{1,p_i}{1}{1:end-1});
        else
            y = pdf(Priors_given_m{2,p_i}{1},Prop_prm.(p_i),Priors_given_m{1,p_i}{1}{:});
        end
        if ~any(y<=0)          
            Prop_input.(Priors_given_m.Properties.VariableNames{p_i}) = Prop_prm.(p_i);
            p_i=p_i+1;
        end
    count = count+1;
    end
end
Ans = {Prop_input,Prop_prm};
end

function [mass] = WeightCalc(cell_args)
% Function that calculates the weight of a particle at generation t>1 from
%
% Assumes that the priors of each parameter are independent. 

[Priors,Prev_wghts,Prop_prm,Prev_gen,PKernel_type,Kernel_sig] = cell_args{:};
                
Allowed_Kerns = {'covarianceMVN','MVN','Piecewise';
                     @MVNKernel,@MVNKernel,@Piecewise};
mass = 1;
for i = 1:width(Prop_prm)
    for j = 1:length(Prop_prm{1,i})
        Param_set=Prop_prm{1,i};
        mass = mass*pdf(Priors{2,i}{1},Param_set(j),Priors{1,i}{1}{:});
    end
end

[found2, where2] = ismember(PKernel_type, Allowed_Kerns(1, :));  %search 1st row of lookup table
if ~found2
    error('Function %s is not a valid function');
end
mass = mass/Allowed_Kerns{2,where2}({Prev_wghts,Prop_prm{1,:},Prev_gen,Kernel_sig});
end

function [Kmatrix] = PerturbationKernMatrix(cell_args)
% Calculates a matrix that will be used for the perturbation kernel in the
% SampleFromLastGen function to draw the next particles.
% Currently only feature for covariance calculations, but may need
% different kernels in the future
PKernel_type = cell_args{1}{1};

switch PKernel_type
    case 'covarianceMVN'
        Kmatrix = 2*cov(cell_args{2}.Variables);
    case 'MVN'
        Kmatrix = cell_args{1}{2};
    case 'Piecewise'
        Kmatrix = cell(2,width(cell_args{2}));
        for i = 1:width(cell_args{2})
            if strcmp('Discrete Uniform',cell_args{1}{2}{i})   
                Kmatrix{1,i} = ceil(sqrt(cov(cell_args{2}{:,i})));
            else                
                Kmatrix{1,i} = cov(cell_args{2}{:,i});
            end
            Kmatrix{2,i} = cell_args{1}{2}{i};
        end
end
end

function [alpha] = MCMC_alpha(cell_args)
% Function that calculates the probability of accepting a particle for the MCMC.
% Assumes that the priors of each parameter are independent. 
[Priors,Prop_prm,Prev_part,PKernel_type,Kernel_sig] = cell_args{:};
            
Allowed_Kerns = {'covarianceMVN','MVN','Piecewise';
                     @MVNKernel,@MVNKernel,@Piecewise};
alpha = 1;
for i = 1:width(Prop_prm)
    alpha = alpha*pdf(Priors{2,i}{1},Prop_prm{1,i},Priors{1,i}{1}{:})...
        /pdf(Priors{2,i}{1},Prev_part{1,i},Priors{1,i}{1}{:});
end

[found2, where2] = ismember(PKernel_type, Allowed_Kerns(1, :));  %search 1st row of lookup table
if ~found2
    error('Function %s is not a valid function');
end
alpha = alpha*Allowed_Kerns{2,where2}({[1],Prev_part{1,:},Prop_prm{1,:},Kernel_sig})...
    /Allowed_Kerns{2,where2}({[1],Prop_prm{1,:},Prev_part{1,:},Kernel_sig});
end

function [All_gen_results_cell] = Unpack(cell_args)
% Unpack results in uncatagorised form and put into a nice, labelled table
if length(cell_args) == 3
    All_gen_results_cell = cell_args{1};
    saved_results = cell_args{2};
    p_set_dim = cell_args{3};
    
    for i_param = 1:p_set_dim(1)
        Ind_svd_r = sum(p_set_dim(2:i_param+1))-p_set_dim(i_param+1); 
        % Find the index that the parameters of this name are in saved_results
        for  k_sub_param = 1:p_set_dim(i_param+1)
            All_gen_results_cell{:,i_param}(:,k_sub_param) = saved_results(:,Ind_svd_r+k_sub_param);
            % For each subparameter of this name, save into the results cell
        end
    end
    All_gen_results_cell{:,end} = saved_results(:,end); % Errors
else
    All_gen_results_cell = cell_args{1};
    saved_results = cell_args{2};
    p_set_dim = cell_args{3};
    No_models = cell_args{4};
    model_vec = cell_args{5};
    for i_model = 1:No_models
        % Unpack the information, due to parallelisation cannot slice arrays
        % and cells and so we need to do it in a slightly clunky way
        row_inds = model_vec == i_model; % Find which parameters belonged to this model
        All_gen_results_cell{i_model} = All_gen_results_cell{i_model}(row_inds,:);
        for i_param = 1:p_set_dim(i_model,1)
            Ind_svd_r = sum(p_set_dim(i_model,2:i_param+1))-p_set_dim(i_model,i_param+1); 
            % Find the index that the parameters of this name are in saved_results
            for  k_sub_param = 1:p_set_dim(i_model,i_param+1)
                All_gen_results_cell{i_model}{:,i_param}(:,k_sub_param) = saved_results(row_inds,Ind_svd_r+k_sub_param);
                % For each subparameter of this name, save into the results cell
            end
        end
        All_gen_results_cell{i_model}{:,end} = saved_results(row_inds,end); % Errors
    end
    
end

end
    %% allowed error calculators

function [err] = SquareDiff(cell_args)
% subfunction of Error. If the method relies on a square difference method,
% this directly calculates it
[Dates,data,Output_T,data_names,Bool_data] = cell_args{:};
w = length(data_names);
err = zeros(1,w);

if Dates(end)+1>length(Output_T.t)
    err = inf;
else
    for i = 1:w 
        err(i) = err(i)+sum((data(:,i)-Output_T.(data_names{i})(Dates(:)+1)).^2);
        err(i) = sqrt(err(i));
    end
    err = sum(err);
end

end

function [err] = SquareDiffMeta(cell_args)
% subfunction of Error. If the method relies on a square difference method,
% this directly calculates it
[Dates,data,Output_T,data_names,Bool_data] = cell_args{:};
w = length(data_names);
err = zeros(1,w);

if Dates(end)+1>length(Output_T.t)
    err = inf;
else
    for i = 1:w 
        err(i) = err(i)+sum(sum(((data(:,i:w:end)-Output_T.(data_names{i})(Dates(:)+1,:).*Bool_data)).^2));
        err(i) = sqrt(err(i));
    end
    err = sum(err);
end

end
%% Allowed weight calculators

function [weightdenom] = MVNKernel(cell_args)
% Function that calculates the denominator of the weighting of a particle
% if the kernel of the smc perturbation step is a Multi-Variate Normal
% distribution
[Prev_wghts,Prop_prm,Prev_gen,Kernel_sig] = cell_args{:};

weightdenom = sum(Prev_wghts.*mvnpdf(Prev_gen,Prop_prm,Kernel_sig));
end

function [weightdenom] = Piecewise(cell_args)

[Prev_wghts,Prop_prm,Prev_gen,Kernel_sig] = cell_args{:};
weightv = NaN(length(Prev_wghts),length(Prop_prm));
for i = 1:length(Prop_prm)
    if strcmp(Kernel_sig{2,i},'Discrete Uniform')
        weightv(:,i) = Prev_wghts.*(pdf('Discrete Uniform',Prev_gen(:,i)-Prop_prm(i)+floor(Kernel_sig{1,i}/2),Kernel_sig{1,i}));
    else        
        weightv(:,i) = Prev_wghts.*mvnpdf(Prev_gen(:,i),Prop_prm(i),Kernel_sig{1,i});
    end
end
weightdenom = sum(sum(weightv));
end

function [q_g] = qt_adaptive(cell_args)
% Function that is as faithful as I can be to Sugiyama et al's (2007) paper
% on an KLIEP for an an adaptive ABC SMC alg.
[NoPart,b,Curr_gen,Prev_gen,sigma,p_set_dim,Priors_given_m] = cell_args{:};

% As per the paper (page 4 end of section 2), we take a sample of the test particles (current generation)
% to be the means of the kernels which make up the set of basis function for the ratio/importance w.
C_template_sample_ind = randperm(NoPart);
C_template_sample_ind = C_template_sample_ind(1:b); 
C_template_sample = Curr_gen(C_template_sample_ind,:); 

% Calculate the matrix A in step 1 of KLIEP main code alg
n_te = length(Curr_gen);
A = NaN(n_te,b);
for i = 1:b
    for j = 1:n_te
        A(j,i) = exp(-norm(Curr_gen(j,:)-C_template_sample(i,:))^2/(2*sigma^2));
    end
end
  
% Calculate B_vec in step 2
n_tr = length(Prev_gen);
B_Matrix = NaN(n_tr,b);
for i = 1:b
    for j = 1:n_tr
        B_Matrix(j,i) = exp(-norm(Prev_gen(j,:)-C_template_sample(i,:))^2/(2*sigma^2));
    end
end
B_vec = transpose(sum(B_Matrix))/n_tr;

% Initialise for main loop
alpha_vec = ones(b,1);
alpha3 = zeros(b,1);
Tries = 0;
epsilon = rand(1)/100;
while Tries < 100 % Optimisation problem for the alpha_vec
    alpha1 = alpha_vec+epsilon*transpose(A)*(ones(NoPart,1)./(A*alpha_vec));
    alpha2 = alpha1+(1-transpose(B_vec)*alpha1)*B_vec/(transpose(B_vec)*B_vec);
    for k = 1:b
        alpha3(k) = max(0,alpha2(k));
    end
    alpha4 = alpha3/(transpose(B_vec)*alpha3);
    alpha_vec = alpha4;
    Tries = Tries + 1;
end

ii_max = max(400,min(20^p_set_dim,1.6e4));
W_hat_vec = zeros(1,ii_max); %Find maximum of the ratio
for ii = 1:ii_max % With a lack of imagination, I have used sampling to find the max
    Prop_prm = AlgorithmModule('SampleFromPrior_Adaptive',{Priors_given_m,p_set_dim});
    W_hat_vec(ii) = AlgorithmModule('W_hat',{alpha_vec,C_template_sample,Prop_prm,sigma});
end
c_g = max(W_hat_vec);
q_g = 1/c_g;
end

function [wx] = W_hat(cell_args)
[alpha_vec,C_sample,x,sigma] = cell_args{:};
b = length(alpha_vec);
K_sig_vec = zeros(b,1);
for i = 1:b
    K_sig_vec(i) = exp(-norm(x-C_sample(i,:))^2/(2*sigma^2));
end
wx = transpose(K_sig_vec)*alpha_vec;
end