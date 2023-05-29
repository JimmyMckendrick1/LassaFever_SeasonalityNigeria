function resp = AlgorithmModule(str_function, cell_args)
% This function is a python-like module, containing sub functions of the
% Approximate Bayesian Computation algorithms that we have for fitting
% epidemiological models.
allowed_functions = {'SampleFromPrior','SampleFromLastGen','Error','WeightCalc',...
    'PerturbationKernMatrix','Unpack','SquareDiff','qt_adaptive';
                     @SampleFromPrior,@SampleFromLastGen,@Error,@WeightCalc,...
    @PerturbationKernMatrix,@Unpack,@SquareDiff,@qt_adaptive};

[found, where] = ismember(str_function, allowed_functions(1, :));  
if ~found
    error('Function %s is not a valid function');
end
resp = allowed_functions{2, where}(cell_args);  
end
%%

function [Ans] = SampleFromPrior(cell_args)
% This function produces a cell containing a particle (a realisation of the set of
% parameters which are being fitted) in both a full input for the model and a results table.
% This is done by iterating through the piecewise independent priors inputted
% into the ABC algorithm, drawing a value for each parameter.

[Priors_given_m,p_set_dim,input] = cell_args{:};
Prop_prm = table(NaN(1,p_set_dim(1)));
Prop_input = input;
Off_set_dists = {'Gamma','Logistic'}; % Can add an offset to these distributions so that they don't start at 0

% For each parameter in the prior distribution, draw from the distribution 
for i = 1:p_set_dim
    Pp_entry = NaN;       
    if ismember(Priors_given_m{2,i}{1},Off_set_dists)
        % For distributions that may be offset:
        Pp_entry = Priors_given_m{1,i}{1}{end}+random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{1:end-1});
    else
        Pp_entry = random(Priors_given_m{2,i}{1},Priors_given_m{1,i}{1}{:});       
    end
    Prop_prm.(i)= Pp_entry; % Shortened results output
    Prop_input.(Priors_given_m.Properties.VariableNames{i}) = Prop_prm.(i); % Full input for model

end
Ans = {Prop_input,Prop_prm};
end

function [Ans] = SampleFromLastGen(cell_args)
% This function will return a perturbed particle that was chosen by the weighting given
% from the previous generation of particles. This can be done by one several
% ways. The first is a multi-variate Gaussian centered at the previous
% "Sampled_Part", perturbed either with the covariance of the previous
% generation or by a manually inputted "MVN" variance matrix. The second
% choice is by a piece-wise chosen kernel. Each individual parameter is
% perturbed with its own kernel, which is manually inputted.

[Sampled_Part,Priors_given_m,p_set_dim,...
    Input,Perturbation,PKernel_type] = cell_args{:};

Prop_prm = table(NaN(1,p_set_dim(1)));
Prop_input = Input;

Off_set_dists = {'Gamma','Logistic'};
SampleFromPrior_Adaptive
Param_Ind = 0;
count = 0;
if strcmp(PKernel_type,"covarianceMVN") || strcmp(PKernel_type,"MVN")
    % For multi-variate normal kernels, do:
    while Param_Ind < sum(p_set_dim(2:end))
        % Samples a new particle by perturbing the sampled particle from
        % the previous generation with a multivariate normal 
        Param_Ind = 0;
        Prop_prm_arr = mvnrnd(Sampled_Part.Variables, Perturbation);

        for i = 1:p_set_dim(1)
            % For each parameter in the particle, check if pdf>0 for the
            % prior of that particular parameter
            if ismember(Priors_given_m{2,i}{1},Off_set_dists)
                y = pdf(Priors_given_m{2,i}{1},Prop_prm_arr(Param_Ind+1)-Priors_given_m{1,i}{1}{end},Priors_given_m{1,i}{1}{1:end-1});
            else
                y = pdf(Priors_given_m{2,i}{1},Prop_prm_arr(Param_Ind+1),Priors_given_m{1,i}{1}{:});
            end
            Pp_entry = Prop_prm_arr(Param_Ind+1);
            % If the pdf <=0 then break the for loop and resample.
            if any(y<=0)
                break
            end
            Param_Ind= Param_Ind+1;
            Prop_prm.(i)= Pp_entry;
            Prop_input.(Priors_given_m.Properties.VariableNames{i}) = Prop_prm.(i);
        end
        count = count+1;
    end
elseif strcmp(PKernel_type,"Piecewise")
    % For piecewise independent kernels, do:
    p_i = 1;
    while p_i < p_set_dim(1)+1
        % For each parameter, perturb:
        if strcmp(Perturbation{2,p_i},'Discrete Uniform')
            % Add the uniform kernel with a mean of the parameter from the
            % sampled particle
            Prop_prm.(p_i) = Sampled_Part.(p_i) + random(Perturbation{2,p_i}, Perturbation{1,p_i})-floor(Perturbation{1,p_i}/2);
        else
            % Add the Gaussian kernel with a mean of the parameter from the
            % sampled particle
            Prop_prm.(p_i) = random(Perturbation{2,p_i}, Sampled_Part.(p_i), Perturbation{1,p_i});
        end
        % Check if the pdf of the prior for that parameter has density > 0
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
% This function calculates the weight of a particle at generation t>1.
% Assumes that the priors of each parameter are independent. 

[Priors,Prev_wghts,Prop_prm,Prev_gen,PKernel_type,Kernel_sig] = cell_args{:};

% Currently supported kernel options
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

PKernel_type = cell_args{1}{1};

switch PKernel_type
    case 'covarianceMVN'
        % Calculate the covariance of the previous generation of particles
        Kmatrix = 2*cov(cell_args{2}.Variables);
    case 'MVN'
        % Manually inputted kernel options for multivariate normal
        Kmatrix = cell_args{1}{2};
    case 'Piecewise'
        % Go through each perturbation kernel of the parameters for the
        % piecewise independent option
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

function [All_gen_results_cell] = Unpack(cell_args)
% Unpack results in uncatagorised form and put into a nice, labelled table

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

end
    %% allowed error calculators

function [err] = SquareDiff(cell_args)
% Subfunction for Error. If the method relies on a square difference method,
% this directly calculates it
[Dates,data,Output_T,data_names] = cell_args{:};
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
%% Allowed weight calculators

function [weightdenom] = MVNKernel(cell_args)
% Function that calculates the denominator of the weighting of a particle
% if the kernel of the smc perturbation step is a Multi-Variate Normal
% distribution
[Prev_wghts,Prop_prm,Prev_gen,Kernel_sig] = cell_args{:};

weightdenom = sum(Prev_wghts.*mvnpdf(Prev_gen,Prop_prm,Kernel_sig));
end

function [weightdenom] = Piecewise(cell_args)
% Function that calculates the denominator of the weighting of a particle
% if the kernel of the smc perturbation step is a Piecewise independent
% distribution
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