function resp = GraphModule(str_function, cell_args)
% This function is a python-like module, containing plotting sub functions 
% that are used in the paper Modelling seasonality in Lassa fever in
% Nigeria
allowed_functions = {'VarToWords', 'TitleCheck','HistogramParameters','hist_layout';
                     @VarToWords, @TitleCheck,@HistogramParameters,@hist_layout};

[found, where] = ismember(str_function, allowed_functions(1, :));  %search 1st row of lookup table
if ~found
    error('Function %s is not a valid function');
end
%if isequal(where,3)
%    allowed_functions{2, where}(cell_args);
%else
%    resp = allowed_functions{2, where}(cell_args); %use matching 2nd row of lookup table
%end 
resp = allowed_functions{2, where}(cell_args);
end

function [char_output] = VarToWords(cell_args)
% This function takes various inputs such as variable names of model
% parameters and function handles, and converts them to latex characters
% for graph titles etc

char_input = cell_args{1};

char_output = char_input;
Variables_to_words_dict = {'beta', 'gamma','phi','beta_rr','beta_rh','beta_hh','I_h','I_r_0','gamma_r';
  "\beta", "\gamma","\phi","\beta_{rr}","\beta_{rh}","\beta_{hh}","Fitted to Cumulative Cases","I_r(0)","\gamma_r"};

if ~ischar(char_input)
    if isequal(char_input,@Model_SEAIR_vSIR)
        char_output = "";
    elseif isequal(char_input,@ModelProp_SEAIR_vSIR)
        char_output = "";
    elseif isequal(char_input,@ModelVector_SIR)
        char_output = "";
    end
else
    [found, where] = ismember(char_input, Variables_to_words_dict(1, :));
    if found
        char_output = Variables_to_words_dict{2, where};
    end    
end

end

function [suptitle_var] = TitleCheck(cell_args)
% Function that converts a
suptitle_var = '';
Length = length(cell_args);

for i = 1:Length
    suptitle_var = append(suptitle_var,' ',GraphModule('VarToWords',cell_args(i)));
end

end

function [F] = HistogramParameters(cell_args)
% Make a single figure  where eaech subplot is a histogram of a one parameter. The
% prior is overlayed

    Param_tbl = cell_args{1}; % One generation of parameters
    Prior = cell_args{2}; % The priors used to fit 
    StartDate=cell_args{3}; % Datetime of the start of the data
    N = length(Param_tbl.Properties.VariableNames)-1; % Number of parameters
    F = figure;
    sgtitle(cell_args{4})
    
    k = 1.8; %Scale the range of the parameters to set suitable viewing limits for plots
    for i = 1:N %For each parameter...
        var_name = Param_tbl.Properties.VariableNames{i}; %Name of param
        param_set = Param_tbl.(var_name); 

        Layout = GraphModule('hist_layout',{N,i}); % Attempt to put table in a nice layout
        subplot(Layout{:})
        min_prm = inf;
        max_prm = -inf;

        if strcmp(var_name,'phi')  
            % If we have phi, which currently is the only variable that is a time variable,
            % we'll need to make sure it is in a nice format.
            StartdateNum = datenum(StartDate,'dd-mm-yyyy');            
            Year = datenum(StartDate(end-3:end),'yyyy'); % What year our data starts
            Offset = StartdateNum-Year; % How much we need to offset what phi = 0 means in real dates

            dateAdj_p_set = 365*param_set+Offset; 
            dateAdj_p_set = mod(dateAdj_p_set,365)+Year; % Date adjusted parameter set

            min_prm = min([min_prm;dateAdj_p_set]);
            max_prm = max([max_prm;dateAdj_p_set]); % Max and min of the parameters

            a = (1-k)*min_prm/2+(1+k)*max_prm/2; % Upper view limit
            b = (1-k)*max_prm/2+(1+k)*min_prm/2; % Lower view limit
            PriorStep = (a-b)/100;
            x = b:PriorStep:a; 
            LB = 365*Prior{1,i}{1}{1}+Year;
            UB = 365*Prior{1,i}{1}{2}+Year;

            % Prior of the pdf between the upper and lower bounds
            pdf_prior = pdf(Prior{2,i}{1},x,LB,UB); 
            plot(x,pdf_prior,'Color','r','LineWidth',1.5)
            % Plot the pdf of the prior for this param, currently only a
            % uniform dist. is supported
            hold on
            h = histogram(dateAdj_p_set,'Normalization','pdf','FaceColor',[100/255,181/255,230/255]);
            xlim([b,a])
            xtks = linspace(b,a,8);
            xticks(xtks)
            datetick('x','dd-mmm','keepticks')            
        else
            min_prm = min([min_prm;param_set]);
            max_prm = max([max_prm;param_set]); % Max and min of the parameters
            a = (1-k)*min_prm/2+(1+k)*max_prm/2;
            b = (1-k)*max_prm/2+(1+k)*min_prm/2;
            PriorStep = (a-b)/100;
            x = b:PriorStep:a;
            pdf_prior = pdf(Prior{2,i}{1},x,Prior{1,i}{1}{:});            
            plot(x,pdf_prior,'Color','r','LineWidth',1.5)
            hold on
            h = histogram(param_set,'Normalization','pdf','FaceColor',[100/255,181/255,230/255]);
            xlim([b,a])
        end
        grid on
        if length(cell_args) == 5
            True_input = cell_args{5};
            line([True_input{1,i},True_input{1,i}], [0, max(h.Values)*1.2], 'Color', 'r', 'LineWidth', 2);
        end
        
        yt = get(gca, 'YTick');
        set(gca, 'YTick', yt) %, 'YTickLabel', yt/numel(param_set))

        title(GraphModule('VarToWords',{var_name})); 
        set(gca,'fontsize',20)
        set(gcf,'Position',[0 0 1200 800])
    end
    legend('Prior distribution','Posterior distribution')
end

function Layout = hist_layout(cell_args)
[N,p] = cell_args{:};
if N == 1
    a = 1;
    b = 1;
    i_star = 1;
elseif N < 7
    a = ceil(N/2);
    b = 2;
    i_star = p;
else 
    a = ceil(N/3);
    b = 3;
    i_star = p;
    if N == 7 && p == 7
        i_star = 8;
    elseif N == 8 && p > 6 
        i_star = 2*p-7;
    end
end
Layout = {a,b,i_star};
end