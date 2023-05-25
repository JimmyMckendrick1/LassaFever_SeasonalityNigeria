function resp = GraphModule(str_function, cell_args)
% This function is a python-like module, containing sub functions of the
% Approximate Bayesian Computation algorithms that I have for fitting my
% epidemiological models.
allowed_functions = {'VarToWords', 'TitleCheck','HistogramParameters',...
    'hist_layout','HistogramParametersMulti','hist_layoutmulti',...
    'HistogramParametersMeta';
                     @VarToWords, @TitleCheck,@HistogramParameters,...
                     @hist_layout, @HistogramParametersMulti,@hist_layoutmulti,...
                     @HistogramParametersMeta};

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

char_input = cell_args{1};

char_output = char_input;
Variables_to_words_dict = {'beta', 'gamma','eps','eps0','eps1','phi','phi1',...
  'phi2','beta_rr','beta_rh','beta_hh','Cases','Infected','I_r_0','gamma_r';
  "\beta", "\gamma","\epsilon", "\epsilon_0","\epsilon_1","\phi","\phi_1",...
  "\phi_2","\beta_{rr}","\beta_{rh}","\beta_{hh}","Fitted to Cumulative Cases","Fitted to Weekly Cases",...
  "I_r(0)","\gamma_r"};

if ~ischar(char_input)
    if isequal(char_input,@Model_SEAIR_const)
        char_output ="SEAIR model with constant forcing" ;
    elseif isequal(char_input,@Model_SEAIR_const_E) 
        char_output ="SEAIR model with constant forcing";
    elseif isequal(char_input,@Model_SEAIR_perGau)
        char_output ="SEAIR model with periodic Guassian forcing";
    elseif isequal(char_input,@Model_SEAIR_perGau_E)
        char_output ="SEAIR model with periodic Guassian forcing";
    elseif isequal(char_input,@Model_SEAIR_sin)
        char_output ="SEAIR model with sinusoidal forcing";
    elseif isequal(char_input,@Model_SEAIR_sin_E)
        char_output ="SEAIR model with sinusoidal forcing";
    elseif isequal(char_input,@Model_SEAIR_step)
        char_output ="SEAIR model with step function forcing";
    elseif isequal(char_input,@Model_SEAIR_step_E)
        char_output ="SEAIR model with step function forcing";
    elseif isequal(char_input,@Model_SEAIR_stepmonth)
        char_output ="SEAIR model with step function forcing";
    elseif isequal(char_input,@Model_SEAIR_stepmonth_E)
        char_output ="SEAIR model with step function forcing";
    end
else
    [found, where] = ismember(char_input, Variables_to_words_dict(1, :));
    if found
        char_output = Variables_to_words_dict{2, where};
    end    
end

end

function [suptitle_var] = TitleCheck(cell_args)

suptitle_var = '';
Length = length(cell_args);

for i = 1:Length
    suptitle_var = append(suptitle_var,' ',GraphModule('VarToWords',cell_args(i)));
end

end

function [F] = HistogramParameters(cell_args)
% Make a single figure of the histograms from the posteriors, with the
% prior overlayed
    Param_tbl = cell_args{1};
    Prior = cell_args{2};
    StartDate=cell_args{3};
    N = length(Param_tbl.Properties.VariableNames)-1; % Number of parameters

    F = figure;
    sgtitle(cell_args{4})
    
    k = 1.8; % Amount we want to scale the range of the parameters to set
    % suitable viewing limits for our plots
    for i = 1:N
        var_name = Param_tbl.Properties.VariableNames{i}; %Name of param
        param_set = Param_tbl.(var_name);

        Layout = GraphModule('hist_layout',{N,i}); % Attempt to put table in a nice layout
        subplot(Layout{:})
        min_prm = inf;
        max_prm = -inf;

        if strcmp(var_name,'phi')  
            % If we have a phi variable, which currently is the only
            % variable that is a time variable, we'll need to make sure it
            % matches up with our perception of time rather than between 0 and 1
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
            if strcmp(Prior{2,i}{1},'Gamma') 
                pdf_prior = pdf(Prior{2,i}{1},x+Prior{1,i}{1}{3},Prior{1,i}{1}{1:2}); %Get the pdf of the prior for this param                
            else
                pdf_prior = pdf(Prior{2,i}{1},x,Prior{1,i}{1}{:});
            end
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
end

function [F] = HistogramParametersMeta(cell_args)
% Make a single figure of the histograms from the posteriors, with the
% prior overlayed
Fnt_sz = 16;
Param_tbl = cell_args{1};
Prior = cell_args{2};
StartDate=cell_args{3};
No_prms = length(Param_tbl.Properties.VariableNames)-1; % Number of parameters
Order = cell_args{4};

k = 1.8; % Amount we want to scale the range of the parameters to set
% suitable viewing limits for our plots
for i = 1:No_prms
    F = figure;
    var_name = Param_tbl.Properties.VariableNames{i}; %Name of param
    sgtitle(GraphModule('VarToWords',{var_name}))
    param_sets = Param_tbl.(var_name);
    for j = 1:width(param_sets)
        param_set = param_sets(:,j);
        Layout = GraphModule('hist_layout',{width(param_sets),j}); % Attempt to put table in a nice layout
        subplot(Layout{:})
        min_prm = inf;
        max_prm = -inf;

        if strcmp(var_name,'phi')  
            % If we have a phi variable, which currently is the only
            % variable that is a time variable, we'll need to make sure it
            % matches up with our perception of time rather than between 0 and 1
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
            if strcmp(Prior{2,i}{1},'Gamma') 
                pdf_prior = pdf(Prior{2,i}{1},x+Prior{1,i}{1}{3},Prior{1,i}{1}{1:2}); %Get the pdf of the prior for this param
            else
                pdf_prior = pdf(Prior{2,i}{1},x,Prior{1,i}{1}{:});
            end
            plot(x,pdf_prior,'Color','r','LineWidth',1.5)
            hold on
            h = histogram(param_set,'Normalization','pdf','FaceColor',[100/255,181/255,230/255]);
            title(Order{j}); 
            xlim([b,a])
        end
        grid on
        if length(cell_args) == 5
            True_input = cell_args{end};
            line([True_input(i),True_input(i)], [0, max(h.Values)], 'Color', 'r', 'LineWidth', 2);
        end
        yt = get(gca, 'YTick');
        set(gca, 'YTick', yt) %, 'YTickLabel', yt/numel(param_set))
        set(gca,'fontsize',Fnt_sz)
        set(gcf,'Position',[0 0 800 800])
    end
end
end

function [F] = HistogramParametersMulti(cell_args)

    Param_cells = cell_args{1};
    Legend_arr = cell_args{2};
    Prior = cell_args{3};
    Starts=cell_args{4};
    N = length(Param_cells{1,1}.Properties.VariableNames)-1; % Number of parameters
    no_hists = length(Param_cells);
    F = figure;
    sgtitle(cell_args{5})
    
    k = 1.8; % Amount we want to scale the range of the parameters to set
    % suitable viewing limits for our plots
    for i = 1:N
        var_name = Param_cells{1,1}.Properties.VariableNames{i}; %Name of param 
        Layout = GraphModule('hist_layout',{N,i}); % Attempt to put table in a nice layout
        subplot(Layout{:})
        min_prm = inf;
        max_prm = -inf;
        if strcmp(var_name,'phi')  
            % If we have a phi variable, which currently is the only
            % variable that is a time variable, we'll need to make sure it
            % matches up with our perception of time rather than between 0 and 1
            for j = 1:no_hists
                Param_tbl = Param_cells{1,j};
                StartDate  = Starts{j};
                StartdateNum = datenum(StartDate,'dd-mm-yyyy');            
                Year = datenum(StartDate(end-3:end),'yyyy'); % What year our data starts
                Offset = StartdateNum-Year; % How much we need to offset what phi = 0 means in real dates
                param_set = Param_tbl.(Param_tbl.Properties.VariableNames{i});
                dateAdj_p_set = 365*param_set+Offset; 
                dateAdj_p_set = mod(dateAdj_p_set,365)+Year; % Date adjusted parameter set

                min_prm = min([min_prm;dateAdj_p_set]);
                max_prm = max([max_prm;dateAdj_p_set]); % Max and min of the parameters

                h = histogram(dateAdj_p_set,'Normalization','pdf');
           
                % Plot the pdf of the prior for this param, currently only a
                % uniform dist. is supported
                hold on
            end
            a = (1-k)*min_prm/2+(1+k)*max_prm/2; % Upper view limit
            b = (1-k)*max_prm/2+(1+k)*min_prm/2; % Lower view limit
            PriorStep = (a-b)/100;
            x = b:PriorStep:a; 
            LB = 365*Prior{1,i}{1}{1}+Year;
            UB = 365*Prior{1,i}{1}{2}+Year;
            pdf_prior = pdf(Prior{2,i}{1},x,LB,UB);
            plot(x,pdf_prior,'Color','r','LineWidth',1.5)
            xlim([b,a])
            xtks = linspace(b,a,8);
            xticks(xtks)
            datetick('x','dd-mmm','keepticks')            
        else
            for j = 1:no_hists
                Param_tbl = Param_cells{1,j};
                param_set = Param_tbl.(Param_tbl.Properties.VariableNames{i});
                min_prm = min([min_prm;param_set]);
                max_prm = max([max_prm;param_set]);   

                hold on
                h = histogram(param_set,'Normalization','pdf');

            end
            a = (1-k)*min_prm/2+(1+k)*max_prm/2;
            b = (1-k)*max_prm/2+(1+k)*min_prm/2;
            PriorStep = (a-b)/100;
            x = b:PriorStep:a;
            if strcmp(Prior{2,i}{1},'Gamma')                
                pdf_prior = pdf(Prior{2,i}{1},x+Prior{1,i}{1}{3},Prior{1,i}{1}{1:2}); %Get the pdf of the prior for this param
            else
                pdf_prior = pdf(Prior{2,i}{1},x,Prior{1,i}{1}{:});
            end
        xlim([b,a])       
        end
                        
        plot(x,pdf_prior,'Color','r','LineWidth',1.5)
        grid on
        if length(cell_args) == 6
            True_input = cell_args{6};
            line([True_input(i),True_input(i)], [0, max(h.Values)], 'Color', 'r', 'LineWidth', 2);
        end
        yt = get(gca, 'YTick');
        set(gca, 'YTick', yt) %, 'YTickLabel', yt/numel(param_set))
           
        title(GraphModule('VarToWords',{var_name})); 
        set(gca,'fontsize',20)
        set(gcf,'Position',[0 0 800 800])
        legend(Legend_arr)
        clear min_prm
        clear max_prm
    end
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

function Layout = hist_layoutmulti(cell_args)
[N,p] = cell_args{:};
switch N 
    case 1
        a = 1;
        b = 1;
        i_star = 1;
    case {2,3,4}
        a = ceil(N/2);
        b = 2;
        i_star = p;
    case {5,6,7,8,9}
        a = ceil(N/3);
        b = 3;
        if p > (a-1)*b && p ~=6
            i_star = 2*p-N+1;
        else
            i_star = p;
        end
    otherwise
        a = ceil(N/4);
        b = 4;
        Pos_vec =  [14,14,15,13,14,15,13,14,15,16];
        i_star = Pos_vec(p-12+(N-12)*(N-11)/2);
end
Layout = {a,b,i_star};
end