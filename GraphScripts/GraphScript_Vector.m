%%              %%%%%%%%%%%%%%Vector%%%%%%%%%%%%%%%%% Single model
% Desktop
addpath('/home/jmckendrick/Documents/MA931/Matlab')
addpath('G:\Work\LHF\LHF-work\GraphnStatsScripts')
%%
addpath('/home/jmckendrick/Documents/LHF-work/GraphnStatsScripts')
addpath('/home/jmckendrick/Documents/MA931/Matlab')
%%
% Load results file
clear
mydir  = pwd;
idcs   = strfind(mydir,filesep);    
newdir = mydir(1:idcs(end)-1); 
load(strcat(newdir,'/TestFits/July2022/FivePrm18_20_res.mat'))
%%
addpath(genpath([pwd,filesep,'matsplit']));
vec_model = @Model_rat_SIR_Freq;
T = length(All_gen_results_cell);

input_v = removevars(Input,{'beta_rh','beta_hh',...
'B_h','p','nu','gamma_h','mu_h_I','mu_h','sigma',...
'N_h_0','S_h_0','E_h_0','A_h_0','I_h_0','C_h_0','D_h_0'});
Test_prms_v = removevars(Priors_given_m,{'beta_rh','beta_hh'});
p_set_dim = zeros(1,width(Test_prms_v)+1);
p_set_dim(1) = width(Test_prms_v); 
tmp_input = input_v;

DataStart = StartWeek;

title_char = {'Vector dynamics of parameter fit', 'to 2018 to 2020 data'};
%% Numbers
X= 18;
dates = DataStart:DataStart+data(end,1);

figure 
sgtitle(title_char)
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 800 800])
ylabel('Number of individuals')
xlabel('Date')
hold on
H = height(All_gen_results_cell{T});

if H>100
    SampleRange = randi([1,H-1],1,100);
else
    SampleRange = 1:H-1;
end
plot(nan,nan,'r')
plot(nan,nan,'g')
plot(nan,nan,'c')
plot(nan,nan,'k')
% Plot a sample of trajectories 
for i_particle = SampleRange
    for kj = 1:p_set_dim(1)
        tmp_input.(Test_prms_v.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{H+1-i_particle,Test_prms_v.Properties.VariableNames{kj}};
    end
    Output_T=vec_model(table2cell(tmp_input));

    hold on
    plot(dates, Output_T{:,2},'Color',[1 0 0 0.1],'LineWidth',1.5);
    plot(dates, Output_T{:,3},'Color',[0 1 0 0.1],'LineWidth',1.5);
    plot(dates, Output_T{:,4},'Color',[0 1 1 0.1],'LineWidth',1.5);

end
% Plot the best fit last and in black
for kj = 1:p_set_dim(1)
    tmp_input.(Test_prms_v.Properties.VariableNames{kj}) =...
        All_gen_results_cell{T}{1,Test_prms_v.Properties.VariableNames{kj}};
end

Output_T=vec_model(table2cell(tmp_input));
hold on
plot(dates,  Output_T{:,2},'Color',[0 0 0],'LineWidth',1.5)
plot(dates,  Output_T{:,3},'Color',[0 0 0],'LineWidth',1.5)
plot(dates,  Output_T{:,4},'Color',[0 0 0],'LineWidth',1.5)
legend('Susceptible','Infected','Recovered','Best Fit','Location','Best')
enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')

%% Quantiles : Separate subdivisions plotted on different graphs
z_values = [1];
Rat_cmpts = {'S_r ','I_r ','R_r '};
NR_cmpts = length(Rat_cmpts);
legend_char = cell(1,NR_cmpts*(length(z_values)+1));
for j = 1:NR_cmpts
    for i = 1:length(z_values)
        legend_char(1,(j-1)*(length(z_values)+1)+i) = strcat(Rat_cmpts(j),sprintf('%d-%d%%',int64((1-z_values(i))*100),int64(z_values(i)*100)));
    end    
    legend_char(1,(j-1)*(length(z_values)+1)+i+1) = strcat(Rat_cmpts(j),'Median');
end
H = height(All_gen_results_cell{T});
SampleRange = 1:H;
Simulated_data_model = NaN(NR_cmpts,length(SampleRange),Input.MaxTime+1);
dates = 0+DataStart:Input.MaxTime+DataStart;

% Run through the sample range to simulate the data
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Test_prms_v.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Test_prms_v.Properties.VariableNames{kj}};
    end
    Output_T = vec_model(table2cell(tmp_input));  
    for i_ratc = 1:NR_cmpts
        Simulated_data_model(i_ratc,i_particle,:) = Output_T{:,i_ratc+1}';
    end
end
% Calculate the median and other credible intervals set in z_values at each
% time point.
CIntervals = NaN(NR_cmpts,2*length(z_values)+1,Input.MaxTime+1);
for i_ratc = 1:NR_cmpts
    CIntervals(i_ratc,end,:) = median(Simulated_data_model(i_ratc,:,:));
end
for i_time =1:Input.MaxTime+1
    for i_z = 1:length(z_values)  
        for i_ratc = 1:NR_cmpts
            tmp_data = sort(Simulated_data_model(i_ratc,:,i_time));
            CIntervals(i_ratc,2*i_z-1,i_time) = tmp_data(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
            CIntervals(i_ratc,2*i_z,i_time) = tmp_data(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
        end
    end
end

% Plots for selected intervals of the simulated data from the fitted model.
figure
X= 18;
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 800 800])
infill_colours = [1 0 0; 0 1 0; 0 1 1];
%
for i_ratc = 1:NR_cmpts
    Curve1 = CIntervals(i_ratc,end,:);
    for i_z = 1:length(z_values)
        Curve2 = CIntervals(i_ratc,2*i_z-1,:);
        Curve3 = CIntervals(i_ratc,2*i_z,:);
        hold on
        plot(dates, Curve2(:),'k','HandleVisibility','off')
        plot(dates, Curve3(:),'k','HandleVisibility','off')
        x2 = [dates, fliplr(dates)];
        inBetween = [Curve2(:)', fliplr(Curve3(:)')];
        fill(x2, inBetween, infill_colours(i_ratc,:), 'FaceAlpha', 0.25*i_z);
    end    
    plot(dates, Curve1(:),'Color',infill_colours(i_ratc,:)/2,'LineWidth',2,'LineStyle','--')
end
set(gcf,'position',[0 0 800 1000])
enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')
[~,h_legend]  = legend(legend_char);    
PatchInLegend = findobj(h_legend, 'type', 'patch');
for i_patch = 1:length(PatchInLegend)
    set(PatchInLegend(i_patch), 'FaceAlpha', 0.25*i_patch);
end

sgtitle(title_char)
ylabel('Number of individuals')
xlabel('Time (days)')
%% Proportions
X= 18;
dates = DataStart:DataStart+data(end,1);

figure 
sgtitle(title_char)
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 800 800])
ylabel('Proportion of population')
xlabel('Date')
hold on
H = height(All_gen_results_cell{T});

if H>100
    SampleRange = randi([1,H-1],1,100);
else
    SampleRange = 1:H-1;
end
plot(nan,nan,'r')
plot(nan,nan,'g')
plot(nan,nan,'c')
plot(nan,nan,'k')
% Plot a sample of trajectories 
for i_particle = SampleRange
    for kj = 1:p_set_dim(1)
        tmp_input.(Test_prms_v.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{H+1-i_particle,Test_prms_v.Properties.VariableNames{kj}};
    end
    Output_T=vec_model(table2cell(tmp_input));
    N = Output_T.S+Output_T.I+Output_T.R;
    hold on
    plot(dates, Output_T{:,2}./N,'Color',[1 0 0 0.1],'LineWidth',1.5);
    plot(dates, Output_T{:,3}./N,'Color',[0 1 0 0.1],'LineWidth',1.5);
    plot(dates, Output_T{:,4}./N,'Color',[0 1 1 0.1],'LineWidth',1.5);

end
% Plot the best fit last and in black
for kj = 1:p_set_dim(1)
    tmp_input.(Test_prms_v.Properties.VariableNames{kj}) =...
        All_gen_results_cell{T}{1,Test_prms_v.Properties.VariableNames{kj}};
end

Output_T=vec_model(table2cell(tmp_input));
hold on
plot(dates,  Output_T{:,2}./N,'Color',[0 0 0],'LineWidth',1.5)
plot(dates,  Output_T{:,3}./N,'Color',[0 0 0],'LineWidth',1.5)
plot(dates,  Output_T{:,4}./N,'Color',[0 0 0],'LineWidth',1.5)
legend('Susceptible','Infected','Recovered','Best Fit','Location','Best')
enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')
%% Quantiles : Separate subdivisions plotted on different graphs
z_values = [1];
Rat_cmpts = {'S_r ','I_r ','R_r '};
NR_cmpts = length(Rat_cmpts);
legend_char = cell(1,NR_cmpts*(length(z_values)+1));
for j = 1:NR_cmpts
    for i = 1:length(z_values)
        legend_char(1,(j-1)*(length(z_values)+1)+i) = strcat(Rat_cmpts(j),sprintf('%d-%d%%',int64((1-z_values(i))*100),int64(z_values(i)*100)));
    end    
    legend_char(1,(j-1)*(length(z_values)+1)+i+1) = strcat(Rat_cmpts(j),'Median');
end
H = height(All_gen_results_cell{T});
SampleRange = 1:H;
Simulated_data_model = NaN(NR_cmpts,length(SampleRange),Input.MaxTime+1);
dates = 0+DataStart:Input.MaxTime+DataStart;

% Run through the sample range to simulate the data
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Test_prms_v.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Test_prms_v.Properties.VariableNames{kj}};
    end
    Output_T = vec_model(table2cell(tmp_input));  
    for i_ratc = 1:NR_cmpts
        Simulated_data_model(i_ratc,i_particle,:) = Output_T{:,i_ratc+1}';
    end
end
% Calculate the median and other credible intervals set in z_values at each
% time point.
CIntervals = NaN(NR_cmpts,2*length(z_values)+1,Input.MaxTime+1);
for i_ratc = 1:NR_cmpts
    CIntervals(i_ratc,end,:) = median(Simulated_data_model(i_ratc,:,:));
end
for i_time =1:Input.MaxTime+1
    for i_z = 1:length(z_values)  
        for i_ratc = 1:NR_cmpts
            tmp_data = sort(Simulated_data_model(i_ratc,:,i_time));
            CIntervals(i_ratc,2*i_z-1,i_time) = tmp_data(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
            CIntervals(i_ratc,2*i_z,i_time) = tmp_data(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
        end
    end
end

% Plots for selected intervals of the simulated data from the fitted model.
figure
X= 18;
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 800 800])
infill_colours = [1 0 0; 0 1 0; 0 1 1];
%
for i_ratc = 1:NR_cmpts
    Curve1 = CIntervals(i_ratc,end,:);
    for i_z = 1:length(z_values)
        Curve2 = CIntervals(i_ratc,2*i_z-1,:);
        Curve3 = CIntervals(i_ratc,2*i_z,:);
        hold on
        plot(dates, Curve2(:),'k','HandleVisibility','off')
        plot(dates, Curve3(:),'k','HandleVisibility','off')
        x2 = [dates, fliplr(dates)];
        inBetween = [Curve2(:)', fliplr(Curve3(:)')];
        fill(x2, inBetween, infill_colours(i_ratc,:), 'FaceAlpha', 0.25*i_z);
    end    
    plot(dates, Curve1(:),'Color',infill_colours(i_ratc,:)/2,'LineWidth',2,'LineStyle','--')
end
set(gcf,'position',[0 0 800 1000])
enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')
[~,h_legend]  = legend(legend_char);    
PatchInLegend = findobj(h_legend, 'type', 'patch');
for i_patch = 1:length(PatchInLegend)
    set(PatchInLegend(i_patch), 'FaceAlpha', 0.25*i_patch);
end

sgtitle(title_char)
ylabel('Number of individuals')
xlabel('Time (days)')
%%

for kj = 1:p_set_dim(i_model)
    tmp_input.(Priors_given_m{i_model}.Properties.VariableNames{kj}) = All_gen_results_cell{cell_ind}{H+1-i_particle,kj};
end
Output_T = M_star(table2cell(tmp_input));
Simulated_data_model_i = Output_T.(error_method{2});    
%% 5-95 rangee for params and median
T = 10;
No_Params = width(Priors_given_m);
R_x_M = NaN(3,No_Params);
for i = 1:No_Params
 
    Tmp_vec = sort(Accepted_params{1,T}{:,i});
    R_x_M(1,i) = Tmp_vec(floor(NoPart/20));
    R_x_M(2,i) = Tmp_vec(floor(NoPart/2));
    R_x_M(3,i) = Tmp_vec(floor(0.95*NoPart));
end
R_x_M = array2table(R_x_M,'VariableNames',Priors_given_m.Properties.VariableNames);

%%              %%%%%%%%%%%%%%Proportion%%%%%%%%%%%%%%%%%
addpath(genpath([pwd,filesep,'matsplit']));
%prop_Models = {@ModelProp_SEAIR_vSIR_pg_f,@Model_SEAIR_vector_SIR_sig_prop};
prop_model = @ModelProp_SEAIR_vSIR_pg_f;
T = 10;%length(tolVec)-1;

%{ 
s = 989.716296086384;  mu_r=1/500; beta_rr=1.90909122362604;
gamma_r=1/90; phi = 0.479938650962018;
N_r_0=1e5; I_r_0=0.5*N_r_0; 
% Humans
B_h = 1/(53.5*365)+1e-4; beta_hh = 0.00130994248902666; beta_rh = 0.00746899555324877;
p = 0.8; nu = 1/14; gamma_h = 1/14; mu_h_I = 0.017857142857143; 
mu_h = 1/(53.5*365);
N_h_0 = 1e6; S_h_0 = 1e5; E_hh_0 = 10; E_hr_0 = 0;
A_h_0 = 20; I_h_0 = 2; 
C_hh_0 = 0; C_hr_0 = 2; D_h_0 = 0;
MaxTime=576;
%}

Input = table(s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
    beta_hh,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,N_h_0,S_h_0,E_hh_0,E_hr_0,...
    A_h_0,I_h_0,C_hh_0,C_hr_0,D_h_0,MaxTime);


if exist('Models','var')    
    No_models = length(Models); 
    p_set_dim = NaN(1,No_models);
    input_size = NaN(1,No_models);
    MaxTime = 576;%data(1,end);%True_input.MaxTime;
    Max_params = 0; % Find the maximum number of parameters tested for a model
    for j_model = 1:No_models
        if width(Priors_given_m{1,j_model}) >= Max_params
            Max_params = width(Priors_given_m{1,j_model});
        end
    end
    prm_index = NaN(Max_params,No_models);
    for i_model = 1:No_models
        input_v = removervars(Inputs{1,i_model},{'beta_rh','beta_hh',...
            'B_h','p','nu','gamma_h','mu_h_I','mu_h',...
            'N_h_0','S_h_0','E_h_0','A_h_0','I_h_0','C_h_0','D_h_0'});

        Test_prms = Priors_given_m{i_model};
        Test_prms = removevars(Test_prms,{'beta_rh','beta_hh'});
        p_set_dim(i_model) = width(Test_prms);
        input_size(i_model) = width(Inputs{1,i_model}); 
        Names_testedvars_model_i = cell(1,p_set_dim(i_model)+1); % Getting the names of tested parameters

        % Get the name and position of the parameters to be tested so that 
        % 1) The input for the model can be easily changed when sampling from 
        % their prior distribution in the main body of the ABC. 2) Label the results
        for i2 = 1:p_set_dim(i_model)
            var_name = Test_prms.Properties.VariableNames{i2};
            Names_testedvars_model_i{i2} = var_name;
            for j = 1:input_size(i_model)
                if strcmp(var_name,input.Properties.VariableNames{j})
                    prm_index(i2,i_model) = j;
                    break
                end
            end
        end 
    end
else
    p_set_dim = width(Priors_given_m);
    input_size = width(Input); 
    Names_testedvars_model_i = cell(1,p_set_dim+1); % Getting the names of tested parameters
    prm_index = NaN(p_set_dim,1);
    MaxTime = data(1,end);

    
    % Get the name and position of the parameters to be tested so that 
    % 1) The input for the model can be easily changed when sampling from 
    % their prior distribution in the main body of the ABC. 2) Label the results
    for i2 = 1:p_set_dim
        var_name = Priors_given_m.Properties.VariableNames{i2};
        Names_testedvars_model_i{i2} = var_name;
        for j = 1:input_size
            if strcmp(var_name,Input.Properties.VariableNames{j})
                prm_index(i2) = j;
                break
            end
        end
    end
end

%%
addpath(genpath([pwd,filesep,'matsplit']));
prop_Models = {@Model_SEAIR_vector_SIR_prop,@Model_SEAIR_vector_SIR_sig_prop};
T = 10;

E_hh_0 = 10; E_hr_0 = 0; I_hr_0 = 0; I_hh_0 = 2; C_hh_0 = 0; C_hr_0 = 2; D_h_0 = 0;
tmp_inputs = {table([1]),table([1])};

input = table2cell(Inputs{1});
    
[s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
    beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,...
    N_h_0,S_h_0,E_h_0,A_h_0,I_h_0,C_h_0,D_h_0,MaxTime]=input{:};

    
tmp_inputs{1} = table(s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
        beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,N_h_0,S_h_0,E_hh_0,E_hr_0,...
        A_h_0,I_hh_0,I_hr_0,C_hh_0,C_hr_0,D_h_0,MaxTime);

input = table2cell(Inputs{4});
[s,phi,beta_rr,a_1,a_2,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
    beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,...
     N_h_0,S_h_0,E_h_0,A_h_0,I_h_0,C_h_0,D_h_0,MaxTime]=input{:};
 
tmp_inputs{2}=table(s,phi,beta_rr,a_1,a_2,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
        beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,N_h_0,S_h_0,E_hh_0,E_hr_0,...
        A_h_0,I_hh_0,I_hr_0,C_hh_0,C_hr_0,D_h_0,MaxTime);
%%
Indd=[1,4];
%Prop = cell(1,2);

for i_model =1
     cell_ind = (T-1)*No_models+Indd(i_model);
    if ~isempty(All_gen_results_cell{cell_ind})
        Prop{i_model} = NaN(2,height(All_gen_results_cell{cell_ind}));
        M_star = prop_Models{i_model};
        tmp_input = tmp_inputs{i_model};
        H = height(All_gen_results_cell{cell_ind});        
        tmp_input = tmp_inputs{i_model};
        % Plot 100, or however much you can, trajectories to see the final
        % generation's fit at a glance
        for i_particle = 1:H
            for kj = 1:p_set_dim(Indd(i_model))
                tmp_input.(Priors_given_m{Indd(i_model)}.Properties.VariableNames{kj}) = All_gen_results_cell{cell_ind}{H+1-i_particle,kj};
            end
            Output_T = M_star(table2cell(tmp_input));
            Prop{i_model}(1,i_particle) = Output_T.C_h_r(end);
            Prop{i_model}(2,i_particle) = Output_T.C_h_h(end);
        end
    end
    Prop{i_model}= rmmissing(Prop{i_model});
end