%%                           Single model Graph production
% A script run in sections that will produce the plots as seen in 

mydir  = pwd;
idcs   = strfind(mydir,filesep);    
upperdir = mydir(1:idcs(end-1)-1); 
addpath(strcat(upperdir,'/Fits_Models'))
addpath(strcat(upperdir,'/GraphScripts'))
%%
% Load results file
clear
mydir  = pwd;
idcs   = strfind(mydir,filesep);    
newdir = mydir(1:idcs(end-1)-1); 
load(strcat(newdir,'/MA931/TestFits/ABCvalidation/AQSMC_val_res_5.mat'))
load(strcat(newdir,'/MA931/TestFits/ABCvalidation/AQSMC_val.mat'))
%%
clear
mydir  = pwd;
idcs   = strfind(mydir,filesep);    
newdir = mydir(1:idcs(end)-1); 
load(strcat(newdir,'/TestFits/FirstPaperResults/FivePrm18_20xN2I3_res.mat'))
%newdir = mydir(1:idcs(end-1)-1); 
%load(strcat(newdir,'/LHF-work/TestFits/July2022/SixPrm18_20_res.mat'))
%% Finding important stuff for plotting before we plot
addpath(genpath([pwd,filesep,'matsplit']));
T = 12;%length(All_gen_results_cell);
p_set_dim = zeros(1,width(Priors_given_m)+1);
p_set_dim(1) = width(Priors_given_m);

tmp_input = Input;
%startdate = datenum('14-12-2018','dd-mm-yyyy');
%DataStart = StartWeek;
DataStart = 737061;
title_char = {'Parameter fit for Nigeria to data from','2018 to 2020 epidemics'};
%model = Model;
%% Quick plots to data that was fitted for 
X= 18;
figure 
sgtitle(title_char)
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 800 800])
ylabel('Number of individuals')
xlabel('Time (days)')
hold on
H = height(All_gen_results_cell{T});

if H>100
    SampleRange = randi([1,H-1],1,100);
else
    SampleRange = 1:H-1;
end
for i_particle = SampleRange(1:end-1)
    for kj = 1:p_set_dim(1)
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{H+1-i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end
    Output_T = model(table2cell(tmp_input));
    Simulated_data_model_i = Output_T.(error_method{2});    
    dates = 0+DataStart:length(Simulated_data_model_i)-1+DataStart; 

    plot(dates,Simulated_data_model_i,'Color',[0.85 0.325 0.098 0.1],'LineWidth',1.5,'HandleVisibility','off');
end
% One plot for legend purposes, the sample simulation
i_particle = SampleRange(end);
for kj = 1:p_set_dim(1)
    tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
        All_gen_results_cell{T}{H+1-i_particle,Priors_given_m.Properties.VariableNames{kj}};
end
Output_T = model(table2cell(tmp_input));
Simulated_data_model_i = Output_T.(error_method{2});    
plot(dates,Simulated_data_model_i,'Color',[0.85 0.325 0.098],'LineWidth',1.5);

% Best fit
i_particle = H;
for kj = 1:p_set_dim(1)
    tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
        All_gen_results_cell{T}{H+1-i_particle,Priors_given_m.Properties.VariableNames{kj}};
end
Output_T=model(table2cell(tmp_input));
Simulated_data_model_i = Output_T.(error_method{2});
ylabel('Number of Cases')
xlabel('Date')

dates = 0+DataStart:length(Simulated_data_model_i)-1+DataStart;
plot(dates,Simulated_data_model_i,'Color',[0 1 1],'LineWidth',1.5)
plot(data(:,1)+DataStart,data(:,2),'Color','k','Marker','+')

enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')
legend('Sampled Trajectories','Best Fit Trajectory','Case Data','Location','Best')
%% Quantiles : Separate subdivisions plotted on different graphs
z_values = [1];
legend_char = cell(1,length(z_values)+2);
for i = 1:length(z_values)
    legend_char{1,i} = sprintf('%d-%d%%',int64((1-z_values(i))*100),int64(z_values(i)*100));
end
legend_char{1,end-1} = 'data';
legend_char{1,end} = 'Median';
H = height(All_gen_results_cell{T});
SampleRange = 1:H;
Simulated_data_model = NaN(length(SampleRange),Input.MaxTime+1);
dates = 0+DataStart:Input.MaxTime+DataStart;

% Run through the sample range to simulate the data
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end
    Output_T = model(table2cell(tmp_input));  
    Simulated_data_model(i_particle,:) = Output_T.(error_method{2}{1})(:)';            
end
% Calculate the median and other credible intervals set in z_values at each
% time point.
CIntervals = NaN(2*length(z_values)+1,Input.MaxTime+1);
CIntervals(end,:) = median(Simulated_data_model(:,:));
for i_time =1:Input.MaxTime+1
    for i_z = 1:length(z_values)  
        tmp_data = sort(Simulated_data_model(:,i_time));
        CIntervals(2*i_z-1,i_time) = tmp_data(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
        CIntervals(2*i_z,i_time) = tmp_data(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
    end
end

% Plots for selected intervals of the simulated data from the fitted model.
figure
X= 18;
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 800 1000])
Curve1 = CIntervals(end,:);
for i_z = 1:length(z_values)
    Curve2 = CIntervals(2*i_z-1,:);
    Curve3 = CIntervals(2*i_z,:);
    hold on
    plot(dates, Curve2,'k','HandleVisibility','off')
    plot(dates, Curve3,'k','HandleVisibility','off')
    x2 = [dates, fliplr(dates)];
    inBetween = [Curve2, fliplr(Curve3)];
    fill(x2, inBetween, [1 14/25 0], 'FaceAlpha', 0.25*i_z);
end
plot(data(:,1)+DataStart,data(:,2),'Color',[0,0,0,1],'Marker','+','LineWidth',2)
plot(dates, Curve1,'b','LineWidth',2)
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
%sgtitle({'SMC fit with qauntile-calculated tolerances', 'Infected'})
ylabel('Number of individuals')
xlabel('Time (days)')
%% Quantiles : ABCrej, MCMC Separate subdivisions plotted on different graphs
p_set_dim = zeros(1,width(Priors_given_m)+1);
p_set_dim(1) = width(Priors_given_m);

tmp_input = Input;
H = height(All_gen_results_cell);
SampleRange = 1:H;
Simulated_data_model = NaN(length(SampleRange),Input.MaxTime+1);
z_values = [0.9,0.75];
% Run through the sample range to simulate the data
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end
    Output_T = model(table2cell(tmp_input));  
    Simulated_data_model(i_particle,:) = Output_T.(error_method{2}{1})(:)';            
end
% Calculate the median and other credible intervals set in z_values at each
% time point.
CIntervals = NaN(2*length(z_values)+1,Input.MaxTime+1);
CIntervals(end,:) = median(Simulated_data_model(:,:));
for i_time =1:Input.MaxTime+1
    for i_z = 1:length(z_values)  
        tmp_data = sort(Simulated_data_model(:,i_time));
        CIntervals(2*i_z-1,i_time) = tmp_data(ceil(z_values(i_z)*length(SampleRange)));
        CIntervals(2*i_z,i_time) = tmp_data(ceil((1-z_values(i_z))*length(SampleRange)));
    end
end
%% Plots for selected intervals of the simulated data from the fitted model.
figure
Curve1 = CIntervals(end,:);
for i_z = 1:length(z_values)
    Curve2 = CIntervals(2*i_z-1,:);
    Curve3 = CIntervals(2*i_z,:);
    hold on
    plot(0:17, Curve2,'k','HandleVisibility','off')
    plot(0:17, Curve3,'k','HandleVisibility','off')
    x2 = [0:17, fliplr(0:17)];
    inBetween = [Curve2, fliplr(Curve3)];
    fill(x2, inBetween, [1 14/25 0], 'FaceAlpha', 0.25*i_z);
end
set(gcf,'position',[0 0 600 600])
plot(data(:,1),data(:,2),'Color',[0,0,0,1],'Marker','+','LineWidth',2)
plot(0:17, Curve1,'b','LineWidth',2)
enddate = 17;
xlim([0,17])
[~,h_legend]  = legend(sprintf('%d-%d%%',int64((1-z_values(1))*100),int64(z_values(1)*100)),...
    sprintf('%d-%d%%',int64((1-z_values(2))*100),int64(z_values(2)*100)),...
    'data','Median');    
PatchInLegend = findobj(h_legend, 'type', 'patch');
set(PatchInLegend(1), 'FaceAlpha', 0.25);
set(PatchInLegend(2), 'FaceAlpha', 0.5);

sgtitle('MCMC rejection fit, Infected')
ylabel('Number of individuals')
xlabel('Time (days)')
%%                  Histograms
%{
if isequal(error_method{1,2},"Infected")
    title_char = {'Vector Model fit to Weekly Data:',' Marginal Posterior Distributions'};
else
    title_char = {'Vector Model fit to Cumulative Data:',' Marginal Posterior Distributions'};               
end
%}
%title_char = {' Marginal posterior distributions ','SMC with quantile-calculated tolerances'};
%title_char = {' Marginal posterior distributions','parameter fit to 2018 - 2020 Nigeria data'};
title_char = {' Marginal posterior distributions','N_r(0) = 5\times 10^6, I_r(0) = 2.5\times 10^5'};
tmp_input = Priors_given_m;
cell_args = cell(1,4); 
cell_args(1:4) = {All_gen_results_cell{T},Priors_given_m,...
    convertStringsToChars(datestr(DataStart)),title_char};

if exist('True_input','var') 
    cell_args{5} = True_input;
    F = GraphModule('HistogramParameters',cell_args);
else
    F = GraphModule('HistogramParameters',cell_args); 
end

set(gcf,'Position',[0 0 800 1000]) % Perhaps change the total size of the
%figure incase things get smushed
%%                    Two Histograms
%%%%%% Outdated
Params_data = {All_gen_results_cell{1,end},All_gen_results_cellO{1,end}};

addpath(genpath([pwd,filesep,'matsplit']));

%Legend_arr = ["Weekly","Cumulative","Prior"];
Legend_arr = ["Edo","Ondo"];%["1","2","3"];
title_char = {'5 parameter fit for Edo and Ondo:','Comparison of Marginal Posterior Distributions'};
Starts = {convertStringsToChars(datestr(StartDate)),...
    convertStringsToChars(datestr(StartWeekO))};

cell_args = {Params_data,Legend_arr,Priors_given_m,...
    Starts,title_char};

F = GraphModule('HistogramParametersMulti',cell_args);
set(gcf,'Position',[0 0 800 1000]) % Perhaps change the total size of the
%figure incase things get smushed

%% Corrplot

[R,PValue] = corrplot(All_gen_results_cell{1,T})

