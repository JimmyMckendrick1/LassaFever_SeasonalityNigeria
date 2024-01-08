%%              %%%%%%%%%%%%%%Vector%%%%%%%%%%%%%%%%% Single model

mydir  = pwd;
idcs   = strfind(mydir,filesep);    
upperdir = mydir(1:idcs(end)-1); 
addpath(strcat(upperdir,'/Fits_Models'))
addpath(strcat(upperdir,'/GraphScripts'))
%%
% Load results file
load(strcat(upperdir,'/Inputs/Manuscript_input.mat'))
load(strcat(upperdir,'/Results/Manuscript_results.mat'))
%%
T =length(All_gen_results_cell);
p_set_dim = zeros(1,width(Priors_given_m)+1);
p_set_dim(1) = width(Priors_given_m);
tmp_input = Input; 
StartWeek = 737061;
DataStart = StartWeek;

title_char = {'Vector dynamics of parameter fit', 'to 2018 to 2020 data'};
%% Numbers
X= 18;
dates = DataStart:DataStart+data(end,1);

figure 
sgtitle(title_char)
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 1600 2000])
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
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end   
    Output_T=model(table2cell(tmp_input));

    hold on
    plot(dates, Output_T.S_r,'Color',[1 0 0 0.1],'LineWidth',1.5);
    plot(dates, Output_T.I_r,'Color',[0 1 0 0.1],'LineWidth',1.5);
    plot(dates, Output_T.R_r,'Color',[0 1 1 0.1],'LineWidth',1.5);

end
% Plot the best fit last and in black
for kj = 1:p_set_dim(1)
    tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
        All_gen_results_cell{T}{1,Priors_given_m.Properties.VariableNames{kj}};
end

Output_T=model(table2cell(tmp_input));
hold on
plot(dates,  Output_T.S_r,'Color',[0 0 0],'LineWidth',1.5)
plot(dates,  Output_T.I_r,'Color',[0 0 0],'LineWidth',1.5)
plot(dates,  Output_T.R_r,'Color',[0 0 0],'LineWidth',1.5)
legend('Susceptible','Infected','Recovered','Best Fit','Location','Best')
enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')

%% Quantiles : Separate subdivisions plotted on different graphs
z_values = [0.95];

Rat_cmpts = {'S_r ','I_r ','R_r '};
tmp_input = Input;
NR_cmpts = length(Rat_cmpts);
legend_char = cell(1,NR_cmpts*(length(z_values)+1));
for j = 1:NR_cmpts
    for i = 1:length(z_values)
        legend_char(1,(j-1)*(length(z_values)+1)+i) = strcat(Rat_cmpts(j),sprintf('%d%%',int64((2*z_values(i)-1)*100)));
    end    
    legend_char(1,(j-1)*(length(z_values)+1)+i+1) = strcat(Rat_cmpts(j),'Median');
end
H = height(All_gen_results_cell{T});
SampleRange = 1:H;
Simulated_data_model = NaN(NR_cmpts,length(SampleRange),tmp_input.MaxTime+1);
dates = 0+DataStart:tmp_input.MaxTime+DataStart;

% Run through the sample range to simulate the data
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end
    Output_T = model(table2cell(tmp_input));  
    for i_ratc = 1:NR_cmpts
        Simulated_data_model(i_ratc,i_particle,:) = Output_T{:,i_ratc+1}';
    end
end
% Calculate the median and other credible intervals set in z_values at each
% time point.
CIntervals = NaN(NR_cmpts,2*length(z_values)+1,tmp_input.MaxTime+1);
for i_ratc = 1:NR_cmpts
    CIntervals(i_ratc,end,:) = median(Simulated_data_model(i_ratc,:,:));
end
for i_time =1:tmp_input.MaxTime+1
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
z_values = [0.95];
Rat_cmpts = {'S_r ','I_r ','R_r '};
tmp_input = Input;
NR_cmpts = length(Rat_cmpts)+1;
legend_char = cell(1,NR_cmpts*(length(z_values)+1));
for j = 1:NR_cmpts-1
    for i = 1:length(z_values)
        legend_char(1,(j-1)*(length(z_values)+1)+i) = strcat(Rat_cmpts(j),sprintf('%d%%',int64((2*z_values(i)-1)*100)),'Credible Interval');
    end    
    legend_char(1,(j-1)*(length(z_values)+1)+i+1) = strcat(Rat_cmpts(j),'Median');
end

for i = 1:length(z_values)
    legend_char(1,(NR_cmpts-1)*(length(z_values)+1)+i) = strcat({'N_r '},sprintf('%d%%',int64((2*z_values(i))*100-1)),'Credible Interval');
end    
legend_char(1,(NR_cmpts-1)*(length(z_values)+1)+i+1) = strcat({'N_r '},'Median');

H = height(All_gen_results_cell{T});
SampleRange = 1:H;
Simulated_data_model = NaN(NR_cmpts,length(SampleRange),tmp_input.MaxTime+1);
dates = 0+DataStart:tmp_input.MaxTime+DataStart;

% Run through the sample range to simulate the data
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end
    Output_T = model(table2cell(tmp_input));  
    N_r = zeros(1,tmp_input.MaxTime+1);
    for i_ratc = 1:NR_cmpts-1
        Simulated_data_model(i_ratc,i_particle,:) = Output_T{:,i_ratc+1}';
        N_r= N_r + Output_T{:,i_ratc+1}';
    end
    Simulated_data_model(NR_cmpts,i_particle,:) = reshape(N_r,[],1,length(N_r));
    Simulated_data_model(1:NR_cmpts-1,i_particle,:) = Simulated_data_model(1:NR_cmpts-1,i_particle,:)./Simulated_data_model(NR_cmpts,i_particle,:);
end
% Calculate the median and other credible intervals set in z_values at each
% time point.
CIntervals = NaN(NR_cmpts,2*length(z_values)+1,tmp_input.MaxTime+1);
for i_ratc = 1:NR_cmpts
    CIntervals(i_ratc,end,:) = median(Simulated_data_model(i_ratc,:,:));
end
for i_time =1:tmp_input.MaxTime+1
    for i_z = 1:length(z_values)  
        for i_ratc = 1:NR_cmpts
            tmp_data = sort(Simulated_data_model(i_ratc,:,i_time));
            CIntervals(i_ratc,2*i_z-1,i_time) = tmp_data(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
            CIntervals(i_ratc,2*i_z,i_time) = tmp_data(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
        end
    end
end

%% Plots for selected intervals of the simulated data from the fitted model.
figure

X= 18;
set(gca,'fontsize',X)
set(gcf,'Position',[0 0 1600 1000])
infill_colours = [1 0 0; 0 1 0; 0 1 1; 0 0 0];
%
for i_ratc = 1:NR_cmpts-1
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
enddate = dates(end);
xlim([DataStart,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')
[~,h_legend]  = legend(legend_char);    
PatchInLegend = findobj(h_legend, 'type', 'patch');
for i_patch = 1:length(PatchInLegend)-1
    set(PatchInLegend(i_patch), 'FaceAlpha', 0.25*i_patch);
end

ylabel('Proportion of individuals')
xlabel('Time (days)')

%% Plot total number of rats
subplot(2,1,1)
set(gca,'fontsize',X)
for i_ratc = NR_cmpts
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
xticks([])
[~,h_legend]  = legend(legend_char);    
PatchInLegend = findobj(h_legend, 'type', 'patch');
for i_patch = length(PatchInLegend)
    set(PatchInLegend(i_patch), 'FaceAlpha', 0.25*i_patch);
end

%sgtitle(title_char)
ylabel('Proportion of individuals') 

%%              %%%%%%%%%%%%%%Proportion%%%%%%%%%%%%%%%%%%
prop_model = @ModelProp_SEAIR_vSIR;

tmp_input = Input;
tmp_input.E_hr_0 = 0;
tmp_input.I_hr_0 = 0;
tmp_input.C_hr_0 = 0;

Prop = NaN(2,H);

% Get the name and position of the parameters to be tested so that 
% 1) The input for the model can be easily changed when sampling from 
% their prior distribution in the main body of the ABC. 2) Label the results
for i_particle = 1:H
    for kj = 1:p_set_dim(1)
        tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
            All_gen_results_cell{T}{i_particle,Priors_given_m.Properties.VariableNames{kj}};
    end
    Output_T = prop_model(table2cell(tmp_input));  
    Prop(1,i_particle) = Output_T.C_hr(end);
    Prop(2,i_particle) = Output_T.C_hh(end);
    Prop(3,i_particle) = 100*Output_T.C_hr(end)/(Output_T.C_hr(end)+Output_T.C_hh(end));

end
Proportions_sorted = sort(Prop(3,:));
Proportions_sorted([int64(0.05*NoPart),int64(0.5*NoPart),int64(0.95*NoPart)])
