
mydir  = pwd;
idcs   = strfind(mydir,filesep);    
upperdir = mydir(1:idcs(end)-1); 
addpath(strcat(upperdir,'/Fits_Models'))
addpath(strcat(upperdir,'/GraphScripts'))
%%
load(strcat(upperdir,'/Inputs/Manuscript_input.mat'))
load(strcat(upperdir,'/Results/Manuscript_results.mat'))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up %%%%%%%%%%%%%%%%%%%%%%%
T = length(Total_iterations);

p_set_dim = zeros(1,width(Priors_given_m)+1);
p_set_dim(1) = width(Priors_given_m);
% Get the name of the parameters to be tested to label the results
Result_prms_gen_T = All_gen_results_cell{T};
for i = 1:p_set_dim(1)
    p_set_dim(i+1) = width(Result_prms_gen_T(1,i));
end

R_0_indicator = NaN(NoPart,ceil(2*data(end,1)/365));
R_rr_0vec = NaN(1,NoPart);
R_rh_0vec = NaN(1,NoPart);
R_hh_0vec = NaN(1,NoPart);
R_rrvec = NaN(NoPart,Input.MaxTime+1);
R_rhvec = NaN(NoPart,Input.MaxTime+1);
R_hhvec = NaN(NoPart,Input.MaxTime+1);
SampleRange = 1:NoPart;%randi([1,NoPart],1,100);
StartDate = StartWeek;
dates = 0+StartDate:Input.MaxTime+StartDate;
%SampleRange = 1:NoPart;

%% %%%%%%%%%%%%%%%%%%%%Effective reproduction rates%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% SIR + SEAIR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = 18;
z_values = [1];
CIntervals_rr = NaN(2*length(z_values)+1,Input.MaxTime+1);
CIntervals_rh = NaN(2*length(z_values)+1,Input.MaxTime+1);
CIntervals_hh = NaN(2*length(z_values)+1,Input.MaxTime+1);
%figure
%set(gca,'fontsize',X)
%set(gcf,'Position',[0 0 1000 1000])

Result_prms_gen_T = All_gen_results_cell{1,T};
if ~isempty(Result_prms_gen_T)
    tmp_input = Input;
    for i_particle = 1:height(Result_prms_gen_T)
        for kj = 1:p_set_dim(1)      
            tmp_input.(Priors_given_m.Properties.VariableNames{kj}) =...
                Result_prms_gen_T{end+1-i_particle,Priors_given_m.Properties.VariableNames{kj}};
        end

        Output_T=model(table2cell(tmp_input));
        Resp = R_0_calc({tmp_input,Output_T});
        %Calculated reproduction ratios for each inter-compartment exchange
        R_rrvec(i_particle,:) = Resp(1:end-1,1);
        R_rhvec(i_particle,:) = Resp(1:end-1,2);
        R_hhvec(i_particle,:) = Resp(1:end-1,3);
        
        A = find(Resp(1:end-1,1)>1); % When do the rats get an epidemic      
        A = A(A>100);
        R_0_indicator(i_particle, 1) = A(1);
        kk = 2;
        for k = 2:length(A)
            if A(k) > A(k-1)+1 
                R_0_indicator(i_particle, kk) = A(k-1);
                R_0_indicator(i_particle, kk+1) = A(k);
                kk = kk+2;
            end
        end
        R_0_indicator(i_particle, 4) = A(end);
        R_rr_0vec(i_particle) = Resp(end,1);
        R_rh_0vec(i_particle) = Resp(end,2);
        R_hh_0vec(i_particle) = Resp(end,3);
        
        %{
        if ismember(i_particle,SampleRange)
            plot(dates, R_rr,'Color',[1 0 0 0.1],'LineWidth',1.5,'HandleVisibility','off');
            hold on
        end
        %}
    end
    CIntervals_rr(end,:) = median(R_rrvec);
    CIntervals_rh(end,:) = median(R_rhvec);
    CIntervals_hh(end,:) = median(R_hhvec);
    for i_time =1:Input.MaxTime+1
        for i_z = 1:length(z_values)  
            tmp_data_rr = sort(R_rrvec(:,i_time));
            tmp_data_rh = sort(R_rhvec(:,i_time));
            tmp_data_hh = sort(R_hhvec(:,i_time));
            CIntervals_rr(2*i_z-1,i_time) = tmp_data_rr(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
            CIntervals_rr(2*i_z,i_time) = tmp_data_rr(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
            CIntervals_rh(2*i_z-1,i_time) = tmp_data_rh(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
            CIntervals_rh(2*i_z,i_time) = tmp_data_rh(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
            CIntervals_hh(2*i_z-1,i_time) = tmp_data_hh(ceil(z_values(i_z)*(length(SampleRange)-1)+1));
            CIntervals_hh(2*i_z,i_time) = tmp_data_hh(ceil((1-z_values(i_z))*(length(SampleRange)-1)+1));
        end
    end
    %{
    plot(NaN,NaN,'Color',[1 0 0])
    sgtitle('Effective Reproductive Rate between rats')
    xlim([StartDate,StartDate+Input.MaxTime])
    xtks = [datenum(2018,12,1),datenum(2019,2:2:12,1),datenum(2020,2:2:7,1)];
    xticks(xtks)    
    datetick('x',1,'keepticks')
    set(gca, 'YScale', 'log')
    legend('R_{rr}, rat effective reproductive rate','Location','Best')%,'R_{hh}, human effective reproductive rate')
    %}
end    
%% Cinterval plots
figure
%subplot(3,1,3)
legend_char = cell(1,length(z_values)+1);
for i = 1:length(z_values)
    legend_char{1,i} = sprintf('%d-%d%%',int64((1-z_values(i))*100),int64(z_values(i)*100));
end
legend_char{1,end} = 'Median';

X= 18;
set(gca,'fontsize',X)
set(gcf,'position',[0 0 800 1000])
Curve1 = CIntervals_rr(end,:);
for i_z = 1:length(z_values)
    Curve2 = CIntervals_rr(2*i_z-1,:);
    Curve3 = CIntervals_rr(2*i_z,:);
    hold on
    plot(dates, Curve2,'k','HandleVisibility','off')
    plot(dates, Curve3,'k','HandleVisibility','off')
    x2 = [dates, fliplr(dates)];
    inBetween = [Curve2, fliplr(Curve3)];
    fill(x2, inBetween, [1 14/25 0], 'FaceAlpha', 0.25*i_z);
end
%plot(data(:,1)+DataStart,data(:,2),'Color',[0,0,0,1],'Marker','+','LineWidth',2)
plot(dates, Curve1,'k','LineWidth',2)

[~,h_legend]  = legend(legend_char);    
PatchInLegend = findobj(h_legend, 'type', 'patch');
for i_patch = 1:length(PatchInLegend)
    set(PatchInLegend(i_patch), 'FaceAlpha', 0.25*i_patch);
end

%sgtitle({'SMC fit with qauntile-calculated tolerances', 'Infected'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
ylabel('Effective reproduction value')
xlabel('Time (date)')
enddate = dates(end);
xlim([StartWeek,enddate])
xtks = [datenum(2018,1:2:11,1),datenum(2019,1:2:11,1),datenum(2020,1:2:7,1)];
xticks(xtks)
datetick('x',1,'keepticks')
sgtitle('Effective Reproduction value for rat-to-rat transmission')
%% R_0 above 1 
%StartdateNum = datenum(StarDate,'dd-mm-yyyy');
Year = datenum('01-01-2018'); % What year our data starts
Offset = StartWeek-Year; % How much we need to offset what phi = 0 means in real dates

figure
set(gcf,'Position',[0 0 1000 1000])

R_0ind_adj = cell(1,1);
R_0z_m = cell(1,1);
RiM2 = NaN(NoPart,1);
Result_prms_gen_T = All_gen_results_cell{1,T};
if ~isempty(Result_prms_gen_T)
    ht = height(Result_prms_gen_T);
    RiM = NaN(2*ht,2);
    for i_particle = 1:ht            
        RiM(2*i_particle-1,1) = R_0_indicator(i_particle,1)+Offset;  
        RiM(2*i_particle,1) = R_0_indicator(i_particle,3)+Offset-365; 
        RiM(2*i_particle-1,2) = R_0_indicator(i_particle,2)+Offset;  
        RiM(2*i_particle,2) = R_0_indicator(i_particle,4)+Offset-365;
    end

    R_0ind_adj =RiM;
    R_0z_m = NaN(3,2);
    for p =1:2
        Res_Sortedj = sort(RiM(:,p));
        R_0z_m(1,p) = Res_Sortedj(floor(0.05*ht));
        R_0z_m(2,p) = Res_Sortedj(floor(0.5*ht));
        R_0z_m(3,p) = Res_Sortedj(floor(0.95*ht));
    end
end 
bx = boxplot(R_0ind_adj+StartWeek,'orientation','horizontal');
set(bx,{'linew'},{2})
set(0, 'DefaultAxesLineWidth', 0.5)
%bx = findobj('Tag','boxplot');
%set(bx.Children,'LineWidth',3)
datetick('x','dd-mmmm')
xlabel('Date')
yticklabels({'Start of rat epidemic','End of rat epidemic'})
xticks([datenum('10-01-2018'),datenum('10-15-2018'),datenum('11-01-2018'),datenum('11-15-2018'),datenum('12-01-2018'),datenum('12-15-2018'),datenum('01-01-2019'),datenum('01-15-2019')])
xticklabels({'01-October','15-October','01-November','15-November','01-December','15-December','01-January','15-January'})
set(gca,'fontsize',14)
title({'Credible Intervals for start and end of period', 'which R^{rr}(t) \geq 1'})
%% Confidence/Credible intervals


z = 0.95;
range_arr = [floor(z*NoPart),floor(0.5*NoPart),NoPart-floor(z*NoPart)];
z_Interval_median = Result_prms_gen_T(1:length(range_arr)+1,1:end-1);
z_Interval_median.Properties.RowNames ={'z','Median','1-z','Best Fit'};

for i = 1:p_set_dim(1)
    % Parameter for Best fit, sorted by error
    var_i = Priors_given_m.Properties.VariableNames{i};
    z_Interval_median(4,var_i) = Result_prms_gen_T(1,var_i);
    % Parameters for range and median, sorted by that parameter
    for j = 1:length(range_arr)
        Res_Sortedj = sort(Result_prms_gen_T{:,var_i});
        z_Interval_median{j,var_i} = Res_Sortedj(range_arr(j));
    end
end

%% If we have a parameter related to period in the year
StartDateStr = datestr(StartDate);
Year = datenum(StartDateStr(end-3:end),'yyyy'); % What year our data starts
Offset = StartDate-Year; % How much we need to offset what phi = 0 means in real dates
dateAdj_p_set = 365*z_Interval_median.phi+Offset;
dateAdj_p_set = mod(dateAdj_p_set,365)+Year;
datestr(dateAdj_p_set)
dateAdjMax = dateAdj_p_set+182.5; % In the case of my models, predominantly phi is 
% the date of the minimum so here is the maximum
datestr(dateAdjMax)
%% Interquantile range, median and best fit for R_0 values
R_0_table = z_Interval_median;
R_0_table.R_rr = NaN(length(range_arr)+1,1);
R_0_table.R_hh = NaN(length(range_arr)+1,1);
R_0_table.R_rh = NaN(length(range_arr)+1,1);
R_0_table = R_0_table(:,end-2:end);
for i = 1:length(range_arr)
    R_0_table.R_rh(end) = R_rh_0vec(1);
    R_0_table.R_rr(end) = R_rr_0vec(1);
    R_0_table.R_hh(end) = R_hh_0vec(1);
    sort_R_rr = sort(R_rr_0vec);
    sort_R_hh = sort(R_hh_0vec);
    sort_R_rh = sort(R_rh_0vec);
    R_0_table.R_rr(1:end-1) = sort_R_rr(range_arr);
    R_0_table.R_hh(1:end-1) = sort_R_hh(range_arr);
    R_0_table.R_rh(1:end-1) = sort_R_rh(range_arr);
    
end
