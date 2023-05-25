function Output_Cell = Model_rat_SIR_Freq(input)

if nargin == 0
% parameters
    s = 10;%18.583901319727340;
    phi = 0.479210428417469;%0.499156424319793;
    beta=0.1;%1.782978817928539;%1.831043660491785;
    gamma=1/90;  mu=1/500; 
    N0=1e6; S0=N0*(gamma+mu)/beta;  
    MaxTime=900;
    I0=NaN;
    input = {s,phi,beta,gamma,mu,N0,I0,MaxTime};
end
% For the purposes of modelling lassa fever this is a choice of rat subsystem
% The model employs a simple SIR system for the disease dynamics, with a 
% Periodic Gaussian function forcing the birth rate, with parameters k,s
% and phi. Births are assumed to be susceptible. For the purposes of LF,
% immunity is forever, deaths are constant
% Transmission rate is frequency dependent

[s,phi,beta,gamma,mu,N0,I0,MaxTime] = deal(input{:});
k = mu/PeriodicGaussian_normalisation(s/2);
R_rr = beta/(gamma+mu);


if isequal(N0,1)
    S0 = max(min(N0/R_rr,N0-1e-4),0);       
else
    S0 = max(min(N0/R_rr,N0-1),0);
end

if isnan(I0)
    if isequal(N0,1)   
        I0 = min(N0-S0,max(N0*mu*(R_rr-1)/beta,1e-4));
    else
        I0 = min(N0-S0,max(N0*mu*(R_rr-1)/beta,1));
    end
end

% Checks all the parameters are valid
if S0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S0);
end

if I0<0 
    error('Initial level of infecteds (%g) is less than zero',I0);
end

if I0>N0 
  error('Initial level of infecteds (%g) is more than the population',I0);
end

if beta<0 
    error('Infection rate (%g) is less than or equal to zero',beta);
end
   
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end

if mu<0 
    error('Birth / Death rate gamma (%g) is less than zero',mu);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end

if s <0
    error('Synchrony (%g) is less than zero',s);
end

if k < 0
    error('Birthing magnitude (%g) is less than zero',k);
end

    
if S0+I0>N0
    warning('Initial level of susceptibles+exposed+infecteds (%g+%g+%g=%g) is greater than total population',S0,I0,N0);
end


S=S0; I=I0; R=N0-S-I;
%n = ceil(MaxTime/365);
% The main iteration 
options = odeset('AbsTol', 1e-5);
[t, pop]=ode23(@Diff_2_2,0:1:MaxTime,[S I R],options,[k s phi beta gamma mu]);


S = pop(:,1); I=pop(:,2); R=pop(:,3); 
Output_Cell = table(t,S,I,R);
%{
t = Output_Cell.t;
S = Output_Cell.S;
I = Output_Cell.I;
R = Output_Cell.R;
N = S+I+R;
t = t+datenum('14-12-2018','dd-mm-yyyy');
figure
startdate = datenum('14-12-2018','dd-mm-yyyy');
enddate = datenum('14-12-2018','dd-mm-yyyy')+900;%datenum('19-07-2020','dd-mm-yyyy');
xtks = [datenum(2018,12,1),datenum(2019,2:2:12,1),datenum(2020,2:2:7,1)];
subplot(3,1,1)
plot(t,S,t,I,t,R,t,N)
xlim([startdate,enddate])
xticks(xtks)
datetick('x','mmm-yyyy','keepticks')
legend("Susceptible","Infected","Recovered","Total",'Location','Best')
subplot(3,1,2)
plot(t,R./(S+I+R))
xticks(xtks);
datetick('x','mmm-yyyy','keepticks')
subplot(3,1,3)
plot(t,I./(S+I+R));
xticks(xtks);
datetick('x','mmm-yyyy','keepticks')
%}
end

function dPop=Diff_2_2(t,pop, parameter)
k = parameter(1); s=parameter(2);
phi=parameter(3); beta=parameter(4); 
gamma=parameter(5); mu=parameter(6);

    
B=k*exp(-s*cos(pi*(t/365-phi))^2);
S=pop(1); I=pop(2); R=pop(3);
N=S+I+R;

dPop=zeros(3,1);
dPop(1)= B*N - beta*S*I/N - mu*S; %S
dPop(2)= beta*S*I/N - (gamma+mu)*I; %I
dPop(3)= gamma*I - mu*R; %R
end