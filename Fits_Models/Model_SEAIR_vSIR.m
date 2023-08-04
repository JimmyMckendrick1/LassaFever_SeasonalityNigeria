function  Output_T = Model_SEAIR_vSIR(input)
% Main model developed in Modelling seasonality in Lassa fever in Nigeria.
% It is titled following the convention of <human system>_v<vector system>,
% which is an susceptible-exposed-asymptomatic/infectious-recovered model
% for humans and an susceptible-infectious-recovered model for the vector,
% Mastomys natalensis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function input: input is a cell containing the following ordered parameters:
% s: Shape parameter for the vector birth/recruitment function 
% phi: Time at which the rate of reproduction is at its minimum (0<=phi<=1)
% beta_rr: Transmission rate, rat-to-rats
% gamma_r: Recovery rate of rats
% mu_r: Natural mortality rate of rats
% N_r_0: Initial number of total rats
% I_r_0: Initial number of infected rats, if left as NaN then this is
% calculated
% B_h: Rate of births for humans
% beta_hh: Transmission rate, human-to-human
% sigma: Relative transmission rate ratio of asymptomatic humans (A_h) to
% symptomatic humans (I_h). Left as 1 for paper
% beta_rh: Transmission rate, rat-to-humans
% p: Proportion of humans that are asymptomatic
% nu: Rate of progression of humans from exposed to an infectious compartment
% gamma_h: Recovery rate of humans
% mu_h_I: Infection induced mortality rate for humans 
% mu_h: Natural mortality rate for humans
% N_h_0: Initial number of total humans
% S_h_0: Initial number of susceptible humans
% E_h_0: Initial number of exposed humans
% A_h_0: Initial number of asymptomatic humans
% I_h_0: Initial number of infected humans
% C_h_0: Initial number of cases
% D_h_0: Initial number of deaths
% MaxTime: Time that the model will run to in days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function output:
% Output_T:
if nargin == 0
% parameters
    % Rats    
    s =100;
    beta_rr = 4;
    beta_rh = 0.0003;
    phi = 0.4;  mu_r=1/500; 
    gamma_r=1/90; N_r_0=10^6; I_r_0=NaN; 
    % Humans
    B_h = 1/(53.5*365)+1e-4;  
    beta_hh = 0.01;
    p = 0.8; nu = 1/14; gamma_h = 1/14; mu_h_I = 0.018; 
    mu_h = 1/(53.5*365); sigma = 1; 
    N_h_0 = 2e8; S_h_0 = 199999960; E_h_0 = 30; A_h_0 = 8; I_h_0 = 2; 
    C_h_0 = 2; D_h_0 = 0;
    MaxTime=917;
    input = {s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
        beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,N_h_0,...
        S_h_0,E_h_0,A_h_0,I_h_0,C_h_0,D_h_0,MaxTime};
end

[s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
    beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,...
    N_h_0,S_h_0,E_h_0,A_h_0,I_h_0,C_h_0,D_h_0,MaxTime]=deal(input{:});

k = mu_r/PeriodicGaussian_normalisation(s/2);
R_rr = beta_rr/(gamma_r+mu_r);

fN = N_r_0;

if isnan(I_r_0)
    if isequal(N_r_0,1) 
        S_r_0 = max(min(fN/R_rr,N_r_0-1e-4),0); 
        I_r_0 = min(N_r_0-S_r_0,max(N_r_0*mu_r*(R_rr-1)/beta_rr,1e-4));
    else   
        S_r_0 = max(min(fN/R_rr,N_r_0-1),0);  
        I_r_0 = min(N_r_0-S_r_0,max(N_r_0*mu_r*(R_rr-1)/beta_rr,1));
    end    
else
    S_r_0 = max(min(fN/R_rr,N_r_0-I_r_0),0);    
end
R_r_0 = N_r_0-S_r_0-I_r_0;
R_h_0 = N_h_0-S_h_0-E_h_0-A_h_0-I_h_0;

% Checks all the parameters are valid
if S_r_0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S_r_0);
end

if I_r_0<=0 
    error('Initial level of infecteds (%g) is less than or equal to zero',I_r_0);
end

if I_r_0>=N_r_0 
  error('Initial level of susceptibles (%g) is more than or equal to the population',I_r_0);
end

if beta_rr<0 
    error('Rate of infection from exposure nu (%g) is less than or equal to zero',beta_rr);
end
   
if gamma_r<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma_r);
end

if mu_r<0 
    error('Birth / Death rate gamma (%g) is less than zero',mu_r);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end
    
if S_r_0+I_r_0>N_r_0
    warning('Initial level of susceptibles+exposed+infecteds (%g+%g+%g=%g) is greater than total population',S_r_0,I_r_0,N_r_0);
end

% The main iteration 
options = odeset('AbsTol', 1e-5);
input_vec = [k s phi beta_rr gamma_r mu_r B_h beta_hh sigma beta_rh p nu gamma_h mu_h_I mu_h];
[t, pop]=ode45(@Diff_2_2,0:MaxTime,[S_r_0 I_r_0 R_r_0 S_h_0 E_h_0 A_h_0 I_h_0 R_h_0 C_h_0 D_h_0],options,input_vec);




A_cell = num2cell(pop,1);
[S_r,I_r,R_r,S_h,E_h,A_h,I_h,R_h,C_h,D_h] = A_cell{:};
Output_T = table(t,S_r,I_r,R_r,S_h,E_h,A_h,I_h,R_h,C_h,D_h);
%{

figure
t = Output_T.t;
S_r = Output_T.S_r;
I_r = Output_T.I_r;
R_r = Output_T.R_r;
subplot(3,1,1)
plot(t,S_r)
subplot(3,1,2)
plot(t,I_r)
subplot(3,1,3)
plot( t, (gamma_r+mu_r)*(S_r+I_r+R_r)/beta_rr)
title("Rat dynamics")
legend("Susceptible","Theoretic",'Location','Best')

figure
S_h = Output_T.S_h;
E_h = Output_T.E_h;
A_h = Output_T.A_h;
I_h = Output_T.I_h;
R_h = Output_T.R_h;
plot(t,S_h,t,E_h,t,A_h,t,I_h,t,R_h)
title("Human dynamics")
subplot(3,1,1)
plot(t,I_h)
legend("Infected")
subplot(3,1,2)
plot(t,S_h)
legend("Susceptible")
subplot(3,1,3)
plot(t,R_h)
legend("Recovered")
%legend("Susceptible","Exposed","Asymptamatics","Infected","Recovered",'Location','Best')
%}
end

function dPop=Diff_2_2(t,pop, parameter)
vart2 = num2cell(parameter);
[k,s,phi,beta_rr, gamma_r, mu_r, B_h, beta_hh, sigma, beta_rh, p, nu, gamma_h,mu_h_I,mu_h] = deal(vart2{:});

dPop=zeros(10,1);

B = k*exp(-s*cos(pi*(t/365-phi))^2);
S_r=pop(1); I_r=pop(2); R_r=pop(3);
N_r=S_r+I_r+R_r;

S_h=pop(4); E_h=pop(5); A_h=pop(6); I_h=pop(7); R_h=pop(8);
N_h = S_h+E_h+A_h+I_h+R_h;

dPop(1)= B*N_r - beta_rr*S_r*I_r/N_r - mu_r*S_r; %S_r
dPop(2)= beta_rr*S_r*I_r/N_r - (gamma_r+ mu_r)*I_r; %I_r
dPop(3)= gamma_r*I_r - mu_r*R_r; %R_r

dPop(4) = B_h*N_h - (beta_rh*I_r + beta_hh*(sigma*A_h+I_h))*S_h/N_h - mu_h*S_h; %S_h
dPop(5) = (beta_rh*I_r + beta_hh*(sigma*A_h+I_h))*S_h/N_h - (nu+mu_h)*E_h; %E_h
dPop(6) = p*nu*E_h - (gamma_h+mu_h)*A_h; %A_h
dPop(7) = (1-p)*nu*E_h - (gamma_h+mu_h_I+mu_h)*I_h; %I_h
dPop(8) = gamma_h*(A_h+I_h) - mu_h*R_h; %R_h
dPop(9) = (1-p)*nu*E_h; %C_h
dPop(10) = mu_h_I*I_h; %D_h
end