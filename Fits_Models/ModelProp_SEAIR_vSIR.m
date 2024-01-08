function  Output_Cell = ModelProp_SEAIR_vSIR(input)
% Function to see the proportions of humans that are infected by other
% humans and those that were infected by rats when running Model_SEAIR_vSIR
% This function takes additional inputs compared to Modle_SEAIR_vSIR as it
% differentiates between humans infected by rats and those by humans.
if nargin == 0
% parameters
    % Rats 
    s =100;
    beta_rr = 4;
    beta_rh = 11;
    phi = 0.4;  mu_r=0.0038; 
    gamma_r=1/90; N_r_0=1; S_r_0 = NaN; I_r_0=NaN; 
    % Humans
    B_h = 1/(53.5*365)+1e-4;  
    beta_hh = 0.01;
    p = 0.8; nu = 1/14; gamma_h = 1/14; mu_h_I = 0.018; 
    mu_h = 1/(53.5*365); sigma = 1; 
    N_h_0 = 2e8; S_h_0 = 0.7*N_h_0; E_hh_0 = 30; E_hr_0 = 0; A_h_0 = 8; 
    I_hh_0 = 2; I_hr_0 = 0; C_hh_0 = 2; C_hr_0 = 0; D_h_0 = 0;
    d_start = 305/365; d_end = 90/365; d_mult = 2;
    MaxTime=917;
    input = {s,phi,beta_rr,gamma_r,mu_r,N_r_0,S_r_0,I_r_0,B_h,...
    beta_hh,sigma,beta_rh,beta_hr,d_start,d_end,d_mult,p,nu,gamma_h,mu_h_I,mu_h,...
    N_h_0,S_h_0,E_hh_0,A_h_0,I_hh_0,C_hh_0,D_h_0,MaxTime,E_hr_0,I_hr_0,C_hr_0};
end

[s,phi,beta_rr,gamma_r,mu_r,N_r_0,S_r_0,I_r_0,B_h,...
    beta_hh,sigma,beta_rh,beta_hr,d_start,d_end,d_mult,p,nu,gamma_h,mu_h_I,mu_h,...
    N_h_0,S_h_0,E_hh_0,A_h_0,I_hh_0,C_hh_0,D_h_0,MaxTime,E_hr_0,I_hr_0,C_hr_0]=deal(input{:});

k = mu_r/PeriodicGaussian_normalisation(s/2);

if d_start > d_end
    d_start = d_start-1;
end
if isnan(S_r_0)
    R_rr = beta_rr/(gamma_r+mu_r);
    Bt = k*exp(-s*cos(pi*(-phi))^2);
    %Bt = mu_r;
    if isnan(I_r_0)
        if isequal(N_r_0,1) 
            I_r_0 = min(N_r_0-1e-4,max(N_r_0*(R_rr*Bt-mu_r)/beta_rr,1e-4));
            S_r_0 = max(min(N_r_0/R_rr,N_r_0-I_r_0),0); 
        else 
            I_r_0 = min(N_r_0-1,max(N_r_0*(R_rr*Bt-mu_r)/beta_rr,1));
            S_r_0 = max(min(N_r_0/R_rr,N_r_0-I_r_0),0); 
        end    
    else
        S_r_0 = max(min(N_r_0/R_rr,N_r_0-I_r_0),0);    
    end
elseif isnan(I_r_0)
    R_rr = beta_rr/(gamma_r+mu_r); 
    Bt = k*exp(-s*cos(pi*(-phi))^2);
    if isequal(N_r_0,1) 
        I_r_0 = min(N_r_0-S_r_0,max(N_r_0*(R_rr*Bt-mu_r)/beta_rr,1e-4));
    else 
        I_r_0 = min(N_r_0-S_r_0,max(N_r_0*(R_rr*Bt-mu_r)/beta_rr,1));
    end   
end
R_r_0 = N_r_0-S_r_0-I_r_0;
R_h_0 = N_h_0-S_h_0-E_hh_0-E_hr_0-A_h_0-I_hh_0-I_hr_0;


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
input_vec = [k s phi beta_rr gamma_r mu_r B_h beta_hh sigma beta_rh p nu gamma_h mu_h_I mu_h d_start d_end d_mult];
[t, pop]=ode45(@Diff_2_2,[0:MaxTime],[S_r_0 I_r_0 R_r_0 S_h_0 E_hr_0 E_hh_0 A_h_0 I_hr_0 I_hh_0 R_h_0 C_hr_0 C_hh_0 D_h_0],...
    options,input_vec);

A_cell = num2cell(pop,1);
[S_r,I_r,R_r,S_h,E_h_r,E_h_h,A_h,I_hr,I_hh,R_h,C_hr,C_hh,D_h]=A_cell{:};
Output_Cell = table(t,S_r,I_r,R_r,S_h,E_h_r,E_h_h,A_h,I_hr,I_hh,R_h,C_hr,C_hh,D_h );

%{
figure
t = Output_Cell{1,1};
plot(t,S_r,t,I_r,t,R_r)
title("Rat dynamics")
legend("Susceptible","Infected","Recovered",'Location','Best')
%}
%{
figure
plot(t,E_hh,t,E_hr)
legend('human-to-human','rat-to-human')
figure
plot(t,C_hh,t,C_hr)
legend('human-to-human','rat-to-human')
%}
end

function dPop=Diff_2_2(t,pop, parameter)
vart2 = num2cell(parameter);
[k,s,phi,beta_rr, gamma_r, mu_r, B_h, beta_hh, sigma, beta_rh, p, nu,...
    gamma_h,mu_h_I,mu_h,d_start, d_end, d_mult] = deal(vart2{:});


B = k*exp(-s*cos(pi*(t/365-phi))^2);
S_r=pop(1); I_r=pop(2); R_r=pop(3);
N_r=S_r+I_r+R_r;

S_h=pop(4); E_hr=pop(5); E_hh=pop(6); A_h=pop(7); I_hr=pop(8);
I_hh = pop(9); R_h=pop(10);
N_h = S_h+E_hr+E_hh+A_h+I_hr+I_hh+R_h;

dPop=zeros(13,1);
if (mod(t,365) > 365*(d_start) && mod(t,365) <= 365*(d_end)) || (mod(t,365) > 365*(d_start+1) && mod(t,365) <= 365*(d_end+1)) 
    beta_rh = d_mult*beta_rh;
end

dPop(1)= B*N_r - beta_rr*S_r*I_r/N_r - mu_r*S_r; %S_r
dPop(2)= beta_rr*S_r*I_r/N_r - gamma_r*I_r - mu_r*I_r; %I_r
dPop(3)= gamma_r*I_r - mu_r*R_r; %R_r

dPop(4) = B_h*N_h - (beta_rh*I_r + beta_hh*(sigma*A_h+I_hr+I_hh))*S_h/N_h - mu_h*S_h; %S_h
dPop(5) = (beta_rh*I_r )*S_h/N_h - (nu+mu_h)*E_hr; %E_hr humans infected by rats, incubation
dPop(6) = beta_hh*(sigma*A_h+I_hr+I_hh)*S_h/N_h - (nu+mu_h)*E_hh; %E_hh humans infected by humans, incubation
dPop(7) = p*nu*(E_hr+E_hh) - (gamma_h+mu_h)*A_h; %A_h
dPop(8) = (1-p)*nu*(E_hr) - (gamma_h+mu_h_I+mu_h)*I_hr; %I_hr humans infected by rats, infectious
dPop(9) = (1-p)*nu*(E_hh) - (gamma_h+mu_h_I+mu_h)*I_hh; %I_hh humans infected by humans, infectious
dPop(10) = gamma_h*(A_h+I_hr +I_hh) - mu_h*R_h; %R_h
dPop(11) = (1-p)*nu*E_hr; %C_h_r
dPop(12) = (1-p)*nu*E_hh; %C_h_h
dPop(13) = mu_h_I*(I_hr+I_hh); %D_h
end