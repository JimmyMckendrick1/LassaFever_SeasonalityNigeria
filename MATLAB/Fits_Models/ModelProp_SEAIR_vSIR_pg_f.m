function  Output_Cell = ModelProp_SEAIR_vSIR_pg_f(input)
% function to see the proportions of humans that are infected by other
% humans and those that were infected by rats when running Model_SEAIR_vSIR_pg_f
if nargin == 0
% parameters
    % Rats 
    s = 989.716296086384;  mu_r=1/500; beta_rr=1.90909122362604;
    gamma_r=1/90; phi = 0.479938650962018;
    N_r_0=1e5; I_r_0=0.5*N_r_0; 
    % Humans
    B_h = 1/(53.5*365)+1e-4; beta_hh = 0.00130994248902666; sigma=1;
    beta_rh = 0.00746899555324877;
    p = 0.8; nu = 1/14; gamma_h = 1/14; mu_h_I = 0.017857142857143; 
    mu_h = 1/(53.5*365);
    N_h_0 = 1e6; S_h_0 = 1e5; E_hh_0 = 10; E_hr_0 = 0;
    A_h_0 = 20; I_hr_0 = 0; I_hh_0 = 2; 
    C_hh_0 = 0; C_hr_0 = 2; D_h_0 = 0;
    MaxTime=576;
    input = {s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
        beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,N_h_0,S_h_0,E_hh_0,E_hr_0,...
        A_h_0,I_hh_0,I_hr_0,C_hh_0,C_hr_0,D_h_0,MaxTime};
end

[s,phi,beta_rr,gamma_r,mu_r,N_r_0,I_r_0,B_h,...
    beta_hh,sigma,beta_rh,p,nu,gamma_h,mu_h_I,mu_h,...
    N_h_0,S_h_0,E_hr_0,A_h_0,I_hr_0,C_hr_0,D_h_0,MaxTime,C_hh_0,E_hh_0,I_hh_0]=deal(input{:});

k = mu_r/PeriodicGaussian_normalisation(s/2);
R_rr = beta_rr/(gamma_r+mu_r);


if isequal(N_r_0,1)
    S_r_0 = max(min(N_r_0/R_rr,N_r_0-1e-4),0);       
else
    S_r_0 = max(min(N_r_0/R_rr,N_r_0-1),0);    
end


if isnan(I_r_0)
    if isequal(N_r_0,1)   
        I_r_0 = min(N_r_0-S_r_0,max(N_r_0*mu_r*(R_rr-1)/beta_rr,1e-4));
    else   
        I_r_0 = min(N_r_0-S_r_0,max(N_r_0*mu_r*(R_rr-1)/beta_rr,1));
    end
end

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


S_r=S_r_0; I_r=I_r_0; R_r=N_r_0-S_r-I_r;
S_h = S_h_0; E_h_r = E_hr_0; E_h_h = E_hh_0; A_h = A_h_0;
I_hr=I_hr_0; I_hh=I_hh_0; R_h=N_h_0-S_h-E_h_r-A_h-I_hh;
C_h_r = C_hr_0; C_h_h = C_hh_0; D_h = D_h_0; 
% The main iteration 
options = odeset('AbsTol', 1e-5);
input_vec = [k s phi beta_rr gamma_r mu_r B_h beta_hh sigma beta_rh p nu gamma_h mu_h_I mu_h];
[t, pop]=ode45(@Diff_2_2,[0:MaxTime],[S_r I_r R_r S_h E_h_r E_h_h A_h I_hr I_hh R_h C_h_r C_h_h D_h],...
    options,input_vec);

[S_r,I_r,R_r,S_h,E_h_r,E_h_h,A_h,I_hr,I_hh,R_h,C_h_r,C_h_h,D_h]=matsplit(pop,1);
Output_Cell = table(t,S_r,I_r,R_r,S_h,E_h_r,E_h_h,A_h,I_hr,I_hh,R_h,C_h_r,C_h_h,D_h );

%{
figure
t = Output_Cell{1,1};
plot(t,S_r,t,I_r,t,R_r)
title("Rat dynamics")
legend("Susceptible","Infected","Recovered",'Location','Best')

figure
plot(t,S_h,t,E_h_h+E_h_r,t,I_h,t,R_h)
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
legend("Susceptible","Exposed","Asymptamatics","Infected","Recovered",'Location','Best')
%}
    %{
figure
plot(t,E_h_h,t,E_h_r)
legend('human-to-human','rat-to-human')
figure
plot(t,C_h_h,t,C_h_r)
legend('human-to-human','rat-to-human')
%}
end

function dPop=Diff_2_2(t,pop, parameter)
vart2 = num2cell(parameter);
[k,s,phi,beta_rr, gamma_r, mu_r, B_h, beta_hh, sigma, beta_rh, p, nu, gamma_h,mu_h_I,mu_h] = deal(vart2{:});


B = k*exp(-s*cos(pi*(t/365-phi))^2);
S_r=pop(1); I_r=pop(2); R_r=pop(3);
N_r=S_r+I_r+R_r;

S_h=pop(4); E_h_r=pop(5); E_h_h=pop(6); A_h=pop(7); I_h_r=pop(8);
I_h_h = pop(9); R_h=pop(10);
N_h = S_h+E_h_r+E_h_h+A_h+I_h_r+I_h_h+R_h;
dPop=zeros(13,1);

dPop(1)= B*N_r - beta_rr*S_r*I_r/N_r - mu_r*S_r; %S_r
dPop(2)= beta_rr*S_r*I_r/N_r - gamma_r*I_r - mu_r*I_r; %I_r
dPop(3)= gamma_r*I_r - mu_r*R_r; %R_r

dPop(4) = B_h*N_h - (beta_rh*I_r + beta_hh*(sigma*A_h+I_h_r+I_h_h))*S_h/N_h - mu_h*S_h; %S_h
dPop(5) = (beta_rh*I_r )*S_h/N_h - (nu+mu_h)*E_h_r; %E_h_r
dPop(6) = beta_hh*(sigma*A_h+I_h_r+I_h_h)*S_h/N_h - (nu+mu_h)*E_h_h; %E_h_h 
dPop(7) = p*nu*(E_h_r+E_h_h) - (gamma_h+mu_h)*A_h; %A_h
dPop(8) = (1-p)*nu*(E_h_r) - (gamma_h+mu_h_I+mu_h)*I_h_r; %I_h_r
dPop(9) = (1-p)*nu*(E_h_h) - (gamma_h+mu_h_I+mu_h)*I_h_h; %I_h_h
dPop(10) = gamma_h*(A_h+I_h_r +I_h_h) - mu_h*R_h; %R_h
dPop(11) = (1-p)*nu*E_h_r; %C_h_r
dPop(12) = (1-p)*nu*E_h_h; %C_h_h
dPop(13) = mu_h_I*(I_h_r+I_h_h); %D_h
end