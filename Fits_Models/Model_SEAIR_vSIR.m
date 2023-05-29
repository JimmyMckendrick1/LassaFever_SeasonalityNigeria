function  Output_T = Model_SEAIR_vSIR_Cfreq_Bpg_f(input)

if nargin == 0
% parameters
    % Rats
    s = 1105.90318393293;  mu_r=1/500; beta_rr=18.0399391317874;%0.30830836377397;
    gamma_r=1/90; phi = 0.438661320162671;
    N_r_0=10^6; I_r_0=NaN; 
    % Humans
    B_h = 1/(53.5*365)+1e-4; beta_hh = 0.000106113392996808; 
    beta_rh = 6.27209803061645e-05;
    p = 0.8; nu = 1/14; gamma_h = 1/14; mu_h_I = 0.017857142857143; 
    mu_h = 1/(53.5*365); sigma = 1;
    N_h_0 = 2e8; S_h_0 = 199999960; E_h_0 = 30; A_h_0 = 8; I_h_0 = 2; 
    C_h_0 = 2; D_h_0 = 0;
    MaxTime=917;
    s =31.2437999646218;
    beta_rr = 4.18653227717189;
    beta_hh = 0.0147131268850148;
    beta_rh = 0.000367378880121168;
    phi = 0.679520808465325;
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
S_h = S_h_0; E_h = E_h_0; A_h = A_h_0; I_h=I_h_0; R_h=N_h_0-S_h-E_h-A_h-I_h;
C_h = C_h_0; D_h = D_h_0;
% The main iteration 
options = odeset('AbsTol', 1e-5);
input_vec = [k s phi beta_rr gamma_r mu_r B_h beta_hh sigma beta_rh p nu gamma_h mu_h_I mu_h];
[t, pop]=ode45(@Diff_2_2,0:MaxTime,[S_r I_r R_r S_h E_h A_h I_h R_h C_h D_h],options,input_vec);



[S_r,I_r,R_r,S_h,E_h,A_h,I_h,R_h,C_h,D_h]=matsplit(pop,1);
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