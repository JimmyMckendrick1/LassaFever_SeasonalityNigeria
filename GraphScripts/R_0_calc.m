function Resp = R_0_calc(cell_args)
tmp_input = cell_args{1};
Output_T = cell_args{2};
sigma = 1;
S_r = Output_T.S_r;
N_r = Output_T.S_r+Output_T.I_r+Output_T.R_r;
S_h = Output_T.S_h;
N_h = Output_T.S_h+Output_T.E_h+Output_T.A_h+Output_T.I_h+Output_T.R_h; 
beta_rh = zeros(1,tmp_input.MaxTime+1);
if tmp_input.d_start > tmp_input.d_end
    tmp_input.d_start = tmp_input.d_start-1;
end
for t = 0:tmp_input.MaxTime
    if (mod(t,365) > 365*(tmp_input.d_start) && mod(t,365) <= 365*(tmp_input.d_end))...
            || (mod(t,365) > 365*(tmp_input.d_start+1) && mod(t,365) <= 365*(tmp_input.d_end+1)) 
        beta_rh(t+1) = tmp_input.d_mult*tmp_input.beta_rh;
    else
        beta_rh(t+1) = tmp_input.beta_rh;
    end
end
R_rr = tmp_input.beta_rr*S_r./((tmp_input.mu_r+tmp_input.gamma_r)*N_r);

R_rh = beta_rh'.*S_h.*Output_T.I_r./((tmp_input.mu_r+tmp_input.gamma_r)*N_h);

R_hh = tmp_input.beta_hh*(tmp_input.p*tmp_input.nu*sigma/((tmp_input.mu_h+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu))+...
(1-tmp_input.p)*tmp_input.nu/((tmp_input.mu_h+tmp_input.mu_h_I+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu)))*S_h./N_h;

R_rr0 = tmp_input.beta_rr/(tmp_input.mu_r+tmp_input.gamma_r);

R_rh0 = tmp_input.beta_rh/(tmp_input.mu_r+tmp_input.gamma_r);

R_hh0 = tmp_input.beta_hh*(tmp_input.p*tmp_input.nu*sigma/((tmp_input.mu_h+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu))+...
(1-tmp_input.p)*tmp_input.nu/((tmp_input.mu_h+tmp_input.mu_h_I+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu)));

Resp =[R_rr,R_rh,R_hh];
Resp = [Resp;[R_rr0,R_rh0,R_hh0]];
end