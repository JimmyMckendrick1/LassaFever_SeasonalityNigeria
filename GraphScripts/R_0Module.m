function resp = R_0Module(str_function, cell_args)
% This function is a python-like module, containing sub functions of the
% Approximate Bayesian Computation algorithms that I have for fitting my
% epidemiological models.
allowed_functions = {'R_0freq', 'R_0dens', 'R_0pwr', 'R_0sigm';
                     @R_0freq, @R_0dens,@R_0pwr,@R_0sigm};

[found, where] = ismember(str_function, allowed_functions(1, :));  %search 1st row of lookup table
if ~found
    error('Function %s is not a valid function');
end
%if isequal(where,3)
%    allowed_functions{2, where}(cell_args);
%else
%    resp = allowed_functions{2, where}(cell_args); %use matching 2nd row of lookup table
%end 
resp = allowed_functions{2, where}(cell_args);
end

function Resp = R_0freq(cell_args)
tmp_input = cell_args{1};
Output_T = cell_args{2};
sigma = 1;
S_r = Output_T.S_r;
N_r = Output_T.S_r+Output_T.I_r+Output_T.R_r;
S_h = Output_T.S_h;
N_h = Output_T.S_h+Output_T.E_h+Output_T.A_h+Output_T.I_h+Output_T.R_h; 

R_rr = tmp_input.beta_rr*S_r./((tmp_input.mu_r+tmp_input.gamma_r)*N_r);

R_rh = tmp_input.beta_rh*S_h./((tmp_input.mu_r+tmp_input.gamma_r)*N_h);

R_hh = tmp_input.beta_hh*(tmp_input.p*tmp_input.nu*sigma/((tmp_input.mu_h+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu))+...
(1-tmp_input.p)*tmp_input.nu/((tmp_input.mu_h+tmp_input.mu_h_I+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu)))*S_h./N_h;

R_rr = [R_rr;tmp_input.beta_rr/(tmp_input.mu_r+tmp_input.gamma_r)];

R_rh = [R_rh;tmp_input.beta_rh/(tmp_input.mu_r+tmp_input.gamma_r)];

R_hh = [R_hh;tmp_input.beta_hh*(tmp_input.p*tmp_input.nu*sigma/((tmp_input.mu_h+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu))+...
(1-tmp_input.p)*tmp_input.nu/((tmp_input.mu_h+tmp_input.mu_h_I+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu)))];

Resp =[R_rr,R_rh,R_hh];
end

function Resp= R_0dens(cell_args)
tmp_input = cell_args{1};
Output_T = cell_args{2};


end

function Resp = R_0pwr(cell_args)
tmp_input = cell_args{1};
Output_T = cell_args{2};


end

function Resp = R_0sigm(cell_args)
tmp_input = cell_args{1};
Output_T = cell_args{2};
sigma=1;
S_r = Output_T.S_r;
N_r = Output_T.S_r+Output_T.I_r+Output_T.R_r;
S_h = Output_T.S_h;
N_h = Output_T.S_h+Output_T.E_h+Output_T.A_h+Output_T.I_h+Output_T.R_h; 

R_rr = tmp_input.beta_rr*S_r./((tmp_input.mu_r+tmp_input.gamma_r)*(N_r.*(1+exp(-tmp_input.a_1*(N_r-tmp_input.a_2)))));

R_rh = tmp_input.beta_rh*S_h./((tmp_input.mu_r+tmp_input.gamma_r)*N_h);

R_hh = tmp_input.beta_hh*(tmp_input.p*tmp_input.nu*sigma/((tmp_input.mu_h+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu))+...
(1-tmp_input.p)*tmp_input.nu/((tmp_input.mu_h+tmp_input.mu_h_I+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu)))*S_h./N_h;

R_rr = [R_rr;tmp_input.beta_rr/(tmp_input.mu_r+tmp_input.gamma_r)];

R_rh = [R_rh;tmp_input.beta_rh/(tmp_input.mu_r+tmp_input.gamma_r)];

R_hh = [R_hh;tmp_input.beta_hh*(tmp_input.p*tmp_input.nu*sigma/((tmp_input.mu_h+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu))+...
(1-tmp_input.p)*tmp_input.nu/((tmp_input.mu_h+tmp_input.mu_h_I+tmp_input.gamma_h)*(tmp_input.mu_h+tmp_input.nu)))];

Resp =[R_rr,R_rh,R_hh];
end
