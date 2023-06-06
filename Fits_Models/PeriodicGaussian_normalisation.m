function [K_s] = PeriodicGaussian_normalisation(s)

fun = @(X) exp(s*cos(X));
I0_S = exp(-s)*integral(fun,0,pi)/pi;

K_s = I0_S;
end