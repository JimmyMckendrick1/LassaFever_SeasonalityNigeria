function [K_s] = PeriodicGaussian_normalisation(s)
% Function to produce normalising constant for recruitment function for
% rats so that year-on-year population is stable.

fun = @(X) exp(s*cos(X));
I0_S = exp(-s)*integral(fun,0,pi)/pi;

K_s = I0_S;
end