clear all;
%Bisect_func( @(x)exp(x)-2-cos(exp(x)-2), 0.5, 1.5, 10^-5, 25);
Bisect_func(@(x)x^4-2*x^3-4*x^2+4*x+4, 0, 2, 10^-2, 100);
