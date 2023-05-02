function [Ahat,Bhat,C,D] = implicit_euler(A,B,C,D,Ts)
% Converts system from continuous to discrete with sampling frequency Ts

n = length(A);
I = eye(n);

Ahat = (I-Ts*A)\I;
Bhat = Ts*((I-Ts*A)\B);