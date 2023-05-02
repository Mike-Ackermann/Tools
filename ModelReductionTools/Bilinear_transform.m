function [Ad,Bd,Cd,Dd] = Bilinear_transform(A,B,C,D)

n = length(A);
I = eye(n);

inv_A = (I-A)\I;
Ad = (A+I)*inv_A;
Bd = sqrt(2)*(inv_A*B);
Cd = sqrt(2)*C*inv_A;
Dd = D + C*(inv_A*B);