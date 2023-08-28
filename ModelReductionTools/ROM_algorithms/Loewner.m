function [Ap,Bp,Cp,Ep] = Loewner(s1,s2,Hs1,Hs2,nmax,epsilon)
% s1 are left sample points
% s2 are right sample points
% Hs1 are transfer function values at s1
% Hs2 are transfer function values at s2

% Assumes length(s1) = length(s2)
[~,idx1] = sort(imag(s1),'ComparisonMethod','real');
s1 = s1(idx1);
Hs1 = Hs1(idx1);
[~,idx2] = sort(imag(s2),'ComparisonMethod','real');
s2 = s2(idx2);
Hs2 = Hs2(idx2);
r = length(s1);

% Solve via Sylvester equations
Z_R = [];
Y_R = [];
Sigma_R = [];
Theta_R = [];
B_R = [];
C_R = [];
index = 1;
while length(Sigma_R) < r
    if isreal(s1(index))
        Sigma_R = blkdiag(Sigma_R,s1(index));
        Y_R = [Y_R,Hs1(index)];
        B_R = [B_R,1];
    else
        T = [real(s1(index)), imag(s1(index)); -imag(s1(index)), real(s1(index))];
        Sigma_R = blkdiag(Sigma_R,T);
        Y_R = [Y_R,[real(Hs1(index)),imag(Hs1(index))]];
        B_R = [B_R,[1,0]];
    end
    index = index + 1;
end
index = 1;
while length(Theta_R) < r
    if isreal(s2(index))
        Theta_R = blkdiag(Theta_R,s2(index));
        Z_R = [Z_R,Hs2(index)];
        C_R = [C_R,1];
    else
        T = [real(s2(index)), imag(s2(index)); -imag(s2(index)), real(s2(index))];
        Theta_R = blkdiag(Theta_R,T);
        Z_R = [Z_R,[real(Hs2(index)),-imag(Hs2(index))]];
        C_R = [C_R,[1,0]];
    end
    index = index + 1;
end

X_L = (C_R.'*Y_R-Z_R.'*B_R).';
X_M = (C_R.'*Y_R*Sigma_R-Theta_R*Z_R.'*B_R).';

L_R = sylvester(Sigma_R.',-Theta_R.',X_L);
L_R = L_R.';
M_R = sylvester(Sigma_R.',-Theta_R.',X_M);
M_R = M_R.';



%%%%%%%%%%%%%%%%%%%%%
B = reshape(Z_R,r,1);
C = reshape(Y_R,1,r);
E = -L_R;
A = -M_R;



%%%%%%%%%%%%%%%%%%%%%
%% remove redundencies
[Y,theta1,X1] = svd([E A]);
[Y2,theta2,X] = svd([E;A]);

sVal = diag(theta2);
sValScaled = sVal/sVal(1);
%semilogy(1:r,sValScaled);
%hold on
%semilogy(1:r,10^(-10)*ones(r,1));

%epsilon was 10^(-13)
pHat = find(sValScaled < epsilon,1);
p = min([pHat,nmax]);%truncate at max allowed degree or at tolerence level

%semilogy(p,sValScaled(p));

%only remove reduncencies if they are present
if all(size(p)>0)
    Ep = Y(:,1:p)'*E*X(:,1:p);
    Ap = Y(:,1:p)'*A*X(:,1:p);
    Bp = Y(:,1:p)'*B;
    Cp = C*X(:,1:p);
else
    Ep = E;
    Ap = A;
    Bp = B;
    Cp = C;
end

%Hr = @(z) Cp*((z*Ep-Ap)\Bp);