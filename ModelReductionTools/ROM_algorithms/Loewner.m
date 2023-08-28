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

Mi = zeros(r); Li = zeros(r);
B = reshape(Hs1,r,1);
C = reshape(Hs2,1,r);

for i = 1:r
    for j = 1:r
        Li(i,j) = (Hs1(i)-Hs2(j))/(s1(i)-s2(j));
        Mi(i,j) = (s1(i)*Hs1(i)-s2(j)*Hs2(j))/(s1(i)-s2(j));
    end
end

E = -Li;
A = -Mi;

%% Keep real
%%%%%%%%%%%%%%%%%%%%%%%
%check for repeated shifts to not multiply by transformation matrix
% repeated_shift = NaN(r,1);
% count = 1;
% for i = 1:2:r
%    if s(i) == s(i+1)
%        repeated_shift(count) = 1;
%    else
%        repeated_shift(count) = 0;
%    end
%    count = count + 1;
% end
% idx = 1:r/2;
% idx = idx(logical(repeated_shift));

% % repeated_shift = imag(s) == 0;
% % repeated_idx = zeros(ceil(r/2));
% % isolated_idx = zeros(ceil(r/2));
% % for i = 1:2:r
% %     if repeated_shift(i) && repeated_shift(i+1)
% %         repeated_idx(floor(i/2)) = 1;
% %     elseif repeated_shift(i)
% %         isolated_idx(floor(i/2)) = 1;
% %     end
% % end

%put the real values in s first
% real_idx = imag(s) == 0;
% indicies = 1:r;
% real_idx = indicies(real_idx);
% s_imag = s(~real_idx);
% s_real = s(real_idx);
% I = eye(2);
% T1 = [1 1; 1i -1i];
% T1c = repmat({T1},1,r/2);
% 
% T1inv = (1/2)*[1 -1i; 1 1i];
% count = 1;
% T = [];
% Tinv = [];
% while count <= r
%     if ismember(count,real_idx)
%         T =blkdiag(T,1);
%         Tinv = blkdiag(Tinv,1);
%         count = count + 1;
%     else
%         T = blkdiag(T,T1);
%         Tinv = blkdiag(Tinv,T1inv);
%         count = count + 2;
%     end
% end
% T1invc = repmat({T1inv},1,r/2);
% for i = 1:length(idx)
%    T1c{idx(i)} = I;
%    T1invc{idx(i)} = I;
% end
% 
% Tinv = blkdiag(T1c{:});
% T = blkdiag(T1invc{:});
% 
% B = Tinv*Bi;
% C = Ci*T;
% A = Tinv*Ai*T;
% E = Tinv*Ei*T;
% L = Tinv*Li*T;
% M = Tinv*Mi*T;
% 
% E = real(E);
% A = real(A);
% B = real(B);
% C = real(C);
% L = real(L);
% M = real(M);

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
Bp = reshape(Z_R,r,1);
Cp = reshape(Y_R,1,r);
Ep = -L_R;
Ap = -M_R;
return


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