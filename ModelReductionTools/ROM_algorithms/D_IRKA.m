function [Ar,Br,Cr,Er,Stilde_vec,count,res_vec]= D_IRKA(A,B,C,E,s,opts)
%A, B,C,E are system matrices
%s is initial guess of poles

%% Place default values in opts
def.converge_tol = 1e-2;
def.max_iter = 50;
if exist('opts','var') == 0
    opts=def;
else
    %Merge default values into opts
    A=fieldnames(def);
    for m=1:length(A)
        if ~isfield(opts,A(m))
            dum=char(A(m));
            dum2=getfield(def,dum);
            opts=setfield(opts,dum,dum2);
        end
    end
end
n = length(A);
r = length(s);
Vr = zeros(n,r);
Wr = zeros(n,r);
%get intial Wr,Vr
for i =1:r
    Vr(:,i) = (s(i)*E-A)\B;
    Wr(:,i) = (s(i)*E.'-A.')\C.';
end
[Vr,~] = qr(Vr,0);
[Wr,~] = qr(Wr,0);

count = 0;
nmax = opts.max_iter;
Stilde_vec = zeros(nmax,r);
res_vec = zeros(nmax,1);
tol = opts.converge_tol;
converged = 0;
diverged = 0;
while ~converged && ~diverged
    count = count + 1;
    Er = Wr'*E*Vr; Ar = Wr'*A*Vr;
    Sold = s;
    %reflect IRKA points outside unit circle
    s = 1./eig(Ar,Er);
    s = sort(s,'ComparisonMethod','real');
%     if count > 2
%         res = min(norm(Sold-s),norm(Stilde_vec(count-2,:).'-s));
%     else
    res = norm(Sold-s)/norm(Sold);
%     end
    res_vec(count) = res;
    Stilde_vec(count,:) = s.';
    for i =1:r
        Vr(:,i) = (s(i)*E-A)\B;
        Wr(:,i) = (s(i)*E'-A')\C';
    end
    [Vr,~] = qr(Vr,0);
    [Wr,~] = qr(Wr,0);
    converged = res < tol;
    diverged = count >= nmax;
end

Er = Wr'*E*Vr; Ar = Wr'*A*Vr;
Br = Wr'*B; Cr = C*Vr;

s = 1./eig(Ar,Er);
Stilde_vec(count+1,:) = s.';
%Stilde = s;

res_vec = res_vec(~(res_vec == 0));
Stilde_vec = Stilde_vec(~(Stilde_vec(:,1)==0),:);