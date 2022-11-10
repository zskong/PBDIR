function [Z,S,P] = DDL(X,F_ini,Z_ini,c,lambda1,lambda2,lambda3,pp,max_iter)
% This code is written by Zisen Kong,2021-12-14.,
% Revised version of the second draft,2022-1-8.
% Revised version of the third draft,2022-3-8.
[m,n] = size(X);
% ---------- Initilization -------- %
miu = 1e-1;%μ
rho = 1.1;%ρ
max_miu = 1e8;%max
tol  = 1e-5;
tol2 = 1e-2;
C1 = zeros(n,n);
%P = ones(m,m);
P=ones(pp,m);
%YY= zeros(n);
Z = Z_ini;
S = Z_ini;  
W = F_ini*F_ini';
clear Z_ini，F_ini;
for iter = 1:max_iter
    % 初始化
%     if iter == 1
%         Z = Z_ini;
%         S = Z_ini;
%         W = F_ini*F_ini';
%         %YY = YY_ini;
%         clear Z_ini，F_ini;
   % end
    S_old = S;
    Z_old = Z;
    P_old = P;
    %% -------- Update Z --------- 
    M_1=X'*P'*P*X;
    e=eye(n);
    %zhong=(2*e+miu*e+2*lambda2*M_1);
    zhong=(2*lambda1*e+miu*e+2*lambda2*M_1);
    vvv=pinv(zhong);
    Z=vvv*(2*lambda2*M_1+miu*S-C1);
    Z = Z- diag(diag(Z));%z_ii=0 
%     for ic = 1:n
%         idx    = 1:n;
%         idx(ic) = [];
%         Z(ic,idx) = EProjSimplex_new(Z(ic,idx));          % 
%     end
    %% -------- Update S --------- 
    D_W= diag(W)*ones(1,n)-W;
    distX = L2_distance_1(P*X,P*X);
    S= -lambda3/(2*miu)*(D_W+D_W')-distX/miu+Z+C1/miu; 
    S= S - diag(diag(S));
    for ic = 1:n
        idx    = 1:n;
        idx(ic) = [];
        S(ic,idx) = EProjSimplex_new(S(ic,idx));          % 
    end
    %% ---------- Update B ----------- %
    B = (S+S')/2;
    L = diag(sum(B)) - B;
    [F, ~, ev] = eig1(L, c, 0);
    W = F*F';
    %% -------- Update P --------- 
    Wz = (S+S')/2;
    Dz = diag(sum(Wz));
    Lz = Dz-Wz;
    Lz1=2*X*Lz*X'+lambda2*(X-X*Z)*(X-X*Z)';
    %[u,~,v]=svd(Lz1');
   % P=v*u';
    [P, ~, ~] = eig1(Lz1, pp, 0);
    %Pkm1=select_vector(Lz1);
    %P=Pkm1(:,1:pp);
    P=P';
    
    %% -------- Update C1 miu -------- %
    A1 = Z-S;
    C1 = C1+miu*A1;
    LL1 = norm(Z-Z_old,'fro');
    LL2 = norm(S-S_old,'fro');
    LL3 = norm(P-P_old,'fro');
    %SLSL(iter) = max(max(LL1,LL2))/norm(X,'fro');
    SLSL = max(max(LL1,LL2),LL3)/norm(X,'fro');
    if miu*SLSL < tol2
        miu = min(rho*miu,max_miu);
    end
   
    
   %% --------- obj ---------- %
    leq1 = max(max(abs(A1(:))));
    stopC = max(leq1);
       
    %obj(iter) = (2*trace(P*X*Lz*X'*P')+lambda1*(norm(Z,'fro')^2)+lambda2*(norm(P*X-P*X*Z,'fro')^2)+lambda3*trace(W'*L))/norm(X,'fro')^2;
    %if iter > 2
        %if abs(obj(iter)-obj(iter-1)) < 10^-7
    if stopC < tol 
            iter
            break;
        %end
    end
    %
    %obj(iter) = (sum(sum(abs(distX.*Z))) + lambda1*Z_obj+ lambda2*sum(sum(abs(E))) + lambda3*(D_Y+D_Y'))/norm(X,'fro');
end
end