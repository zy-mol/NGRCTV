function [DenoisedHSI,U,V] = RCTV_IRNN(Nhsi, beta,lambda, tau, r, rho)
%function[output1,output2,…]=functionname(intput1,intput2,…)%左边输出右边输入
%% Input Varaible
%   Nhsi   --- Polluted HSI data
%   beta   --- trade-off parameter of Gaussian noise
%   lambda --- trade-off parameter of sparse noise
%   rk     --- rank of needed-to-be recovered Tensor along mode-3
%% Output Variable
%   DenoisedHSI --- Denoised HSI
%   U           --- Representative Coefficient Image(RCI)
%   V           --- Orthogonal Basis

%% Initializing admm variables  初始化ADMM参数
tol     = 1e-5;    %收敛容忍度
maxIter = 300;    %最大迭代次数
max_mu  = 1e6;    %最大惩罚参数
[M,N,p] = size(Nhsi);    %输入数据维数

%% 数据预处理--将三维高光谱数据(Nhsi)重塑为一个二维矩阵
D       = zeros(M*N,p) ;    %预分配一个M*N行P列的全零矩阵
for i=1:p
    bandp = Nhsi(:,:,i);
    D(:,i)= bandp(:);    %将二维图像展开成一列(:)，存入D的第i列
end

%% 参数默认值处理--检查输入参数数量，并为缺失参数赋默认值
if nargin < 2    %表示实际输入参数的数量，如果输入参数少于2个（即缺失beta）
    beta = 50;
end
if nargin < 3    % 如果输入参数少于3个（即缺失lambda）
    lambda = 0.3;
end
if nargin < 4
    tau = [0.8,0.8];
end
if nargin < 5
    r = 6;
end

%% % choose penalty in IRNN
%  fun1 = 'lp' ;        gamma1 = 0.7;
% fun = 'scad' ;      gamma =0.1 ;
% fun2 = 'logarithm' ; gamma2 = 0.5; % or 1
 fun1 = 'mcp' ;       gamma1 =0.1;    %0.9\1.5
%  fun1 = 'etp' ;       gamma1 = 0.01;
 fun2 = 'logarithm'  ;      gamma2 = 1;

%% initialize--通过归一化将数据调整到标准尺度
Y = D;
norm_two = lansvd(Y, 1, 'L');     % 计算Y的谱范数（最大奇异值）
%	lansvd语法：s = lansvd(Y, k, 'L')	对待分解矩阵Y计算最大(L表示最大/S表示最小)的第k个奇异值
norm_inf = norm( Y(:), inf);	 % 计算Y的无穷范数（最大绝对值）
%	norm语法：norm_p= norm(Y(:), p)	将Y展开为列向量，计算其p范数（inf指定计算无穷范数）
dual_norm = max(norm_two, norm_inf);	%取最大值可确保归一化后，矩阵的整体尺度和极端值均被合理控制??
Y = Y / dual_norm;
mu = 1/norm_two;		% 初始化ADMM的惩罚参数
normD   = norm(D,'fro');	% 计算D的Frobenius范数



%% FFT setting 
h               = M;
w               = N;
d               = r;	% r控制分解后的成分数量
sizeU           = [h,w,d];

%%  
%	函数otf = psf2otf(psf, output_size);	将 ??点扩散函数（PSF）?? 转换为 ??光学传递函数（OTF）
%		输入：psf--点扩散函数；output_size--输出otf的尺寸；输出：otf--频域表示		
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
determ  =  Eny_x + Eny_y;



%% Initializing optimization variables			**
[u,s,v]= svd(D,'econ');
U              = u(:,1:r)*s(1:r,1:r);
V              = v(:,1:r);
X              = U*V';
S              = zeros(M*N,p);%Sparse
E              = zeros(M*N,p);%Gaussian

M1 = zeros(M*N*r,1);  % multiplier for Dx_U=G1
W1 = ones(M*N*r,1);
M2 = zeros(M*N*r,1);  % multiplier for Dy_U=G2
W2 = ones(M*N*r,1);
M3 = zeros([M*N,p]);  % multiplier for D-UV^T-E
W3 = ones(M*N,p);
% main loop
iter = 0;
eps =1e-4;
tic

%% 依次更新G1/G2-U-V-E-S-乘子-检查停止条件
while iter<maxIter
    iter          = iter + 1;  
    
    %% -Update G1,G2		**
        % Compute gradients
        % IRNN proximal operator for G1
        R1= diff_x(U,sizeU)+M1/mu;
        H1 = reshape(R1, sizeU);
        G1 = irnn_row(H1, tau(1)/mu, gamma1, fun1, tau(1));    
        G1 = reshape(G1, size(R1));        

        % IRNN proximal operator for G2 
        R2 = diff_y(U,sizeU)+M2/mu;
        H2 = reshape(R2, sizeU);
        G2 = irnn_col(H2, tau(2)/mu, gamma1 , fun1, tau(2));
        G2 = reshape(G2, size(R2));
 
 

    %% -Update U
    diffT_p  = diff_xT(G1-M1/mu,sizeU)+diff_yT(G2-M2/mu,sizeU);
    temp     = (D-E-S+M3/mu)*V;
    numer1   = reshape( diffT_p +temp(:), sizeU);
    x        = real( ifftn( fftn(numer1) ./ (determ + 1+eps) ) );
    U        = reshape(x,[M*N,r]);

    %% -Update V
    [u,~,v]     = svd((D-E-S+M3/mu)'*U,'econ');
    V           = u*v';

    %% -Update E
    E = mu*(D-U*V'-S+M3/mu)/(2*beta+mu);

    %% -Update S			**
    if lambda >100
        % 忽略掉高斯噪音影响
        S  = 0;
    else
        T = D - U*V' - E + M3/mu;
        S = nst_prox(T, lambda/mu, W3);
    end

    %% -Update Multiplier
    leq1 = diff_x(U,sizeU)-G1;
    leq2 = diff_y(U,sizeU)-G2;
    leq3 = D - U*V' - E - S;
    stopC1 = norm(leq1,'fro')/normD;
    stopC2 = norm(leq2,'fro')/normD;
    stopC3 = norm(leq3,'fro')/normD;

    if mod(iter,10)==0
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
                ',U_x rele = ' num2str(stopC1,'%2.3e') ',U_y rele = ' num2str(stopC2,'%2.3e')...
                ',X-UV = ' num2str(stopC3,'%2.3e')]);
    end

    if stopC1<tol && stopC2<tol && stopC3 <tol
        break;
    else
        M1 = M1 + mu*leq1;
        M2 = M2 + mu*leq2;
        M3 = M3 + mu*leq3;
        mu = min(max_mu,mu*rho);
%         epsilon = 1e-6 ;
%         W3(S~=0) = 1 ./(abs(S(S~=0))+ epsilon) ;  
        hfun_sg = str2func([fun2 '_sg']);
        W3 = hfun_sg(S, gamma2 ,lambda);
     end 
end
DenoisedHSI = reshape(U*V',[M,N,p]);
end

