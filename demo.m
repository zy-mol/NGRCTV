clear all;clc;	

% addpath('D:\RCTV\RCTV-IRNN\PROPACK');
% addpath('D:\RCTV\RCTV-IRNN\nonconvex_funs');
% addpath('D:\RCTV\RCTV-IRNN\prox');
% addpath('D:\RCTV\RCTV-IRNN\quality_assess');
% addpath('D:\RCTV\RCTV-IRNN\prox_operators');
% addpath('D:\RCTV\RCTV-IRNN\TV_operator');

load('hsi_PaC.mat')
% Ori_H = Ori_H(1:2:200, 1:2:200, :);
Ori_H = data;
Ohsi = Normalize(Ori_H);	%数据归一化
Nhsi      = Ohsi;
[M,N,p]   = size(Ohsi);


rng(42);  % 可以任意指定一个数字
sparselevel = 0.2;
sparsesigma  = sparselevel* ones(p,1); 
% noiselevel = 0.1;
% gausssigma = noiselevel *ones(p, 1) ; 

%% Gaussian noise	添加高斯噪声
% for ii = 1:p
%     Nhsi(:,:,ii) = Ohsi(:,:,ii) + gausssigma(ii) * randn(M, N);
% end
for ii = 1:p
    Nhsi(:,:,ii) = imnoise(Nhsi(:,:,ii), 'salt & pepper', sparsesigma(ii));
end

%%调用去噪算法RCTV
r = 13;
beta = 50;%50
lambda = 0.3;%1, 5,0.5,0.6,0.7
tau = [0.8, 0.8];% 0.8,need to fine tune,
rho     = 1.3;    %惩罚参数增长系数1.075;1.1;1.5
tic	;%计算算法运行时间
output_image = RCTV_IRNN(Nhsi, beta,lambda, tau, r, rho);
time = toc;
[mpsnr,mssim,ergas] = msqia(Ohsi, output_image)	%计算去噪质量指标
fprintf('运行时间: %.2f 秒\n', time);
