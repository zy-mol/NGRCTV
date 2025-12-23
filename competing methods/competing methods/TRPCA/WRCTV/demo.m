clear all;clc;	

% addpath('D:\RCTV\RCTV-IRNN\PROPACK');
% addpath('D:\RCTV\RCTV-IRNN\nonconvex_funs');
% addpath('D:\RCTV\RCTV-IRNN\prox');
% addpath('D:\RCTV\RCTV-IRNN\quality_assess');
% addpath('D:\RCTV\RCTV-IRNN\prox_operators');
% addpath('D:\RCTV\RCTV-IRNN\TV_operator');

load('Simu_indian.mat')	%加载名为Simu_indian.mat的模拟高光谱数据集
%load('pure_DCmall.mat')
Ohsi = Normalize(Ori_H);	%数据归一化
Nhsi      = Ohsi;
[M,N,p]   = size(Ohsi);

%设置噪声参数
noiselevel = 0.075*rand(p,1);		%高斯噪声强度
ratio = 0.15*rand(p,1);		%椒盐噪声密度
%% Gaussian noise	添加高斯噪声
for i = 1:p
     Nhsi(:,:,i)=Ohsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
%添加椒盐噪声
for i = 1:p
     Nhsi(:,:,i)=imnoise(Nhsi(:,:,i),'salt & pepper',ratio(i));
end

%%调用去噪算法RCTV
r = 13;
beta = 50;%50
lambda = 0.70;%1, 5,0.5,0.6
tau = [0.7, 0.7];% 0.8,need to fine tune,
rho     = 1.4;    %惩罚参数增长系数1.075;1.11;1.10;1.584;1.53,1.528
tic	;%计算算法运行时间
output_image = RCTV_IRNN(Nhsi, beta,lambda, tau, r, rho);
time = toc;
[mpsnr,mssim,ergas] = msqia(Ohsi, output_image)	%计算去噪质量指标
fprintf('运行时间: %.2f 秒\n', time);
