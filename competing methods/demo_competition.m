%% 
clc;
clear;close all;
rng('default');rng(1997);   %随机数生成器,便于结果复现
addpath(genpath('lib'));
addpath(genpath('data'));
dataName = 'pure_DCmall';   % pure_DCmall
%  Please make sure the HSI is a cubic of size [height, width, band] and in range [0, 1].
%  You can use other tensor data such as RGB Image, Video, CT/MRI for test. 
%  Note some parameter might need reset for other methods. 
%  But the proposed TCTV is parameter free, where the trade-off paramter 'lambda' is determined theoretically 
dataRoad = ['data/', dataName];
resultsRoad = ['results/TRPCA/results_for_', dataName];
if ~exist(resultsRoad); mkdir(resultsRoad); end
%% Set enable bits
% 
Run_RPCA      = 1;% 1 or 0 <=> run or not
Run_SNN       = 1;  
Run_KBR       = 1;
Run_TCTV      = 1; 
Run_E3DTV    = 1;
Run_LRTV      = 1;
Run_RCTV       =1;
Run_WRCTV    =1;
% 
% Run_RPCA      = 0;% 1 or 0 <=> run or not
% Run_SNN       = 0;  
% Run_KBR       = 0;
% Run_TCTV      = 0; 
% Run_E3DTV    = 0;
% Run_LRTV      = 0;
% Run_RCTV       =0;
% Run_WRCTV    =1;

getGIF         = 1;  % if save gif result or not
getBand        = 1;  % if save one slected band of the denoised HSI as a gray image result or not
selected_band  = 25;
getPseudoImage = 1;  % if save three slected bands of the denoised HSI as a pseudo color image result or not
selected_bands = [49, 27, 7];
%% Load Data    
methodName = {'Noisy', 'RPCA',  'SNN', 'KBR', 'TCTV', 'E3DTV', 'LRTV', 'RCTV', 'WRCTV' };
Mnum = length(methodName);
load(dataRoad);  % load data
% Ori_H = data;
Ohsi = Ori_H;
[height, width, band] = size(Ohsi);
dim = [height, width, band];

%% Add impluse noise
i = 1;
sparselevel = 0.1;
sparsesigma  = sparselevel* ones(band,1);
%保存路径
disp(['=== the sp noise level is ', num2str(sparselevel), ' ===']);
saveRoad = ['results/TRPCA/results_for_', dataName, '/', 's', erase(num2str(mean(sparselevel)),'.')];
% % % % ↑打印当前噪声级别并保存

% % 添加高斯噪声
% noiselevel = 0.1;
% gausssigma = noiselevel *ones(band, 1) ;  
% disp(['=== the gauss noise level is ', num2str(noiselevel), ' ===']);
% saveRoad = ['results/TRPCA/results_for_', dataName, '/', 'g', erase(num2str(noiselevel),'.')];

% % 添加稀疏+高斯噪声
% sparselevel = 0.1;
% sparsesigma  = sparselevel* ones(band,1); 
% noiselevel = 0.1;
% gausssigma = noiselevel *ones(band, 1) ;  
% disp(['=== the sp noise level is ', num2str(sparselevel), ' ===']);
% disp(['=== the gauss noise level is ', num2str(noiselevel), ' ===']);
% saveRoad = ['results/TRPCA/results_for_', dataName, '/', 'g', erase(num2str(noiselevel),'.'), '_s', erase(num2str(mean(sparselevel)),'.')];

if ~exist(saveRoad); mkdir(saveRoad); end
if exist([saveRoad, '/QA_Results.mat']); load([saveRoad, '/QA_Results.mat']); end
if exist([saveRoad, '/Results.mat']); load([saveRoad, '/Results.mat']); end
if getGIF
    if ~exist([saveRoad, '/GIF']); mkdir([saveRoad, '/GIF']); end
    togetGif(Ohsi, [saveRoad, '/GIF/Ohsi']); 
end
if getBand
    if ~exist([saveRoad,'/Band']); mkdir([saveRoad,'/Band']); end
    imwrite(Ohsi(:,:,selected_band), [saveRoad, '/Band/Ohsi.jpg']); 
end
if getPseudoImage
    if ~exist([saveRoad,'/PseudoImage']); mkdir([saveRoad,'/PseudoImage']); end
    imwrite(PseudoImage(Ohsi, selected_bands), [saveRoad, '/PseudoImage/Ohsi.jpg']);
end

rng(42);  % 可以任意指定一个数字
% 添加椒盐噪声
for ii = 1:band
    Nhsi(:,:,ii) = imnoise(Ohsi(:,:,ii), 'salt & pepper', sparsesigma(ii));
end

% %添加高斯噪声
% for ii = 1:band
%     Nhsi(:,:,ii) = Ohsi(:,:,ii) +gausssigma(ii) * randn(height, width);
% end

% % 添加高斯+稀疏噪声
% for ii = 1:band
%     Nhsi(:,:,ii) = Ohsi(:,:,ii) + gausssigma(ii) * randn(height, width);
% end
% for ii = 1:band
%     Nhsi(:,:,ii) = imnoise(Nhsi(:,:,ii), 'salt & pepper', sparsesigma(ii));
% end

% 评估噪声图像质量
Results{i} = Nhsi;
[MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
enList = 1;

%% Run RPCA
i = i+1;
if Run_RPCA
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the rpca_m.m in RPCA file
    D = zeros(height* width, band);
     for j=1:band
        bandp = Nhsi(:,:,j);
        D(:,j)= bandp(:);
    end
    tic;
    A_hat = rpca_m(D);
    Results{i} = reshape(A_hat, [height, width,band]);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run SNN
i = i+1;
if Run_SNN
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']); 
    % see the trpca_snn.m in SNN file
    % the used code is from Lu Canyi, A unified alternating direction method of multipliers by majorization minimization, TPAMI, 2018 
    opts = [];
    opts.DEBUG = 1;
    tic
    alpha=[1, 1, 200]; 
    Results{i} = trpca_snn(Nhsi, alpha, opts);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run KBR
i = i+1;
if Run_KBR
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the KBR_RPCA.m or Demo_TC_MSI.m in KBR file
    
    beta           = 2.5*sqrt(max(dim));
    gamma          = beta*100;
    Par.maxIter    = 1000;
    Par.lambda     = 0.1;
    Par.mu         = 10;
    Par.tol        = 1e-5;
    Par.rho        = 1.1;
    
    tic
    Results{i} =   KBR_RPCA(Nhsi,beta,gamma,Par);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run TCTV
i = i+1;
if Run_TCTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the TCTV_TRPCA.m in TCTV file
    opts = [];
    opts.rho = 1.25; % larger rho makes the algorithm faster while lose the accuracy
    opts.directions = [1,2,3]; % consider the lowrankness and smoothness both along the spatial and spectral directions
    tic
    Results{i} = TCTV_TRPCA(Nhsi, opts);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath('TCTV'));
    enList = [enList, i];
end

%% Run E3DTV
i = i+1;
if Run_E3DTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the EnhancedTV.m in Enhanced3DTV file
%     r = [6, 6, 6];
    r = [13, 13, 13];
    tau = 0.004 *sqrt(height * width);
    tic;
    Results{i} = EnhancedTV(Nhsi, tau, r);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run LRTV
i = i+1;
if Run_LRTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the LRTV_accelerate.m in LRTV file
    tau = 0.005;
    lambda = 5/sqrt(height*width);
%     rank = 6;
    rank = 13;
    
    tic
    Results{i} = LRTV_accelerate(Nhsi, tau, lambda, rank);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run RCTV
i = i+1;
if Run_RCTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the RCTV_TRPCA.m in RCTV file
    r = 13;%HSI_13
    beta =50;
    lambda = 1;
    tau = [0.8, 0.8];
    tic;
    Results{i} = RCTV(Nhsi, beta,lambda, tau, r);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run WRCTV
i = i+1;
if Run_WRCTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the WRCTV_TRPCA.m in WRCTV file
%     r = 6;
    r = 13;
    beta = 50;%50，30
    lambda = 0.30;%1, 5,0.5,0.6
    tau = [0.8, 0.8];% 0.8,need to fine tune,0.7-48.62
    rho = 1.05;%1.05
    Results{i} = RCTV_IRNN(Nhsi, beta,lambda, tau, r, rho);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Show result
fprintf('\n');

fprintf('================== QA Results =====================\n');
fprintf(' %8.8s    %5.5s    %5.5s    %5.5s    %5.5s    %5.5s    %5.5s  \n',...
    'Method', 'MPSNR', 'MSSIM', 'MFSIM', 'ERGAS', 'MSAM', 'Time');
for i = 1:length(enList)
    fprintf(' %8.8s   %5.3f    %5.3f    %5.3f    %5.3f    %5.3f    %5.3f   \n',...
        methodName{enList(i)}, MPSNR(enList(i)), MSSIM(enList(i)), MFSIM(enList(i)), ERGAS(enList(i)), MSAM(enList(i)), Time(enList(i)));
end

fprintf('================== Show Results =====================\n');
close all;
showHSIResult(Results,Ohsi,methodName,enList,selected_band,band);


%% save results
All = [MPSNR; MSSIM; MFSIM; ERGAS; MSAM; Time];
save([saveRoad,'/QAResults'], 'All', 'MPSNR', 'MSSIM', 'MFSIM', 'ERGAS', 'MSAM', 'Time');
save([saveRoad,'/Results'], 'Results');