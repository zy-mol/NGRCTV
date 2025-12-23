%% IRNN Proximal Operator (For Gi updates)
function X = irnn_row(R, tau, gamma, fun, lambda)
% Input:
%   R     : Input matrix (vectorized)
%   tau   : Regularization parameter
%   gamma : Non-convexity parameter (0<gamma<1)
%   sizeX : Original size [M,N,r]

%%
 for k = 1:size(R,3)
     for i = 1:size(R,1)
         R_row = squeeze(R(i, :, k));
         R_norm = norm(R_row,2); 
         hfun_sg = str2func([fun '_sg']);
         w = hfun_sg(R_row,gamma,lambda);
         if R_norm >  tau*w
            % 情况1: 进行缩放
            X(i, :, k) = (R_norm - tau*w)/R_norm .* R_row;
         else
            % 情况2: 置零
            X(i, :, k) = zeros(size(R_row));
         end
     end
        
end



%%
% for i = 1:size(R,1)   % 遍历R的每一行
%     R_row = R(i,:);   % 提取第i行（对应[R]_{i,:}）
%     R_norm = norm(R_row,2);  % 计算行向量的l2范数（对应‖[R]_{i,:}‖2）
%     hfun_sg = str2func([fun '_sg']);
%     w = hfun_sg(R_row,gamma,lambda);
%     if R_norm >  tau*w
%         % 情况1: 进行缩放
%         X(i,:) = (R_norm - tau*w)/R_norm .* R_row;
%     else
%         % 情况2: 置零
%         X(i,:) = zeros(size(R_row));
%     end
% end

%%
%         hfun_sg = str2func([fun '_sg']);
%         w = hfun_sg(R(:),gamma,lambda);
%         X = sign(R).* max( abs(R) - tau*w, 0);
%         
        


%%  X = sign(R).* max( abs(R) - tau*w, 0);
end