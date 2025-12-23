function X = irnn_col(R, tau, gamma, fun, lambda)
% Input:
%   R     : Input matrix (vectorized)
%   tau   : Regularization parameter
%   gamma : Non-convexity parameter (0<gamma<1)
%   sizeX : Original size [M,N,r]

%%
 for k = 1:size(R,3)
     for j = 1:size(R,2)
         R_col = squeeze(R(:, j, k));
         R_norm = norm(R_col,2); 
         hfun_sg = str2func([fun '_sg']);
         w = hfun_sg(R_col,gamma,lambda);
         if R_norm >  tau*w
            % 情况1: 进行缩放
            X(:, j, k) = (R_norm - tau*w)/R_norm .* R_col;
         else
            % 情况2: 置零
            X(:, j, k) = zeros(size(R_col));
         end
     end
        
end
end