%% NST Proximal Operator (For S updates)
function S = nst_prox(R, tau ,w)
S = sign(R) .* max(abs(R) - tau.*w, 0);
end

% function out=softthre_s(a,tau,w)
% out = sign(a).* max( abs(a) - tau*w, 0);
% end
