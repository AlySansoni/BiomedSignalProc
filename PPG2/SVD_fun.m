%% Singular Decomposition

function[U5,S5,V5] = SVD_fun(sig)

%% Embedding
N = length(sig);
L = 1000; % fs/fl < L < N/2 --- 65 < L < 3500
K = N-L+1;
Tx = zeros(K,L);
for i=1:K
    Tx(i,:) = sig(1,i:(i+L-1));
end

%% SVD
[U,S,V] = svd(Tx, 'econ'); % U: KxL - S: diagonal matrix - V: LxL
num_top5 = length(S)*0.05;
top_S = max(S); 
% U, S, V of the top 5 values
U5 = U(:,1:num_top5);
S5 = top_S(1:num_top5);
V5 = V(:,1:num_top5);
end