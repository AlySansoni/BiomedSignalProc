%% Singular Spectrum Analysis

function[sigR] = SSA_fun(sig)

% sig = load('tmp').U;
sig = sig';

%% Singular Decomposition
[U,S,V] = SVD_fun(sig);

% Keeping the first 3 largest eigenvalues obtained from SVD
U = U(:,1:3);
S = S(1,1:3);
V = V(:,1:3);

%% Spectral Reconstruction
N = length(sig);
L = 1000;
K = N-L+1;

% Grouping
Tj = zeros(3,K,L);
for i=1:3
    Tj(i,:,:) = (U(:,i) * S(1,i) * V(:,i)');
end

% Diagonal Averaging
X = zeros(3,K,L);
for i=1:3
    T = squeeze(Tj(i,:,:));
    Tr = rot90(T);
    h = zeros(1,K);
    count = 1;
    for j=-(L-1):K-1
        d = diag(Tr,j);
        h(1,count) = sum(d) / length(d);
        count = count+1;
    end
    X(i,:,:) = hankel(h(1:K),h(K:end));
end

X = X(1,:,:) + X(2,:,:) + X(3,:,:);
X = squeeze(X(1,:,:));

for i=1:K
    sigR(1,i:(i+L-1)) = X(i,:);
end

