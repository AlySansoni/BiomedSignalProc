%% Artifact removal

close all;
clear;
clc;

% loading a single signal
% First input SNR = -10 dB
signal = load('Input/Signal_10.mat').zeroPhase;
labels = load('Input/Labels_10.mat').labels;

%Second input SNR = -20 dB
% signal = load('Input/Signal_20.mat').zeroPhase;
% labels = load('Input/Labels_20.mat').Predict_label;

figure();
plot(signal);
title('Corrupted Signal');
axis([20000 70000 -2 2]);
set(gca,'XTick',20000:5000:70000);
set(gca,'XTickLabel',40:10:140);
xlabel('Seconds');
ylabel('Signal amplitude [au]');
% % segments
% xline(21000, '--g', '0');
% xline(24500, '--g', '0');
% xline(28000, '--r', '1');
% xline(31500, '--r', '1');
% xline(35000, '--r', '1');
% xline(38500, '--r', '1');
% xline(42000, '--r', '1');
% xline(45500, '--r', '1');
% xline(49000, '--r', '1');
% xline(52500, '--r', '1');
% xline(56000, '--r', '1');
% xline(59500, '--r', '1');
% xline(63000, '--g', '0');
% xline(66500, '--g', '0');

% finding corrupted segments
flag = find(labels == 1);

% taking a clean segment as reference
if flag > 2
    sigClean = signal((3500*(flag-2))+1:3500*(flag-1));
end

% considering the first corrupted segment
a = 1;
% sigCorr = signal((3500*(flag(a)-1))+1:3500*(flag(a))); % a-th segment
sigCorr = signal((3500*(flag(a)))+1:3500*(flag(a)+1)); % (a+1)-th segment
    
% figure();
% subplot(2,1,1);
% plot(sigCorr);
% title('sigCorr');
% subplot(2,1,2);
% plot(sigClean);
% title('sigClean');
    
%% Step 1-2: computing SVD on corrupted and clean segments and keeping the top 5% of their components
% We called sqrt(lambda) as LB for a better readability
[U_clean,LB_clean,V_clean] = SVD_fun(sigClean); % U: KxL - S: diagonal matrix - V: LxL
[U_corr,LB_corr,V_corr] = SVD_fun(sigCorr); % U: KxL - S: diagonal matrix - V: LxL
    
%% Step 3: replacing corrupted eigenvalues with corresponding clean eigenvalues
LB_corr = LB_clean;

%% Step 4-5: working on frequencies
Uf_clean = abs(fft(U_clean));
Uf_corr = abs(fft(U_corr));
U_match = zeros(length(Uf_clean(1,:)));
i = 1;
    
for j=1:length(Uf_clean(1,:))
    % Step 4: choosing components with HR frequency range 0.66 < Fs < 3 Hz
    Uf_clean_cut = Uf_clean(Uf_clean(:,j) > 0.66,j);
    Uf_clean_cut = Uf_clean_cut(Uf_clean_cut < 3);
        
    Uf_corr_cut = Uf_corr(Uf_corr(:,j) > 0.66,j);
    Uf_corr_cut = Uf_corr_cut(Uf_corr_cut < 3);
        
    % Step 5: frequency matching
    Uf_clean_max = max(Uf_clean_cut);
    Uf_corr_max = max(Uf_corr_cut);

    if abs(Uf_clean_max - Uf_corr_max) < 0.5 % frequency matching error = 0.5
        U_match(i) = j;
        i = i+1;
    end
end
    
U_match = nonzeros(U_match);
U_match = U_match';
        
%% Step 6: applying the basic SSA algorithm iteratively
U_sig_corr = U_corr(:,U_match(:)); % corrupted components obtained from frequency matching
U_sig_clean = U_clean(:,U_match(:)); % clean components obtained from frequency matching
len = length(U_sig_clean(:,1));
U_final = zeros(length(U_sig_corr),length(U_match));

% loop on component obtained from step 5
for z =1:length(U_match)
    iter = 4; % after 4 iterations the result converges
    Ur = zeros(len,iter);
    DM_corr = zeros(iter,1); % discarding metric
    Ur(:,1) = SSA_fun(U_sig_corr(:,z)); % U reconstructed
    DM_corr(1,1) = sum(Ur(:,1)) / len;  % discarding metric for components achieved from SSA iterations
    DM_clean = sum(U_sig_clean(:,z)) / len; % discarding metric for counterpart clean components
    DM_best = DM_corr(1,1); 
    U_final(:,z) = Ur(:,1);
    % selecting components with the closest DM ???
    for i=2:iter
        Ur(:,i) = SSA_fun(Ur(:,i-1));
        DM_corr(i,1) = sum(Ur(:,i)) / len;
        % verifying if the last DM is the best (yes, it is)
        if abs(DM_clean - DM_corr(i,1)) < abs(DM_clean-DM_best)
            DM_best = DM_corr(i,1);
            U_final(:,z) = Ur(:,i);
        end
    end
    
    % plotting an example of clean, corrupted and reconstructed eigenvectors
    if z == 10
        figure();
        subplot(3,1,1);
        plot(U_sig_clean(:,z));
        axis([0 100 -0.05 0.05]);
        title('Clean EigenVector');
        subplot(3,1,2);
        plot(U_sig_corr(:,z));
        axis([0 100 -0.05 0.05]);
        title('Corrupted EigenVector');
        subplot(3,1,3);
        plot(U_final(:,z));
        axis([0 100 -0.05 0.05]);
        title('Reconstructed EigenVector');
    end
end
    
%% Step 7: reconstruct the corrupted PPG segment 
N = 3500;
L = 1000; % fs/fl < L < N/2 --- 65 < L < 3500
K = N-L+1;
Tx = zeros(K,L);
    
for k = 1:length(U_match)
    X = zeros(3,K,L);
    
    % replacing the corrupted eigenvector with the reconstructed one
    T = U_final(:,k)*LB_corr(U_match(k))*V_corr(:,U_match(k))';
    Tr = rot90(T);
    h = zeros(1,K);
    count = 1;
    % computing diagonal averaging
    for j=-(L-1):K-1
        d = diag(Tr,j);
        h(1,count) = sum(d) / length(d);
        count = count+1;
    end
    Tx = Tx + hankel(h(1:K),h(K:end));
  
    %Tx = Tx + (U_final(:,k)*LB_corr(U_match(2,k))*V_corr(:,U_match(2,k))');
end
    
for i=1:K
    sigR(1,i:(i+L-1)) = Tx(i,:);
    %Tx(i,:) = sig(1,i:(i+L-1));
end
    
figure();
subplot(3,1,1);
plot(sigClean);
title('Clean Segment');
subplot(3,1,2);
plot(sigCorr);
title('Corrupted Segment');
subplot(3,1,3);
plot(sigR);
title('Reconstructed Segment');
