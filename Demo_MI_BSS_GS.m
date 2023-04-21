%% This a demo for Robust blind separation of smooth graph signals using minimization of graph regularized mutual information
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Robust blind separation of smooth graph signals using minimization of graph regularized mutual information." Digital Signal Processing 132 (2022): 103792.
clc; close all; clear;
%% Generates the results of Figures 6 and 7 in our paper.
%% Add necessary functions:
addpath('./NeededFunctions')
%%
t0 = tic;
%% setting input parameters:
param.mu = 0.1; % the learning rate in the gradient descent step of the MI-BSS-GS method
param.lambda = 1e1; % balances the weight between minimization of mutual information and graph smoothness terms, e.g., 1e1
param.Tol = 1e-2; % tolerance for reaching convergence, e.g., 1e-2
param.MaxIter = 5000; % maximum iterations for convergence
%% Loading and normalization of two images:
s1 = imresize(imread('lena.png'), 0.25); % reading and resizing the first image
W = AssignGridWeights_8thNeighbour(double(s1)); % A function that assigns grid weights in 8-Neighbour strategy, as discussed in our paper
Adj_cell{1} = W; % A cell for storing adjacencies
L_cell{1} = (diag(W*ones(size(s1, 1)^2,1)) - W)/sum(W*ones(size(s1, 1)^2,1), []); % A cell for storing Laplacians
[N, N, N2] = size(s1); % initial size of the image
s11 = mat2gray(s1(1:N, 1:N)); % converting the image to grayscale
s1 = vec(s11)';
s1 = s1 / sqrt(mean(s1.^2)); % normalization
 
s2 = imresize(imread('man.png'), 0.5); % reading and resizing the second image
W = AssignGridWeights_8thNeighbour(double(s2));% A function that assigns grid weights in 8-Neighbour strategy, as discussed in our paper
Adj_cell{2} = W;% A cell for storing adjacencies
L_cell{2} = (diag(W*ones(size(s2, 1)^2,1)) - W)/sum(W*ones(size(s2, 1)^2,1), []);% A cell for storing Laplacians
s22 = mat2gray(s2(1:N, 1:N));% initial size of the image
s2 = vec(s22)'; %reshape(s22,1,N^2);
s2 = s2 / sqrt(mean(s2.^2));% normalization

%% Perform the related analysis explained in our manuscript:
a_vec = 0.05 : 0.05 : 0.45; % The span of mixing rate (alpha) in our paper
for a = a_vec 
disp(['>>> Mixing rate: ', num2str(a)])
% mixing image sources:
x1 = (1 - a) * s1 + a * s2;
x2 = a * s1 + (1 - a) * s2;
%% Apply MI-BSS method (MI-BSS-GL with lambda=0): 

[y, B1] = MI_BSS([x1;x2], param);
y1_p_eq_2 = y(1,:)/sqrt(mean(y(1,:).^2));
y2_p_eq_2 = y(2,:)/sqrt(mean(y(2,:).^2));


SNR11_p_eq_2(find(a_vec==a)) = 10*log10(mean((s1/norm(s1)).^2)/mean((s1/norm(s1)-y1_p_eq_2/norm(y1_p_eq_2)).^2));
SNR22_p_eq_2(find(a_vec==a)) = 10*log10(mean((s2/norm(s2)).^2)/mean((s2/norm(s2)-y2_p_eq_2/norm(y2_p_eq_2)).^2));


%% MI-BSS-GS:

[y_GS, B2] = MI_BSS_GS([x1;x2], L_cell, param);
y1_p_eq_GS = y_GS(1,:)/sqrt(mean(y_GS(1,:).^2));
y2_p_eq_GS = y_GS(2,:)/sqrt(mean(y_GS(2,:).^2));

SNR11_p_eq_GS(find(a_vec==a)) = 10*log10(mean((s1/norm(s1)).^2)/mean((s1/norm(s1)-y1_p_eq_GS/norm(y1_p_eq_GS)).^2));
SNR22_p_eq_GS(find(a_vec==a)) = 10*log10(mean((s2/norm(s2)).^2)/mean((s2/norm(s2)-y2_p_eq_GS/norm(y2_p_eq_GS)).^2));

%%
if a == a_vec(end-2) %0.35
    
    x11 = reshape(x1, N, N);

    x22 = reshape(x2, N, N);

    y1_MIBSS = reshape(y1_p_eq_2, N, N);
 
    y2_MIBSS = reshape(y2_p_eq_2, N, N);
    
    y1_MIBSSGS = reshape(y1_p_eq_GS, N, N);
 
    y2_MIBSSGS = reshape(y2_p_eq_GS, N, N);
    
end

end
%%

figure; 
plot(SNR11_p_eq_2, '-->', 'LineWidth', 2); hold on; plot(SNR22_p_eq_2, '--o', 'LineWidth', 2);  
plot(SNR11_p_eq_GS, '->', 'LineWidth', 2); plot(SNR22_p_eq_GS, '-o', 'LineWidth', 2);
legend('MI-BSS: fig1', 'MI-BSS: fig2', 'MI-BSS-GS: fig1', 'MI-BSS-GS: fig2', 'GraphJADE: fig1', 'GraphJADE: fig2')
xticks(1:length(a_vec)); xticklabels(a_vec); xlabel('\alpha'); ylabel('SNR Output (db)')
%%
figure;
subplot(2, 4, 1); imshow(s11); title('true image 1'); subplot(2, 4, 5); imshow(s22); title('true image 2')
subplot(2, 4, 2); imshow(x11, []); title('mixed image 1'); subplot(2, 4, 6); imshow(x22, []); title('mixed image 2');
subplot(2, 4, 3); imshow(y1_MIBSS, []); title('recovered image 1 by MI-BSS'); subplot(2, 4, 7); imshow(y2_MIBSS, []); title('recovered image 2 by MI-BSS');
subplot(2, 4, 4); imshow(y1_MIBSSGS, []); title('recovered image 1 by MI-BSS-GS'); subplot(2, 4, 8); imshow(y2_MIBSSGS, []); title('recovered image 2 by MI-BSS-GS');

%%
t1 = toc(t0);

disp(['>>>>> run-time: ', num2str(round(t1/60,2)), ' minutes']);


