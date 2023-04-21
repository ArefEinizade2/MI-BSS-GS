function [y, B] = MI_BSS_GS(x, L, param)
% Implements the MI-BSS-KG method for blind separation of smooth graph signal sources with known graphs, as proposed in our paper.
% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Robust blind separation of smooth graph signals using minimization of graph regularized mutual information." Digital Signal Processing 132 (2022): 103792.
%%
% Usage: 
%   >> [y, B] = MI_BSS_GS(x, L, param);

%   x: nxT, n is the number of sources, T is the number of observed temporal samples
%   L: a n-length cell containg the original Laplacian matrices of the graph signal sources, where L{i} has the size of T x T
%   param.mu: the learning rate used in the Gradient Descent step, e.g., 0.1
%   param.lambda: balances the weight between minimization of mutual information and graph smoothness terms, e.g., 1e-2
%   param.Tol: tolerance for reaching convergence
%   param.MaxIter: maximum iterations for convergence

%   y: nxT, the estimated sources
%   B: nxn, the estimated unmixing matrix
%% Initialization:

n = size(x,1);

B = eye(size(x,1));

for iiii = 1 : n

    Phi{iiii} = x*L{iiii}*x'/trace(L{iiii});

end

L_temp = L;

lambda_e1_e2 = 0;

mu = param.mu;
lambda = param.lambda;
Tol = param.Tol;
MaxIter = param.MaxIter;

error = zeros(1, MaxIter);
%% First applying MI-BSS in small iterations for an initial estimation:
for iter1 = 1 : 1000

    disp(['MI-BSS-GS iter ', num2str(iter1)])

    y = B * x;

    t_nonGraph = tic;

    Gamma = zeros(size(x));

    for i = 1:n-1

        for j = i+1:n

            Beta_star = SFD([y(i,:); y(j,:)]);

            Beta = Beta_star;

            Gamma(i,:) = Gamma(i,:) + Beta(1,:);

            Gamma(j,:) = Gamma(j,:) + Beta(2,:);

        end

    end

    MI_grad = Gamma*x'/size(x,2);

    new_B = B - mu * MI_grad;

    sigma = std(y, [], 2);

    B = new_B ./ sigma;

    if error(iter1) < 1e-2   
        break;  
    end


end
%% Sequential Graph Assignment step

Set_idx_Smooth = 1 : n;

Set_idx_Sources = randperm(n);

L_temp = {};

for p1 = Set_idx_Sources

   smooth_term = [];

   for p2 = Set_idx_Smooth

       smooth_term = [smooth_term, B(p1, :) * Phi{p2} * B(p1, :)'];                    

   end

    [~, min_idx] = min(smooth_term);                

    L_temp{p1} = L{Set_idx_Smooth(min_idx)};                

    Set_idx_Smooth = setdiff(Set_idx_Smooth, Set_idx_Smooth(min_idx));

end
   % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%% Implementing MI-BSS-GS using a loop till convergence:

for iter = iter1 + 1 : param.MaxIter

    disp(['MI-BSS-GS iter ', num2str(iter)])

    y = B * x;

    Gamma = zeros(size(x));

    % Pairwise estimation of Score Difference Functions:      
    for i = 1:n-1

        for j = i+1:n

            Beta_star = SFD([y(i,:); y(j,:)]);

            Beta = Beta_star;

            Gamma(i,:) = Gamma(i,:) + Beta(1,:);

            Gamma(j,:) = Gamma(j,:) + Beta(2,:);

        end

    end

    MI_grad = Gamma*x'/size(x,2);

    Theta = zeros(size(B));

    for iii = 1:size(x,1)

        Theta(iii,:) = B(iii,:)*x*L_temp{iii}*x';

    end

    MI_grad = MI_grad + ((lambda) * (Theta));

    % Gradient Descent Step:
    new_B = B - mu * MI_grad;

    sigma = std(y, [], 2);

    new_B = new_B ./ sigma;

    error(iter) = norm((new_B - B), 'fro')/norm(B,'fro');      

    B = new_B;

   if error(iter) < Tol  
        break;  
   end


end
%%
y = B * x;

end