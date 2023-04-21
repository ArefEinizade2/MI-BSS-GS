function [y, B] = MI_BSS(x, param)
% Implements the MI-BSS method for blind separation of sources using minimization of mutual information, as proposed in the following paper.
% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% Babaie-Zadeh, Massoud, and Christian Jutten. "A general approach for mutual information minimization and its application to blind source separation." Signal Processing 85, no. 5 (2005): 975-995.
%%
% Usage: 
%   >> [y, B] = MI_BSS(x, param);

%   x: nxT, n is the number of sources, T is the number of observed temporal samples
%   param.mu: the learning rate used in the Gradient Descent step, e.g., 0.1
%   param.Tol: tolerance for reaching convergence
%   param.MaxIter: maximum iterations for convergence

%   y: nxT, the estimated sources
%   B: nxn, the estimated unmixing matrix

%% Initialization:

n = size(x,1);

B = eye(size(x,1));

mu = param.mu;
Tol = param.Tol;
MaxIter = param.MaxIter;
%% Gradient Descent loop:   

for iter = 1 : MaxIter

    disp(['MI-BSS iter ', num2str(iter)])

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