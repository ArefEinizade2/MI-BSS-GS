function Beta = SFD(y)
% Implements the Score_Function_Difference (SFD) for two estimated sources as proposed in the following paper:
% Babaie-Zadeh, Massoud, and Christian Jutten. "A general approach for mutual information minimization and its application to blind source separation." Signal Processing 85, no. 5 (2005): 975-995.
%% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Robust blind separation of smooth graph signals using minimization of graph regularized mutual information." Digital Signal Processing 132 (2022): 103792.
%%
% Usage: 
%   >> Beta = SFD(y);
%   y: nxT, n is the number of sources and must be 2, T is the number of temporal samples
%   Beta: nxT: the Score_Function_Difference
%%
Psi = MSF(y);

Phi = JSF(y);

Beta = Psi - Phi;

end