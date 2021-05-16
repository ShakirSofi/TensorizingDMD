function [M,invS,invN,N] = PINV_TT(tt,l)
% Pinv of (Lth) unfolding using algorithm in tensor_train based dynamic mode decomposition 
%%TT_inverse of reduced matrix
% dimensionality of tensor
d=tt.d; 
% rank of tensor
r=tt.r; 
% tensor core
cr=tt.core; 
% core position in tt.core
ps=tt.ps;
n=tt.n;

% compution for (L-th) tensor core
core_l =  cr(ps(l):ps(l+1)-1); 
core_l = reshape(core_l,[r(l)*n(l),r(l+1)]);
% SVD 

[~,S,~] = svd(core_l,'econ');
invS = diag(1./diag(S)); 
% compute M

M = tt_tensor;
M.d = l; M.r = r(1:l+1); M.n = n(1:l);
M.core = cr(1:ps(l+1)); M.ps = ps(1:l+1);
% compute N and invN using formulate provided in pdf file
 
N = tt_tensor;
N.d = d-l; N.r = r(l+1:end); N.n = n(l+1:end);
N.core = cr(ps(l+1):end);
N.ps = ps(l+1:end) - ps(l+1) + 1;
tmp = pinv(full(N)); 
invN = tt_matrix(tmp) ;
% Return Results M,N, invN, invS
end
