function [Phi,Lambda,omega,F] = DMD_TT(ttX,ttY,dt)
%% This code is according Stefan Klus (paper: tensor-train dynamic mode decomposition)
% Input:
% ttX: input sequence X in TT-format
% ttY: input sequence Y in TT-format
% Result:
% Phi: Exact dmd modes
% Lambda: eigenvalues in dmd
% omega: temporal frequency 

% parameters of ttX
d = ttX.d ; r = ttX.r ; n = ttX.n ; 
cr = ttX.core ; ps = ttX.ps ;

% Compute M, invS, invN from ttX  % See Tensor train _dynamic mode
% decomposition
[M,invS,invN,N] = PINV_TT(ttX,d-1) ;

% Compute P from ttY
P = tt_tensor;
P.core = ttY.core(1:ttY.ps(end-1));
P.ps = ttY.ps(1:end-1);
P.d = ttY.d-1;
P.n = ttY.n(1:end-1);
P.r = ttY.r(1:end-1);

% Compute Q from ttY
Q = tt_tensor;
Q.core = ttY.core(ttY.ps(end-1):ttY.ps(end)-1);
Q.ps = ttY.ps(end-1:end) - ttY.ps(end-1) + 1;
Q.d = 1;
Q.n = ttY.n(end);
Q.r = ttY.r(end-1:end);

% Compute M^T . P
MtP = dot(ctranspose(M),P); % rd*sd 
% MtP = full(M)'*full(P) ; % same

% Compute Q . N^+ - contraction of  core with orthonormized N^T
% QNinv = dot(Q,invN); % sd*rd 
QNinv = full(Q)*full(invN); % sd*rd 

% Compute reduced matrix A 
A = MtP*QNinv*invS ;

% Compute W - eignvector, lambda - eignvalue of reduced A
[W, Lambda] = eig(A);
omega = log(diag(Lambda))/dt/(2*pi) ;
F = abs(imag(omega));

% Compute Exact TT-DMD mode
tt_Phi = ttY;
tmp = QNinv*invS*W*pinv(Lambda); 
tt_Phi.core(tt_Phi.ps(end-1):end) = [];
tt_Phi.core = cat(1,tt_Phi.core, tmp(:)) ;
tt_Phi.ps(end) = ttY.ps(end-1) + length(tmp(:));
tt_Phi.n(end) = size(Lambda,1);
Phi = full(tt_Phi,tt_Phi.n');


