clear;
close all;
clc;
n = 200;
m = 80;
x = linspace(-15,15,n);
t = linspace(0,8*pi,m); 
dt = t(2) - t(1); 
[Xgrid,T] = meshgrid(x,t);

%% Create two spatio-temporal patterns
f1 = 0.5*cos(Xgrid) .* (1+0*T).*sin(Xgrid);  % time-independent!
f2 = (sech(Xgrid).*tanh(Xgrid)) .* (2*exp(1j*2.8*T));

%% Combine signals and make data matrix
X = (f1 + f2)'; % Data Matrix
%X=ctranspose(X);
figure(1);
surfl(real(X')); 
shading interp; colormap('default'); view(-20,60);
title('Mixed signal')
%%  Simple example for first
X1 = tt_tensor(X(:,1:end-1));
X2 = tt_tensor(X(:,2:end));
[Phi,Lambda,omega,F]=DMD_TT(X1,X2,0.05);  % Tensor based dmd
%%
bg = find(abs(omega)<1e-2);
fg = setdiff(1:2, bg);
omega_fg = omega(fg); % foreground
Phi_fg = Phi(:,fg); % DMD foreground modes

omega_bg = omega(bg); % background
Phi_bg = Phi(:,bg); % DMD background mode
%%
%%{
%% Compute DMD Background Solution
b = Phi_bg \ X(:, 1);
X_bg = zeros(numel(omega_bg), length(t));
for tt = 1:length(t),
    X_bg(:, tt) = b .* exp(omega_bg .* t(tt));
end;
X_bg = Phi_bg * X_bg;
X_bg = X_bg(1:n, :);

figure(3);
surfl(real(X_bg')); 
shading interp; colormap('default'); view(-20,60);
title('Background signal')
%% Compute DMD FOreground Solution
b = Phi_fg \ X(:, 1);
X_fg = zeros(numel(omega_fg), length(t));
for tt = 1:length(t),
    X_fg(:, tt) = b .* exp(omega_fg .* t(tt));
end;
X_fg = Phi_fg * X_fg;
X_fg = X_fg(1:n, :);

figure(4);
surfl(real(X_fg')); 
shading interp; colormap('default'); view(-20,60);
title('Foreground signal')
size(b),size(X(:,1)),size(Phi_bg)
%}