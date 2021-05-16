%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close, clc

%%  For test
%spatio-temporal Generated signal  
%x=-5:0.1:5; y=-6:0.1:6; t=0:0.1:10*pi;
%[X,Y,T]=meshgrid(x,y,t);
%data = exp(-(X.^2+0.5*Y.^2)).*(cos(2*T))+(sech(X).*tanh(X).*exp(-0.2*Y.^2)).*sin(T);
%%
n_shots=100; % Number of snapshots 
%% Data for plotting Snapshots from .npyfile
data=readNPY('snapshots.npy'); % for reading npy files in matlab we need add this add this to matlab https://github.com/kwikteam/npy-matlab
xtt=tt_tensor(data(:,:,1:n_shots));
ytt=tt_tensor(data(:,:,2:n_shots+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Phi,evals,~,~] = DMD_TT(xtt,ytt,0.001);

%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% for the example using data from numpy file like in actuall paper plot
%%% different modes as below
% eigen modes corresponding to eigenvalues around 1
fhandle = plotCylinder(real(Phi(:,:,51)));
fhandle = plotCylinder(real(Phi(:,:,52)));
fhandle = plotCylinder(real(Phi(:,:,54)));
fhandle = plotCylinder(real(Phi(:,:,61)));
fhandle = plotCylinder(real(Phi(:,:,70)));
fhandle = plotCylinder(real(Phi(:,:,90)));
%}