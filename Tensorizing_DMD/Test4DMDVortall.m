%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close, clc

%%  For test
%spito-temporal Generated signal  
%x=-5:0.1:5; y=-6:0.1:6; t=0:0.1:10*pi;
%[X,Y,T]=meshgrid(x,y,t);
%data = exp(-(X.^2+0.5*Y.^2)).*(cos(2*T))+(sech(X).*tanh(X).*exp(-0.2*Y.^2)).*sin(T);

load VORTALL
data=VORTALL;

data=reshape(data,numel(data)/(151*199),199,151);
size(data)
%%
n_shots=51; % Number of snapshots 
%% Data Vortall
xtt=tt_tensor(data(:,:,1:n_shots));
ytt=tt_tensor(data(:,:,2:n_shots+1));
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Phi,evals,~,~] = DMD_TT(xtt,ytt,0.005);
toc
%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
size(Phi)
fhandle = plotCylinder(real(reshape(Phi(:,:,1),199,449)));
fhandle = plotCylinder(real(reshape(Phi(:,:,2),199,449)));
fhandle = plotCylinder(real(reshape(Phi(:,:,5),199,449)));
fhandle = plotCylinder(real(reshape(Phi(:,:,7),199,449)));
