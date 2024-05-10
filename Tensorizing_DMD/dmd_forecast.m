function [forecast] = dmd_forecast(data, n_shots, steps)

disp('Starting DMD..')

%data=ncread('\small_df.nc', 'TMAX');
%data1=permute(data, [2,1,3]);

%% Data
data = data(:,:,1:n_shots+1);
[mm,nn, TT] = size(data);
dt = 0.05;
t = (0:TT-1+steps)*dt;
%% Applying DMD
xtt=tt_tensor(data(:,:,1:end-1));
ytt=tt_tensor(data(:,:,2:end));
[ttPhi,Lambda,omega,F] = DMD_TT(xtt,ytt,dt);

%% Forecasting 
r = size(ttPhi,3);
x1 = data(:, :, 1);       % time = 0
[M,invS,invN,N] = PINV_TT(ttPhi, ttPhi.d-1);
b = full(invN)*invS*tensorprod(full(M), x1, [1 2]);

t_dyn = zeros(r, length(t));
for i = 1:length(t)
   t_dyn(:, i) = (b.*exp(omega*t(i))); 
end
forecast = real(tensorprod(full(ttPhi), t_dyn, [3], [1])); % the last 10 mode-3 slices are prediction maps.
forecast = forecast(:,:,n_shots+2:end);
end