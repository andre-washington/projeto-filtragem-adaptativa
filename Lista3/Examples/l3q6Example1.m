clear all
close all

%simulation length
N = 1000;

%channel length
M = 5;

%number of independent trials
T = 100;

cascade_impulse_response = zeros(1,2*M-1);
for j = 1 : T 
%training signal
u = randn(1,N);

%channel to be equalized 
c = randn(M,1);
c = c / norm(c);

%channel output
z = filter(c,1,u);

%additive noise to the channel output
SNR = 30;
var_v = var(z) * 10^(-SNR/10);
v = var_v^0.5 * randn(1,N);

%input to the equalizer
x = z + v;

%NLMS channel equalization
w = zeros(M,1);
x_regressor = zeros(1,M);
step = 0.1;
epsilon = 10^(-6);
for k = 4 : N
x_regressor = [x(k) x_regressor(1:M-1)];
e = u(k-3) - x_regressor * w;
w = w + step * x_regressor' * e / (x_regressor * x_regressor' + epsilon);
end
cascade_impulse_response = cascade_impulse_response + conv(w,c)';
display(j);
end
figure;
stem(cascade_impulse_response/T);
title('cascade channel-equalizer impulse response');
xlabel('taps');