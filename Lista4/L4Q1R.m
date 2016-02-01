clear, clc, close all;

num_ite = 10;
num_taps = 3;
input = cos(pi/3*(1:num_ite+1))';
lambda = 0.98;
delta = 100;
w = [ 1 0 0]';
sd = delta*eye(num_taps);
pd = zeros(num_taps, 1);

x = zeros(num_taps, 1);
err_vec = zeros(num_ite, 1);
y = zeros(num_ite, 1);
for k = 1:num_ite 
   d = input(k+1);
   x = [input(k); x(1:end-1)]; 
   e = d - x'*w;
   phi = sd*x;
   
   sd = 1/lambda*(sd - (phi*phi')/(lambda + phi'*x));
   
   w = w + e*sd*x;
   
   y(k) = w'*x;
   err_vec(k) = d - y(k);
end

plot(input(2:end), 'r');
hold on;
plot(y, 'b--');