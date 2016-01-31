
% clc
clear

num_ite = 3000; % iterations number
num_taps = 12; % quantity of filter taps

input_var = 1; %input variance
input_mean = 0; %input mean

noise_var = 10^-3;% noise variance
noise_mean = 0; %noise mean

% using eq 3.37 of textbook and the fact that the number of coefficents of
% the system and the filter are equivalent
jmin = noise_var;

input = sqrt(input_var)*randn(num_ite, 1) + input_mean; %sinal de entrada
noise = sqrt(noise_var)*randn(num_ite, 1) + noise_mean; %ruído

% finding mu max 
mu_max = 1/trace(input_var.*eye(num_taps));
mu = [mu_max/2, mu_max/10, mu_max/50]; 
idx = [2, 10, 50];

% defining the system according to the given transfer function
h_num = [1 0 0 0 0 0 0 0 0 0 0 0 -1];
h_den = [1 -1];
output_sys = filter(h_num, h_den, input);

% reference signal is equal to system output + measurement noise
ref_sig = output_sys + noise;

%break_point = num_ite;
exp_mse = zeros(1,3);
misadj = zeros(2,3);

for j=1:3 % for each mu
  %error vector
error_vec = zeros(num_ite, 1);

% filter output vector
output_flt = zeros(num_ite, 1);

% initialiazing the filter taps
lms_flt = zeros(num_taps, 1);

    for i=num_taps:num_ite   
       output_flt(i) = dot(input(i-11:i), lms_flt);
       error_vec(i) = ref_sig(i) - output_flt(i);
       lms_flt = lms_flt + 2*mu(j)*input(i-11:i)*error_vec(i);
        
       %if(error_vec(i)^2 < 10^-8)
       %     break_point = i;
       %     break;
       %end
       %if(break_point > 200)
       exp_mse(1,j) = mean(error_vec(num_ite-200:num_ite).^2);
       %else
       %     exp_mse(1,j) = mean(error_vec(1:break_point).^2);
       %end
       
       misadj(1,j) = abs(exp_mse(1,j)-jmin)/jmin;
       misadj(2,j) = (mu(j)*trace(input_var.*eye(num_taps)))/ ...
           (1-mu(j)*trace(input_var.*eye(num_taps)));
    end
    
figure(2*j-1)
semilogy(error_vec.^2);
title(sprintf('Comportamento do Erro Quadrático LMS (1/%d*mu-max)', idx(j)));
xlabel('Num. de Iterações');%axis([0 500 10^-10 10^2]);
hold on;
ylabel('Erro Quadrático');
plot(ones(1, num_ite)*jmin);
legend('LMS', 'Erro Mínimo');
hold off;


figure(2*j)
freqz(lms_flt);
title(sprintf('Resposta em Frequência do Filtro FIR (1/%d*mu-max)', idx(j)));

%plot(output_flt(1:break_point));
%hold on;
%plot(ref_sig(1:break_point), 'r--');

end

figure(7)
freqz(h_num, h_den);
title('Resposta em Frequência do Sistema');