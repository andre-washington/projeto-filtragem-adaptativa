clear, clc, close all

num_taps = 15; % number of filter taps
num_symt = 500; % number of symbols of training mode (tm)
num_symd = 5000; % number of symbols of control-by-decision mode (dm)
const_sizet = 4; % constellation syze of tm
const_sized = 256; % constellation size of dm

num_ite = num_symt + num_symd; % number of total iterations

SNR = 30; % given signal to noise ratio

mu = 0.4; % given step factor
gama = 10^-6; % 
   
%generating the data signal for training (4QAM) and
%control-by-decision(16QAM)
data_t = randi(const_sizet, num_symt, 1) - 1; 
data_d = randi(const_sized, num_symd, 1) - 1;

% modeling the data with 4QAM and 16QAM

signal_t = qammod(data_t, const_sizet);
signal_d = qammod(data_d, const_sized);

% channel input signal
signal = [signal_t; signal_d];

% channel coefficients
H_num = [0.5 1.2 1.5 -1];
H_den = 1;

% computing channel output
ch_out = filter(H_num, H_den, signal);

% Computing the noise variance for each constellation map for 30db SNR
% Normalization factor QAM sqrt( 2/3(M-1) ), M being the constellation size

n_var_4qam = (norm(H_num).^2* 2/3*(const_sizet-1)) * 10^(-SNR/10);
n_var_16qam = (norm(H_num).^2* 2/3*(const_sized-1)) * 10^(-SNR/10);   

%the noise should be complex and the variance divided for each part (real/imag)
%noise for the training part
noise_t = sqrt(n_var_4qam/2)*(randn(num_symt, 1)+1i*(randn(num_symt, 1)));

%noise for the control by decision part
noise_d = sqrt(n_var_16qam/2)*(randn(num_symd, 1)+1i*(randn(num_symd, 1)));

% computing the equalizer input (adding awgn noise to channel output)
x = ch_out + [noise_t ; noise_d]; 

w = zeros(num_taps, 1); % initialiazing the equalizer taps

eq_out = zeros(1, num_ite); % vector for storage of equalizer outputs
err_vec = zeros(1,num_ite); % vector for error storage

init = zeros(num_taps - 1, 1); % vector for convolution sliding window

% part I: trainning mode.

% it's necessary to wait for a number of samples before perform the
% comparation with the input (section 2.10.4, page 57, Diniz)
% Half of filter length
delay = ceil(num_taps/2);
for k = 1:num_symt  
	inp = [x(k); init];
    init = inp(1:end-1);
    eq_out(k) = w'*inp;
    if(k > delay)
    err_vec(k) = signal(k-delay) - eq_out(k); % Compara o sinal de entrada com saida (adiantado). Por exemplo(sinal(s) 1 com saida (eq_out) 8)
    w = w + mu*conj(err_vec(k))*inp/(inp'*inp+gama);
    end
end

% ploting the equalizer output in time
figure(1)
plot3(real(eq_out(1:num_symt)),imag(eq_out(1:num_symt)),1:num_symt,'r.'); 

% part II: control by decision mode. 

for k = num_symt+1:num_symt+num_symd % É nessessario esperar um momento (amostras) até fazer a comparação com o sinal (Seção 2.10.4 pag.57 Diniz) . Nesse caso 7 ou 8(metade do comprimento do filtro), pois tem comprimento 15.
	inp = [x(k); init];
    init = inp(1:end-1);
    
    eq_out(k) = w'*inp;
   
    d = qam_decisor(eq_out(k), const_sized);
   
    err_vec(k) = d - eq_out(k); 
    w = w + mu*conj(err_vec(k))*inp/(inp'*inp+gama);
    
end

figure(3)
plot3(real(eq_out(num_symt+1:end)),imag(eq_out(num_symt+1:end)),1:num_symd,'r.'); % plota da modulação 256 QAM

figure(4)
semilogy(1:num_ite, conj(err_vec).*err_vec, 'red');
xlabel('plot(e.*e)');
% axis([-1 600 -1 1]);

