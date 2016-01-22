clear, clc, close all;

num_ite = 100;
num_taps = [2 3];
signal = sin((1:num_ite+1))';
SNR = [3 1000];

lambda =  [0.90 0.1];
%lambda = 0.5;
%SNR = 30;
%num_taps= 3;
idx = 1;
for i = 1:length(num_taps)
    for j = 1:length(SNR)
        for l = 1:length(lambda)
            n_var = var(signal) * 10^(-SNR(j)/10);
            
            %ruido deve ser complexo e a variancia divida para cada parte (real/imag)
            noise = sqrt(n_var)*(randn(num_ite+1, 1));
            % transposte deve ser usada no lugar do hermitiano
            
            input = signal + noise; % sinal desejado na entrada do equalizador + ruido
            
            delta = 1/var(input);
            
            sd = delta*eye(num_taps(i));
            pd = zeros(num_taps(i), 1);
            
            x = zeros(num_taps(i), 1);
            err_vec = zeros(num_ite, 1);
            y = zeros(num_ite, 1);
            
            for k = 1:num_ite
                d = input(k+1);
                if (k < num_taps(i))
                    x = [input(k); x(1:end-1)];
                else
                    x = input(k:-1:k-num_taps(i)+1);
                end
                sd = 1/lambda(l)*(sd - (sd*x*x'*sd)/(lambda(l) + x'*sd*x));
                
                pd = lambda(l)*pd + d*x;
                
                w = sd*pd;
                
                y(k) = w'*x;
                err_vec(k) = d - y(k);
            end
            
            figure(idx)
            plot(input(2:end), 'r');
            title(['Num de Coef =' num2str(num_taps(i)) ' SNR =' num2str(SNR(j)) ' Fator de esq =' num2str(lambda(l))]);
            hold on;
            plot(y, 'b--');
            hold off
            
            idx = idx+1;
            
            figure(idx)
            plot(err_vec.^2);
            title(['Num de Coef =' num2str(num_taps(i)) ' SNR =' num2str(SNR(j)) ' Fator de esq =' num2str(lambda(l))]);
            idx = idx+1;
            
            clear w;
        end
    end
end