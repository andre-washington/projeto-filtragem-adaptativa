clear all;  
close all;
clc

%% Questão 6 - item B - LMS

%% Canal
 canal=[0.5 1.2 1.5 -1]; % from chapter 3. 

%% Definições inicias
trainningSamples = [150 300 500];
errorNumTranning = [];
for sampleIndex = 1:3
    compEqualizador = 15;   % Comprimento do Equalizador
    % Delta = floor(length(canal) + floor(compEqualizador/2)); % Aproximação do atraso do canal + metade do comp. do equalizador 
    Delta = 8; % Aproximação do atraso do canal + metade do comp. do equalizador 
    mu = 0.001;               % Fator de passo
    epsilon = 1e-2;         

    treinamento = trainningSamples(sampleIndex);                  % comprimento de treinamento com dados BPSK (500 símbolos)
    opDirDec    = 5000;                 % comprimento da operação dirigida por decisão com dados QAM (5000 símbolos)
    N=treinamento + opDirDec;           % número total de símbolos

    numTentativas = 100; % Número de rodadas para se obter o erro médio (Monte carlo)
    errorTrial = zeros(numTentativas,N);


    for trial=1:numTentativas
    

        qam = 16;              % Escolher entre 4/16/64/256
        SNR = 30;                % SNR da saída do canal

        dadosTrein = zeros(1,treinamento+Delta); % Dados de Treinamento (515 dados)
        s=zeros(1,N);

        % Variância do ruído durante a operação dirigida por decisão
        % Símbolos QAM não possuem variância unitária QAM 
        switch qam
        case 4
           sigma_s = sqrt(2);
        case 16
           sigma_s = sqrt(10);
        case 64
           sigma_s = sqrt(42);
        case 256
           sigma_s = sqrt(170);
        end

        sigma_v_dd = sqrt((sigma_s)^2*norm(canal)^2/10^(SNR/10));   % durante a operação dirigida por decisão
        sigma_v_tt = sqrt(norm(canal)^2/10^(SNR/10));               % durante treinameto

        v = zeros(1,N);   % Definição do ruído v durante as duas fases
        v(1:treinamento) = (sigma_v_tt/2)*(randn(1,treinamento)+1j*randn(1,treinamento));
        v(treinamento:N) = (sigma_v_dd/2)*(randn(1,N-treinamento+1)+1j*randn(1,N-treinamento+1));

        % Criando os dados...
        s(1:treinamento) = (sign(randn(1,treinamento))+1j*sign(randn(1,treinamento)))/(sqrt(2)); 
        dadosTrein(Delta+1:treinamento+Delta) = s(1:treinamento); % dados de treinamento 
        % geração dados QPSK. Possuem variância unitária 

        for i = 1:opDirDec
         x     = rand;   
         xint  = round((qam-1)*x);                    % gerando inteiros, entre 0 e 255 
         xcoor = modmap(xint,1,1,'qam',qam);          % dando coordenadas reais e imaginárias 
         s(treinamento+i) = xcoor(1) + 1j*xcoor(2);    % combinando coordenadas em números complexos
        end 

        y = filter(canal,1,s);  % dados passam pelo canal
        r = y+v;                % adiciona ruído 

        %% Inicialização da equalização
        w  = zeros(compEqualizador,1);  %  coeficientes do equalizador 
        u  = zeros(compEqualizador,1);  %  vetor de regressãoregressor vector (row) 
        e  = zeros(1,N);                %  vetor de erros
        num_erros=0;

        %% Equalização adaptativa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modo de treinamento
        for i = 1:treinamento+Delta  % Apenas 515 símbolos
           u  = [r(i);u(1:compEqualizador-1)];
           estimativa_s(i) = w'*u;
           d(i) = dadosTrein(i); 
           e(i) = d(i)- estimativa_s(i);
           w = w + mu*conj(e(i))*u;  %% LMS complexo
           %w = w + mu*conj(e(i))*u/(norm(u)^2+epsilon);  %% LMS normalizado
        end

        % Operação dirigida por decisão
        for i = treinamento+Delta+1:N
           u  = [r(i);u(1:compEqualizador-1)];
           %estimativa_s(i) = u*w;
           estimativa_s(i) = w'*u;
           switch qam 
            case 4
              check_s(i)=slicer4(estimativa_s(i));
            case 16
              check_s(i)=slicer16(estimativa_s(i));
            case 64
              check_s(i)=slicer64(estimativa_s(i));
            case 256
              check_s(i)=slicer256(estimativa_s(i));
           end

           d(i) = check_s(i); % Decisão dirigida por decisão
           e(i) = d(i)- estimativa_s(i);
           w = w + mu*conj(e(i))*u;  %% LMS complexo
           %w = w + mu*conj(e(i))*u/(norm(u)^2+epsilon);  %% LMS normalizado
            %    if (check_s(i) ~= s(i-Delta))
            %      num_erros = num_erros + 1;
            %    end
        end
        errorTrial(trial,:) = e(:);
        estimativa_treinamento{sampleIndex} = estimativa_s;
    end
    %keyboard;
    errorNumTranning{sampleIndex} = mean(errorTrial.*conj(errorTrial),1);
end
%keyboard;
figure(1)
for sampleIndex = 1:3
    plot(log10(errorNumTranning{sampleIndex}));hold on;   
end
 legend(...
 sprintf('%d Amostras',trainningSamples(1)),...
 sprintf('%d Amostras',trainningSamples(2)),...
 sprintf('%d Amostras',trainningSamples(3))...
 );hold on;
title('Evolução do erro médio quadrático');
ylabel('log_{10}{(MSE)}');
xlabel('Iterações');

%%%% Representação da evolução do erro de estimação 
  
%%%% Representação das constelações dos sinais transmitidos, recebidos e
%%%% equalizados
 
figure(2)
subplot(131); 
plot(real(estimativa_treinamento{1}(trainningSamples(1)+Delta+1:end)),imag(estimativa_treinamento{1}(trainningSamples(1)+Delta+1:end)),'rx');grid; 
title('Saída do Equalizador (N=150)');
axis([-4 4 -4 4]);
subplot(132); 
plot(real(estimativa_treinamento{2}(trainningSamples(2)+Delta+1:end)),imag(estimativa_treinamento{2}(trainningSamples(2)+Delta+1:end)),'rx');grid; 
title('Saída do Equalizador (N=300)');
axis([-4 4 -4 4]);
subplot(133); 
plot(real(estimativa_treinamento{3}(trainningSamples(3)+Delta+1:end)),imag(estimativa_treinamento{3}(trainningSamples(3)+Delta+1:end)),'rx');grid; 
title('Saída do Equalizador (N=500)');
axis([-4 4 -4 4]);

%keyboard;


 

