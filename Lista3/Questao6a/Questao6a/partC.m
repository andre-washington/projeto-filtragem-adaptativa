clear all;  
close all;
clc

%% Questão 6 - item B - LMS

%% Canal
 canal=[0.5 1.2 1.5 -1]; % from chapter 3. 

%% Definições inicias
trainningSamples = [500];
errorNumTranning = [];
for sampleIndex = 1:1
    compEqualizador = 15;   % Comprimento do Equalizador
    %Delta = floor(length(canal) + floor(compEqualizador/2)); % Aproximação do atraso do canal + metade do comp. do equalizador 
    Delta = 8; % Aproximação do atraso do canal + metade do comp. do equalizador 
    mu = 0.4;               % Fator de passo
    epsilon = 1e-2;         
    vardelta = 0.1;
    treinamento = trainningSamples(sampleIndex);                  % comprimento de treinamento com dados BPSK (500 símbolos)
    opDirDec    = 5000;                 % comprimento da operação dirigida por decisão com dados QAM (5000 símbolos)
    N=treinamento + opDirDec;           % número total de símbolos
    alpha = 0.8;
    numTentativas = 100; % Número de rodadas para se obter o erro médio (Monte carlo)
    errorTrial = zeros(numTentativas,N);
    errorTrialN = zeros(numTentativas,N);


    for trial=1:numTentativas
    

        qam = 256;              % Escolher entre 4/16/64/256
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
        %debug = [];
        for i = 1:opDirDec
         x     = rand;   
         xint  = round((qam-1)*x);                    % gerando inteiros, entre 0 e 255 
         xcoor = modmap(xint,1,1,'qam',qam);          % dando coordenadas reais e imaginárias 
         s(treinamento+i) = xcoor(1) + 1j*xcoor(2);    % combinando coordenadas em números complexos
        %debug = [debug xcoor];
        end 

        y = filter(canal,1,s);  % dados passam pelo canal
        r = y+v;                % adiciona ruído 

        %% Inicialização da equalização
        w  = zeros(compEqualizador,1);  %  coeficientes do equalizador 
        wN = zeros(compEqualizador,1);  %  coeficientes do equalizador LMS-Newton
        u  = zeros(compEqualizador,1);  %  vetor de regressãoregressor vector (row) 
        e  = zeros(1,N);                %  vetor de erros
        eN = zeros(1,N);
        num_erros=0;
        Rinv = vardelta*eye(compEqualizador);

        %% Equalização adaptativa %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modo de treinamento
        for i = 1:treinamento+Delta  % Apenas 515 símbolos
           u  = [r(i);u(1:compEqualizador-1)];
           estimativa_s(i) = w'*u;
           estimativa_s_newton(i) = wN'*u;
           d(i) = dadosTrein(i); 
           e(i) = d(i)- estimativa_s(i);
           eN(i) = d(i)- estimativa_s_newton(i);
           Rinv = (1/(1-alpha))*(Rinv - (Rinv*u*u'*Rinv)/(((1-alpha)/alpha)+u'*Rinv*u));
           % w = w + mu*conj(e(i))*u;  %% LMS complexo
           w = w + mu*conj(e(i))*u/(norm(u)^2+epsilon);  %% LMS normalizado
           wN = wN + 2*mu*conj(eN(i))*Rinv*u; %% LMS Newton
        end

        % Operação dirigida por decisão
        for i = treinamento+Delta+1:N
           u  = [r(i);u(1:compEqualizador-1)];
           %estimativa_s(i) = u*w;
           estimativa_s(i) = w'*u;
           estimativa_s_newton(i) = wN'*u;
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
           eN(i) = d(i)- estimativa_s_newton(i);
           Rinv = (1/(1-alpha))*(Rinv - (Rinv*u*u'*Rinv)/(((1-alpha)/alpha)+u'*Rinv*u));
           % w = w + mu*conj(e(i))*u;  %% LMS complexo
           w = w + mu*conj(e(i))*u/(norm(u)^2+epsilon);  %% LMS normalizado
           wN = wN + 2*mu*conj(eN(i))*Rinv*u; %% LMS Newton
        end
        errorTrial(trial,:) = e(:);
        errorTrialN(trial,:) = eN(:);
    end
end
figure(1)
semilogy(mean(errorTrial.*conj(errorTrial),1));hold on;
semilogy(mean(errorTrialN.*conj(errorTrialN),1));hold on;
legend('LMS-Norm','LMS-Newton');hold on;
ylabel('MSE');
xlabel('Iterações');
title('Evolução do Erro Médio Quadrático');

%%%% Representação da evolução do erro de estimação 
  
%%%% Representação das constelações dos sinais transmitidos, recebidos e
%%%% equalizados
 
figure(2)
subplot(131);  
plot(real(s(treinamento+Delta+1:N)),imag(s(treinamento+Delta+1:N)),'ro');
grid; 
title('Sequência Transmitida'); 
axis([-16 16 -16 16]);
subplot(132); 
plot(real(estimativa_s(treinamento+Delta+1:N)),imag(estimativa_s(treinamento+Delta+1:N)),'rx');grid; 
title('Saída do Equalizador -- LMS-Norm');
axis([-16 16 -16 16]);
subplot(133); 
plot(real(estimativa_s_newton(treinamento+Delta+1:N)),imag(estimativa_s(treinamento+Delta+1:N)),'rx');grid; 
title('Saída do Equalizador-- LMS-Newton');
axis([-16 16 -16 16]);
% Saída do equalizador a 4-QAM
samples = 49;
showPlots = [50 100 150 200];
figure(3)
for sampleIndex=1:4
   subplot(2,2,sampleIndex); 
   plot(real(estimativa_s(showPlots(sampleIndex)-samples:showPlots(sampleIndex))),imag(estimativa_s(showPlots(sampleIndex)-samples:showPlots(sampleIndex))),'rx');grid; 
    title(sprintf('Saída do Equalizador para N = %d (4-QAM)',showPlots(sampleIndex)));
end

% Representação da constelação no treinamento versus iteração
figure(4)
plot3(real(estimativa_s_newton(Delta+1:treinamento)),imag(estimativa_s_newton(Delta+1:treinamento)),Delta+1:treinamento,'rx');grid; 
zlabel('Número de iterações');
title('LMS-Newton - Saída do Equalizador - Temporal');

figure(5)
plot3(real(estimativa_s(Delta+1:treinamento)),imag(estimativa_s(Delta+1:treinamento)),Delta+1:treinamento,'rx');grid; 
zlabel('Número de iterações');
title('LMS-Norm - Saída do Equalizador - Temporal');


 

