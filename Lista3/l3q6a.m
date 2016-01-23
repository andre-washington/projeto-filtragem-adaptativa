clear all

% clf (figure(1))
% clf (figure(2))
% clf (figure(3))
% clf (figure(4))

rep = 500; % numero de repetições
SNR = 30;

% sinal 4 qam
constelation = [1+1i,1-1i,-1+1i,-1-1i];
data  = ceil(4*rand(rep,1));
s = transpose(constelation(data));  % sinal de treinamento
% figure(3);
% plot(1:length(s), s);



H_num = [0.5 1.2 1.5 -1];
H_den = 1;
sh = filter(H_num, H_den, s);
n_var_4qam = var(sh) * 10^(-SNR/10);

%ruído deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(1,rep)+1i*(randn(1, rep)));
% transposte deve ser usada no lugar do hermitiano
x=sh+transpose(n); % sinal desejado na entrada do equalizador + ruido

% part I: trainning.
w(1,1:15) = 0;
% xw(1,4) = 0;
xr = zeros(1, 15);
passo = 0.4;
gama = 10^(-6);
erro1 = zeros(1,500);
% figure(4);
% hold on;
for k=15:rep
    
    xr = [x(k) xr(1:15-1)]; % varre os elementos de x
    xw = xr*w';
%     plot(k, xw,'+g');
    e=s(k-14) - xw; %xr * w';
    erro1(k)=e;
    w = w + passo * xr * conj(e) / (xr * xr' + gama);
end

figure(1)
plot(1:500,real(erro1.*erro1));
xlabel('plot(e.*e)');

% part II: decision block included.

rep = 5000;
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
dataIn = randi([0 1],k*rep,1);
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

s2 = qammod(dataSymbolsIn,M,0);

sh2 = filter(H_num, H_den, s2);
n_var_16qam = var(sh2) * 10^(-SNR/10);


%ruído deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_16qam/2)*(randn(1,rep)+1i*(randn(1, rep)));
% transposte deve ser usada no lugar do hermitiano
x2=sh2+transpose(n); % sinal desejado na entrada do equalizador + ruido
xr = zeros(1, 15);
xd =  zeros(1,4995);
xw =  zeros(1,4995);
erro = zeros(1,5000);
for k=16:rep
    xr = [x2(k) xr(1:15-1)]; % varre os elementos de x
    xw(k) = xr*w';
    
    % decision block
%     xd(k) = round(real(xw(k))) + round(imag(xw(k)))*1i ;
    xd(k) = decisor(xw(k), 3);
%     if(real(xd(k))>3)
%         xd(k)=3+imag(xd(k))*1i;
%     end
%     if(real(xd(k))<-3)
%         xd(k)= -3 + imag(xd(k))*1i;
%     end
%     if(imag(xd(k))>3)
%         xd(k)= real(xd(k)) + 3i;
%     end
%     
%     if(imag(xd(k))<-3)
%         xd(k)= real(xd(k)) - 3i;
%     end
    % end of decision block
    
    e= xd(k) - xw(k);
    erro(k)=e;
    w = w + passo * xr * conj(e) / (xr * xr' + gama);
    
    
end
h = comm.ConstellationDiagram('Title','Customized Constellation for QAM', ...
    'XLimits',[-4 4],'YLimits',[-4 4], ...
    'ReferenceConstellation',s2, ...
    'ReferenceMarker','o','ReferenceColor',[0 1 0]);
step(h,xw(k)')
% figure(5)
% plot(1:length(xd), xd);

figure(2)
plot(1:5000,real(erro.*erro),'red');
axis([-1 rep -1 1]);
xlabel('plot(e.*e)');

[num, rate ]= biterr(abs(real(s2))',abs(real(xd)))









