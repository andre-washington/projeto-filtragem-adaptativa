% L. 3 Q. 4. 
clear, clc, close all;

amostras = 500;
H = [1; 1.6];
dados = randn(amostras,1);
x = conv(dados,H);
% filter coefficient inicialization
w_gradiente =  [0;0];
w_newton = w_gradiente;
w_lms = w_gradiente;
w_lms_normalized = w_gradiente;

% auto-correlation matrix
R = zeros(2,2);
data = randn(amostras,1);
xr = conv(data,H);
for i = 2:length(xr)
    v = [xr(i); xr(i - 1)];
    R = R + v * v';
end
R = R/(amostras);

% correlation vector
data = randn(amostras,1);
xp = conv(data,H);

p = zeros(2,1);
for i = 1:length(xp)-1
    v = [xp(i); xp(i + 1)];
    p = p + v * data(i);
end
p = p/(length(x)-2);
p = flipud(p);


w_opt = R \ p; % calculate the optimum filter
Cr = cond(R); % condition number
wf_opt = R\ [corr(x(2:length(x)), x(1:length(x)-1) ); 0] ;%slide 96 filtragem otima
% wf_opt = [ w_opt(1) w_opt(2);w_opt(2) w_opt(1)]*p
%  plots the zeros of w_opt
[hz1,hp,ht]=zplane(roots(w_opt)); hold on
[hz3,hp3,ht3]=zplane(roots(wf_opt)); 
%  plots the zeros of the channel
[hz2,hp2,ht2]=zplane(roots(H));hold off
axis([-2 2 -2 2]);
set(findobj(hz1, 'Type', 'line'), 'Color', 'b'); 
set(findobj(hz2, 'Type', 'line'), 'Color', 'r'); 
set(findobj(hz3, 'Type', 'line'), 'Color', 'g'); 
legend([hz1,hz2, hz3], 'Zero do filtro', 'Zero do canal','Zero do preditor');

% error surface
m = -0.3:0.01:0.8;
n = m;
J = zeros(length(m),length(n));
for i=1:length(m)
    for j=1:length(n)
        aux = [m(i);n(j)];
        J(i,j) = 1 - 2*aux'*p + aux'*R*aux;
    end
end
J = J';

figure(2)
subplot(2,2,1); contour(m,n,J,18);
title('LMS');xlabel('w(1)'); ylabel('w(2)');
subplot(2,2,2); contour(m,n,J,18);
title('NLMS');xlabel('w(1)'); ylabel('w(2)');
subplot(2,2,3); contour(m,n,J,18);
title('Gradiente Determinístico');xlabel('w(1)'); ylabel('w(2)');
subplot(2,2,4); contour(m,n,J,18);
title('Newton');xlabel('w(1)'); ylabel('w(2)');

figure(4)
subplot(2,1,1);
contour(m,n,J,8);
title('RLS');

figure(2)
subplot(2,2,1); hold on;
mse_lms = do_algorithm(amostras, dados, x, w_lms, p, R, 'lms');
hold off;
subplot(2,2,2); hold on;
mse_nlms = do_algorithm(amostras, dados, x, w_lms_normalized, p, R, 'nlms');
hold off;
subplot(2,2,3); hold on;
mse_gradient = do_algorithm(amostras, dados, x, w_gradiente, p, R, 'grad');
hold off;
subplot(2,2,4); hold on;
mse_newton = do_algorithm(amostras, dados, x, w_newton, p, R, 'newton');
hold off;
% rls...
figure(4)
subplot(2,1,1);hold on;
mse_rls = do_algorithm(amostras, dados, x, w_newton, p, R, 'rls');
hold off;

figure(3)
subplot(2,2,1);
semilogy(mse_lms);
title('LMS'); xlabel('amostras'); ylabel('MSE');
subplot(2,2,2);
semilogy(mse_nlms);
title('NLMS'); xlabel('amostras'); ylabel('MSE');
subplot(2,2,3);
semilogy(mse_gradient);
title('Gradiente Determinístico'); xlabel('amostras'); ylabel('MSE');
subplot(2,2,4);
semilogy(mse_newton);
title('Newton'); xlabel('amostras'); ylabel('MSE');

figure(4)
subplot(2,1,2);
semilogy(mse_rls);
title('RLS'); xlabel('amostras'); ylabel('MSE');
