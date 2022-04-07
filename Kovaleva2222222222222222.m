%ФАЗОВАЯ МОДУЛЯЦИЯ(ФМ) ВАРИАНТ 1
clc 
clear all
nfig = 1;
close all 
figure('position',[100,100,700,500]);
f0 = 1200;
Vmod = 600;
Vinf = 600;
T0 = 1 / f0;
T = 1 / Vmod;
q = 2 ^ round(Vinf * T);
A = 1;
L0 = 50;
dt = T0 / L0;   
t = 0 : dt : T;
si = zeros(q, length(t));
E = zeros(q, 1);

for i = 0:q-1
    si(i + 1, :) = A*cos(2*pi*f0*t - 2*pi*i/q);
    figure(nfig);
    subplot(2,1,i+1);
    plot(t, si(i + 1,:),'LineWidth',2)
    xlabel('t');
    ylabel('si(t)');
    title(['Сигнал i=', num2str(i)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    axis([0 T -1.5 1.5]);
    grid on
    
end
for i = 0:q-1
    E(i+1) = sum(si(i + 1, :).^2) *dt;
    fprintf('E%d = %d\n', i+1, E(i+1)); 
end

% Задание номер 2. Спектры
nfig = 2;
dF = 1/(10*T);
F = 0:dF:6500;
S = zeros(q, length(F));
W=1/T;
% fprintf('W%d  ширина спектра = %d\n', W); 
for i = 0:q-1
    S(i+1,:) = A*(T/2)*cos(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(-1j*pi*F*T)...
        +(1/1j)*A*(T/2)*sin(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(-1j*pi*F*T);
    figure(nfig); 
    subplot(2,1,i+1);
    plot(F, abs(S(i + 1, :)),'m','LineWidth', 2);
    hold on;
    grid on;
    xlabel('f');
    ylabel('S(f)');
    title(['Спектр сигнала Si(t) при i = ', num2str(i) + 1]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 

end
    %legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    nfig = 3;
    N=600;
l=round((1)*rand(1,N));
ss=zeros(1,length(F));
%F = 0:dF:6500;
for k=1:N
    ss=ss+S(l(k)+1,:).*exp(-1j*2*pi*(k-1)*F*T);
end
ss1=ss/N;
%figure(nfig);
% subplot(2,1,2)
   % plot(F(length(F)/2:length(F)),abs(ss(length(F)/2:length(F))));
   plot(F,abs(ss));
   title('Спектр последовательностей');


