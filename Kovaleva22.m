%ФАЗОВАЯ МОДУЛЯЦИЯ(ФМ) ВАРИАНТ 1
clc 
clear all
close all 
nfig = 1;
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

i = 0;
    si1(1,:) = A*cos(2*pi*f0*t - 2*pi*i/q);
    ss=si1+Si1(l(k)+1,:)
    figure(nfig);
   % subplot(2,1,i+1);
    plot(t, si(i + 1,:),'LineWidth',2)
    xlabel('t');
    ylabel('si(t)');
    title(['Сигнал i=', num2str(i)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    axis([0 T -1.5 1.5]);
    grid on

for i = 0:q-1
    E(i+1) = sum(si(i + 1, :).^2) *dt;
    fprintf('E%d = %d\n', i+1, E(i+1)); 
end

% Задание номер 2. Спектры
nfig = 2;
dF = 1/(10*T);
F = 0:dF:6500;
S = zeros(q, length(F));

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
    title(['Спектр сигнала Si(t) дискретной ФМ при i = ', num2str(i) + 1]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
end

    
 nfig=4;
ss = zeros(2, length(F));
ss(2, :) = -1*A*(T/2)*((sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*pi*F*T);
ss(1, :) = A*(T/2)*((sinc((F-f0)*T)+sinc((F+f0)*T)));
   for i=0:q-1
   figure(nfig);
   subplot(2,1,i+1)
    plot(F,abs(ss(i+1, :)),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title(['Амплитудный спектр отрезка гармоники сигнала при i=', num2str(i)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
   end

    nfig=6;
I=["0,0,1,1,0,0,0,1,1,0","1,1,1,1,1,1,1,1,1,1"];
ss = zeros(2, length(F));
ss(2, :) = A*(T/2)*(-10*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
ss(1, :) = A*(T/2)*(6*(sinc((F-f0)*T)+sinc((F+f0)*T))-4*(sinc((F-f0)*T)+sinc((F+f0)*T))).*exp(-1j*2*pi*F*T);
   for i=0:q-1
   figure(nfig);
   subplot(2,1,i+1)
    plot(F,abs(ss(i+1, :)),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title(['спектр последовательностей сигнала при N=10, I=',I(i+1)]);
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
   end
   
nfig=7;
N=500;
l=round((q-1)*rand(1,N));
ss=zeros(1,length(F));
for k=1:N
    ss=ss+S(l(k)+1,:).*exp(-1j*2*pi*(k-1)*F*T);
end
figure(nfig);
    plot(F,abs(ss),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title('N=500');
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on

nfig=8;
N=1000;
l=round((q-1)*rand(1,N));
ss=zeros(1,length(F));
for k=1:N
    ss=ss+S(l(k)+1,:).*exp(-1j*2*pi*(k-1)*F*T);
end
figure(nfig);
    plot(F,abs(ss),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title('N=1000');
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
    
        
    nfig=9;
N=90;
l=round((q-1)*rand(1,N));
ss=zeros(1,length(F));
for k=1:N
    ss=ss+S(l(k)+1,:).*exp(-1j*2*pi*(k-1)*F*T);
end
figure(nfig);
    plot(F,abs(ss),'LineWidth',2);
    xlabel('F');
    ylabel('si(F)');
   title('N=9');
    legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
    grid on
    
    nfig=10;
N=8;
figure(nfig);
for i = 0:q-1
    S(i+1,:) = A*(T/2)*cos(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(-1j*pi*F*T)...
        +(1/1j)*A*(T/2)*sin(2*pi*i/q)*(sinc((F - f0)*T)+sinc((F + f0)*T)).*exp(-1j*pi*F*T);
   % figure(nfig); 
   % subplot(2,1,i+1);
  %plot(F(length(F)/2:length(F)),abs(S(i + 1,(length(F)/2:length(F))),'m');
 % plot(F(length(F)/2:length(F)),abs(S(i+1,length(F)/2:length(F))));
    plot(F,abs(S(i+1,:)));
    hold on;
    grid on;
end
 xlabel('f');
    ylabel('S(f)');
   
    %legend('f0 = 1200,Vmod = 600,Vinf = 600'); 
l=round((q-1)*rand(1,N));
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
 
