%������� ���������(��) ������� 1
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
L0 = 32;
dt = T0 / L0; 
df=f0/L0;
t = 0 : dt : T;
s = zeros(q, length(t));


for i = 0:q-1
    Sn = A*cos(2*pi*f0*t - 2*pi*i/q);
    s(i+1,:) = Sn;
end

baz1 = sqrt(2/T)*cos(2*pi*f0*t);
baz1 = baz1/norm(baz1);
    
baz2 = sqrt(2/T)*sin(2*pi*f0*t);
baz2 = baz2/norm(baz2);    


fprintf('������������ �������� ������� f(1) � f(1) = %1.f \n',round(baz1*baz1'));
fprintf('������������ �������� ������� f(1) � f(2) = %1.f \n',-1*round(baz1*baz2'));
fprintf('������������ �������� ������� f(2) � f(2) = %1.f \n',round(baz2*baz2'));

%s1 = zeros(1, length(s(i+1,:)));
%s2 = zeros(1, length(s(i+1,:)));  
s1 = zeros(q, 1);
s2 = zeros(q, 1);

fprintf('����� f(1) = %1.f \n',round(norm(baz1)));
fprintf('����� f(2) = %1.f \n',round(norm(baz2)));

for i = 1:q
    s1(i) = s(i,:)*baz1';
    s2(i) = s(i,:)*baz2';
    
    plot(s1(i),s2(i),'o');
    hold on;
    grid on;
    
end

plot(0,0,'o');
hold on;
grid on;

axis([-6 6 -6 6],'square');



%��������� ��������

SNRdB = 1:10; %����� �������� ��������� ������ ���
Pe=zeros(1,length(SNRdB));
SNRdBtheor = 1:10;
SNRtheor = 10.^(SNRdBtheor/10);
Petheor = qfunc(sqrt(2*SNRtheor));

nErrMax = 5;
for nSNR = 1:length(SNRdB) %���� �� ��������� ��������� ������/��� � �� (y)
    nErr = 0;%��������� �������� �������� ����� ������ 
    nRun = 0;%��������� �������� �������� ����� ���������
    SNR = 10^(SNRdB(nSNR)/10); %Y
    sigma2=(sum(sum(s1.^2)))/(q*2*SNR);
    sigma=sqrt(sigma2);
    
    while nErr < nErrMax %���� ������������� ��� ����� �������� ��������� ������/���
       % �������� ������� i � ��������� 0,1,q-1 ; ��������� ������ �� ������ ������  
        i = floor(q*rand)+1;
        r = s(i,:)+sigma*randn(1,length(t)); %������������� ����������� � ������
        r1 = sum(r.*baz1); %������������� ���������

              for index = 1:q
                  if(r1<0)
                      resi=2 ;
                      if 20 < nRun && nRun < 100
                            plot(r1,'red.','MarkerSize',5); 
                             title('������ �����������');
                      end
                  else
                     resi=1; 
                     if 20 < nRun && nRun < 100
                           plot(r1,'b.','MarkerSize',5); 
                           title('������ �����������');
                      end
                   end
              end
             if resi~=i
                    nErr = nErr + 1; 
             end   
             nRun = nRun + 1;% ���������� �������� ����� ���������
    end
    disp([SNRdB(nSNR),nErr,nErrMax,nRun]);
    Pe(nSNR)=nErr/nRun; % ���������� ����������������� ������ ����������� ������
      %������� ����� ������������ ������ �� ����� ���������
end
figure(2);
axis('square');
semilogy(SNRdBtheor,Petheor,'black',SNRdB,Pe,'red.-','MarkerSize',10);
title('nErrMax=5');
xlabel('SNRdB')
ylabel('Pe')
hold on;
grid on;


