clc; close all; clear all;
%Modelovanje i generisanje analognih signala i analiza prenosa ovih signala
%kroz sisteme prenosa koriscenjem frekvencijske modulacije
%GENERISANJE SIGNALA NA ULAZU U SISTEM I SIGNALA ODBIRAKA

Fs=512000;
Nsample=65536;
t=(0:1/Fs:(Nsample-1)/Fs);
Nfft=1024; % menja se sa 1024 2048 4096
fosa2=(0:Fs/Nfft:Fs-Fs/Nfft);
snr = [5,14,26, 15];
f0=Fs/4;
k=2;

%Periodicni signal na ulazu u sistem
x=0.5*cos(2*pi*t*0.5*1000)+0.25*cos(2*pi*t*1*1000)+0.15*cos(2*pi*t*1.5*1000)+0.1*cos(2*pi*t*2*1000);
x_integral= 0.5*sin(2*pi*t*0.5*1000)/(2*pi*0.5*1000)+0.25*sin(2*pi*t*1*1000)/(2*pi*1*1000)+0.15*sin(2*pi*t*1.5*1000)/(2*pi*1.5*1000)+0.1*sin(2*pi*t*2*1000)/(2*pi*2*1000);
figure, plot(t,x); xlabel('t [ms]'); title('vremenski oblik ulaznog signala');
%spektar ulaznog signala
figure, stem( fosa2, abs(fft(x,Nfft))); xlim([0 Fs]);xlabel('f [Hz]');title('frekvencijski spektar ulaznog signala');
x_f = zeros((65536/Nfft), Nfft);
x_f(1, :) = abs(fft(x(1:Nfft),Nfft));
for br =2:(65536/Nfft)
    x_f(br,:) = abs(fft(x((br-1)*Nfft:br*Nfft),Nfft));
end
x_zbir_2 = sum((x_f.^2)./(Fs*Nfft));
x_usrednjeno_2 = x_zbir_2./(65536/Nfft); 
figure, stem( fosa2, x_usrednjeno_2);  xlim([0 Fs]);xlabel('f [Hz]');ylabel('sgss');title('SGSS ulaznog signala');

%MODULACIJA
Ufm= cos(2*pi*f0*t + 2*pi*k*x_integral);
figure, plot(t,Ufm); xlabel('t[ms]'); title('vremenski oblik modulisanog signala');
%spektar modulisanog signala
figure, stem( fosa2, abs(fft(Ufm,Nfft))); xlim([0 Fs]);xlabel('f [Hz]');title('frekvencijski spektar  modulisanog signala');
Ufm_f = zeros((65536/Nfft), Nfft);
Ufm_f(1, :) = abs(fft(Ufm(1:Nfft),Nfft));
for br =2:(65536/Nfft)
    Ufm_f(br,:) = abs(fft(Ufm((br-1)*Nfft:br*Nfft),Nfft));
end
Ufm_zbir_2 = sum((Ufm_f.^2)./(Fs*Nfft));
Ufm_usrednjeno_2 = Ufm_zbir_2./(65536/Nfft); 
figure, stem( fosa2, Ufm_usrednjeno_2);  xlim([0 Fs]);xlabel('f [Hz]');ylabel('sgss');title('SGSS  modulisanog signala');

%dodavanje suma
Ufm_awgn = awgn(Ufm, snr(4)); %menjati 1 2 3 4

%DEMODULACIJA
 idealna_on=true; %true false
 greska_sinhronizacije_faze=pi/4;
h_nf=fir1(20,2200/Fs);
[H_nf,w_nf]=freqz(h_nf,1);
 if(idealna_on)
     y_i = Ufm_awgn*2.*cos(2*pi*f0*t);
     y_q= Ufm_awgn*2.*sin(2*pi*f0*t);
     y_i_nf=filter(h_nf,1,y_i);
     y_q_nf=filter(h_nf,1,y_q);
 else
    y_i = Ufm_awgn*2.*cos(2*pi*f0*t +greska_sinhronizacije_faze);
    y_q= Ufm_awgn*2.*sin(2*pi*f0*t +greska_sinhronizacije_faze);
    y_i_nf=filter(h_nf,1,y_i);
    y_q_nf=filter(h_nf,1,y_q);
 end

 %Ffm=atan(y_q_nf./y_i_nf);
 Ffm=angle(y_q_nf*1i+y_i_nf);
 y_diff=diff(Ffm)/(2*pi*k);
 y=zeros(1,65536);
 y(1:65535)=y_diff(1:65535);
figure, plot((0:1/Fs:(Nsample-1)/Fs),y); xlabel('t [ms]');title('vremenski oblik demodulisanog signala');
%spektar demodulisanog signala
figure, stem( fosa2, abs(fft(y,Nfft))); xlim([0 Fs]);xlabel('f [Hz]');title('frekvencijski spektar demodulisanog signala');
y_f = zeros((65536/Nfft), Nfft);
y_f(1, :) = abs(fft(y(1:Nfft),Nfft));
for br =2:(65536/Nfft)
    y_f(br,:) = abs(fft(y((br-1)*Nfft:br*Nfft),Nfft));
end
y_zbir_2 = sum((y_f.^2)./(Fs*Nfft));
y_usrednjeno_2 = y_zbir_2./(65536/Nfft); 
figure, stem( fosa2, y_usrednjeno_2);  xlim([0 Fs]);xlabel('f [Hz]');ylabel('sgss');title('SGSS demodulisanog signala');
