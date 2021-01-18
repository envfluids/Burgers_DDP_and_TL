clc
clear all
close all

L=100.0;
nu=0.02;
A=sqrt(2)*1e-2;

N=1024;
dt=0.01;
s=20; %ratio of LES and DNS time steps

% number of time steps
M=10000000;

% time steps between samples
P=1;

x=[0:N-1]'*L/N;

kx=[0:N/2 -N/2+1:-1]'*2.0*pi/L;

u_old=sin(2*pi*2*x/L+randn*2*pi);

un_old=fft(u_old);
Fn_old=1i*kx.*fft(0.5*(u_old).^2);

U_DNS=zeros(N,M/P);
f_store = zeros(size(U_DNS));
U_DNS(:,1)=u_old;
z=0;

u=u_old;
un=zeros(N,1);

f=zeros(N,1);
for kk=1:3
    C1=randn;
    C2=randn;
    f=f+C1*A/sqrt(kk*s*dt)*cos(2*pi*kk*x/L+2*pi*C2);
end
fn=fft(f);

for m=2:M
    Fn=1i*kx.*fft(0.5*u.^2);

    if(mod(m,s)==0)  
        f= zeros(size(f));

        for kk=1:3
            C1=randn;
            C2=randn;
            f=f+C1*A/sqrt(kk*s*dt)*cos(2*pi*kk*x/L+2*pi*C2);
        end
        fn=fft(f);
    end
    
    C=nu*dt*0.5*(kx).^2; 
    un=((1.0-C).*un_old-0.5*dt*(3.0*Fn-Fn_old)+dt*fn)./(1.0+C);
    
    un_old=un;
    u=real(ifft(un));
    Fn_old=Fn;
    
    if(mod(m,P)==0)
        z=z+1;
        U_DNS(:,z) = u;
        f_store(:,z+1) = f;
    end
end

f_store = f_store(:,1:s:end);

save('DNS_Burgers_s_20.mat','U_DNS','-v7.3')
save('DNS_Force_LES_s_20.mat','f_store','-v7.3')
