clearvars, close all
% implementazione modello SIR per la simulazione della diffusione di epidemie

global mu beta

N=1;
I0=0.35e-4;
S0=1-I0;
R0=0;
lt=0;
b0=0.35;
beta=(1-lt)*b0;
mu=0.14;
h=0.5  %guess iniziale per h
tol=1e-8;


%%%SENZA LOCKDOWN

y0=[S0;I0;R0]; % condizione iniziale 

[t, u]=Runge_Kutta_Fehlberg(@SIR,[0 200],y0,h,tol);

S=u(:,1);
I=u(:,2);
R=u(:,3);

figure(1),plot(t,S,'b',t,I,'r',t,R,'k')
legend('S','I','R')
title('senza lockdown')
set(gca, 'Fontsize', 14)


%%LOCKDOWN BREVE E STRETTO
y0=[S0;I0;R0]; % condizione iniziale 

[t1, u]=Runge_Kutta_Fehlberg(@SIR,[0 40],y0,h,tol);

S1=u(:,1);
I1=u(:,2);
R1=u(:,3);

lt=0.7;
beta=(1-lt)*b0

y0=[S1(end);I1(end);R1(end)];

[t2, u]=Runge_Kutta_Fehlberg(@SIR,[40 70],y0,h,tol);

S2=u(:,1);
I2=u(:,2);
R2=u(:,3);

lt=0;
beta=(1-lt)*b0

y0=[S2(end);I2(end);R2(end)];

[t3,u]=Runge_Kutta_Fehlberg(@SIR,[70 200],y0,h,tol);

S3=u(:,1);
I3=u(:,2);
R3=u(:,3);


figure(2),plot(t1,S1,'b',t1,I1,'r',t1,R1,'k',t2,S2,'b',t2,I2,'r',t2,R2,'k',t3,S3,'b',t3,I3,'r',t3,R3,'k')
legend('S','I','R')
title('lockdown breve e stretto')
set(gca, 'Fontsize', 14)


%LOCKDOWN BREVE E DEBOLE
y0=[S0;I0;R0]; % condizione iniziale 

[t1, u]=Runge_Kutta_Fehlberg(@SIR,[0 40],y0,h,tol);

S1=u(:,1);
I1=u(:,2);
R1=u(:,3);

lt=0.3;
beta=(1-lt)*b0;


y0=[S1(end);I1(end);R1(end)];

[t2,u]=Runge_Kutta_Fehlberg(@SIR,[40 70],y0,h,tol);

S2=u(:,1);
I2=u(:,2);
R2=u(:,3);

lt=0;
beta=(1-lt)*b0;

y0=[S2(end);I2(end);R2(end)];

[t3, u]=Runge_Kutta_Fehlberg(@SIR,[70 200],y0,h,tol);

S3=u(:,1);
I3=u(:,2);
R3=u(:,3);


figure(3),plot(t1,S1,'b',t1,I1,'r',t1,R1,'k',t2,S2,'b',t2,I2,'r',t2,R2,'k',t3,S3,'b',t3,I3,'r',t3,R3,'k')
legend('S','I','R')
title('lockdown breve e debole')
set(gca, 'Fontsize', 14)