function [t,u] = Runge_Kutta4(f,tvet,A,n)

t0 = tvet(1);
tf = tvet(2);

h = (tf-t0)/n;    %passo di discretizz
t = [t0; zeros(n,1)];

u = [A;zeros(n,1)];
for i = 1:n
    t(i+1) = t(i)+h;
    K1 = f(t(i),u(i)); %passo EE
    K2 = f(t(i)+0.5*h, u(i)+1/2*h*K1);
    K3 = f(t(i)+0.5*h, u(i)+1/2*h*K2);
    K4 = f(t(i)+h, u(i)+h*K3);
    u(i+1) = u(i)+1/6*(K1+2*K2+2*K3+K4)*h;
    
end



return
    