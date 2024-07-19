function [t,u] = formula(f,tvet,y0,n)

t0 = tvet(1);
tf = tvet(2);

h = (tf-t0)/n;    %passo di discretizzazione
t = [t0; zeros(n,1)];

u = [y0; zeros(n,1)];
for i = 1:n
    t(i+1) = t(i) + h;
    K1 = f(t(i), u(i));
    K2 = f(t(i) + 1/3*h, u(i) + 1/3*h*K1);
    K3 = f(t(i) + 2/3*h, u(i) - 1/3*h*K1 + h*K2);
    K4 = f(t(i) + h, u(i) + h*K1 - h*K2 + h*K3);
    
    u(i+1) = u(i) + 1/8*h*(K1 + 3*K2 + 3*K3 + K4);
    
end



return
    