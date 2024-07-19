function [t,u] = Heun(f,tspan,y0,n)

t0 = tspan(1);
tf = tspan(2);

h = (tf-t0)/n;    %passo di discretizz
t = [t0; zeros(n,1)];

u = [y0; zeros(n,1)];
for i = 1:n
    t(i+1) = t(i)+h;
    K1 = f(t(i), u(i)); %passo EE
    K2 = f(t(i)+h, u(i)+h*K1);
    
    u(i+1) = u(i)+1/2*h*(K1+K2);
    
end



return
    