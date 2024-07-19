function [t, u] = Runge_Kutta_Fehlberg_non_ad(f, tvet, y0, n)

t0 = tvet(1);
tf = tvet(2);

h = (tf-t0)/n;    %passo di discretizz
t = [t0; zeros(n, 1)];

u = [y0;zeros(n, 1)];

[c, A, b, bhat, be]  =  costanti_RK45();

for i = 1:n
    t(i+1) = t(i)+h;
    K1 = f(t, u(i)); % Eulero esplicito 
    K2 = f(t+c(2)*h, u(i)+A(2, 1)*h*K1);
    K3 = f(t+c(3)*h, u(i)+A(3, 1)*h*K1+A(3, 2)*h*K2);
    K4 = f(t+c(4)*h, u(i)+A(4, 1)*h*K1+A(4, 2)*h*K2+A(4, 3)*h*K3);
    K5 = f(t+c(5)*h, u(i)+A(5, 1)*h*K1+A(5, 2)*h*K2+A(5, 3)*h*K3+...
      A(5, 4)*h*K4);
    K6 = f(t+c(6)*h, u(i)+A(6, 1)*h*K1+A(6, 2)*h*K2+A(6, 3)*h*K3+...
      A(6, 4)*h*K4+A(6, 5)*h*K5);
   
    u(i+1) = u(i)+h*(bhat(1)*K1+bhat(2)*K2+bhat(3)*K3+bhat(4)*K4+bhat(5)*K5+bhat(6)*K6);
    
end



return
    