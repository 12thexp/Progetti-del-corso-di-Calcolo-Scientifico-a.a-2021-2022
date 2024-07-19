clearvars, close all

T = 3000;
tspan = [0 T];
y0 = 1200;
ya = 300;
k = 2.20e-12;

% sol esatta
h = @(t,y) -(log(y+ya)+2*atan(y/ya)-log(y-ya)-log(y0+ya)+log(y0-ya)-2*atan(y0/ya))/(4*ya^3) + k*t;
% fimplicit(h, [0 T 300 1200], 'bo-');
% hold on

% sol esatta t = t(y)
g = @(y) 1/k*(log(y+ya)+2*atan(y/ya)-log(y-ya)-log(y0+ya)+log(y0-ya)-2*atan(y0/ya))/(4*ya^3);
% plot(g(yap), yap,'ko-');

y_end = fzero(@(y) g(y)- tspan(end), 400);  %valore di y in t = 3000


nn = [20 40 80 160 320 640 1280];

% Runge-Kutta 4 -----------------------------------------------------------
for i=1:numel(nn)
    n = nn(i);

    [t4,u4] = Runge_Kutta4(f, tspan, y0,n);
    err_rk4(i) = abs(u4(end)-y_end);
    
    %ordine di convergenza:
    if i>1
        p_rk4(i-1) = log(err_rk4(i-1)/err_rk4(i))/log(2);
    end
    % grafici:
%         figure
%         fimplicit(h, [0 T 300 1200], 'bo-');
%         hold on
%         plot(t4, u4, 'r-');
%         title("RK4 |n = " + n);
%         hold off
end


% Heun --------------------------------------------------------------------
for i=1:numel(nn)
    n = nn(i);

    [t4,u4] = Heun(f,tspan, y0, n);
    err_heun(i) = abs(u4(end)-y_end);
    if i>1
        %ordine di convergenza:
        p_heun(i-1) = log(err_heun(i-1)/err_heun(i))/log(2);
    end
    %grafici:
%         figure
%         fimplicit(h, [0 T 300 1200], 'bo-');
%         hold on
%         plot(t4, u4, 'g*-');
%         title("Heun | n = " + n);
%         hold off
end

% Runge-Kutta-Fehlberg non ad -----------------------------------------------------------
for i=1:numel(nn)
    n = nn(i);

    [t4,u4] = Runge_Kutta_Fehlberg_non_ad(f, tspan, y0,n);
    err_rkf45(i) = abs(u4(end)-y_end);
    
    %ordine di convergenza:
    if i>1
        p_rkf45(i-1) = log(err_rkf45(i-1)/err_rkf45(i))/log(2);
    end
    % grafici:
%         figure
%         fimplicit(h, [0 T 300 1200], 'bo-');
%         hold on
%         plot(t4, u4, 'r-');
%         title("RKF45 |n = " + n);
%         hold off
end


% formula -----------------------------------------------------------------
for i=1:numel(nn)
    n = nn(i);

    [t4,u4] = formula(f, tspan, y0,n);
    err_form(i) = abs(u4(end)-y_end);
    
    if i>1
        %ordine di convergenza:
        p_form(i-1) = log(err_form(i-1)/err_form(i))/log(2);
    end  
    % grafici:
%         figure
%         fimplicit(h, [0 T 300 1200], 'bo-');
%         hold on
%         plot(t4, u4, 'r-');
%         title("formula |n = " + n);
%         hold off
    
end


%grafico errore
figure
loglog(nn, err_rk4, 'o-', nn, err_heun, 'o-', nn, err_rkf45, 'o-', nn, err_form, 'o-','LineWidth',2);
xlabel('nr. passi', 'FontSize',12);
ylabel('errore', 'FontSize',12)
leg = legend('RK4', 'Heun', 'RKF45', 'formula', 'FontSize',12);
set(leg,'position',[0 0 0.3 0.3])
title("Andamento dell'errore")
grid on

% p_rk4
% p_heun
% p_rkf45
% p_form
% 

