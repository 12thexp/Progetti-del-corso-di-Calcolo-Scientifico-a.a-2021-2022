function [tt,y]=Runge_Kutta_Fehlberg(f,tvet,y0,h,tol)
    
    [c,A,b,bhat,be]=costanti_RK45();
    atol=1e-13;  %tolleranza assoluta
    rtol=tol;  %tolleranza relativa
    alpha=0.8;
    
    k=0;  %contatore delle iterazioni
    t0=tvet(1);
    tf=tvet(2);
    
    %ICs
    t(1)=t0;
    tt(1)=t0;
    y(1,:)=y0;
    wi=y0;
    
    i=1;
    lastit=0;
    while ~lastit  %fino a che lastit=0 rimane dendro al ciclo
        %stretch step se t+1.1*h eccede tf
        if t+1.1*h > tf
            h= tf-t;
            lastit=1;
        end
        
        %calcola uno step per entrambi i metodi RK4, RK5
        K1=f(t,wi);  %passo di Eulero Esplicito
        K2=f(t+c(2)*h, wi+A(2,1)*h*K1);
        K3=f(t+c(3)*h, wi+A(3,1)*h*K1+A(3,2)*h*K2);
        K4=f(t+c(4)*h, wi+A(4,1)*h*K1+A(4,2)*h*K2+A(4,3)*h*K3);
        K5=f(t+c(5)*h, wi+A(5,1)*h*K1+A(5,2)*h*K2+A(5,3)*h*K3+A(5,4)*h*K4);
        K6=f(t+c(6)*h, wi+A(6,1)*h*K1+A(6,2)*h*K2+A(6,3)*h*K3+A(6,4)*h*K4+A(6,5)*h*K5);
        
        %soluzione RKF4  (non serve esplicitamente)
        w=wi+h*(b(1)*K1+b(2)*K2+b(3)*K3+b(4)*K4+b(5)*K5+b(6)*K6);
        
        %soluzione RKF5
        z=wi+h*(bhat(1)*K1+bhat(2)*K2+bhat(3)*K3+bhat(4)*K4+bhat(5)*K5+bhat(6)*K6);
        
        %calcolo stima dell'errore 
        e=h*norm(be(1)*K1+be(2)*K2+be(3)*K3+be(4)*K4+be(5)*K5+be(6)*K6);
        
        
        % gestione passo di integrazione
        tollT=rtol*norm(wi)+atol;  %mix fra tolleranza relativa e assolutata   (garantisce un po di più di robustezza)
        if e<=tollT  % la stima è ok!
            t=t+h;
            h=alpha*h*(tollT/e)^0.2;   %formula empirica
            i=i+1;
            tt(i)=t;
            wi=z;
            y(i, :)=z;
            k=0;
        elseif k==0   %tolleranza non raggiunta la prima volta
            h=alpha*h*(tollT/e)^0.2;
            k=k+1;
            lastit=0;      
        else %tolleranza non raggiunta più di una volta
            h=h/2;
            lastit=0;       
        end    
    
    end

    


return