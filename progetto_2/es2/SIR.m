function F=SIR(t,y)
    global mu beta 
    S=y(1);
    I=y(2);
    R=y(3);
    
    F=zeros(3,1);
    F(1)=-beta*S*I;
    F(2)=beta*S*I-mu*I;
    F(3)=mu*I;
    
return