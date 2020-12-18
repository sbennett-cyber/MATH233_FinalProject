function x = biconjgrad(A, b, x)
%%Algorithm from Wikipedia
    r = b - A * x;
    rHat=r;
    rhoOld=1;
    alpha=1;
    omega=1;
    v=x*0;
    p=x*0;
   
    for i = 1:length(b)
        rhoNew=rHat'*r;
        beta=(rhoNew/rhoOld)*(alpha/omega);
        p= r+beta*(p-omega*v);
        v=A*p;
        alpha=rhoNew/(rHat'*v);
        h=x+alpha*p;
        if sqrt(r'*r) < 1e-10
            x=h;
            break
        end
        s=r-alpha*v;
        t=A*s;
        omega=(t'*s)/(t'*t);
        x=h+omega*s;
        if sqrt(s'*s) < 1e-10
            break
        end
        rhoOld=rhoNew;
        r=s-omega*t;  
    end
end