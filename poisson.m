function Y = poisson(lambda)

P0 = exp(-lambda);
Pp(1) = P0;
ipp_up = round(10*lambda); 

    if ipp_up<10;
    ipp_up=10;
    end

    for i=2:ipp_up;
        Pp(i) = lambda*Pp(i-1)/(i-1);
    end
   
    QQ=0;
    x=rand;
    kout = ipp_up;
    
    for m=1:ipp_up;
    QQ = QQ + Pp(m);
        if x<QQ;
        kout = m-1;
        break
        end
    end
    Y = kout;
    
