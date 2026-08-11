% Ecuacion de Burgers, shock - 1D
tt = linspace(0,4);
xx = linspace(0,4);
for i = 1:numel(tt)
    t = tt(i);
    u = nan(numel(xx),1);
    for k = 1:numel(xx)
        u(k) = sol(xx(k),t);
    end
    plot(xx,u), axis([0 4 0 1])
    pause(.5)
end

function u = sol(x,t)
if(t<1)
    if(x<=t)
        u = 1;
    elseif(x<=1)
        u = 1-(x-t)/(1-t);
    else
        u=0;
    end
else
    if(x>(1+t)/2)
        u = 0;
    else
        u = 1;
    end
end
end