% performs euler method on function f
% on range [xinit, xend]
% y(xinit) = yinit
% h - step
function [x,y] = EulerMethod(f,xinit,xend,yinit,h)
    N = floor((xend-xinit)/h);

    x = [xinit zeros(1, N)];
    y = [yinit zeros(1, N)];

    for i=1:N    
        x(i+1) = x(i)+h;    
        y(i+1) = y(i) + h*feval(f,x(i),y(i));
    end

end