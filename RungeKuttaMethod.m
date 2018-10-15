% performs runge kutta method on function f
% on range [xinit, xend]
% y(xinit) = yinit
% h - step
function [ x, y ] = RungeKuttaMethod(f,xinit,xend,yinit,h)
    N = floor((xend-xinit)/h);

    x = xinit:h:xend; 
    y = [yinit zeros(1, N)];
    
    for i=1:(length(x)-1) 
        k1 = feval(f, x(i),y(i)); 
        k2 = feval(f,x(i)+0.5*h,y(i)+0.5*h*k1); 
        k3 = feval(f,(x(i)+0.5*h),(y(i)+0.5*h*k2)); 
        k4 = feval(f,(x(i)+h),(y(i)+k3*h)); 

        y(i+1) = y(i) + (1/6)*(k1+2*k2+2*k3+k4)*h;
    end
    
end