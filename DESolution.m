function [ a ] = DESolution( x, c )
    a = -(3*exp(3*c+3/2*x^2))/(exp(3*c+3/2*x^2)-1);
end