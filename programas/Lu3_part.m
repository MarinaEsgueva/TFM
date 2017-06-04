function Lu=Lu3_part(x,y)
Lu=(-2*(y*(1-y)+x*(1-x))+(2+(1-x)*y*(1-y)-x*y*(1-y))^2+(2+(1-y)*x*(1-x)-x*y*(1-x))^2)*exp(2*x+2*y+x*y*(1-x)*(1-y));
