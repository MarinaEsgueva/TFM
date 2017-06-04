function Lphi=Lgaussiana(r,epsilon)
%Operador de la Laplace de la gaussiana.
Lphi=4*epsilon^4.*r.^2.*exp(-(epsilon*r).^2)-4*epsilon^2*exp(-(epsilon*r).^2);