function Lphi=Limq(r,epsilon)
%Operador de Laplace de la IMQ
Lphi=epsilon^2*((epsilon*r).^2-2)./(1+(epsilon*r).^2).^(5/2);