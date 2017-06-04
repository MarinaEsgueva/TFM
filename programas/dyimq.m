function dphi=dyimq(r,dy,epsilon)
%Derivada parcial con respecto a y de la RBF IMQ.
dphi=-epsilon*dy/(1+(epsilon*r)^2)^(5/2);