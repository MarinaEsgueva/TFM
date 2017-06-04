function f=aprox2(rbf,ctr_int,ctr_bd,coef,eps_int,eps_front,x)
%Calcula el valor de la solución aproximada en un punto distinguiendo el
%parámetro de forma en la frontera y en el interior.
%Input:
%    rbf: función rbf elegida
%    ctr_int: centros en el interior
%    ctr_front: centros en la frontera
%    coef: coeficentes del aproximante
%    eps_int: parámetro de forma en el interior
%    eps_front: parámetro de forma en la frontera
%    x: punto en el que se calcula la aproximación
ni=size(ctr_int,2);
nb=size(ctr_bd,2);
for i=1:ni
    rbf_eval(i)=feval(rbf,norm(x-ctr_int(:,i)),eps_int);
end
for i=1:nb
    rbf_eval(ni+i)=feval(rbf,norm(x-ctr_bd(:,i)),eps_front);
end
f=rbf_eval*coef; 