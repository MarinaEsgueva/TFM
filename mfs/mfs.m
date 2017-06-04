clear all
%Método de las soluciones fundamentales. 
%Problemas homogeneos.
r=100;
%Función de las condiciones de frontera
frontera='u3';
particular='Lu3_part';
sol='u3';
%Fijamos el número de puntos de colocación
N=50; 
h=4/N;
v=linspace(0,4-h,N);

for i=1:N
    if v(i)<1
        x(:,i)=[v(i);0];
    elseif 1<=v(i)&&v(i)<2
        x(:,i)=[1;v(i)-1];
    elseif 2<=v(i)&&v(i)<3
        x(:,i)=[3-v(i);1];
    else
        x(:,i)=[0;4-v(i)];
    end
end

centro=[.5;.5]; 
%radio de los centros2

alpha=2*pi/N;
%centros
for i=0:(N-1)
    y(1,i+1)=centro(1)+r*sin(i*alpha);
    y(2,i+1)=centro(2)+r*cos(i*alpha);
end

%matriz de colocación

for i=1:N
    for j=1:N
        M(i,j)=log(norm(x(:,i)-y(:,j))^2);
    end
end

%vector de términos independientes

for i=1:N
    f(i,1)=feval(frontera,x(1,i),x(2,i))-feval(particular,x(1,i),x(2,i));
end
coef=M\f;

%Representacion de la solucion
nval=40;
v=linspace(0,1,nval);
[X,Y]=meshgrid(v,v);
test=[X(:)';Y(:)'];

for i=1:size(test,2)
    ap(i)=aprox(coef,y,test(:,i))+feval(particular,test(1,i),test(2,i));
    err(i)=ap(i)-feval(sol,test(1,i),test(2,i));
end

Z=reshape(ap,nval,nval);
ZZ=reshape(err,nval,nval);
figure(1)
surf(X,Y,Z)
figure(2)
surf(X,Y,ZZ)
ECM=norm(err)/nval
cond(M)



