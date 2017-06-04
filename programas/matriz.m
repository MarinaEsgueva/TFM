function L=matriz(rbf,Lrbf,xi,xb,ctr,eps1)
N=size(ctr,2);
ni=size(xi,2);
nb=size(xb,2);
%Matriz de distancias.
for i=1:N
    for j=1:ni
        d(j,i)=norm(xi(:,j)-ctr(:,i));
    end
    for j=1:nb
        d(ni+j,i)=norm(xb(:,j)-ctr(:,i));
    end
end
%Matriz de colocación
for i=1:ni
    L(i,1:N)=feval(Lrbf,d(i,1:N),eps1);
end
for i=1:nb
    L(ni+i,1:N)=feval(rbf,d(ni+i,1:N),eps1);
end
