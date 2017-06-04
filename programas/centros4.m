function [ECM,N,ctr,xi,xb,coef,condic,eps1]=centros4(simul_Lu,simul_frontera,rbf,Lrbf,xi,b1,b2,b3,b4,dmin,sol,Nmax)
%Incialización de valores. 
ni=size(xi,2); %número de puntos en el interior.
xb=unique([b1 b2 b3 b4]','rows')'; %Se eliminan los posibles puntos repetidos en las esquinas. 
nb=size(xb,2); %Número de puntos en la frontera.
N=ni+nb; %Número de centros
ctr=[xi xb];

%Parámetros para el cálculo del parámetro de forma óptimo. 
minim=0.1;
maxim=10;
paso=0.2;


%Generamos la rejilla de posibles nuevos centros:
v=[0:dmin:1];
[X,Y]=meshgrid(v,v);
new_ctr=[X(:)';Y(:)'];
M=size(new_ctr,2);
plot(new_ctr(1,:),new_ctr(2,:),'*')
%Eliminamos de la rejilla los centros
ind=zeros(1,M);
for i=1:N
    for j=1:M
        if norm(new_ctr(:,j)-ctr(:,i))<0.05
            ind(j)=1;
        end
    end
end
 new_ctr(:,ind==1)=[];
 M=size(new_ctr,2);
 
    
while N<Nmax
    
    %Parámetro de forma óptimo. 
    [eps1,~]=eps_optimo(simul_Lu,simul_frontera,xi,xb,ctr,rbf,Lrbf,minim,maxim,paso);
    
    %Desplazamiento de los centros y puntos interiores. 
    coord_x=sort(unique(xi(1,:))); 
    coord_y=sort(unique(xi(2,:)));
    dx=diff([0 coord_x 1]);
    dy=diff([0 coord_y 1]);
    d1=min([dx dy]);
    yi=xi+0.3*d1*ones(2,ni);
    
    %Puntos de la frontera. Hay que actualizar lado a lado 
    bb1=[b1(1,1:end-1)+0.3*diff(b1(1,:));zeros(1,size(b1,2)-1)];
    bb2=[b2(1,2:end)-0.3*diff(b2(1,:));ones(1,size(b2,2)-1)];
    bb3=[zeros(1,size(b3,2)-1);b3(2,2:end)-0.3*diff(b3(2,:))];
    bb4=[ones(1,size(b4,2)-1);b4(2,1:end-1)+0.3*diff(b4(2,:))];
    yb=[bb1 bb2 bb3 bb4];
    ctr2=[yi yb];
    
    %Matriz de distancias para x e y
    for j=1:N
        for i=1:ni
            di(i,j)=norm(xi(:,i)-ctr(:,j));
            di2(i,j)=norm(yi(:,i)-ctr2(:,j));
            
        end
        for i=1:nb
            db(i,j)=norm(xb(:,i)-ctr(:,j));
            db2(i,j)=norm(yb(:,i)-ctr2(:,j));
        end

    end
%Matrices de colocación. 
    for i=1:ni
        Li(i,1:N)=feval(Lrbf,di(i,1:N),eps1); 
        Mi(i,1:N)=feval(Lrbf,di2(i,1:N),eps1);
        
    end
    for i=1:nb
        Lb(i,1:N)=feval(rbf,db(i,1:N),eps1); 
        Mb(i,1:N)=feval(rbf,db2(i,1:N),eps1);
    end
    %Vector de términos independientes. 
    for i=1:ni
        fi(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
        gi(i,1)=feval(simul_Lu,yi(1,i),yi(2,i));
    end
    for i=1:nb
        fb(i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
        gb(i,1)=feval(simul_frontera,yb(1,i),yb(2,i));
    end
    %Calculamos los coeficientes de la solución. 
    coef=[Li;Lb]\[fi;fb];
    coef2=[Mi;Mb]\[gi;gb];
    condic=cond([Li;Lb]);
    
    %Diferencia entre aproximaciones en los posibles nuevos centros
    aprox1=[];
    aprox2=[];
    for i=1:M
        aprox1(i)=feval('aprox',rbf,ctr,coef,eps1,new_ctr(:,i));
        aprox2(i)=feval('aprox',rbf,ctr2,coef2,eps1,new_ctr(:,i));
    end
    error=abs(aprox1-aprox2)./(1+aprox1);
    [m,n]=max(error);
    nuevo=new_ctr(:,n);
    new_ctr(:,n)=[];
    M=M-1;
     
    % Comprobamos a que conjunto de centros iniciales pertenece
    if abs(nuevo(2))<1.e-8
        b1=[b1 nuevo];
        b1=[sort(b1(1,:));zeros(1,size(b1,2))];
         nb=nb+1;
        xb=[xb nuevo];
    elseif abs(nuevo(2)-1)<1.e-8
        b2=[b2 nuevo];
        b2=[sort(b2(1,:));ones(1,size(b2,2))];
         nb=nb+1;
        xb=[xb nuevo];
    elseif abs(nuevo(1))<1.e-8
        b3=[b3 nuevo];
        b3=[zeros(1,size(b3,2));sort(b3(2,:))];
         nb=nb+1;
        xb=[xb nuevo];
    elseif abs(nuevo(1)-1)<1.e-8
        b4=[b4 nuevo];
        b4=[ones(1,size(b4,2));sort(b4(2,:))];
        nb=nb+1;
        xb=[xb nuevo];
    else
        xi=[xi nuevo];
        ni=ni+1;
    end
    N=N+1;
    ctr=[ctr nuevo];
    
    
end

%Rejilla fija para calcular el ECM. 
neval=40;
v=linspace(0,1,neval);
[X,Y]=meshgrid(v,v);
ev_points=[X(:)'; Y(:)'];

%Matriz de distancias para x e y
for j=1:N
    for i=1:ni
        di(i,j)=norm(xi(:,i)-ctr(:,j));
    end
    for i=1:nb
        db(i,j)=norm(xb(:,i)-ctr(:,j));
    end

end

%Calcular la matriz de colocación:
for i=1:ni
    Li(i,1:N)=feval(Lrbf,di(i,1:N),eps1); 
end
for i=1:nb
    Lb(i,1:N)=feval(rbf,db(i,1:N),eps1); 
end
    %Vector de términos independientes. 
for i=1:ni
    fi(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
end
for i=1:nb
    fb(i,1)=feval(simul_frontera,xb(1,i),xb(2,i));    
end
%Calculamos los coeficientes de la solución. 
coef=[Li;Lb]\[fi;fb];
condic=cond([Li;Lb]);


for i=1:size(ev_points,2)
    aproxim(i)=aprox(rbf,ctr,coef,eps1,ev_points(:,i));
    fun(i)=feval(sol,ev_points(1,i),ev_points(2,i));
    error(i)=abs(aproxim(i)-fun(i))/max([1,abs(fun(i))]);
    diferencias(i)=abs(aproxim(i)-aprox(rbf,ctr2,coef2,eps1,ev_points(:,i)));
end
Z=reshape(error,40,40);
ZZ=reshape(aproxim,40,40);
ZZZ=reshape(diferencias,40,40);
figure(1)
surf(X,Y,Z);
title('Error')
figure(2)
surf(X,Y,ZZ);
figure(3)
plot(ctr(1,:),ctr(2,:),'*')
figure(4)
surf(X,Y,ZZZ);
title('Diferencias')
ECM=1/neval*norm(error);

    
    
    
        
        
    
    
    
    
    