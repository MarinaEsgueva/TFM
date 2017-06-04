function [ECM,N,ctr,xi,xb,dx,coef]=centros2(simul_Lu,simul_frontera,rbf,Lrbf,eps1,tol,sol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elección de la distribución de centros de las RBF.
%Comienza con N0 nodos distribuidos en una rejilla y va añadiendo centros
%en puntos donde el error sea mayor que una tolerancia dada. 
%Válida en dominios [0,1]x[0,1]
%EN ENTRADA: 
%   simul_Lu: simulador de la función f(x,y) de la EDP con Lu(x,y)=f(x,y)
%   simul_frontera simulador de la función g(x,y) de las condiciones de
%   frontera: u(x,y)=g(x,y)
%   rbf: base de  rbf elegida
%   Lrbf: operador diferencial aplicado a la RBF.
%   epsilon: parámetro de forma elegido para las RBF
%   Nmax: número máximo de centros 
%   tol1: tolerancia para el test de parada sobre el error comentido en la
%   aproximación. 
%   tol2: tolerancia para añadir un check-point al conjunto de centros
%   xb: puntos de colocación en la frontera.


% EN SALIDA:
%   error: estimación del error cometido en la resolución de la edp
%   N: número de centros utilizado en la aproximación
%   ctr: conjunto de centros
%   x: conjunto de puntos de colocación
%   dx: menor separación entre centros. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Puntos de la frontera. Definimos cada parte del borde. 
v=linspace(0,1,10);
b1=[v(1:10);zeros(1,10)]; %y=0
b2=[v(1:10);ones(1,10)]; %y=1
b3=[zeros(1,10);v(1:10)]; %x=0
b4=[ones(1,10);v(1:10)]; % x=1

xb=unique([b1 b2 b3 b4]','rows')'; %Eliminar esquinas. 


%Rejilla inicial de puntos interiores.  
v=v(2:9);
[X,Y]=meshgrid(v,v);
xi=[X(:)';Y(:)'];

ctr=[xb xi]; 
ni=size(xi,2);
nb=size(xb,2);
N=ni+nb;
ind=true; %Inicializamos para entrar en el bucle while. 
k=1;
while any (ind)
    %Desplazamos las puntos iniciales interiores
    coord_x=unique(xi(1,:));
    coord_y=unique(xi(2,:));
    dx=diff([0 coord_x 1]);
    dy=diff([0 coord_y 1]);
    d=min([dx dy]);
    yi=xi+0.3*d*ones(2,ni);
   
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
%Matriz de colocación. 
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

    %Rejilla de puntos en los que se evalua la diferencia entre soluciones 
    zi=[xi,[zeros(1,size(b3,2)-2);b3(2,2:end-1)],[b1(1,1:end-1);zeros(1,size(b1,2)-1)]]+0.5*d*ones(2,ni+size(b3,2)+size(b1,2)-3);
    
    z1=[b1(1,1:end-1)+0.5*diff(b1(1,:));zeros(1,size(b1,2)-1)];
    z2=[b2(1,2:end)-0.5*diff(b2(1,:));ones(1,size(b2,2)-1)];
    z3=[zeros(1,size(b3,2)-1);b3(2,2:end)-0.5*diff(b3(2,:))];
    z4=[ones(1,size(b4,2)-1);b4(2,1:end-1)+0.5*diff(b4(2,:))];
    figure(k)
    plot(zi(1,:),zi(2,:),'*')
    %Evaluamos la diferencia entre soluciones  en la rejilla z. 
    ind_int=zeros(1,ni);
    for i=1:ni
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(zi(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(zi(:,i)-ctr2(:,j)),eps1);
        end
         
        if abs(rbf_eval1*coef-rbf_eval2*coef2)>tol
            ind_int(i)=1; %Añadimos el punto al conjunto de puntos de colocación. 
        end
    end
    %Actualización del conjunto de puntos del interior. 
    xi=[xi zi(:,ind_int==1)];
    ni=size(xi,2);
    
    
    ind_b1=zeros(1,size(z1,2));
    for i=1:length(z1);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z1(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z1(:,i)-ctr2(:,j)),eps1);
        end
        if abs(rbf_eval1*coef-rbf_eval2*coef2)>tol
            ind_b1(i)=1;
        end
    end
    b1=[b1 z1(:,ind_b1==1)];
    
    ind_b2=zeros(1,size(z2,2));
    for i=1:length(z2);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z2(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z2(:,i)-ctr2(:,j)),eps1);
        end 
        if abs(rbf_eval1*coef-rbf_eval2*coef2)>tol
            ind_b2(i)=1;
        end
    end
    b2=[b2 z2(:,ind_b2==1)];
    
    ind_b3=zeros(1,size(z3,2));
    for i=1:length(z3);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z3(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z3(:,i)-ctr2(:,j)),eps1);
        end 
        if abs(rbf_eval1*coef-rbf_eval2*coef2)>tol
            ind_b3(i)=1;
        end
    end
    
    b3=[b3 z3(:,ind_b3==1)];
    
    ind_b4=zeros(1,size(z4,2));
    for i=1:length(z4);
        for j=1:N
            rbf_eval1(j)=feval(rbf,norm(z4(:,i)-ctr(:,j)),eps1);
            rbf_eval2(j)=feval(rbf,norm(z4(:,i)-ctr2(:,j)),eps1);
        end 
        if abs(rbf_eval1*coef-rbf_eval2*coef2)>tol
            ind_b4(i)=1;
        end
    end
    b4=[b4 z4(:,ind_b4==1)];
    %Actualizamos puntos de la frontera
    xb=unique([b1 b2 b3 b4]','rows')';
    nb=size(xb,2);
    %Actualizamos centros.
    ctr=[xi xb];
    N=ni+nb;
    
    ind=[ind_int ind_b1 ind_b2 ind_b3 ind_b4]; 
    dx=0.5*d; %Minima distancia entre centros (interior)
    k=k+1;
end

for j=1:N
        for i=1:ni
            di(i,j)=norm(xi(:,i)-ctr(:,j));
            
        end
        for i=1:nb
            db(i,j)=norm(xb(:,i)-ctr(:,j));    
        end
end
%Matriz de colocación. 
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
%Rejilla fija para calcular el ECM. 
neval=40;
v=linspace(0,1,neval);
[X,Y]=meshgrid(v,v);
ev_points=[X(:)'; Y(:)'];
for i=1:size(ev_points,2)
    for j=1:N
        rbf_eval(j)=feval(rbf,norm(ev_points(:,i)-ctr(:,j)),eps1);
    end
    aprox(i)=rbf_eval*coef;
    error(i)=abs(feval(sol,ev_points(1,i),ev_points(2,i))-rbf_eval*coef);
end
Z=reshape(error,40,40);
ZZ=reshape(aprox,40,40);
figure(1)
surf(X,Y,Z);
figure(2)
surf(X,Y,ZZ);
ECM=1/neval*norm(error);
        
  