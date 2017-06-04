function [est,N,ctr,x,dx,coef]=centros2(N,simul_Lu,simul_frontera,rbf,Lrbf,eps1,Nmax,tol1,tol2.xb,xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elección de la distribución de centros de las RBF.
%Comienza con N0 nodos distribuidos en una rejilla y va añadiendo centros
%en puntos donde el error sea mayor que una tolerancia dada. 
%Válida en dominios [0,1]x[0,1]
%EN ENTRADA: 
%   N0: número inicial de centros.
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

% EN SALIDA:
%   error: estimación del error cometido en la resolución de la edp
%   N: número de centros utilizado en la aproximación
%   ctr: conjunto de centros
%   x: conjunto de puntos de colocación
%   dx: menor separación entre centros. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Comenzamos generando una rejilla con el número inicial de puntos:

% El número de nodos debe ser un cuadrado:
if rem(sqrt(N),1)~=0 
    disp('El número inicial de nodos debe ser un cuadrado perfecto')
    return 
end
%Puntos para evaluar el error cuadrático medio. 
neval=40;
v=linspace(0,1,neval);
[X,Y]=meshgrid(v,v);
ev_points=[X(:)'; Y(:)'];

ctr=[xb xi];

%Resolvemos la EDP 
ni=size(xi,2); %Número de puntos interiores.
nb=size(xb,2);  %Número de puntos en la frontera. 
%Calculamos las matrices de distancias entre los puntos del interior y los
%centros y entre los puntos de la frontera y los centros. 
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

%Rejilla de checkponits

%Check-points interiores
v1=[0 v]+(dx/2)*ones(1,nv+1);
nv1=nv+1;
v=sort([0 v v1]);
nv=length(v);
[X,Y]=meshgrid(v1,v1);
check_int=[X(:)'; Y(:)'];
nci=size(check_int,2); %Número de puntos check interiores.
%Check-points en la frontera
check_frontera=[[zeros(1,nv1);v1],[ones(1,nv1);v1],[v1;zeros(1,nv1)],[v1;ones(1,nv1)]];
ncb=size(check_frontera,2); %Número de puntos check de la frontera. 
%Estimación del error

for i=1:nci
    for j=1:N
        Lrbf_eval(j)=feval(Lrbf,norm(check_int(:,i)-ctr(:,j)),eps1);
    end
    error(i)=abs(feval(simul_Lu,check_int(1,i),check_int(2,i))-Lrbf_eval*coef);%/max(1,abs(feval(simul_Lu,check_int(1,i),check_int(2,i))));
end

for i=1:ncb
    for j=1:N
        rbf_eval(j)=feval(rbf,norm(check_frontera(:,i)-ctr(:,j)),eps1);
    end
    error(i+nci)=abs(feval(simul_frontera,check_frontera(1,i),check_frontera(2,i))-rbf_eval*coef);%/max(1,abs(feval(simul_frontera,check_frontera(1,i),check_frontera(2,i))));
end
est=1/sqrt(nci+ncb)*norm(error);
N1=0; %Inicializamos valor para que entre en while
while est>tol1 && N<Nmax && N>N1
  ni1=ni; %Guardamos el número de la iteración anterior de puntos interiores
  nb1=nb;
  N1=N; 
  dx=dx/2;
  ind_interior=ones(1,nci);
  ind_frontera=ones(1,ncb);
  for i=1:nci
      if error(i)>tol2
          ctr=[ctr check_int(:,i)];% Actualizamos el conjunto de centros
          N=N+1;
          xi=[xi check_int(:,i)];
          ni=ni+1;
          ind_interior(i)=0; %Este punto se va a eliminar del conjunto de check-points
      end
      
  end
  check_int=check_int(:,ind_interior==1);
  for i=1:ncb
      if error(nci+i)>tol2
          ctr=[ctr check_frontera(:,i)];% Actualizamos el conjunto de centros
          N=N+1;
          xb=[xb check_frontera(:,i)];
          nb=nb+1;
          ind_frontera(i)=0; %Este punto se va a eliminar del conjunto de check-points
      end
      
  end
  check_frontera=check_frontera(:,ind_frontera==1);
  
  %Refinar rejilla
  
  v1=v+(dx/2)*ones(1,nv);
  nv1=nv;
  v=sort([v v1]);
  nv=length(v);
  [X,Y]=meshgrid(v1,v1);
  check_int=[check_int,[X(:)'; Y(:)']];
  nci=size(check_int,2); %Número de puntos check interiores.
  
  %Check-points en la frontera
  check_frontera=[check_frontera, [zeros(1,nv1);v1],[ones(1,nv1);v1],[v1;zeros(1,nv1)],[v1;ones(1,nv1)]];
  ncb=size(check_frontera,2); %Número de puntos check de la frontera.
  
  %Actualización de la matriz de distancias:
  %Para las filas correspondientes a los puntos de colocación anteriores
  %solo calculamos distancias con respecto a nuevos puntos.
  for j=(N1+1):N
    for i=1:ni1
        di(i,j)=norm(xi(:,i)-ctr(:,j));
    end
    for i=1:nb1
        db(i,j)=norm(xb(:,i)-ctr(:,j));
    end
    
  end 
  %Para las filas correspondientes a nuevos puntos de colocación hay que
  %añadir distancias con todos los centros
  
  for j=1:N
      if ni>ni1 %Comprobar si se han añadido puntos interiores
          for i=(ni1+1):ni
            di(i,j)=norm(xi(:,i)-ctr(:,j));
          end
      end
      if nb>nb1 %conporbar si se han añadido puntos de frontera
        for i=(nb1+1):nb
            db(i,j)=norm(xb(:,i)-ctr(:,j));
        end
      end
    
  end 
  %Actualizamos la matriz de colocación
  %Primero actualizamos las filas correspondientes a los puntos anteriores
  for i=1:ni1
    Li(i,(N1+1):N)=feval(Lrbf,di(i,(N1+1):N),eps1);  
  end
  for i=1:nb1
    Lb(i,(N1+1):N)=feval(rbf,db(i,(N1+1):N),eps1);
  end
  %Añadimos las nuevas filas correspondientes a los nuevos puntos
  if ni>ni1
      for i=(ni1+1):ni
        Li(i,1:N)=feval(Lrbf,di(i,1:N),eps1);  
      end
  end
  if nb>nb1
      for i=(nb1+1):nb
        Lb(i,1:N)=feval(rbf,db(i,1:N),eps1);
      end
  end
  
  %Actualización del vector de términos independientes
  if ni>ni1
      for i=(ni1+1):ni
        fi(i,1)=feval(simul_Lu,xi(1,i),xi(2,i));
      end
  end
  if nb>nb1
    for i=(nb1+1):nb
        fb(i,1)=feval(simul_frontera,xb(1,i),xb(2,i));
    end
  end
  
  coef=[Li;Lb]\[fi;fb];

  %Nueva estimación del error. 
  for i=1:nci
    for j=1:N
        Lrbf_eval(j)=feval(Lrbf,norm(check_int(:,i)-ctr(:,j)),eps1);
    end
    error(i)=abs(feval(simul_Lu,check_int(1,i),check_int(2,i))-Lrbf_eval*coef);%/max(1,abs(feval(simul_Lu,check_int(1,i),check_int(2,i))));
  end

   for i=1:ncb
      for j=1:N
          rbf_eval(j)=feval(rbf,norm(check_frontera(:,i)-ctr(:,j)),eps1);
      end
      error(i+nci)=abs(feval(simul_frontera,check_frontera(1,i),check_frontera(2,i))-rbf_eval*coef);%/max(1,abs(feval(simul_frontera,check_frontera(1,i),check_frontera(2,i))));
   end
   est=1/sqrt(nci+ncb)*norm(error);  
end

if N==N1
    disp('El algoritmo no está añadiendo más centros')
end
if N>=Nmax
    disp('Se ha alcanzado el número máximo de centros permitidos')
end
if est<=tol1
    disp('La estimación del error es menor que la tolerancia')
end
x=[xi xb];