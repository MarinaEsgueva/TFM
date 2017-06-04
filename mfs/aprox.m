function f=aprox(coef,centros,x)
N=size(centros,2);
for i=1:N
    s(i)=log(norm(x-centros(:,i))^2);
end
f=s*coef;
