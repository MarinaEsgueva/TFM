function f=aprox(rbf,ctr,coef,epsilon,x)
N=size(ctr,2);
for i=1:N
    rbf_eval(i)=feval(rbf,norm(x-ctr(:,i)),epsilon);
end
f=rbf_eval*coef; 
