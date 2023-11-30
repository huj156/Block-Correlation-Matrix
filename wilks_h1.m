function T=wilksal(sample_m,vec_p,A)
if isrow(vec_p)
    vec_p=vec_p';
end

[p,n]=size(sample_m);
k=length(vec_p);
if max(abs(sum(sample_m,2)))<0.00001
    cnt_sample_m=sample_m;
    fd=n;
else
cnt_sample_m=sample_m-mean(sample_m,2)*ones(1,n);
   fd=n-1;
end
r=size(A,2);
shift=log(det(eye(p)+A*A'));

for j=1:k
         indj=sum(vec_p(1:j));
         Aj=A(indj-vec_p(j)+1:indj,1:r);
         Ij=eye(vec_p(j));
        shift=shift-log(det(Ij+Aj*Aj'));         
end

    
meanjy=-(-log(1-p/fd))*(p-fd+1/2)+sum((-log(1-vec_p/fd)).*(vec_p-fd+1/2))+shift;
varjy=2*(-log(1-p/fd))-2*sum(-log(1-vec_p/fd));

X=mat2cell(cnt_sample_m,vec_p,n);
if p>n-2
    error('Not applicable for Wilks statistic');
else
    T_jy=log(det(1/n*(cnt_sample_m*cnt_sample_m')));
    for i=1:k
        T_jy=T_jy-log(det(1/n*(X{i}*X{i}')));
    end
    T=-(T_jy-meanjy)/sqrt(varjy);
end
end
