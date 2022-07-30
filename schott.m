function T=schott(sample_m,vec_p)
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

meanL2=p+(sum(sum(vec_p*vec_p'))-trace(vec_p*vec_p'))/fd;
varL2=4*(sum(sum((vec_p.*(fd-vec_p))*(vec_p.*(fd-vec_p))'))...
    -trace((vec_p.*(fd-vec_p))*(vec_p.*(fd-vec_p))'))/fd^4;

X=mat2cell(cnt_sample_m,vec_p,n);
if max(vec_p)>n-2
    error('Not applicable for Schott statistic');
else
    T_L2=X{1}'*(X{1}*X{1}')^(-1)*X{1};
    for i=2:k
        T_L2=T_L2+X{i}'*(X{i}*X{i}')^(-1)*X{i};
    end
    L2=trace(T_L2^2);
   T=(L2-meanL2)/sqrt(varL2);
end
end
