function T=jbz(sample_m,vec_p)
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

X=mat2cell(cnt_sample_m,vec_p,n);

if sum(vec_p(1:k-1))>n-2
    error('Not applicable for JBZ statistic');
else
    T_jbz=0;
    e=0;
    v=0;
    for i=2:k
        barpi=sum(vec_p(1:i-1));
        tildeYi=cnt_sample_m(sum(vec_p(1:i-1))+1:p,:);
        Piminus1=X{i-1}'*(X{i-1}*X{i-1}'/n)^(-1)*X{i-1}/n;
        Fi=1/barpi*tildeYi* Piminus1*tildeYi'*(1/(fd-barpi)*tildeYi* (eye(n)-Piminus1)*tildeYi')^(-1);
        r1i=vec_p(i)/barpi;
        r2i=vec_p(i)/(fd-barpi);
        alphani=r1i/r2i;
        Mi=Fi*(Fi+alphani*eye(p-barpi))^(-1);
        T_jbz=T_jbz+trace(Mi);
        egi=r2i/(r1i+r2i);
        e=e+vec_p(i)*egi;
        h=sqrt(r1i+r2i-r1i*r2i);
        vgi=2*h^2*r1i^2*r2i^2/(r1i+r2i)^4;
        v=v+vgi;
    end
    T=(T_jbz-e)/sqrt(v);
end
end
