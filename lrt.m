function T=lrt(sample_m,vec_p)
if isrow(vec_p)
    vec_p=vec_p';
end

[p,n]=size(sample_m);
k=length(vec_p);
meanjy=-(-log(1-p/(n-1)))*(p-n+1+1/2)+sum((-log(1-vec_p/(n-1))).*(vec_p-n+1+1/2));
varjy=2*(-log(1-p/(n-1)))-2*sum(-log(1-vec_p/(n-1)));
cnt_sample_m=sample_m-mean(sample_m,2)*ones(1,n);
X=mat2cell(cnt_sample_m,vec_p,n);
if p>n-2
    error('Not applicable for LRT statistic');
else
    T_jy=log(det(1/n*(cnt_sample_m*cnt_sample_m')));
    for i=1:k
        T_jy=T_jy-log(det(1/n*(X{i}*X{i}')));
    end
    T=-(T_jy-meanjy)/sqrt(varjy);
end
end
