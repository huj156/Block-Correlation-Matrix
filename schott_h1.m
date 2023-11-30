function T=schottal(sample_m,vec_p,A,r)
if isrow(vec_p)
    vec_p=vec_p';
end
[p,n]=size(sample_m);
k=length(vec_p);

% if max(abs(sum(sample_m,2)))<0.00001
%     cnt_sample_m=sample_m;
%     fd=n;
% else
% cnt_sample_m=sample_m-mean(sample_m,2)*ones(1,n);
%     fd=n-1;
% end
 cnt_sample_m=sample_m;
fd=n;

r=size(A,2);

shift=0;

for i=1:k
    for j=1:k
         indi=sum(vec_p(1:i));
         indj=sum(vec_p(1:j));
         Ai=A(indi-vec_p(i)+1:indi,1:r);
         Aj=A(indj-vec_p(j)+1:indj,1:r);
         Ii=eye(vec_p(i));
         Ij=eye(vec_p(j));
        if i~=j
            shift2=(1-vec_p(i)/fd)*(1-vec_p(j)/fd)*trace(Ai'*(Ii+Ai*Ai')^(-1)*Ai*Aj'*(Ij+Aj*(Aj'))^(-1)*Aj);
        else
            shift2=0;
        end
        shift=shift+shift2;         
    end
end
meanL2=p+(sum(sum(vec_p*vec_p'))-trace(vec_p*vec_p'))/fd+shift;
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