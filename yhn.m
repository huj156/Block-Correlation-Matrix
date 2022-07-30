function T=yhn(sample_m,vec_p)
if isrow(vec_p)
    vec_p=vec_p';
end
[p,n]=size(sample_m);
k=length(vec_p);
cnt_sample_m=sample_m-mean(sample_m,2)*ones(1,n);
X=mat2cell(cnt_sample_m,vec_p,n);

S=1/(n-1)*(cnt_sample_m)*cnt_sample_m';
K=1/(n-1)*trace((cnt_sample_m'*cnt_sample_m).*(cnt_sample_m'*cnt_sample_m));
SigmaF=(n-1)/(n*(n-2)*(n-3))*((n-1)*(n-2)*trace(S^2)+(trace(S))^2-n*K);

vec_Sigmad=zeros(k,1);
for g=1:k
    Sg=1/(n-1)*X{g}*X{g}';
    Kg=1/(n-1)*trace((X{g}'*X{g}).*(X{g}'*X{g}));
    vec_Sigmad(g)=(n-1)/(n*(n-2)*(n-3))*((n-1)*(n-2)*trace(Sg^2)+(trace(Sg))^2-n*Kg);
end
v=4/(n*(n-1))*(sum(sum(vec_Sigmad*vec_Sigmad'))-trace(vec_Sigmad*vec_Sigmad'));
T=(SigmaF-sum(vec_Sigmad))/sqrt(v);
end

