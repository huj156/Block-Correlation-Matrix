clear;
% s=rng(0);
% save('s','s');
load('s.mat');
rng(s);
N=1000;
p=32;
%p=32*5;
rho=[0:0.5:4];
size_rho=length(rho);
for ii=1:size_rho
 if rho(ii)==1
    rho(ii)=0.9;
 end
end
alpha05=norminv(0.95,0,1);

D=3;
    M=1;
    for G=1:3
        for S=1:3
            if S==1
                vec_p=ones(4,1)*p/4;
            elseif S==2
                vec_p=ones(p/2,1)*2;
            elseif S==3
                vec_p=[2,2,p/2-2,p/2-2]';
            end
            k=length(vec_p);
            if G==1
                n=p*2;
            elseif G==2
                n=p+3;
            elseif G==3
                n=max(vec_p)*3;
            end

            powerplot=figure('color',[1 1 1]);
            for jj=1:3
                power_schott=zeros(1,size_rho);
                power_wilks=zeros(1,size_rho);
                power_jbz=zeros(1,size_rho);
                power_yhn=zeros(1,size_rho);
                for ii=1:size_rho
                    Sigma=eye(p);
                    res_schott=zeros(1,N);
                    res_wilks=zeros(1,N);
                    res_jbz=zeros(1,N);
                    res_yhn=zeros(1,N);
                for j=1:N
                        if D==1
                            X=randn(p,n);
                            %Z=randn(r,n);
                        elseif D==2
                            X=(chi2rnd(1,p,n)-1)/sqrt(2);
                            %Z=(chi2rnd(1,r,n)-1)/sqrt(2);
                        elseif D==3
                            X=trnd(5,p,n)/sqrt(5/3);
                            %Z=trnd(5,r,n)/sqrt(5/3);
                        end
                    if jj==1
                        if p>n
                           r=p;
                        else
                           r=n;
                        end
                        if D==1
                            Z=randn(r,n);
                        elseif D==2
                            Z=(chi2rnd(1,r,n)-1)/sqrt(2);
                        elseif D==3
                            Z=trnd(5,r,n)/sqrt(5/3);
                        end
                       for ti=1:k
                           At=zeros(vec_p(ti),r);
                           for tj=1:vec_p(ti)
                               At(tj,tj)=rho(ii)*p^(-1/2)/5;
                           end
                           Samplepart2=At*Z;
                           if ti==1
                              Sample2=Samplepart2;
                           else
                              Sample2=[Sample2;Samplepart2];
                           end
                       end
                       sample_m=sqrtm(Sigma)*X+Sample2;
                    elseif jj==2
                         r=length(vec_p);
                         Z=trnd(5,r,n)/sqrt(5/3);
                        for ti=1:k
                            At=zeros(vec_p(ti),r);
                            At(1,1)=rho(ii);
                            Samplepart2=At*Z;
                           if ti==1
                              Sample2=Samplepart2;
                           else
                              Sample2=[Sample2;Samplepart2];
                           end
                        end
                        sample_m=sqrtm(Sigma)*X+Sample2;
                    elseif jj==3
                         r=length(vec_p);
                        if D==1
                            Z=randn(r,n);
                        elseif D==2
                            Z=(chi2rnd(1,r,n)-1)/sqrt(2);
                        elseif D==3
                            Z=trnd(5,r,n)/sqrt(5/3);
                        end
                         for ti=1:k
                            if ti<3
                            At=zeros(vec_p(ti),r);
                            At(1,1)=rho(ii);
                            else
                                At=zeros(vec_p(ti),r);
                            end
                            Samplepart2=At*Z;
                           if ti==1
                              Sample2=Samplepart2;
                           else
                              Sample2=[Sample2;Samplepart2];
                           end
                        end
                        sample_m=sqrtm(Sigma)*X+Sample2;     
                      end


                        res_schott(j)=schott(sample_m,vec_p);
                        if sum(vec_p(1:k-1))<n-1
                            res_jbz(j)=jbz(sample_m,vec_p);
                        end
                        res_yhn(j)=yhn(sample_m,vec_p);
                        if p<n-1
                            res_wilks(j)=wilks(sample_m,vec_p);
                        end
                end
                    power_schott(ii)=sum(real(res_schott)>alpha05)/N;
                    power_jbz(ii)=sum(real(res_jbz)>alpha05)/N;
                    power_yhn(ii)=sum(real(res_yhn)>alpha05)/N;
                    power_wilks(ii)=sum(real(res_wilks)>alpha05)/N;
                end



                subplot(1,3,jj)
                if jj==1
                    %rho_plot=rho*0.1; %p=32
                    rho_plot=rho/5/p^(-1/2); %p=160
                else
                    rho_plot=rho;
                end
                plot(rho_plot,power_schott,'-ko', 'LineWidth',2,'MarkerSize',10,'DisplayName','Schott')
                hold on
                if all(power_jbz==0)==0
                    plot( rho_plot,power_jbz,'-.r^', 'LineWidth',2,'MarkerSize',10,'DisplayName','JBZ')
                end
                plot( rho_plot,power_yhn, ':bs','LineWidth',2,'MarkerSize',10,'DisplayName','YHN')
                if all(power_wilks==0)==0
                    plot(rho_plot,power_wilks,'--gh','LineWidth',2,'MarkerSize',10,'DisplayName','Wilks')
                end
                subplotname=sprintf('H%d',jj);
                title(subplotname);
                legend('Schott','JBZ','YHN','Wilks');
                set(gca,'FontSize',20);
                set(gca,'Position', [0.05+0.33*(jj-1),0.27,0.27,0.45])
                grid on

                xlabel('\rho','FontSize',20);
                ylabel('power','FontSize',20);
            end
            ztitle=sgtitle(['Empirical power of D',num2str(D), ...
                'G',num2str(G),'S',num2str(S),'p',num2str(p)]);
               ztitle.FontSize = 40;
            saveas(powerplot,['Empirical power of D',num2str(D), ...
                'G',num2str(G),'S',num2str(S),'p',num2str(p)]);
            close(powerplot);
            end
    end