clear;
% s=rng(0);
% save('s','s');
load('s.mat');
rng(s);
N=1000;
p=32;
%p=32*5;
rho=[0:0.2:0.8];
size_rho=length(rho);
for D=1:3
    for G=1:3
        for S=1:3
            for H=1:3
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

          for ii=1:size_rho      
                res_schott=zeros(1,N);
                res_wilks=zeros(1,N);


                if H==1
                   rho=[0:0.2:0.8]*p^(-1/2);
                    if p>n
                       r=p;
                    else
                       r=n;
                    end

                   for ti=1:k
                       At=zeros(vec_p(ti),r);
                       for tj=1:vec_p(ti)
                       At(tj,tj)=rho(ii);
                       end
                       Apart2=At;
                       if ti==1
                          A=Apart2;
                       else
                          A=[A;Apart2];
                       end
                   end
                elseif H==2
                    rho=[0:0.2:0.8];
                    r=length(vec_p);
                    for ti=1:k
                            At=zeros(vec_p(ti),r);
                            At(1,1)=rho(ii);
                            Apart2=At;
                           if ti==1
                              A=Apart2;
                           else
                              A=[A;Apart2];
                           end
                    end
                elseif H==3
                    rho=[0:0.2:0.8];
                    r=length(vec_p);
                    for ti=1:k
                         if ti<3
                            At=zeros(vec_p(ti),r);
                            At(1,1)=rho(ii);
                         else
                            At=zeros(vec_p(ti),r);
                         end
                            Apart2=At;
                           if ti==1
                              A=Apart2;
                           else
                              A=[A;Apart2];
                           end
                    end
                end
      

                for j=1:N
                    if D==1
                        X=randn(p,n);
                        Z=randn(r,n);
                    elseif D==2
                        X=(chi2rnd(1,p,n)-1)/sqrt(2);
                        Z=(chi2rnd(1,r,n)-1)/sqrt(2);
                    elseif D==3
                        X=trnd(5,p,n)/sqrt(5/3);
                        Z=trnd(5,r,n)/sqrt(5/3);
                    end
           
                     for ti=1:k
                       indi=sum(vec_p(1:ti));
                       At=A(indi-vec_p(ti)+1:indi,1:r);
                       Samplepart2=At*Z;
                       if ti==1
                          Sample2=Samplepart2;
                       else
                          Sample2=[Sample2;Samplepart2];
                       end
                    end
                    
                    sample_m=X+Sample2;
                    
               
                    
                    
                    res_schott(j)=schott_h1(sample_m,vec_p,A);
                    if p<n-1
                        res_wilks(j)=wilks_h1(sample_m,vec_p,A);
                    end
                    
                end
                
                x = -4.5:0.01:4.5;
                norm=normpdf(x,0,1);
                
                
                sizeplot=figure('color',[1 1 1]);
                ztitle=sgtitle(['Simulated distribution of D',num2str(D), ...
                    'G',num2str(G),'S',num2str(S),'p',num2str(p),'rho',num2str(rho(ii)),'H',num2str(H)]);
               ztitle.FontSize = 40;
                subplot(1,2,1)
                [kernel ,xi] = ksdensity(res_schott);
                plot(xi,kernel,'--','Color','k','LineWidth',2);
                xlim([-4.2 4.2]);
                ylim([0 0.45]);
                hold on
                plot(x,norm,'-','Color','k','LineWidth',2);
                title('Schott')
                set(gca,'FontSize',20);
                set(gca,'XTick', [-4:2:4])
                set(gca,'YTick', [0:0.1:0.4])
                set(gca,'Position', [0.05,0.4,0.4,0.3])
                legend('Sample', 'Asymptotic')
                grid on
                
                if all(res_wilks==0)==0
                    subplot(1,2,2)
                    [kernel ,xi] = ksdensity(res_wilks);
                    plot(xi,kernel,'--','Color','k','LineWidth',2);
                    xlim([-4.2 4.2]);
                    ylim([0 0.45]);
                    hold on
                    plot(x,norm,'-','Color','k','LineWidth',2);
                    title('Wilks')
                    set(gca,'FontSize',20);
                    set(gca,'XTick', [-4:2:4])
                    set(gca,'YTick', [0:0.1:0.4])
                    set(gca,'Position', [0.52,0.4,0.4,0.3])
                    legend('Sample', 'Asymptotic')
                    grid on
                end

                saveas(sizeplot,['Simulated distribution of D',num2str(D), ...
                    'G',num2str(G),'S',num2str(S),'p',num2str(p),'rho',num2str(ii),'H',num2str(H)]);
                close(sizeplot);
          end
            end
        end
    end
end

