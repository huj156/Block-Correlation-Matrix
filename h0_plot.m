clear;
% s=rng(0);
% save('s','s');
load('s.mat');
rng(s);
N=1000;
p=32;
%p=32*5;
for D=1:3
    for M=1:3
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
                
                if M==1
                    Sigma=eye(p);
                elseif M==2
                    B=ones(1,vec_p(1))*0.5+[1:vec_p(1)]/(vec_p(1)+1);
                    B=diag(sqrt(B));
                    Sij=abs([1:vec_p(1)]'*ones(1,vec_p(1))-([1:vec_p(1)]'*ones(1,vec_p(1)))');
                    Sigma=B*(0.3.^(Sij.^(1/3)))*B;
                    for jj=2:k
                        B=ones(1,vec_p(jj))*0.5+[1:vec_p(jj)]/(vec_p(jj)+1);
                        B=diag(sqrt(B));
                        Sij=abs([1:vec_p(jj)]'*ones(1,vec_p(jj))-([1:vec_p(jj)]'*ones(1,vec_p(jj)))');
                        Sigma=blkdiag(Sigma,B*(0.3.^(Sij.^(1/3)))*B);
                    end
                elseif M==3
                    Srnd=unifrnd(1,5,vec_p(1),2*vec_p(1));
                    Sigma=Srnd*Srnd'/vec_p(1);
                    for jj=2:k
                        Srnd=unifrnd(1,5,vec_p(jj),2*vec_p(jj));
                        Sigma=blkdiag(Sigma,Srnd*Srnd'/vec_p(jj));
                    end
                end
                
                res_schott=zeros(1,N);
                res_wilks=zeros(1,N);
                res_jbz=zeros(1,N);
                res_yhn=zeros(1,N);
                
                for j=1:N
                    if D==1
                        X=randn(p,n);
                    elseif D==2
                        X=(chi2rnd(1,p,n)-1)/sqrt(2);
                    elseif D==3
                        X=trnd(5,p,n)/sqrt(5/3);
                    end
                    sample_m=sqrtm(Sigma)*X;
                    
                    res_schott(j)=schott(sample_m,vec_p);
                    if sum(vec_p(1:k-1))<n-1
                        res_jbz(j)=jbz(sample_m,vec_p);
                    end
                    res_yhn(j)=yhn(sample_m,vec_p);
                    if p<n-1
                        res_wilks(j)=wilks(sample_m,vec_p);
                    end
                    
                end
                
                x = -4.5:0.01:4.5;
                norm=normpdf(x,0,1);
                
                
                sizeplot=figure('color',[1 1 1]);
                
                subplot(1,4,1)
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
                set(gca,'Position', [0.02,0.4,0.225,0.3])
                legend('Sample', 'Asymptotic')
                grid on
                
                if all(res_wilks==0)==0
                    subplot(1,4,2)
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
                    set(gca,'Position', [0.27,0.4,0.225,0.3])
                    legend('Sample', 'Asymptotic')
                    grid on
                end
                
                subplot(1,4,3)
                [kernel ,xi] = ksdensity(res_yhn);
                plot(xi,kernel,'--','Color','k','LineWidth',2);
                xlim([-4.2 4.2]);
                ylim([0 0.45]);
                hold on
                plot(x,norm,'-','Color','k','LineWidth',2);
                title('YHN')
                set(gca,'FontSize',20);
                set(gca,'XTick', [-4:2:4])
                set(gca,'YTick', [0:0.1:0.4])
                set(gca,'Position', [0.52,0.4,0.225,0.3])
                legend('Sample', 'Asymptotic')
                grid on
                
     
                 if all(res_jbz==0)==0
                    subplot(1,4,4)
                    [kernel ,xi] = ksdensity(res_jbz);
                    plot(xi,kernel,'--','Color','k','LineWidth',2);
                    xlim([-4.2 4.2]);
                    ylim([0 0.45]);
                    hold on
                    plot(x,norm,'-','Color','k','LineWidth',2);
                    title('JBZ')
                    set(gca,'FontSize',20);
                    set(gca,'XTick', [-4:2:4])
                    set(gca,'YTick', [0:0.1:0.4])
                    set(gca,'Position', [0.77,0.4,0.225,0.3])
                    legend('Sample', 'Asymptotic')
                    grid on
                 end
                sgtitle(['Simulated distribution of D',num2str(D),'M',num2str(M), ...
                    'G',num2str(G),'S',num2str(S),'p',num2str(p)],'FontSize',40)
                
                saveas(sizeplot,['Simulated distribution of D',num2str(D),'M',num2str(M), ...
                    'G',num2str(G),'S',num2str(S),'p',num2str(p)]);
                close(sizeplot);
            end
        end
    end
end



