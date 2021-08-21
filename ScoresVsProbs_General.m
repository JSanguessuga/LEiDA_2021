%%%%%%%%
%
% Compare Psichological and Physiological markers with probability and
% dwell times of LEiDA networks
%
% Modified by João Sanguessuga, April 2021
%
%
%%%%%%%%
kmeans_file  = 'LEiDA_Kmeans_results';
load(kmeans_file,'Kmeans_results','rangeK');
load Scores.mat
stats_file = 'LEiDA_stats';
load(stats_file,'P','LT')

disp('School_years')

for k=1:length(rangeK)
    
    for c=1:rangeK(k)
        
        a=squeeze(P(:,k,c));
        [cc,p]=corrcoef(a,School_years);
        Pkc_pval_school_years(k,c)=p(2);

        if Pkc_pval_school_years(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end  
        
        b=squeeze(LT(:,k,c));
        [cc,p]=corrcoef(b,School_years);
        LTkc_pval_school_years(k,c)=p(2);

        if LTkc_pval_school_years(k,c)<(0.05)
            disp(['signif correlation with LT K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end    
    end
end

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    
    for c=1:rangeK(k)        
        if Pkc_pval_school_years(k,c)>0.05
            semilogy(rangeK(k),Pkc_pval_school_years(k,c),'*k');
            text(rangeK(k),Pkc_pval_school_years(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_school_years(k,c)<0.05 && Pkc_pval_school_years(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Pkc_pval_school_years(k,c),'*r');
            text(rangeK(k),Pkc_pval_school_years(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_school_years(k,c)<(0.05/rangeK(k)) && Pkc_pval_school_years(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_school_years(k,c),'*g');
            text(rangeK(k),Pkc_pval_school_years(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_school_years(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_school_years(k,c),'*b');
            text(rangeK(k),Pkc_pval_school_years(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs School Years')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    
    for c=1:rangeK(k)        
        if LTkc_pval_school_years(k,c)>0.05
            semilogy(rangeK(k),LTkc_pval_school_years(k,c),'*k');
            text(rangeK(k),LTkc_pval_school_years(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_school_years(k,c)<0.05 && LTkc_pval_school_years(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LTkc_pval_school_years(k,c),'*r');
            text(rangeK(k),LTkc_pval_school_years(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_school_years(k,c)<(0.05/rangeK(k)) && LTkc_pval_school_years(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_school_years(k,c),'*g');
            text(rangeK(k),LTkc_pval_school_years(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_school_years(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_school_years(k,c),'*b');
            text(rangeK(k),LTkc_pval_school_years(k,c),[' ' num2str(c)])
        end
    end
end

title('Lifetimes Networks vs School Years')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

disp('Age')

for k=1:length(rangeK)
    
    for c=1:rangeK(k)
        
        a=squeeze(P(:,k,c));
        [cc,p]=corrcoef(a,Age);
        Pkc_pval_age(k,c)=p(2);

        if Pkc_pval_age(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end
        
        b=squeeze(LT(:,k,c));
        [cc,p]=corrcoef(b,Age);
        LTkc_pval_age(k,c)=p(2);

        if LTkc_pval_age(k,c)<(0.05)
            disp(['signif correlation with LT K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end    
    end
end

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK)
    
    for c=1:rangeK(k)        
        if Pkc_pval_age(k,c)>0.05
            semilogy(rangeK(k),Pkc_pval_age(k,c),'*k');
            text(rangeK(k),Pkc_pval_age(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_age(k,c)<0.05 && Pkc_pval_age(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Pkc_pval_age(k,c),'*r');
            text(rangeK(k),Pkc_pval_age(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_age(k,c)<(0.05/rangeK(k)) && Pkc_pval_age(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_age(k,c),'*g');
            text(rangeK(k),Pkc_pval_age(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_age(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_age(k,c),'*b');
            text(rangeK(k),Pkc_pval_age(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs Age')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    
    for c=1:rangeK(k)        
        if LTkc_pval_age(k,c)>0.05
            semilogy(rangeK(k),LTkc_pval_age(k,c),'*k');
            text(rangeK(k),LTkc_pval_age(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_age(k,c)<0.05 && LTkc_pval_age(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LTkc_pval_age(k,c),'*r');
            text(rangeK(k),LTkc_pval_age(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_age(k,c)<(0.05/rangeK(k)) && LTkc_pval_age(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_age(k,c),'*g');
            text(rangeK(k),LTkc_pval_age(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_age(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_age(k,c),'*b');
            text(rangeK(k),LTkc_pval_age(k,c),[' ' num2str(c)])
        end
    end
end

title('Lifetimes Networks vs Age')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

disp('Cortisol')

for k=1:length(rangeK)
    
    for c=1:rangeK(k)
         
        a=squeeze(P(:,k,c));
        [cc,p]=corrcoef(a,Cortisol);
        Pkc_pval_cortisol(k,c)=p(2);

        if Pkc_pval_cortisol(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end
        
        b=squeeze(LT(:,k,c));
        [cc,p]=corrcoef(b,Cortisol);
        LTkc_pval_cortisol(k,c)=p(2);

        if LTkc_pval_cortisol(k,c)<(0.05)
            disp(['signif correlation with LT K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end    
    end
end

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK)
    
    for c=1:rangeK(k)        
        if Pkc_pval_cortisol(k,c)>0.05
            semilogy(rangeK(k),Pkc_pval_cortisol(k,c),'*k');
            text(rangeK(k),Pkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_cortisol(k,c)<0.05 && Pkc_pval_cortisol(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Pkc_pval_cortisol(k,c),'*r');
            text(rangeK(k),Pkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_cortisol(k,c)<(0.05/rangeK(k)) && Pkc_pval_cortisol(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_cortisol(k,c),'*g');
            text(rangeK(k),Pkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_cortisol(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_cortisol(k,c),'*b');
            text(rangeK(k),Pkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs Cortisol')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK)
    
    for c=1:rangeK(k)        
        if LTkc_pval_cortisol(k,c)>0.05
            semilogy(rangeK(k),LTkc_pval_cortisol(k,c),'*k');
            text(rangeK(k),LTkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_cortisol(k,c)<0.05 && LTkc_pval_cortisol(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LTkc_pval_cortisol(k,c),'*r');
            text(rangeK(k),LTkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_cortisol(k,c)<(0.05/rangeK(k)) && LTkc_pval_cortisol(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_cortisol(k,c),'*g');
            text(rangeK(k),LTkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_cortisol(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_cortisol(k,c),'*b');
            text(rangeK(k),LTkc_pval_cortisol(k,c),[' ' num2str(c)])
        end
    end
end

title('Lifetimes Networks vs Cortisol')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off


disp('MMSE')

for k=1:length(rangeK)
    
    for c=1:rangeK(k)
        
        a=squeeze(P(:,k,c));
        [cc,p]=corrcoef(a,MMSE);
        Pkc_pval_MMSE(k,c)=p(2);

        if Pkc_pval_MMSE(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end
        
        b=squeeze(LT(:,k,c));
        [cc,p]=corrcoef(b,MMSE);
        LTkc_pval_MMSE(k,c)=p(2);

        if LTkc_pval_MMSE(k,c)<(0.05)
            disp(['signif correlation with LT K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end    
    end
end

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    
    for c=1:rangeK(k)        
        if Pkc_pval_MMSE(k,c)>0.05
            semilogy(rangeK(k),Pkc_pval_MMSE(k,c),'*k');
            text(rangeK(k),Pkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_MMSE(k,c)<0.05 && Pkc_pval_MMSE(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Pkc_pval_MMSE(k,c),'*r');
            text(rangeK(k),Pkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_MMSE(k,c)<(0.05/rangeK(k)) && Pkc_pval_MMSE(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_MMSE(k,c),'*g');
            text(rangeK(k),Pkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_MMSE(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_MMSE(k,c),'*b');
            text(rangeK(k),Pkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs MMSE')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    
    for c=1:rangeK(k)        
        if LTkc_pval_MMSE(k,c)>0.05
            semilogy(rangeK(k),LTkc_pval_MMSE(k,c),'*k');
            text(rangeK(k),LTkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_MMSE(k,c)<0.05 && LTkc_pval_MMSE(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LTkc_pval_MMSE(k,c),'*r');
            text(rangeK(k),LTkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_MMSE(k,c)<(0.05/rangeK(k)) && LTkc_pval_MMSE(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_MMSE(k,c),'*g');
            text(rangeK(k),LTkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_MMSE(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_MMSE(k,c),'*b');
            text(rangeK(k),LTkc_pval_MMSE(k,c),[' ' num2str(c)])
        end
    end
end

title('Lifetimes Networks vs MMSE')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

disp('GDS')

for k=1:length(rangeK)
    
    for c=1:rangeK(k)
        
        a=squeeze(P(:,k,c));
        [cc,p]=corrcoef(a,GDS);
        Pkc_pval_GDS(k,c)=p(2);

        if Pkc_pval_GDS(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end
        
        b=squeeze(LT(:,k,c));
        [cc,p]=corrcoef(b,GDS);
        LTkc_pval_GDS(k,c)=p(2);

        if LTkc_pval_GDS(k,c)<(0.05)
            disp(['signif correlation with LT K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc(2)) ' p=' num2str(p(2))])
        end    
    end
end

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK)
    
    for c=1:rangeK(k)        
        if Pkc_pval_GDS(k,c)>0.05
            semilogy(rangeK(k),Pkc_pval_GDS(k,c),'*k');
            text(rangeK(k),Pkc_pval_GDS(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_GDS(k,c)<0.05 && Pkc_pval_GDS(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Pkc_pval_GDS(k,c),'*r');
            text(rangeK(k),Pkc_pval_GDS(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_GDS(k,c)<(0.05/rangeK(k)) && Pkc_pval_GDS(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_GDS(k,c),'*g');
            text(rangeK(k),Pkc_pval_GDS(k,c),[' ' num2str(c)])
        end
        if Pkc_pval_GDS(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Pkc_pval_GDS(k,c),'*b');
            text(rangeK(k),Pkc_pval_GDS(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs GDS')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK)
    
    for c=1:rangeK(k)        
        if LTkc_pval_GDS(k,c)>0.05
            semilogy(rangeK(k),LTkc_pval_GDS(k,c),'*k');
            text(rangeK(k),LTkc_pval_GDS(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_GDS(k,c)<0.05 && LTkc_pval_GDS(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LTkc_pval_GDS(k,c),'*r');
            text(rangeK(k),LTkc_pval_GDS(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_GDS(k,c)<(0.05/rangeK(k)) && LTkc_pval_GDS(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_GDS(k,c),'*g');
            text(rangeK(k),LTkc_pval_GDS(k,c),[' ' num2str(c)])
        end
        if LTkc_pval_GDS(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LTkc_pval_GDS(k,c),'*b');
            text(rangeK(k),LTkc_pval_GDS(k,c),[' ' num2str(c)])
        end
    end
end

title('Lifetimes Networks vs GDS')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off
