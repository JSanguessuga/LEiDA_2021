%%%%%%%%
%
% Compare Psichological and Physiological markers with probability and 
% dwell times of LEiDA networks
%
%
%  
%
%%%%%%%%

load LEiDA_OCD_for_stats

% Remove Patient 
Index_Patients=Index_Patients([1:32 34:end]);

%% Patient specific scores
Patient_Scores=[6 7 8 13 14 15];

OCD_P_pval=zeros(length(Patient_Scores),length(rangeK),max(rangeK));
OCD_LT_pval=zeros(length(Patient_Scores),length(rangeK),max(rangeK));

for score=1:length(Patient_Scores)
    
    disp(' ')
    disp([BDOCD.Properties.VariableNames{Patient_Scores(score)}  ' - Score ' num2str(Patient_Scores(score))])

    Scores=table2array(BDOCD(:,Patient_Scores(score)));
    
    for k=1:length(rangeK)
        
        for c=1:rangeK(end)
            
            P90_OCD=squeeze(P90(Index_Patients,k,c));
            [~, p]=corrcoef(P90_OCD,Scores(Index_Patients));
            
            OCD_P_pval(score,k,c)=p(2);
          
            LT90_OCD=squeeze(LT90(Index_Patients,k,c));
            [~, p]=corrcoef(LT90_OCD,Scores(Index_Patients));
            
            OCD_LT_pval(score,k,c)=p(2);
            
        end
    end
end

%% Global scores
Global_Scores=[4 5 9 10 11 12];

ALL_P_pval=zeros(length(Global_Scores),length(rangeK),max(rangeK));
ALL_LT_pval=zeros(length(Global_Scores),length(rangeK),max(rangeK));

for score=1:length(Global_Scores)
    
    disp(' ')
    disp([BDOCD.Properties.VariableNames{Global_Scores(score)}  ' - Score ' num2str(Global_Scores(score))])

    Scores=table2array(BDOCD(:,Global_Scores(score)));
    s_include=find(~isnan(Scores));
    Scores=Scores(s_include);
    
    for k=1:length(rangeK)
        
        for c=1:rangeK(end)
            
            P90_all=squeeze(P90(s_include,k,c));
            [~, p]=corrcoef(P90_all,Scores);
            
            ALL_P_pval(score,k,c)=p(2);
            
            LT90_all=squeeze(LT90(s_include,k,c));
            [~, p]=corrcoef(LT90_all,Scores);
            
            ALL_LT_pval(score,k,c)=p(2);
            
        end
    end
end

%% Plot Patient scores 

figure('Name','Probability and Lifetimes vs Patient Scores')

for score=1:length(Patient_Scores)
    
    subplot(2,length(Patient_Scores),score)
    
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'--','Color',[0 0.7 .3],'LineWidth',1)

    for k=1:length(rangeK)
        for c=1:rangeK(k)
            if OCD_P_pval(score,k,c)>0.05
                semilogy(rangeK(k),OCD_P_pval(score,k,c),'.k');
                %text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
            end
            if OCD_P_pval(score,k,c)<0.05 && OCD_P_pval(score,k,c)>(0.05/rangeK(k))
                semilogy(rangeK(k),OCD_P_pval(score,k,c),'.r');
                %text(rangeK(k),OCD_P_pval(score,k,c),[' ' num2str(c)])
            end
            if OCD_P_pval(score,k,c)<(0.05/rangeK(k))
                semilogy(rangeK(k),OCD_P_pval(score,k,c),'*','Color',[0 0.7 .3]);
                text(rangeK(k),OCD_P_pval(score,k,c),[' ' num2str(c)])
                
            end
        end
    end
    title(BDOCD.Properties.VariableNames{Patient_Scores(score)},'Interpreter', 'none');
    xlim([0 21])
    ylim([0.0005 1])
    box off
    if score==1
            ylabel({'Probabilities vs score','Pearson  p-value'})
    end

    subplot(2,length(Patient_Scores),score+length(Patient_Scores))
    
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'--','Color',[0 0.7 .3],'LineWidth',1)

    for k=1:length(rangeK)
        for c=1:rangeK(k)
            if OCD_LT_pval(score,k,c)>0.05
                semilogy(rangeK(k),OCD_LT_pval(score,k,c),'.k');
                %text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
            end
            if OCD_LT_pval(score,k,c)<0.05 && OCD_LT_pval(score,k,c)>(0.05/rangeK(k))
                semilogy(rangeK(k),OCD_LT_pval(score,k,c),'.r');
                %text(rangeK(k),OCD_LT_pval(score,k,c),[' ' num2str(c)])
            end
            if OCD_LT_pval(score,k,c)<(0.05/rangeK(k))
                semilogy(rangeK(k),OCD_LT_pval(score,k,c),'*','Color',[0 0.7 .3]);
                text(rangeK(k),OCD_LT_pval(score,k,c),[' ' num2str(c)])
                
            end
        end
    end
    xlim([0 21])
    ylim([0.0005 1])
    xlabel('Number of states K')
    if score==1
            ylabel({'Duration vs score','Pearson  p-value'})
    end
    box off
end

%% Plots networks that relate significantly with each Patient score

load LEiDA90_results.mat Kmeans_results 
load AAL_labels.mat label116

for score=1:length(Patient_Scores)
    % Bonferroni correction
    OCD_LT_pcorr=squeeze(OCD_LT_pval(score,:,:)).*rangeK';  % Correct p-values by the number of states;
    [k,c]=ind2sub([length(rangeK),max(rangeK)],find(OCD_LT_pcorr<0.05 & OCD_LT_pcorr>0));

    if k
        
        Scores=table2array(BDOCD(Index_Patients,Patient_Scores(score)));
        
        [~,b]=sort(k);
        
        for h=b'
            figure
            subplot(4,2,1)
            V=Kmeans_results{k(h)}.C(c(h),:);
            %V=round(V*1e2)/1e2;
            V=V/max(abs(V));
            
            Render_Net(V,'AAL116','g')
            title(['K=' num2str(rangeK(k(h))) ' C=' num2str(c(h))])
            
            subplot(4,2,3) 
            Render_Net(V,'AAL116','g')
            view(0,0)            
            
            subplot(4,2,5)
            plot(squeeze(LT90(Index_Patients,k(h),c(h))),Scores,'.')
            xlabel('Duration (s)')
            [cc, p]=corrcoef(squeeze(LT90(Index_Patients,k(h),c(h))),Scores);
            title({['r=' num2str(cc(2))],['p=' num2str(p(2),'%1.2e')]})
            ylabel( BDOCD.Properties.VariableNames{Patient_Scores(score)})
            
            subplot(4,2,[2 4 6 8])
            hold on
            barh(find(V<0),V(V<0),'FaceColor','k','EdgeColor','none','Barwidth',.5)
            barh(find(V>=0),V(V>=0),'FaceColor','green','EdgeColor','none','Barwidth',.5)
            set(gca,'YTick',1:90,'Fontsize',7)
            set(gca,'YTickLabel',label90(:,:))
            set(gca,'Ydir','reverse')
            ylim([0 117])
            xlim([-1 1])
            title({'BOLD phase','projection into V1'},'Fontsize',12)
            grid on

            
            
        end
    end
    
        % Bonferroni correction
    OCD_P_pcorr=squeeze(OCD_P_pval(score,:,:)).*rangeK';  % Correct p-values by the number of states;
    [k,c]=ind2sub([length(rangeK),max(rangeK)],find(OCD_P_pcorr<0.05 & OCD_P_pcorr>0));

    if k
        
        figure('Name',['Probabilitites vs ' BDOCD.Properties.VariableNames{Patient_Scores(score)}])
        Scores=table2array(BDOCD(:,Patient_Scores(score)));
        
        [~,b]=sort(k);
        
        u=0;
        for h=b'
            u=u+1;
            subplot(numel(b),2,2*u-1)
            V=Kmeans_results{k(h)}.C(c(h),:);
            Render_Net(V,'AAL116','g')
            title(['K=' num2str(rangeK(k(h))) ' C=' num2str(c(h))])
            
            subplot(numel(b),2,2*u)
            plot(squeeze(P90(Index_Patients,k(h),c(h))),Scores,'.')
            xlabel('Probability')
            ylabel( BDOCD.Properties.VariableNames{Patient_Scores(score)})
            title(['p=' num2str(OCD_P_pcorr(k(h),c(h)),'%1.2e') ])
            
        end
    end
    
end



%% Plot Global scores 
figure('Name','Probability and Lifetimes vs Global Scores')

for score=1:length(Global_Scores)
    
    subplot(3,2,score)
    
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05/length(Global_Scores)*ones(1,length(rangeK)+2),'b--','LineWidth',1)  
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05/length(Global_Scores)./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'--','Color',[0 0.7 .3],'LineWidth',1)

    for k=1:length(rangeK)
        for c=1:rangeK(k)
            if ALL_P_pval(score,k,c)>0.05
                semilogy(rangeK(k),ALL_P_pval(score,k,c),'*k');
                %text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
            end
            if ALL_P_pval(score,k,c)<0.05 && ALL_P_pval(score,k,c)>(0.05/length(Global_Scores)/rangeK(k))
                semilogy(rangeK(k),ALL_P_pval(score,k,c),'*r');
                text(rangeK(k),ALL_P_pval(score,k,c),[' ' num2str(c)])
            end
            if ALL_P_pval(score,k,c)<(0.05/length(Global_Scores)/rangeK(k))
                semilogy(rangeK(k),ALL_P_pval(score,k,c),'*','Color',[0 0.7 .3]);
                text(rangeK(k),ALL_P_pval(score,k,c),[' ' num2str(c)])
            end
        end
    end
    title(BDOCD.Properties.VariableNames{Global_Scores(score)},'Interpreter', 'none');
end

figure('Name','Lifetimes')

for score=1:length(Global_Scores)
    
    subplot(3,2,score)
    
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05/length(Global_Scores)*ones(1,length(rangeK)+2),'b--','LineWidth',1)  
    semilogy(rangeK(1)-1:rangeK(end)+1,0.05/length(Global_Scores)./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'--','Color',[0 0.7 .3],'LineWidth',1)

    for k=1:length(rangeK)
        for c=1:rangeK(k)
            if ALL_LT_pval(score,k,c)>0.05
                semilogy(rangeK(k),ALL_LT_pval(score,k,c),'*k');
                %text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
            end
            if ALL_LT_pval(score,k,c)<0.05 && ALL_LT_pval(score,k,c)>(0.05/length(Global_Scores)/rangeK(k))
                semilogy(rangeK(k),ALL_LT_pval(score,k,c),'*r');
                text(rangeK(k),ALL_LT_pval(score,k,c),[' ' num2str(c)])
            end
            if ALL_LT_pval(score,k,c)<(0.05/length(Global_Scores)/rangeK(k))
                semilogy(rangeK(k),ALL_LT_pval(score,k,c),'*','Color',[0 0.7 .3]);
                text(rangeK(k),ALL_LT_pval(score,k,c),[' ' num2str(c)])
            end
        end
    end
    title(BDOCD.Properties.VariableNames{Global_Scores(score)},'Interpreter', 'none');
end


%% Plots networks that relate significantly with each Global score

for score=1:length(Global_Scores)
    % Bonferroni correction
    ALL_LT_pcorr=squeeze(ALL_LT_pval(score,:,:)).*rangeK'*length(Global_Scores);  % Correct p-values by the number of states;
    [k,c]=ind2sub([length(rangeK),max(rangeK)],find(ALL_LT_pcorr<0.05 & ALL_LT_pcorr>0));

    if k
        
        figure('Name',['State Duration vs ' BDOCD.Properties.VariableNames{Global_Scores(score)}])
        Scores=table2array(BDOCD(:,Global_Scores(score)));
        
        [~,b]=sort(k);
        
        u=0;
        for h=b'
            u=u+1;
            subplot(numel(b),2,2*u-1)
            V=Kmeans_results{k(h)}.C(c(h),:);
            Render_Net(V,'AAL116','r')
            title(['K=' num2str(rangeK(k(h))) ' C=' num2str(c(h))])
            
            subplot(numel(b),2,2*u)
            plot(squeeze(LT90(:,k(h),c(h))),Scores,'.')
            xlabel('Duration (s)')
            title(['p=' num2str(ALL_LT_pcorr(k(h),c(h)),'%1.2e') ])
            ylabel( BDOCD.Properties.VariableNames{Global_Scores(score)})
        end
    end
    
        % Bonferroni correction
    ALL_P_pcorr=squeeze(ALL_P_pval(score,:,:)).*rangeK'*length(Global_Scores);  % Correct p-values by the number of states;
    [k,c]=ind2sub([length(rangeK),max(rangeK)],find(ALL_P_pcorr<0.05 & ALL_P_pcorr>0));

    if k
        
        figure('Name',['Probabilitites vs ' BDOCD.Properties.VariableNames{Global_Scores(score)}])
        Scores=table2array(BDOCD(:,Global_Scores(score)));
        
        [~,b]=sort(k);
        
        u=0;
        for h=b'
            u=u+1;
            subplot(numel(b),2,2*u-1)
            V=Kmeans_results{k(h)}.C(c(h),:);
            Render_Net(V,'AAL116','r')
            title(['K=' num2str(rangeK(k(h))) ' C=' num2str(c(h))])
            
            subplot(numel(b),2,2*u)
            plot(squeeze(P90(:,k(h),c(h))),Scores,'.')
            xlabel('Probability')
            ylabel( BDOCD.Properties.VariableNames{Global_Scores(score)})
            title(['p=' num2str(ALL_P_pcorr(k(h),c(h)),'%1.2e') ])
            
        end
    end
    
end
