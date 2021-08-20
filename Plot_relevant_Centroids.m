function Plot_relevant_Centroids

%%%%%%%
%
% Code to make detailed Figure of any given network
%
%%%%%%

prompt_1='Qual é o K que pretende avaliar? K= ';
K=input(prompt_1);
prompt_2='Qual é o c que pretende avaliar? c= ';
c=input(prompt_2);

N=90;

Parcellation='AAL116';

stats_file = 'LEiDA_stats';
load(stats_file,'Kmeans_results','rangeK','P','P_pval','LT','LT_pval','Index_Good','Index_Poor')
V_1=Kmeans_results{rangeK==K}.C(c,:);
V_2=Kmeans_results{rangeK==K}.C;
V_3=V_2'./max(abs(V_2'));
Vc=V_1'./max(abs(V_1));
load AAL_labels.mat label90
load Scores.mat

% Same color code from the Yeo paper (CHANGE HERE THE COLORS)
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
Yeo_names= {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

[cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo(Kmeans_results,rangeK);

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation]));
Brain_Mask=niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex=smooth3(Brain_Mask>0);

figure
disp(['k=' num2str(K) ' c=' num2str(c)])

%% First Plot the arrows view from top

[~, net]= max(cc_V_yeo7(rangeK==K,c,:));

subplot_tight(4,3,1)
hold on
cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
reducepatch(cortexpatch,0.01);
isonormals(scortex,cortexpatch);

% Color all areas with positive V with transparency proportional to
% contribution
n_pos=find(V_3(:,c)>0);
if numel(n_pos)>0
    for region=1:length(n_pos)
        region_Vol=Volume==n_pos(region); % 3D volume with 1 only in the voxels of that area
        sregion=smooth3(region_Vol>0);

        if p_V_yeo7(rangeK==K,c,net)<0.05/K
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
        else
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
        end
    end
else
    region_Vol=Volume>0 & Volume<=90; % 3D volume with 1 only in the voxels of that area
    sregion=smooth3(region_Vol>0);
    patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
end
if p_V_yeo7(rangeK==K,c,net)<0.05/K
    title({['K=' num2str(K) ' C=' num2str(c), newline, Yeo_names{net}]} ) %,Yeo_names{c}})
elseif c==1
    title({['K=' num2str(K) ' C=' num2str(c), newline, 'Global Mode']} ) %,Yeo_names{c}})
else
    title({['K=' num2str(K) ' C=' num2str(c), newline, '']} ) %,Yeo_names{c}})
end
material dull; lighting gouraud
daspect([1 1 1])
view(-90,90)  
camlight;
xlim([5 105]); ylim([5 85]); zlim([0 80])
axis off

%% Same but view from the side

subplot_tight(4,3,2)
hold on
cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
reducepatch(cortexpatch,0.01);
isonormals(scortex,cortexpatch);

n_pos=find(V_3(:,c)>0);
if n_pos
    for region=1:length(n_pos)
        region_Vol=Volume==n_pos(region); % 3D volume with 1 only in the voxels of that area
        sregion=smooth3(region_Vol>0);

        if p_V_yeo7(rangeK==K,c,net)<0.05/K
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
        else
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
        end
    end
else
    region_Vol=Volume>0 & Volume<=90; % 3D volume with 1 only in the voxels of that area
    sregion=smooth3(region_Vol>0);
    patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))

end
material dull; lighting gouraud
daspect([1 1 1])
view(0,0)
camlight;
xlim([5 105]); ylim([5 85]); zlim([0 80])
axis off

% Probabilities
subplot(4,3,[4 5])
hold on
Patients=squeeze(P(Index_Poor,rangeK==K,c));  
Controls=squeeze(P(Index_Good,rangeK==K,c));
boxplot([Controls;Patients]',[ones(1,numel(Controls)) 2*ones(1,numel(Patients))]);
hold on
set(gca,'XTickLabel',{'Good', 'Poor'})
%ylim([0 .62])
ylabel('Probability','fontsize',12)
box off

% Lifetimes
subplot(4,3,[7 8])
hold on
Patients=squeeze(LT(Index_Poor,rangeK==K,c));  
Controls=squeeze(LT(Index_Good,rangeK==K,c));
boxplot([Controls;Patients]',[ones(1,numel(Controls)) 2*ones(1,numel(Patients))]);
%ylim([0 30])
set(gca,'XTickLabel',{'Good', 'Poor'})
ylabel('Duration (s)','fontsize',12)
box off

subplot(4,3,[10 11])
hold on
plot(squeeze(P(Index_Good,rangeK==K,c)),School_years(Index_Good),'r*');
plot(squeeze(P(Index_Poor,rangeK==K,c)),School_years(Index_Poor),'b*');
xlabel('Probability','fontsize',12)
ylabel('School Years','fontsize',12)

%Order=[1:2:N N:-2:2];
%Vc=Vc(Order);
cmap=jet(8);
subplot(4,3,[3 6 9 12])
hold on
for n=1:N
    if Vc(n)~=-1
        bin=find(-1:0.25:1>=Vc(n),1)-1;
    else
        bin=1;
    end
    barh(n,Vc(n),'FaceColor',cmap(bin,:),'EdgeColor','none','Barwidth',.5)
end
%barh(find(Vc>=0),Vc(Vc>=0),'FaceColor','red','EdgeColor','none','Barwidth',.5)
set(gca,'YTick',1:N,'Fontsize',7)
set(gca,'YTickLabel',label90(:,:))
set(gca,'Ydir','reverse')
set(gca,'XTick',-1:.25:1,'Fontsize',10)
ylim([0 N+1])
xlim([-1 1])
title({'BOLD phase','projection into V1'},'Fontsize',12)
grid on

