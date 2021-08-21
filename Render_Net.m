function Render_Net(V,Parcellation,color)

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation]));
Brain_Mask=niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex=smooth3(Brain_Mask>0);

Centroid_Vol=zeros(size(Volume));
% To color all areas above zero
for n=find(V>=0)
    Centroid_Vol(Volume==n)=1;
end
sregion=smooth3(Centroid_Vol>0);

hold on
% First plot a transparent cortex
cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
reducepatch(cortexpatch,0.01);
isonormals(scortex,cortexpatch);

% Color in Black if no significant overlap with any Yeo RSN
patch(isosurface(sregion,0.3), 'FaceColor', color, 'EdgeColor', 'none') %'FaceAlpha', V(n))

material dull
lighting gouraud
view(-90,90) % Top view
daspect([1 1 1])
camlight;
xlim([5 105])
ylim([5 85])
zlim([0 80])
axis off