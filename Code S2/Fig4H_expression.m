%% Code for Fig. 4H 


clear
close all
load('nk_data_240503.mat')
load('colormaps')
% load('origin_list.mat')
% load('tumor_cat_unique')
load('validation_cell_line_ligand_v2')
% scatter(b7h6_all(2:13,1),pNK_05_1(2:13),150,'MarkerFaceColor',[218,41,28]/255,'MarkerEdgeAlpha',0);


valid_ind=[1:35];
A=mfi_all(valid_ind,[1:5]);

B=mfi_all(valid_ind,[1:5]);

valid_A=validation_cell_line_ligand([1:5],[1:5]);
% outliers1=A<prctile(A,10);
% outliers2=A>prctile(A,90);
% for i=1:6
%     B(outliers1(:,i),i)=prctile(A(:,i),10);
%     B(outliers2(:,i),i)=prctile(A(:,i),90);
% 
% end

% 
% 
BB=B(:,1:5);
% BBB=B(:,6);
% B=(B-min(B))./(max(B)-min(B));
% B=(B-min(B))./(max(B)-min(B));
B(:,1:5)=(B(:,1:5)-min(min(BB)))./(max(max(BB))-min(min(BB)));
% B(:,6)=(B(:,6)-min(B(:,6)))./(max(B(:,6))-min(B(:,6)));

valid_A(:,1:5)=(valid_A(:,1:5)-min(min(BB)))./(max(max(BB))-min(min(BB)));
% new_A(:,6)=(new_A(:,6)-min(BBB))./(max(BBB)-min(BBB));


C=[validation_cell_line_cyto_1_1/100];

%%
figure
ax1=axes;

im1=imagesc(ax1,1:5,1:5,valid_A');

colormap_magma=magma(100);
colormap_magma=[colormap_magma(15:end,:);repmat([1,1,0.8],3,1)];

coll=1;
hold on




ax2=axes;

im2=imagesc(ax2,1:5,6,C);
caxis([0.1,0.6])
% xticks([])
% yticks([])
linkaxes([ax1,ax2])
% ax1.Visible='off';


tt=text(0*ones(1,6),(1:6)+0.,["MIC","BP1","BP2","BP3","B7","NK1"],'FontSize',11,'FontWeight','bold','HorizontalAlignment', 'right');
set(tt,'Rotation',0)

% text(-3.3,-1.5,"Red: Liquid Tumor",'FontSize',12,'FontWeight','bold','Color',[218,41,28]/255)
% text(-3.3,-0.5,"Blue: Solid Tumor",'FontSize',12,'FontWeight','bold','Color',[0,61,165]/255)
for i=0:6
    line([0.5,5.5],[0.5+i,0.5+i],'Color',[coll,coll,coll],'LineWidth',2)
end
% 
for i=0:5
    line([i+0.5,i+0.5],[0.5, 6.5],'Color',[coll,coll,coll],'LineWidth',2)
end
names=text(1:5,6.8*ones(1,5),["CCRF-SB","Jiyoye","H9","RPMI6666","RPMI8226"],'Color',[0,0,0]/255,'FontSize',13,'FontWeight','bold','HorizontalAlignment','right');
set(names,'Rotation',90)
 
ax2.Visible='off';
ax2.XTick = []; 
ax2.YTick = []; 


colormap(ax1,colormap_magma(end:-1:1,:));
colormap(ax2,colormap_waves(1:200,:));

set([ax1,ax2],'Position',[.05 .10 .90 .80]); 
set([ax1,ax2],'Box','off'); 
set([ax1,ax2],'Visible','off'); 
set([ax1,ax2],'XLim',[-0.5,5.5]); 
set([ax1,ax2],'YLim',[0,7.5]);
%%
