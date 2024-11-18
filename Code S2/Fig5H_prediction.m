clear
% close all
load('nk_data_240503.mat')
% load('k562_hla_202107')
% load pastel_rainbow
% load colormap_202107
% load kd_oe_202107
 load('serial_KD_OE.mat')

%%
% Screening data : 13 liquid tumor cell lines
x_data=mfi_all(2:13,5:6);
y_data=pNK_05_1(2:13);


x_data=[x_data];
x_data=repmat(x_data,3,1);
% 
% y_data=[pNK_05_1(1:13);kd_pNK_05_1_k562_copy;k562_hla_pNK_05_1;...
%     pNK_1_1(1:13); kd_pNK_1_1_k562_copy;k562_hla_pNK_1_1;];
y_data=[pNK_05_1(2:13);pNK_05_1(2:13);pNK_1_1(2:13);];
%%

% c0_b7h6_hla=2*rand*ones(1,8);
% c0_b7h6_hla=2*rand*ones(1,8);
c0_b7h6_hla=2*rand(1,8);

 g= @(const, x_data) comp_hill_b7h6_hla(const, x_data);
[coef_value2,res2] = lsqcurvefit(g,c0_b7h6_hla,x_data,y_data,[zeros(1,8)],[100,100,100,100*ones(1,5)]);
res2


% [X,Y]=meshgrid(xx,yy);
% coef_value2

% coef_value2=[45.5737    7.2999    0.9897    2.9155  100.0000    2.3746    1.5377];
c=coef_value2(1);
k_b0=coef_value2(2);
h_b=coef_value2(3);
k_h=coef_value2(4);
h_h=coef_value2(5);
scaling1=coef_value2(6);
scaling2=coef_value2(7);
base1=coef_value2(8);

% h_b=5;
x_b7=1:0.01:9;

Y=3;
k_b=k_b0;
Z= c* (base1+((x_b7-1).^h_b)./((x_b7-1).^h_b+k_b.^h_b));
Z=scaling2*Z;
Z=min(100,Z);   
Z=max(0,Z);

Y2=1.5;
k_b2=k_b0;
Z2= c* (base1+((x_b7-1).^h_b)./((x_b7-1).^h_b+k_b2.^h_b));
Z2=scaling2*Z2;
Z2=min(100,Z2);   
Z2=max(0,Z2);



%% gradient data loading
% 
% % load('serial_KD_OE.mat')
% % 
% % 
% % figure
% % hold on
% % scatter([U937_OE(:,1)],[U937_OE(:,3)],80,[255,149,0]/255,'MarkerFaceColor',[1,1,1])
% % scatter([U937_OE(:,1)],[U937_OE(:,4)],80,[255,149,0]/255,'MarkerFaceColor',[1,1,1])
% % scatter([Raji_OE(:,1);Raji_OE(:,1)],[Raji_OE(:,3);Raji_OE(:,4)],80,[0,80,255]/255,'MarkerFaceColor',[1,1,1])
% % scatter([K562_KD(:,1)],[K562_KD(:,2)],80,[255,45,85]/255,'MarkerFaceColor',[1,1,1])
% 
% % load('hla_kd_full')
% % scatter(log10(k562_hla_x_full(:,1)),k562_hla_x_full(:,3)-k562_hla_x_full(:,4),120,'MarkerFaceColor',[255,59,48]/255,'MarkerEdgeAlpha',0)
% % scatter(log10(raji_hla_x_full(:,1)),raji_hla_x_full(:,3),120,'MarkerFaceColor',[0,61,185]/255,'MarkerEdgeAlpha',0)
% 
% load('sgb2m')
% figure
% hold on
% for i=1
% scatter(ligand2(1:2,1),mean(cytotoxicity(1:2,1:2),2),80,[255,59,48]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand2(3:4,1),mean(cytotoxicity(3:4,1:2),2),80,[255,149,0]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand2(5:10,1),mean(cytotoxicity(5:10,1:2),2),80,[0,122,255]/255,'MarkerFaceColor',[1,1,1])
% end
% errorbar(ligand2(1:2,1),mean(cytotoxicity(1:2,1:2),2),std(cytotoxicity(1:2,1:2),0,2),'Color',[255,59,48]/255,'LineStyle','none')
% errorbar(ligand2(3:4,1),mean(cytotoxicity(3:4,1:2),2),std(cytotoxicity(3:4,1:2),0,2),'Color',[255,149,0]/255,'LineStyle','none')
% errorbar(ligand2(5:10,1),mean(cytotoxicity(5:10,1:2),2),std(cytotoxicity(5:10,1:2),0,2),'Color',[0,122,255]/255,'LineStyle','none')
% 
% for i=3
% scatter(ligand(1:2,1),mean(cytotoxicity(1:2,3:6),2),80,[255,59,48]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand(3:4,1),mean(cytotoxicity(3:4,3:6),2),80,[255,149,0]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand(5:10,1),mean(cytotoxicity(5:10,3:6),2),80,[0,122,255]/255,'MarkerFaceColor',[1,1,1])
% end
% errorbar(ligand(1:2,1),mean(cytotoxicity(1:2,3:6),2),std(cytotoxicity(1:2,3:6),0,2),'Color',[255,59,48]/255,'LineStyle','none')
% errorbar(ligand(3:4,1),mean(cytotoxicity(3:4,3:6),2),std(cytotoxicity(3:4,3:6),0,2),'Color',[255,149,0]/255,'LineStyle','none')
% errorbar(ligand(5:10,1),mean(cytotoxicity(5:10,3:6),2),std(cytotoxicity(5:10,3:6),0,2),'Color',[0,122,255]/255,'LineStyle','none')
hold on

plot(x_b7,Z,'LineWidth',1,'Color',[255,59,48]/255)

plot(x_b7,Z2,'LineWidth',1,'Color',[0,122,255]/255)

scatter(Raji_OE(:,1),Raji_OE(:,3))
scatter(Raji_OE(:,1),Raji_OE(:,4))
scatter(U937_OE(:,1),U937_OE(:,3))
scatter(U937_OE(:,1),U937_OE(:,4))

xticks([1,9])
xlim([1,9])
ylim([0,110])
yticks(0:20:100)
set(gca,'TickDir','out'); 
xlabel('B7H6 MFI')
ylabel('cytotoxicity')

%%
% figure
% hold on
% for i=1:2
% scatter(ligand2(1:2,1),cytotoxicity(1:2,i),80,[255,59,48]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand2(3:4,1),cytotoxicity(3:4,i),80,[255,149,0]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand2(5:10,1),cytotoxicity(5:10,i),80,[0,122,255]/255,'MarkerFaceColor',[1,1,1])
% end
% 
% 
% for i=3:6
% scatter(ligand(1:2,1),cytotoxicity(1:2,i),80,[255,59,48]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand(3:4,1),cytotoxicity(3:4,i),80,[255,149,0]/255,'MarkerFaceColor',[1,1,1])
% scatter(ligand(5:10,1),cytotoxicity(5:10,i),80,[0,122,255]/255,'MarkerFaceColor',[1,1,1])
% end

% plot(x_b7,Z,'LineWidth',1,'Color',[255,59,48]/255)
% 
% plot(x_b7,Z2,'LineWidth',1,'Color',[0,122,255]/255)
% xticks([1,9])
% xlim([1,9])
% ylim([0,110])
% yticks(0:20:100)
% set(gca,'TickDir','out'); 
% xlabel('B7H6 MFI')
% ylabel('cytotoxicity')
%%

function y=comp_hill_b7h6_hla(const, x_data) 
% [m,~]=size(x_data);
c=const(1);
k_b0=const(2);
h_b=const(3);
k_h=const(4);
% k_h=1.3;
h_h=const(5);
% base0=const(6);
scaling1=const(6);
scaling2=const(7);

base1=const(8);
% y=zeros(m,1); 
b7=x_data(:,1)-1;
hla=x_data(:,2)-1;
y= c* (base1+(b7.^h_b)./(b7.^h_b+k_b0.^h_b));
% y=y.*(base0+(k_h^h_h)./(x_data(:,2).^h_h+k_h^h_h));
kd_cut=length(x_data)/2;

% y(kd_cut+1:end)=scaling*y(kd_cut+1:end);

y(13:24)=scaling1*y(13:24);
y(25:36)=scaling2*y(25:36);

y=min(100,y);
y=max(0,y);
% y=log(y);
end