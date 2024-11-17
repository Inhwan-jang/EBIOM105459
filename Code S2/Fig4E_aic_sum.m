%% Code for Fig. 4C and 4E solid: Run Fig.4C heatmap after this

% Last updated: 240502
close 
clear
load('nk_data_230622.mat')
leuk_ind=2:13;
x_data=mfi_all(leuk_ind,:)-1;
x_data(x_data<0)=0;
nkg2d_sum=sum(x_data(:,1:4),2);

corr(nkg2d_sum,pNK_025_1(2:13))
corr(nkg2d_sum,pNK_05_1(2:13))
corr(nkg2d_sum,pNK_1_1(2:13))

%%

x_intercept=min(sum(x_data(:,1:4),2));
x_slope=(max(sum(x_data(:,1:4),2))-min(sum(x_data(:,1:4),2)))/10;
%% Data normalization
x_data=[x_data;x_data;x_data];
x_data=sum(x_data(:,1:4),2);

x_data=10*(x_data-min(x_data))./(max(x_data)-min(x_data));

y_data=[pNK_025_1(leuk_ind,1);pNK_05_1(leuk_ind,1);pNK_1_1(leuk_ind,1)];
%%

% leuk_names=["K562","Raji","Ramos", "Jurkat", "U937", "CCRF-CEM", "MOLT-4", "THP-1","Kasumi-1","Kasumi-3", "KG-1", "MV4-11","HL-60"];
leuk_names=["Raji","Ramos", "Jurkat", "U937", "CCRF-CEM", "MOLT-4", "THP-1","Kasumi-1","Kasumi-3", "KG-1", "MV4-11","HL-60"];

%%
c07=ones(1,17);

f= @(const, x_data) comp_hill(const, x_data);
[coef,res7] = lsqcurvefit(f,c07,x_data,y_data,[zeros(1,17)],[100*ones(1,17)]);

[~,cell_sort]=sort(y_data,'descend');
% hold on
% scatter(1:length(y_data),y_data(cell_sort),150,'MarkerFaceColor',[255,255,255]*0.4/255,'MarkerEdgeColor',[1,1,1])
% scatter(1:length(y_data),comp_hill(coef_list(240,:),x_data(cell_sort,:),para_list(240,:)),150,'MarkerFaceColor',[255,45,85]/255,'MarkerEdgeColor',[1,1,1])


 %%
%%


figure(3)
hold on


xx=(0:0.1:25)';
xxx=[zeros(length(xx),4),xx];

yy=comp_hill([coef(1:15),1,1], xxx);

yy=yy*coef(16);

scatter(sum(mfi_all(2:13,1:4),2),pNK_05_1(2:13),50,'MarkerFaceColor',[218,41,28]/255,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.8);

plot(1+xx*x_slope+x_intercept, yy,'Color',[75,75,75,100]/255,'LineWidth',1)
plot(1+(-5:0.1:0)*x_slope+x_intercept, coef(16)*coef(10)*ones(size((-5:0.1:0))),'Color',[75,75,75,100]/255,'LineWidth',1)


xlim([1,10])
xticks(1:10)
xticklabels(1:10)
xtl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xtl,'fontsize',8,'FontWeight','bold')
xlabel('B7H6 level')

ylim([0,100])
yticks(0:50:100)
ytl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',ytl,'fontsize',8,'FontWeight','bold')
ylabel('Specific cell lysis(%)')


set(gca, 'LineWidth',1)
set(gca,'TickDir','out');
% title('Leukemia','FontSize',20)
set(gcf, 'Position',[0,0,200,160])
[h1,icons]=legend({'Experiment','Fitted'},'fontsize',8);
% icons2=findobj(icons,'type','Group');
icons(3).Children.MarkerSize=6;
icons(3).Children.MarkerEdgeColor=[1,1,1];
% icons(5).Children.FaceAlpha=0.8;


set(h1, 'position',[0.35,0.7,0.2,0.2])
legend boxoff



function y=comp_hill(const, x_data) 
[m,~]=size(x_data);
cc=const(1);
kk=const(6:10);
c0=const(10);
K=kk(1);
y=zeros(m,1); 
h=2;
hh=const(11:15);


y=y+cc*x_data.^(h)./(x_data.^(h)+K.^(h));

y=y+c0;
y(13:24)=const(16)*y(13:24);
y(25:36)=const(17)*y(25:36);
% y=min(100,y);
% y=max(0,y);
end
