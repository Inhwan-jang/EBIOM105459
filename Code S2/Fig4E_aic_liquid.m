%% Code for Fig. 4C and 4E solid: Run Fig.4C heatmap after this

% Last updated: 240502
close 
clear
load('nk_data_230622.mat')
leuk_ind=2:13;
x_data=mfi_all(leuk_ind,:)-1;
x_data(x_data<0)=0;

b7_intercept=min(x_data(:,5));
b7_slope=(max(x_data(:,5))-min(x_data(:,5)))/10;
%% Data normalization
x_data=[x_data;x_data;x_data];
x_data=10*(x_data-min(x_data))./(max(x_data)-min(x_data));

y_data=[pNK_025_1(leuk_ind,1);pNK_05_1(leuk_ind,1);pNK_1_1(leuk_ind,1)];
%%

% leuk_names=["K562","Raji","Ramos", "Jurkat", "U937", "CCRF-CEM", "MOLT-4", "THP-1","Kasumi-1","Kasumi-3", "KG-1", "MV4-11","HL-60"];
leuk_names=["Raji","Ramos", "Jurkat", "U937", "CCRF-CEM", "MOLT-4", "THP-1","Kasumi-1","Kasumi-3", "KG-1", "MV4-11","HL-60"];

cost_list=zeros(32,1);

c07=ones(1,17);

para_list=dec2base(0:31,2)-'0';
coef_list=zeros(32,17);

for ind=1:32
    para=para_list(ind,:);
    f= @(const, x_data) comp_hill(const, x_data,para);
    [yuri,res7] = lsqcurvefit(f,c07,x_data,y_data,[zeros(1,17)],[100*ones(1,17)]);
    cost_list(ind )=res7;
    coef_list(ind,:)=yuri;
end
[~,cell_sort]=sort(y_data,'descend');
% hold on
% scatter(1:length(y_data),y_data(cell_sort),150,'MarkerFaceColor',[255,255,255]*0.4/255,'MarkerEdgeColor',[1,1,1])
% scatter(1:length(y_data),comp_hill(coef_list(240,:),x_data(cell_sort,:),para_list(240,:)),150,'MarkerFaceColor',[255,45,85]/255,'MarkerEdgeColor',[1,1,1])


% aic0=13*log(cost_list/13)+2*(10-2*sum(para_list==2,2));
v=sum(para_list==1,2)+4;
v(1)=3;

w=length(y_data);
aic0=w*log(cost_list/w)+2*v+2*v.*(v+1)./(w-v-1);

[aic,b]=sort(aic0);

%%
aic_diff=aic-aic(1);

aic_w=exp(-0.5*aic_diff)/sum(exp(-0.5*aic_diff));
% histogram(aic_w)
aic_w_sum=cumsum(aic_w);
para_list2=para_list(b,:);
para_list_neg=para_list2==0;
para_list_none=para_list2==2;
para_list_acti=para_list2==1;


 %%

num_models=(sum(aic_w_sum<0.95))+1;
figure(2)
para_list3=1.3*(para_list_neg*2+para_list_none); % order change

ligand_list=["MICAB","ULBP1","ULBP2","ULBP3","B7H6"];
for subplot_i=1:5
    subplot(2,3,subplot_i)
    hold on
    histogram(para_list3(1:num_models,subplot_i),0:1,'FaceColor',[0,61,165]/255,'EdgeAlpha',0.6,'LineWidth',2, 'EdgeColor',[1,1,1]*0.)    
    histogram(para_list3(1:num_models,subplot_i),1.1:2.1,'FaceColor',[150,150,150]/255,'EdgeAlpha',0.6,'LineWidth',2, 'EdgeColor',[1,1,1]*0.)        
    histogram(para_list3(1:num_models,subplot_i),2.2:3.2,'FaceColor',[216,41,28]/255,'EdgeAlpha',0.6,'LineWidth',2, 'EdgeColor',[1,1,1]*0.)        

   
    xlim([-0.2,3.3])
    xticks([0.5,1.6,2.7])
    xticklabels({'Acti.','None','Inhi.'})
    xtickangle(45)
    ylim([0,num_models])
    yticks([0,num_models])
    
    ylabel('Number of models')
    
    xtl = get(gca,'XTickLabel');  
    set(gca, 'LineWidth',1.5)
    set(gca,'XTickLabel',xtl,'fontsize',15,'FontWeight','bold')
    title(ligand_list(subplot_i),'FontSize',15,'FontWeight','bold') 
end
set(gcf, 'position',[0,0,1200,800])
%%


figure(3)
hold on


xx=(0:0.1:10)';
xxx=[zeros(length(xx),4),xx];

yy=comp_hill([coef_list(2,1:15),1,1], xxx, para_list(2,:));

yy=yy*coef_list(2,16);

scatter(mfi_all(2:13,5),pNK_05_1(2:13),50,'MarkerFaceColor',[218,41,28]/255,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.8);

plot(1+xx*b7_slope+b7_intercept, yy,'Color',[75,75,75,100]/255,'LineWidth',1)
plot(1+(-5:0.1:0)*b7_slope+b7_intercept, coef_list(2,16)*coef_list(2,10)*ones(size((-5:0.1:0))),'Color',[75,75,75,100]/255,'LineWidth',1)


xlim([1,5])
xticks(1:5)
xticklabels(1:5)
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



function y=comp_hill(const, x_data, para) 
[m,~]=size(x_data);
cc=const(1:5);
kk=const(6:10);
c0=const(10);
K=kk(1);
y=zeros(m,1); 
h=2;
hh=const(11:15);
% for i=1:5
%     if para(i)==1
%             y=y+cc(i)*x_data(:,i).^(hh(i))./(x_data(:,i).^(hh(i))+kk(i).^(hh(i)));
% %             y=y+cc(i)*x_data(:,i);
%     elseif para(i)==0
%             y=y+cc(i)*kk(:,i).^(hh(i))./(x_data(:,i).^(hh(i))+kk(i).^(hh(i))) ;
%     end
% end
for i=1:5
    if para(i)==1
            y=y+cc(i)*x_data(:,i).^(h)./(x_data(:,i).^(h)+K.^(h));
%             y=y+cc(i)*x_data(:,i);
%     elseif para(i)==0
%             y=y+cc(i)*kk(:,i).^(h)./(x_data(:,i).^(h)+kk(i).^(h)) ;
    end
end
y=y+c0;
y(13:24)=const(16)*y(13:24);
y(25:36)=const(17)*y(25:36);
% y=min(100,y);
% y=max(0,y);
end
