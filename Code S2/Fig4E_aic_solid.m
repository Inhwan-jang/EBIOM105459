%% Code for Fig. 4C and 4E solid: Run Fig.4C heatmap after this

% Last updated: 240502

close 
clear
load('nk_data_230622.mat')

solid_ind=[15:35];
x_data=mfi_all(solid_ind,:)-1;
x_data(x_data<0)=0;

b7_intercept=min(x_data(:,5));
b7_slope=(max(x_data(:,5))-min(x_data(:,5)))/10;
x_data=10*(x_data-min(x_data))./(max(x_data)-min(x_data));

x_data=[x_data;x_data];

red=[219,41,28]/255;
blue=[0,61,165]/255;

y_data=[pNK_05_1(solid_ind,1);pNK_2_1(solid_ind,1)];

%%
% y_data=y_data(2:13);
% x_data=mfi_all(logical(1-leuk_ind),:);
% x_data=x_data(:,1:5);
% y_data=cyto_pnk(logical(1-leuk_ind),1);

cost_list=zeros(32,1);

c07=ones(1,16);

para_list=dec2base(0:31,2)-'0';
coef_list=zeros(32,16);

for ind=1:32
    para=para_list(ind,:);
    f= @(const, x_data) comp_hill(const, x_data,para);
    [yuri,res7] = lsqcurvefit(f,c07,x_data,y_data,[zeros(1,16)],[100*ones(1,16)]);
    cost_list(ind )=res7;
    coef_list(ind,:)=yuri;
end
[~,cell_sort]=sort(y_data,'descend');


v=sum(para_list==1,2)+3;
v(1)=2;
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


% ligand_list=["MICAB","ULBP1","ULBP2","ULBP3","B7H6"];

%%

figure(3)
hold on



xx=(0:0.1:10)';
xxx=[zeros(length(xx),4),xx];
% model 2 (predicting cytotoxicity with using only B7H6) 

yy=comp_hill([coef_list(2,1:10),1,1], xxx, para_list(2,:));

yy=yy;

scatter(mfi_all(solid_ind,5),pNK_05_1(solid_ind),50,'^','MarkerFaceColor',blue,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.8);

plot(1+xx*b7_slope+b7_intercept, yy,'Color',[75,75,75,100]/255,'LineWidth',1)
plot(1+(-5:0.1:0)*b7_slope+b7_intercept, coef_list(2,11)*coef_list(2,10)*ones(size((-5:0.1:0))),'Color',[75,75,75,100]/255,'LineWidth',1)



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
[h1,icons]=legend({'Experiment','Fitted'},'FontSize',8);
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

for i=1:5
    if para(i)==1
            y=y+cc(i)*x_data(:,i).^(h)./(x_data(:,i).^(h)+K.^(h));

    end
end
y=y+c0;
y(22:42)=const(11)*y(22:42);

end