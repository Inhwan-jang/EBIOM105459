%% Code for Fig. 4H : Run Fig.4E_aic_liquid before this

figure(30)
hold on
% b(1)=18;
% load('nk_data_230622')
load('validation_cell_line_ligand_v2.mat')
validation_cell_line_b7h6=validation_cell_line_ligand(:,5);

validation_ind=[1:5];
xx=(0:0.1:5)';
xxx=[zeros(length(xx),4),xx];

% normalizing B7H6 
b7_intercept=min(x_data(:,5));
b7_slope=(max(x_data(:,5))-min(x_data(:,5)))/10;

% Predicting cytotoxicity at E:T=1:1
prediction=comp_hill([coef_list(2,1:15),1,1],[zeros(5,4), (validation_cell_line_b7h6-b7_intercept)/b7_slope], para_list(2,:))*coef_list(2,17);
experiment_cyto=validation_cell_line_cyto_1_1(validation_ind);
scatter(prediction(validation_ind),experiment_cyto,50,'MarkerFaceColor',[255,255,255]/255,'MarkerEdgeAlpha',1,'MarkerFaceAlpha',0,'MarkerEdgeColor',[0,0,0]/255);


plot([0,100], [0,100],'Color',[75,75,75,100]/255,'LineWidth',1)


xlim([0,100])
xticks(0:50:100)
xticklabels(0:50:100)
xtl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xtl,'fontsize',8,'FontWeight','bold')
xlabel('Prediction')

ylim([0,100])
yticks(0:50:100)
ytl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',ytl,'fontsize',8,'FontWeight','bold')
ylabel('Specific cell lysis(%)')

set(gca, 'LineWidth',1)
set(gca,'TickDir','out');
% title('Leukemia','FontSize',20)
set(gcf, 'Position',[0,0,200,160])

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

end
