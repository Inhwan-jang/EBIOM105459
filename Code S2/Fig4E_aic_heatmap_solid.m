%% Code for Fig. 4E : Run Fig4E_aic_solid.m before this

load('colormaps')
load('nk_data_230622.mat')

close all
nrow=32;

A=zeros(nrow,5);
solid_ind=14:35;
max_mfi=10;
x_data=mfi_all(solid_ind,:)-1;

for row_i=1:nrow
    for col_i=1:5
        temp_max=10;
        temp_val=para_list2(row_i,col_i);
        if temp_val==1
            A(row_i,col_i)=1;
            
            temp_val=coef_list(b(row_i),col_i)*temp_max.^2./(temp_max.^2+coef_list(b(row_i),6)^2);
            temp_alpha=temp_val/70;
        else
            A(row_i,col_i)=0;
            temp_alpha=0;
        end
       
        A(row_i,col_i)=temp_alpha*A(row_i,col_i);
    end
end
figure(215)
imagesc(A(:,[1:5]))
colormap(colormap_rwb2(256:-1:1,:))

yticks([])

caxis([-1,1])
xticks([1:5])

set(gca,'FontSize',8,'FontWeight','bold')
set(gca,'xcolor',[0,61,165]/255)
set(gca,'ycolor',[0,61,165]/255)
ylabel('model rank','Color',[0,0,0])
xticklabels(["\color{black} MIC","\color{black} BP1","\color{black} BP2","\color{black} BP3","\color{black} B7"]);

set(gca,'Linewidth',2)
set(gca,'xaxisLocation','top')

set(gcf,'Position',[100,100,150,200])