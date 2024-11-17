clear; close all;
load("B7_MIC_OE_240321.mat")
coef_0=rand(1,5);
g= @(const, x_data) cyto_fitting(const, x_data);


B7H6_dat=[B7H6_OE];
MICAB_dat=[MICAB_OE];
[coef_b7,res_b7] = lsqcurvefit(g,coef_0,B7H6_dat(:,1),B7H6_dat(:,2),[zeros(1,5)],[100,100,100,100,100]);

[coef_mic,res_mic] = lsqcurvefit(g,coef_0,MICAB_dat(:,1),MICAB_dat(:,2),[zeros(1,5)],[100,100,100,100,100]);



c_b7=coef_b7(1);
k_b7=coef_b7(2);
h_b7=coef_b7(3);
base_b7=coef_b7(5);


c_mic=coef_mic(1);
k_mic=coef_mic(2);
h_mic=coef_mic(3);
base_mic=coef_mic(5);

% h_b=5;

X=1:0.01:5;

y_b7= c_b7* (base_b7+((X-1).^h_b7)./((X-1).^h_b7+k_b7.^h_b7));
y_b7=min(100,y_b7);

y_mic= c_mic* (base_mic+((X-1).^h_mic)./((X-1).^h_mic+k_mic.^h_mic));
y_mic=min(100,y_mic);
hold on
% scatter(log10(X_self(1:31)), Y_self(1:31))
plot(X,y_b7,'Color',[0,122,255]/255)
plot(X,y_mic,'Color',[255,149,0]/255)

scatter(B7H6_dat(:,1),B7H6_dat(:,2),'MarkerEdgeAlpha',0,'MarkerFaceColor',[0,122,255]/255)
scatter(MICAB_dat(:,1),MICAB_dat(:,2),'MarkerEdgeAlpha',0,'MarkerFaceColor',[255,149,0]/255)

xlabel("MFI")
ylabel("cyto (0.5:1)")
xticks([1,5])
xlim([1,5])
set(gca,'TickDir','out')
% xlim([0,1])
ylim([0,30])
yticks([0:15:30])

function y=cyto_fitting(const, x_data) 
c=const(1);
k_b0=const(2);
h_b=const(3);
% scaling1=const(4);
base1=const(5);
b7=x_data(:,1)-1; b7=max(b7,0);
% hla=x_data(:,2)-1; hla=max(hla,0);
% k_b=k_b0*(hla.^h_h)./(hla.^h_h+k_h^h_h);
k_b=k_b0;
y= c* (base1+(b7.^h_b)./(b7.^h_b+k_b.^h_b));

% y(1:32)=scaling1*y(1:32);
% y(50:end)=scaling2*y(50:end);

end