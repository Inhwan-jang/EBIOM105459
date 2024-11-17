clear
% close all

load('draw_beta_liquid_230622.mat')
% load('draw_beta_20220415.mat')

%%
figure
subplot(3,2,1)
histogram(draw_beta(2e5:1e6,1),0:10:100,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])
    
subplot(3,2,2)
histogram(draw_beta(2e5:1e6,2).*draw_gamma(2e5:1e6,2),0:10:100,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,3)
histogram(draw_beta(2e5:1e6,3).*draw_gamma(2e5:1e6,3),0:10:100,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,4)
histogram(draw_beta(2e5:1e6,4).*draw_gamma(2e5:1e6,4),0:10:100,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,5)
histogram(draw_beta(2e5:1e6,5).*draw_gamma(2e5:1e6,5),0:10:100,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,6)
histogram(draw_beta(2e5:1e6,6).*draw_gamma(2e5:1e6,6),0:20:200,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,100,200])
set(gcf,'Position',[100,100,1000,400])
%%
subplot(3,2,1)

histogram(draw_gamma(2e5:1e6,2),0:0.5:1,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,2)
histogram(draw_gamma(2e5:1e6,3),0:0.5:1,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,3)
histogram(draw_gamma(2e5:1e6,4),0:0.5:1,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,4)
histogram(draw_gamma(2e5:1e6,5),0:0.5:1,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,50,100])

subplot(3,2,5)
histogram(draw_gamma(2e5:1e6,6),0:0.5:1,'Normalization','probability', 'FaceColor',[255,55,64]/255,'FaceAlpha',0.3)
ylim([0,1])
yticks([0,1])
xticks([0,100,200])
set(gcf,'Position',[100,100,1000,400])

%%
figure
subplot(3,2,1)
plot(draw_beta(1:1000:1e6,1),'Color',[250,100,100,80]/255,'LineWidth',0.75)
yticks([0,50,100])
ylim([0,100])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])


subplot(3,2,2)
plot(draw_beta(1:1000:1e6,2).*draw_gamma(1:1000:1e6,2),'Color',[250,100,100,80]/255,'LineWidth',0.75)
yticks([0,200,400])
ylim([0,400])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])

subplot(3,2,3)
plot(draw_beta(1:1000:1e6,3).*draw_gamma(1:1000:1e6,3),'Color',[250,100,100,80]/255,'LineWidth',0.75)
yticks([0,100,200])
ylim([0,200])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])


subplot(3,2,4)
plot(draw_beta(1:1000:1e6,4).*draw_gamma(1:1000:1e6,4),'Color',[250,100,100,80]/255,'LineWidth',0.75)
yticks([0,100,200])
ylim([0,200])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])


subplot(3,2,5)
plot(draw_beta(1:1000:1e6,5).*draw_gamma(1:1000:1e6,5),'Color',[250,100,100,80]/255,'LineWidth',0.75)
yticks([0,100,200])
ylim([0,200])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])


subplot(3,2,6)
plot(draw_beta(1:1000:1e6,6).*draw_gamma(1:1000:1e6,6),'Color',[250,100,100,80]/255,'LineWidth',0.75)
yticks([0,200,400])
ylim([0,400])
xticks([0,1e3])
xlim([0,1e3])

xticklabels(["0","1 \times 10 ^6"])

set(gcf,'Position',[100,100,1000,400])

%%
figure
subplot(2,2,1)
histogram(draw_c_vec(2e5:1e6,1),0:0.5:5,'Normalization','probability', 'FaceColor',[51,171,95]/255,'FaceAlpha',0.3)
ylim([0,0.5])
yticks([0,0.5])

subplot(2,2,2)
histogram(draw_c_vec(2e5:1e6,2),0:0.5:5,'Normalization','probability', 'FaceColor',[51,171,95]/255,'FaceAlpha',0.3)
ylim([0,0.5])
yticks([0,0.5])

subplot(2,2,3)
histogram(draw_sigma(2e5:1e6),0:4:40,'Normalization','probability', 'FaceColor',[51,171,95]/255,'FaceAlpha',0.3)
ylim([0,0.5])
yticks([0,0.5])

subplot(2,2,4)
histogram(draw_K(2e5:1e6),0:10:100,'Normalization','probability', 'FaceColor',[51,171,95]/255,'FaceAlpha',0.3)
ylim([0,0.5])
yticks([0,0.5])

%%
figure
subplot(2,2,1)
plot(draw_c_vec(1:1000:1e6,1),'Color',[51,171,95,100]/255)
ylim([0,10])
yticks([0,5,10])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])

subplot(2,2,2)
plot(draw_c_vec(1:1000:1e6,2),'Color',[51,171,95,100]/255)
ylim([0,10])
yticks([0,5,10])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])

subplot(2,2,3)
plot(draw_sigma(1:1000:1e6),'Color',[51,171,95,100]/255)
ylim([0,40])
yticks([0,40])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])

subplot(2,2,4)
plot(draw_K(1:1000:1e6),'Color',[51,171,95,100]/255)
ylim([0,200])
yticks([0,200])
xticks([0,1e3])
xlim([0,1e3])
xticklabels(["0","1 \times 10 ^6"])
%%
figure
start_ind=2e5;
end_ind=1e6;
bar(mean(draw_gamma(start_ind:end_ind,[2,3,4,5,6]),1))
yticks([0,0.5,1])
ylim([0,1])
set(gca,'TickDir','out')