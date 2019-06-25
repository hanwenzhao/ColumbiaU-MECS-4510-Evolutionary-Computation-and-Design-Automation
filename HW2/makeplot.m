%% Read File
[eval1,~,RM_1] = textread('Random_1.txt', '%d %f %f');
[~,~,RM_2] = textread('Random_2.txt', '%d %f %f');
[~,~,RM_3] = textread('Random_3.txt', '%d %f %f');
[~,~,RM_4] = textread('Random_4.txt', '%d %f %f');

[~,~,RMHC_1] = textread('RMHC_1.txt', '%d %f %f');
[~,~,RMHC_2] = textread('RMHC_2.txt', '%d %f %f');
[~,~,RMHC_3] = textread('RMHC_3.txt', '%d %f %f');
[~,~,RMHC_4] = textread('RMHC_4.txt', '%d %f %f');

[eval2,~,GP1_1] = textread('GP1_1.txt', '%d %f %f');
[~,~,GP1_2] = textread('GP1_2.txt', '%d %f %f');
[~,~,GP1_3] = textread('GP1_3.txt', '%d %f %f');
[~,~,GP1_4] = textread('GP1_4.txt', '%d %f %f');

[~,~,GP2_1] = textread('GP2_1.txt', '%d %f %f');
[~,~,GP2_2] = textread('GP2_2.txt', '%d %f %f');
[~,~,GP2_3] = textread('GP2_3.txt', '%d %f %f');
[~,~,GP2_4] = textread('GP2_4.txt', '%d %f %f');

[eval3,~,GP2_LP_1] = textread('GP2_LP_1.txt', '%d %f %f');
[~,~,GP2_LP_2] = textread('GP2_LP_2.txt', '%d %f %f');
[~,~,GP2_LP_3] = textread('GP2_LP_3.txt', '%d %f %f');
[~,~,GP2_LP_4] = textread('GP2_LP_4.txt', '%d %f %f');

%% Calculation
% RM
RM = mean([RM_1,RM_2,RM_3,RM_4],2);
% define x coordinates for error bar
ErrorX = 0:60000:300000;
ErrorX = ErrorX(2:end);
% calcuate y coordinates for error bar
RM_ErrorY = RM(ErrorX);
% calculate error
RM_Error = [RM_1(ErrorX),RM_2(ErrorX),RM_3(ErrorX),RM_4(ErrorX)];
RM_Error = std(RM_Error,0,2)/2;

% RMHC
RMHC = mean([RMHC_1,RMHC_2,RMHC_3,RMHC_4],2);
% calcuate y coordinates for error bar
RMHC_ErrorY = RMHC(ErrorX);
% calculate error
RMHC_Error = [RMHC_1(ErrorX),RMHC_2(ErrorX),RMHC_3(ErrorX),RMHC_4(ErrorX)];
RMHC_Error = std(RMHC_Error,0,2)/2;

% GP1
GP1 = mean([GP1_1,GP1_2,GP1_3,GP1_4],2);
% define x coordinates for error bar
ErrorX2 = ErrorX/100;
% calcuate y coordinates for error bar
GP1_ErrorY = GP1(ErrorX2);
% calculate error
GP1_Error = [GP1_1(ErrorX2),GP1_2(ErrorX2),GP1_3(ErrorX2),GP1_4(ErrorX2)];
GP1_Error = std(GP1_Error,0,2)/2;

% GP2
GP2 = mean([GP2_1,GP2_2,GP2_3,GP2_4],2);
% calcuate y coordinates for error bar
GP2_ErrorY = GP2(ErrorX2);
% calculate error
GP2_Error = [GP2_1(ErrorX2),GP2_2(ErrorX2),GP2_3(ErrorX2),GP2_4(ErrorX2)];
GP2_Error = std(GP2_Error,0,2)/2;

% GP2_LP
GP2_LP = mean([GP2_LP_1,GP2_LP_2,GP2_LP_3,GP2_LP_4],2);
% calcuate y coordinates for error bar
ErrorX3 = ErrorX/200;
GP2_LP_ErrorY = GP2_LP(ErrorX3);
% calculate error
GP2_LP_Error = [GP2_LP_1(ErrorX3),GP2_LP_2(ErrorX3),GP2_LP_3(ErrorX3),GP2_LP_4(ErrorX3)];
GP2_LP_Error = std(GP2_LP_Error,0,2)/2;

%%
%%%%%%%%%%%%%%%define color############
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.9290, 0.6940, 0.1250];
color4 = [0.4940, 0.1840, 0.5560];
color5 = [0.4660, 0.6740, 0.1880];

figure;hold on
plot(eval1,RM,'color',color1,'LineWidth',1.5)
plot(eval1,RMHC,'color',color2,'LineWidth',1.5)
plot(eval2,GP1,'color',color3,'LineWidth',1.5)
plot(eval2,GP2,'color',color4,'LineWidth',1.5)
plot(eval3,GP2_LP,'color',color5,'LineWidth',1.5)


errorbar(ErrorX,RM_ErrorY,RM_Error,'.','color',color1)
errorbar(ErrorX,RMHC_ErrorY,RMHC_Error,'.','color',color2)
errorbar(ErrorX2*100,GP1_ErrorY,GP1_Error,'.','color',color3)
errorbar(ErrorX2*100,GP2_ErrorY,GP2_Error,'.','color',color4)
errorbar(ErrorX3*200,GP2_LP_ErrorY,GP2_LP_Error,'.','color',color5)
set(gca,'YScale','log');
xlabel('Evaluations')
ylabel('Mean Absolute Error')
xlim([0 300000])
legend('Random Search','Hill Climber','GP (Deterministic Crowding)','GP (Convention Selection)', 'GP (Convention Selection with Large Population)','Location','best')
saveas(gcf,'performancePlot.png')