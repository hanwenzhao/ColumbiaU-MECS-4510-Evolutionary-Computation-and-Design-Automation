%% Dot Plot
data = dotData;
eval = (1:3000);
figure
hold on
for i = 1:3000
   plot(eval(i)*ones(1,200),data(i,:),'.'); 
end
set(gca,'YScale','log');
xlabel('Generations')
ylabel('Population Errors')
title('Dot Plot for GP(Determinstic Crowding)')
saveas(gcf,'dotPlot.png')

%% Diversity Plot
figure
hold on
diversity = zeros(1,3000);
for i = 1:3000
   diversity(i) = length(unique(data(i,:)))*100/200; 
end
plot(1:3000,diversity)
xlabel('Generations')
ylabel('Diversity (%)')
title('Diversity Plot for GP(Determinstic Crowding)')
saveas(gcf,'diversityPlot.png')

%% Convergence Plot
figure
hold on
convergence = zeros(1,3000);
threshold = 10^-2;
for i = 1:3000
   convergence(i) = length(find(data(i,:)<threshold))*100/200; 
end
plot(1:3000,convergence)
xlabel('Generations')
ylabel('Convergence (%)')
title('Convergence Plot for GP(Determinstic Crowding) Threshold = 0.01')
saveas(gcf,'convergencePlot.png')

%% Validation
figure
hold on
plot(1:1000,errorAndValidation(:,1),'LineWidth',2)
plot(1:1000,errorAndValidation(:,1))
set(gca,'YScale','log');
xlabel('Generations')
ylabel('Mean Absolute Error')
legend('Trainning Data','Testing Data')
title('Validation Curve for GP(Determinstic Crowding)')
saveas(gcf,'validationPlot.png')
