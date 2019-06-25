data = textread('populationDistance.txt', '', 'delimiter', ' ');
data = data(:,1:end-1);
bestDis = max(data,[],2);
%% learning curve
figure
plot(bestDis)
%xlim([0 4000])
xlabel('Generation')
ylabel('Distance (m)')
title('Learning Curve')
% %% Dot plot
% figure
% plot(data,'.k')
% %xlim([0 2606])
% xlabel('Generation')
% ylabel('Distance (m)')
% title('Dot Plot')
% %% Diversity chart
% figure
% uni = zeros(209,1);
% for i = 1:209
%     uni(i) = 100 * length(unique(data(i,:)))/128;
% end
% plot(uni)
% xlabel('Generation')
% ylabel('Diversity (%)')
% title('Diversity Plot')