[actualTime, frameTime, PE, KE, TE] = textread('breathing.txt', '%f %f %f %f %f');

color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.9290, 0.6940, 0.1250];

figure
hold on
plot(frameTime, PE, 'color', color1)
plot(frameTime, KE, 'color', color2)
plot(frameTime, TE, 'color', color3)
legend('Potential Energy', 'Kinetic Energy', 'Total Energy')
xlabel("Frame Time (s)")
ylabel("Energy (J)")
title('Breathing Cube')
saveas(gcf,'BreathingCubeEnergy.png')