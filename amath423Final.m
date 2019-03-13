% AMATH 423
% Levi Davis and Oliver Speltz
close all; clc;

% Default parameter values
j1 = 80; %80; default
j2 = 0.0005;
j3 = 0.5;
j4 = 0.16;
j5 = 0.05;
ta1 = 1500;
ta2 = 0.5;
ta3 = 60;
ta4 = 100;
ta5 = 0;

stimFun = @(t) (t>=0).*125 - 125.*(t>30);

initialGuess = [0.003, 0.003, 0.003];

[T,Y] = ode45(@(t,y) neuronFireODE(t,y,stimFun,j1,j2,j3,j4,j5,ta1,ta2,...
    ta3,ta4,ta5), [0,1000], initialGuess);

plot(T, Y(:,1))

% 
% % make sure to set this to the parameter being varied. For plotting
% % purposes
% varying = j1;
% 
% % compute values and determine stability
% [pkm, actin, rna, hs, stability, indexes] = computeSS(j1,j2,j3,j4,j5);
% 
% stable = [];
% stableIndex = [];
% unstable = [];
% unstableIndex = [];
% for j = 1:length(pkm)
%     for i = 1:length(pkm{j})
%         if stability{j}(i) > 0
%             stable = [stable, pkm{j}(i)];
%             stableIndex = [stableIndex, varying(indexes{j}(i))];
%         else
%             unstable = [unstable, pkm{j}(i)];
%             unstableIndex = [unstableIndex, varying(indexes{j}(i))];
%         end
%     end
% end
% 
% 
% figure(1)
% scatter(stableIndex, stable, '.'); hold on;
% scatter(unstableIndex, unstable, '.');
% legend('Stable', 'Unstable');
% ylabel('[PKM]'); xlabel('j1');
% 
