
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')
%% basic settings

calciumDir = [baseDir, 'Calcium\'];

colors = [248 171 66; 40 170 225] / 255;







%% Figure 7A
tuneFn = 'angle_tuning_summary_preAnswer_perTouch_NC_PTC';
tune = load([calciumDir, tuneFn]);
learnerInd = [1:4,7,9];
tune.naive = tune.naive(learnerInd);
numMice = length(learnerInd);

jk027length = length(tune.naive(2).touchID);
jk027inds = find(tune.naive(2).touchID > 5000);
fn = fieldnames(tune.naive);
for i = 1 : length(fn)
    if size(tune.naive(2).(fn{i}), 1) == jk027length
        tune.naive(2).(fn{i}) = tune.naive(2).(fn{i})(jk027inds,:);
    elseif size(tune.naive(2).(fn{i}), 2) == jk027length
        tune.naive(2).(fn{i}) = tune.naive(2).(fn{i})(:,jk027inds);
    end
end

maxResponseTuned = cell(numMice,2);
meanOtherResponseTuned = cell(numMice,2);
angleSelectivityTuned = cell(numMice,2);
angleSelectivityNormTuned = cell(numMice,2);
for mi = 1 : numMice
    maxResponseTuned{mi,1} = cellfun(@(x) max(cellfun(@mean, x)), tune.naive(mi).val);
    maxResponseTuned{mi,2} = cellfun(@(x) max(cellfun(@mean, x)), tune.expert(mi).val);
    
    meanOtherResponseTuned{mi,1} = cellfun(@(x) (sum(cellfun(@mean, x)) - max(cellfun(@mean, x)))/6, tune.naive(mi).val);
    meanOtherResponseTuned{mi,2} = cellfun(@(x) (sum(cellfun(@mean, x)) - max(cellfun(@mean, x)))/6, tune.expert(mi).val);
    
    angleSelectivityTuned{mi,1} = maxResponseTuned{mi,1} - meanOtherResponseTuned{mi,1};
    angleSelectivityTuned{mi,2} = maxResponseTuned{mi,2} - meanOtherResponseTuned{mi,2};
    
    angleSelectivityNormTuned{mi,1} = (maxResponseTuned{mi,1} - meanOtherResponseTuned{mi,1}) ./ maxResponseTuned{mi,1};
    angleSelectivityNormTuned{mi,2} = (maxResponseTuned{mi,2} - meanOtherResponseTuned{mi,2}) ./ maxResponseTuned{mi,2};
end

%%
figure, hold on
tempMat = cellfun(@mean, angleSelectivityTuned);
[~,p] = ttest(tempMat(:,1), tempMat(:,2));
errorbar([1,2], mean(tempMat), sem(tempMat), 'k', 'lines', 'no')
for i = 1 : 2
    bar(i, mean(tempMat(:,i)), 'facecolor', colors(i,:))
end

for mi = 1 : numMice
    plot([1,2], tempMat(mi,:), 'ko-')
end
errorbar([1,2], mean(tempMat), sem(tempMat), 'k', 'linewidth', 2, 'lines', 'no')

xlim([0.5 2.5]), xticks([1, 2])
xticklabels({'Naive', 'Expert'}), xtickangle(45)
ylim([0 0.6]), yticks(0:0.2:0.6)
ylabel('Mean angle selectivity')
title(sprintf('p = %f',p))






%%
%% 7B
%% Not-selective neurons among touch neurons

tuneFn = 'angle_tuning_summary_preAnswer_perTouch_NC_PTC';
tune = load([calciumDir, tuneFn]);
learnerInd = [1:4,7,9];
tune.naive = tune.naive(learnerInd);
numMice = length(learnerInd);
jk027length = length(tune.naive(2).touchID);
jk027inds = find(tune.naive(2).touchID > 5000);
fn = fieldnames(tune.naive);
for i = 1 : length(fn)
    if size(tune.naive(2).(fn{i}), 1) == jk027length
        tune.naive(2).(fn{i}) = tune.naive(2).(fn{i})(jk027inds,:);
    elseif size(tune.naive(2).(fn{i}), 2) == jk027length
        tune.naive(2).(fn{i}) = tune.naive(2).(fn{i})(:,jk027inds);
    end
end


propNsTuned = zeros(numMice,2); % 1: naive 2: expert
for mi = 1 : numMice
    if mi == 2
        tempInd = find(tune.naive(mi).touchID > 5000);
        propNsTuned(mi,1) = 1 - sum(tune.naive(mi).tuned(tempInd)) / length(tempInd);
    else
        propNsTuned(mi,1) = 1 - sum(tune.naive(mi).tuned) / length(tune.naive(mi).tuned);
    end
    propNsTuned(mi,2) = 1 - sum(tune.expert(mi).tuned) / length(tune.expert(mi).tuned);
end

figure, hold on
for i = 1 : size(propNsTuned,1)
    plot(propNsTuned(i,:), 'o-', 'color', [0.6 0.6 0.6])
end

[~,p] = ttest(propNsTuned(:,1) - propNsTuned(:,2));

errorbar(mean(propNsTuned), sem(propNsTuned), 'ro', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
xlim([0.5 2.5])
ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel('Proportion (/touch cells)')
title({'Not-selective'; sprintf('p = %f', p)})











%%
%% 7C
%%
% Quantification of angle selectivity increase
% compared across preferred angles


angles = 45:15:135;
atShapes = cell(numMice,length(angles),2);
for mi = 1 : numMice
    % naive
    indTunedAngle = tune.naive(mi).tunedAngle;
    
    for ai = 1 : length(angles)
        indAngle = find(indTunedAngle == angles(ai));
        atShapes{mi,ai,1} = cell2mat(cellfun(@(x) cell2mat(cellfun(@mean, x', 'un', 0)),tune.naive(mi).val(indAngle), 'un', 0));
    end
    
    % epxert
    indTunedAngle = tune.expert(mi).tunedAngle;
    
    for ai = 1 : length(angles)
        indAngle = find(indTunedAngle == angles(ai));
        atShapes{mi,ai,2} = cell2mat(cellfun(@(x) cell2mat(cellfun(@mean, x', 'un', 0)),tune.expert(mi).val(indAngle), 'un', 0));
    end
    
end


asNaive = cell(1,length(angles));
asExpert = cell(1,length(angles));
for ai = 1 : length(angles)
    asNaive{ai} = cell2mat(cellfun(@(x) max(x,[],2) - (sum(x,2) - max(x,[],2))/6, atShapes(:,ai,1), 'un', 0));
    asExpert{ai} = cell2mat(cellfun(@(x) max(x,[],2) - (sum(x,2) - max(x,[],2))/6, atShapes(:,ai,2), 'un', 0));
end

%
barOffset = 0.2*15;
barWidth = 0.4*15;
p = zeros(1,length(angles));
figure, hold on
for ai = 1 : length(angles)
    bar(angles(ai)-barOffset, mean(asNaive{ai}), barWidth, 'facecolor', colorsTransient(1,:))
    bar(angles(ai)+barOffset, mean(asExpert{ai}), barWidth, 'facecolor', colorsTransient(2,:))
    if ai == 1
        legend({'Naive', 'Expert'}, 'box', 'off', 'autoupdate', false)
    end
    errorbar(angles(ai)-barOffset, mean(asNaive{ai}), sem(asNaive{ai}), 'k')
    errorbar(angles(ai)+barOffset, mean(asExpert{ai}), sem(asExpert{ai}), 'k')
    
    [~, p(ai)] = ttest2(asNaive{ai}, asExpert{ai});
end

ylabel('Angle selectivity')
xlim([35 145])
xticks(angles)
xlabel('Preferred angle (\circ)')



%%
%% 7D - Tuning curves (45 and 135 degrees)
%%

figure, subplot(211), hold on
angleInd = 1;
naiveShape = cell2mat(atShapes(:,angleInd,1));
expertShape = cell2mat(atShapes(:,angleInd,2));
plot(angles, mean(naiveShape), 'color', colorsTransient(1,:))
plot(angles, mean(expertShape), 'color', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'box', 'off', 'autoupdate', false)
boundedline(angles, mean(naiveShape), sem(naiveShape), 'cmap', colorsTransient(1,:))
boundedline(angles, mean(expertShape), sem(expertShape), 'cmap', colorsTransient(2,:))
xlim([40, 140]), xticks(angles), xlabel('Object angle (\circ)')
ylim([0 0.8]), yticks(0:0.2:0.8), ylabel('Response')
title('Preferring 45\circ')

subplot(212), hold on
angleInd = 7;
naiveShape = cell2mat(atShapes(:,angleInd,1));
expertShape = cell2mat(atShapes(:,angleInd,2));
plot(angles, mean(naiveShape), 'color', colorsTransient(1,:))
plot(angles, mean(expertShape), 'color', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'box', 'off', 'autoupdate', false)
boundedline(angles, mean(naiveShape), sem(naiveShape), 'cmap', colorsTransient(1,:))
boundedline(angles, mean(expertShape), sem(expertShape), 'cmap', colorsTransient(2,:))
xlim([40, 140]), xticks(angles), xlabel('Object angle ()')
ylim([0 0.8]), yticks(0:0.2:0.8), ylabel('Response')
title('Preferring 135\circ')




%%
%% S7G (intermediate angles) 
%%
figure
for i = 1:5
    subplot(1,5,i), hold on
    angleInd = i+1;
    naiveShape = cell2mat(atShapes(:,angleInd,1));
    expertShape = cell2mat(atShapes(:,angleInd,2));
    plot(angles, mean(naiveShape), 'color', colorsTransient(1,:))
    plot(angles, mean(expertShape), 'color', colorsTransient(2,:))
    if i == 1
        legend({'Naive', 'Expert'}, 'box', 'off', 'autoupdate', false)
    end
    boundedline(angles, mean(naiveShape), sem(naiveShape), 'cmap', colorsTransient(1,:))
    boundedline(angles, mean(expertShape), sem(expertShape), 'cmap', colorsTransient(2,:))
    xlim([40, 140]), xticks(angles), xlabel('Object angle (\circ)')
    ylim([0 0.8]), yticks(0:0.2:0.8), ylabel('Response')
    title(['Preferring ', num2str(angles(angleInd)), '\circ'])
end

