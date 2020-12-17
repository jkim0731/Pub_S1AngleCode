function [p, m] = graph_whisker_feature_naive_expert(naiveWhisker, expertWhisker, ylabelStr, titleStr, statTest, colors4face, colors4edge)
% Specialized format for plotting whisker features across training (input
% should remain same as old)
% For Angle tuning in S1 paper 
% 2020/07/08 JK

naiveWhisker(find(~isfinite(naiveWhisker(:)))) = nan;
expertWhisker(find(~isfinite(expertWhisker(:)))) = nan;

barOffset = 0.2;
barWidth = 0.4;
xpos = [8,9,10,11,13,12,1,2,3,4,5,6];

figure, hold on
bar(xpos - barOffset, nanmean(naiveWhisker), barWidth, 'facecolor', colors4face(1,:), 'edgecolor', colors4edge(1,:))
bar(xpos + barOffset, nanmean(expertWhisker), barWidth, 'facecolor', colors4face(2,:), 'edgecolor', colors4edge(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, nanmean(naiveWhisker), sem(naiveWhisker), 'k', 'lines', 'no')
errorbar(xpos + barOffset, nanmean(expertWhisker), sem(expertWhisker), 'k', 'lines', 'no')
xlabel('Whisker features')
ylabel(ylabelStr)
title(titleStr)
xlim([0.5 13.5])


p = nan(1,13);
m = cell(1,13);
if strcmpi(statTest, 'paired')
    for i = 1 : 12
        [~, p(xpos(i)), m{xpos(i)}] = paired_test(naiveWhisker(:,i), expertWhisker(:,i));
    end
elseif strcmpi(statTest, '2sample')
    for i = 1 : 12
        [~, p(xpos(i))] = ttest2(naiveWhisker(:,i), expertWhisker(:,i));
    end
end

sigInd = find(p<0.05);

title( {titleStr; sprintf('Significant features: %s',num2str(sigInd))} )