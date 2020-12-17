
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')
%% basic settings

calciumDir = [baseDir, 'Calcium\'];
load([calciumDir, 'piezo_awake_matching'], 'matchedID', 'tformAll', 'cellAwakeAll', 'cellPiezoAll')
load([calciumDir, 'cellmaps_and_tuning_active'], 'naive')
load([calciumDir, 'summaryPiezoVsAwake_New'], 'isTouch', 'isTuned', 'piezoTuningCurve', 'whiskerEncoding', 'whiskerIot')

awakeTune = load([calciumDir, 'angle_tuning_summary_preAnswer_perTouch_NC_PTC'], 'naive');

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = [4,3,3,1,7,2,1,3,3,3,3,3];
piezoInd = [1:7,9,11];

for pi = 1 : length(piezoInd)
    piezo(pi) = load(sprintf('%spiezo\\newAnimalNaive%03d',calciumDir, mice(piezoInd(pi))), 'newAnimal');
end

matchedMiceInd = [1,1,2,3,3,4,ones(1,8)*5,ones(1,7)*6,ones(1,4)*7, ones(1,4)*9, ones(1,8)*11];

daysApart = [7,6,7,2,17,3,2,5,5,4];
piezoAngles = 0:45:315;
polarplotRad = deg2rad([piezoAngles,0]);
objectAngles = 45:15:135;

%%
%% Figure S5A
%%
% Mean tuned angle and length of all passive deflection-responsive neurons
% and angle-tuned neurons
meanTunedAngleAll = cell(9,1);
ampTunedAngleAll = cell(9,1);
indTouch = cell(9,1);
indTuned = cell(9,1);
for mi = 1 : 9
    tempMean = cell(8,1);
    tempAmp = cell(8,1);
    tempTouch = cell(8,1);
    tempTuned = cell(8,1);
    for pi = 1 : 8
        tempMean{pi} = cellfun(@(x) rad2deg(circ_mean(deg2rad(piezoAngles)', mean(x)')), ...
            piezo(mi).newAnimal(pi).meanResponse)';
        tempAmp{pi} = cellfun(@(x) circ_r(deg2rad(piezoAngles)', mean(x)'), ...
            piezo(mi).newAnimal(pi).meanResponse)';
        tempTouch{pi} = piezo(mi).newAnimal(pi).active';
        tempTuned{pi} = piezo(mi).newAnimal(pi).tuned' == 1;
    end
    meanTunedAngleAll{mi} = cell2mat(tempMean);
    ampTunedAngleAll{mi} = cell2mat(tempAmp);
    indTouch{mi} = cell2mat(tempTouch);
    indTuned{mi} = cell2mat(tempTuned);
end

%
tempAngle = cell2mat(meanTunedAngleAll);
tempAmp = cell2mat(ampTunedAngleAll);
tempTouch = cell2mat(indTouch);
tempTuned = cell2mat(indTuned);

figure
polarplot(deg2rad(tempAngle(find(tempTouch==1))), tempAmp(find(tempTouch==1)), 'k.'), hold on
polarplot(deg2rad(tempAngle(find(tempTuned==1))), tempAmp(find(tempTuned==1)), 'k.')
% legend({'All passive deflection-responsive neurons', 'All passive deflection angle-tuned neurons'})
rlim([0 0.8])

ax = gca;
ax.ThetaZeroLocation = 'left';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = piezoAngles;


%%
%% S5B
%%
% Mean tuning curve of all passive deflection responsive neurons and angle-tuned neurons
% Mouse-averaged

% From all piezo responsive neurons

allPiezoTuningCurve = cell(9,8);
for mi = 1 : length(piezo)
    for pi = 1 : length(piezo(mi).newAnimal)
        allPiezoTuningCurve{mi,pi} = cell2mat(cellfun(@mean, piezo(mi).newAnimal(pi).meanResponse(find(piezo(mi).newAnimal(pi).active))', 'un', 0));
    end
end

mousePiezoTuningCurve = cell(9,1);
for mi = 1 : length(piezo)
    mousePiezoTuningCurve{mi} = cell2mat(allPiezoTuningCurve(mi,:)');
end

% From all piezo angle-tuned neurons
allTunedPiezoTuningCurve = cell(9,8);
for mi = 1 : length(piezo)
    for pi = 1 : length(piezo(mi).newAnimal)
        allTunedPiezoTuningCurve{mi,pi} = cell2mat(cellfun(@mean, piezo(mi).newAnimal(pi).meanResponse(find(piezo(mi).newAnimal(pi).tuned==1))', 'un', 0));
    end
end

mousePiezoTunedTuningCurve = cell(9,1);
for mi = 1 : length(piezo)
    mousePiezoTunedTuningCurve{mi} = cell2mat(allTunedPiezoTuningCurve(mi,:)');
end

% draw

tempMatTouch = cell2mat(cellfun(@nanmean, mousePiezoTuningCurve, 'un', 0));
tempMatTouch = [tempMatTouch, tempMatTouch(:,1)];

tempMatTuned = cell2mat(cellfun(@nanmean, mousePiezoTunedTuningCurve, 'un', 0));
tempMatTuned = [tempMatTuned, tempMatTuned(:,1)];

figure,
polarplot(polarplotRad, mean(tempMatTouch), 'k-', 'linewidth', 1), hold on
polarplot(polarplotRad, mean(tempMatTuned), 'r-', 'linewidth', 1)
legend({'All touch responsive', 'All deflection angle-tuned'}, 'autoupdate', false)

polarplot(polarplotRad, mean(tempMatTouch) + sem(tempMatTouch) , '-', 'color', [0.6 0.6 0.6], 'linewidth', 0.5)
polarplot(polarplotRad, mean(tempMatTouch) - sem(tempMatTouch) , '-', 'color', [0.6 0.6 0.6], 'linewidth', 0.5)

polarplot(polarplotRad, mean(tempMatTuned) + sem(tempMatTuned) , '-', 'color', [1 0.6 0.6], 'linewidth', 0.5)
polarplot(polarplotRad, mean(tempMatTuned) - sem(tempMatTuned) , '-', 'color', [1 0.6 0.6], 'linewidth', 0.5)


% touchAnovaP = anova1(tempMatTouch, [], 'off');
[~,touchPmean, touchMmean] = paired_test(mean(tempMatTouch(:,[3,7]),2), mean(tempMatTouch(:,[1,5]),2));
% [~,touchPmax, touchMmax] = paired_test(max(tempMatTouch(:,[3,7]),[],2), max(tempMatTouch(:,[1,5]),[],2));
[~,tunedPmean, tunedMmean] = paired_test(mean(tempMatTuned(:,[3,7]),2), mean(tempMatTuned(:,[1,5]),2));
% [~,tunedPmax, tunedMmax] = paired_test(max(tempMatTuned(:,[3,7]),[],2), max(tempMatTuned(:,[1,5]),[],2));

title({'Mouse-averaged piezo tuning curve, all angle-tuned neurons'; ...
    sprintf(' Touch p = %s (%s), Tuned p = %s (%s) (mean horizontal vs mean vertical)',num2str(touchPmean,3), touchMmean, num2str(tunedPmean,3), tunedPmean)})
%     sprintf(' Touch p = %s (%s), Tuned p = %s (%s) (max horizontal vs max vertial)',num2str(touchPmax,3), touchMmax, num2str(tunedPmax,3), tunedMmax);})

ax = gca;
ax.ThetaZeroLocation = 'left';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = piezoAngles;

% mean(tempMatTouch(:,[3,7]),2) / mean(tempMatTouch(:,[1,5]),2)



%%
%% 5A-G
%%
% example co-tuned neuron
% it's map, passive responses, and active responses

%%
%% 5A - Passive deflection under light anesthesia
%%
%%

%%
%% 5B - Passive deflection response
%%
coTunedInd = intersect( find(isfinite(isTuned{ind,1})), find(isfinite(isTuned{ind,2})) );

%
%
%
cti = 3;
%
%
%
piezoCi = mod(matchedID{ind,2}(coTunedInd(cti)),1000);
figure, imagesc(piezo(animalInd).newAnimal(piezoPlaneInd).trialOrdered2p{piezoCi}(:,1:end-1)), colormap gray, colorbar

%%
%% 5C - Vector sum of deflection response
%%
meanResponse = piezo(animalInd).newAnimal(piezoPlaneInd).meanResponse{piezoCi};
tempMat = [meanResponse, meanResponse(:,1)];
vecDeg = isTuned{ind,4}(coTunedInd(cti));
vecLen = isTuned{ind,5}(coTunedInd(cti));

figure,
polarplot(polarplotRad, mean(tempMat), 'k-'), hold on
polarplot(polarplotRad, mean(tempMat)+sem(tempMat), '-', 'color', [0.6 0.6 0.6])
polarplot(polarplotRad, mean(tempMat)-sem(tempMat), '-', 'color', [0.6 0.6 0.6])
polarplot(ones(1,2)*deg2rad(vecDeg), [0 vecLen], 'k-', 'linewidth', 3)
ax = gca;
ax.ThetaZeroLocation = 'left';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = piezoAngles;

%%
%% 5D - Active 7-angle task
%%

%%
%% 5E - Active touch response
%%
mouse = mice(mouseInd);
session = sessions(mouseInd);
load(sprintf('%s%03d\\UberJK%03dS%02d_NC',calciumDir, mouse, mouse, session), 'u')
spkCi = u.trials{u.planeTrialInds{1}(1)}.neuindSession == matchedID{ind,1}(coTunedInd(cti));
numFrames = floor(u.frameRate * 4); % up to 3 sec
spkAwakeCotuned = cell2mat(cellfun(@(x) x.spk(spkCi,1:numFrames), u.trials(u.planeTrialInds{1}), 'un', 0));
awakeCotunedAngles = cellfun(@(x) x.angle, u.trials(u.planeTrialInds{1}));
[sortAngles, sorti] = sort(awakeCotunedAngles);
spkAwakeCotuned = spkAwakeCotuned(sorti,:);

figure, imagesc(spkAwakeCotuned), colormap gray, colorbar

angleBoundaries = find(diff(sortAngles))+1;


%%
%% 5F - Active touch angle-tuning curve
%%
awakeCi = find(awakeTune.naive(mouseInd).touchID == matchedID{ind,1}(coTunedInd(cti)));

figure, boundedline(objectAngles, cellfun(@mean, awakeTune.naive(mouseInd).val{awakeCi}), ...
    cellfun(@sem, awakeTune.naive(mouseInd).val{awakeCi}), 'k-' )
xlim([40 140]), xticks(objectAngles)




%%
%% 5G - cell map
%%

ind = 25; % JK039 plane #4
% awakeXuse = 1:640;
% piezoYuse = 110:end;

mouseInd = matchedMiceInd(ind);
awakePlaneInd = floor(matchedID{ind,1}(1)/1000);

animalInd = find(piezoInd == mouseInd);
piezoPlaneInd = floor(matchedID{ind,2}(1)/1000);

tform = tformAll{ind};
mimgPiezo = piezo(animalInd).newAnimal(piezoPlaneInd).mimgPiezo(110:end,:);
mimgAwake = imwarp(naive(mouseInd).mimg{awakePlaneInd}(:,1:640), tform, 'OutputView', imref2d(size(mimgPiezo)));




%%
cellPiezoBinary = zeros(size(cellPiezoAll{ind}));
cellPiezoBinary(find(cellPiezoAll{ind})) = 1;
cellPiezoColor = zeros([size(cellPiezoAll{ind}),3]);
cellPiezoColor(:,:,1) = cellPiezoBinary;
cellPiezoColor(:,:,3) = cellPiezoBinary;
figure, imshowpair(mimgPiezo(107:474, 61:end), cellPiezoColor(107:474, 61:end,:), 'blend', 'Scaling', 'joint'), title('Passive deflection')




%%
cellAwakeBinary = zeros(size(cellAwakeAll{ind}));
cellAwakeBinary(find(cellAwakeAll{ind})) = 1;
cellAwakeColor = zeros([size(cellAwakeAll{ind}),3]);
cellAwakeColor(:,:,2) = cellAwakeBinary;
figure, imshowpair(adapthisteq(mat2gray(mimgAwake(107:474, 61:end))), cellAwakeColor(107:474, 61:end,:), 'blend', 'Scaling', 'joint'), title('Passive deflection')


%%
figure, imshowpair(cellAwakeAll{ind}(107:474, 61:end), cellPiezoAll{ind}(107:474, 61:end)), title('Merged')



%% Passive deflection example cell map
cellmapPiezo = zeros(size(cellPiezoAll{ind}));
cellmapPiezo(cellPiezoAll{ind} == matchedID{ind,2}(coTunedInd(1))) = 1;
figure, imshowpair(mimgPiezo, cellmapPiezo)

%% Active touch example cell map

cellmapAwake = zeros(size(cellAwakeAll{ind}));
cellmapAwake(cellAwakeAll{ind} == matchedID{ind,1}(coTunedInd(1))) = 1;
figure, imshowpair(mimgAwake, cellmapAwake)





%%
%% 5H - Venn Diagram
%%





%%
%% 5I - Mean tuning curves, with preferred object angle
%%

ptcPiezoTuned = cell(1,3);

allAwakeTunedAngle = cell2mat(isTuned(:,1));
ind45 = find(allAwakeTunedAngle == 45);
indMid = find(allAwakeTunedAngle >= 75 & allAwakeTunedAngle <= 105);
ind135 = find(allAwakeTunedAngle == 135);
ptcAll = cell2mat(piezoTuningCurve);


ptcPiezoTuned{1} = ptcAll(intersect(ind45,indPiezoTuned),:);
ptcPiezoTuned{2} = ptcAll(intersect(indMid,indPiezoTuned),:);
ptcPiezoTuned{3} = ptcAll(intersect(ind135,indPiezoTuned),:);
ptcPiezoTuned{4} = ptcAll(indPiezoTuned,:);

figure,
polarplot(deg2rad([piezoAngles,0]),[nanmean(ptcPiezoTuned{4}), nanmean(ptcPiezoTuned{4}(:,1))], 'k-', 'linewidth',3), hold on
polarplot(deg2rad([piezoAngles,0]),[mean(ptcPiezoTuned{1}), mean(ptcPiezoTuned{1}(:,1))], 'b-', 'linewidth',1.5)
polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{2}), mean(ptcPiezoTuned{2}(:,1))], 'g-', 'linewidth',1.5)
polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{3}), mean(ptcPiezoTuned{3}(:,1))], 'r-', 'linewidth',1.5)

polarplot(deg2rad([piezoAngles,0]), [nanmean(ptcPiezoTuned{4}), nanmean(ptcPiezoTuned{4}(:,1))] ...
    + [sem(ptcPiezoTuned{4}), sem(ptcPiezoTuned{4}(:,1))], '-', 'linewidth',0.5, 'color', [0.6 0.6 0.6]),
polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{4}), mean(ptcPiezoTuned{4}(:,1))] ...
    - [sem(ptcPiezoTuned{4}), sem(ptcPiezoTuned{4}(:,1))], '-', 'linewidth',0.5, 'color', [0.6 0.6 0.6])

polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{1}), mean(ptcPiezoTuned{1}(:,1))] ...
    + [sem(ptcPiezoTuned{1}), sem(ptcPiezoTuned{1}(:,1))], 'c-', 'linewidth',0.5),
polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{1}), mean(ptcPiezoTuned{1}(:,1))] ...
    - [sem(ptcPiezoTuned{1}), sem(ptcPiezoTuned{1}(:,1))], 'c-', 'linewidth',0.5)

polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{2}), mean(ptcPiezoTuned{2}(:,1))]...
    + [sem(ptcPiezoTuned{2}), sem(ptcPiezoTuned{2}(:,1))], '-', 'color', [0.4 1 0.4], 'linewidth',0.5),
polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{2}), mean(ptcPiezoTuned{2}(:,1))]...
    - [sem(ptcPiezoTuned{2}), sem(ptcPiezoTuned{2}(:,1))], '-', 'color', [0.4 1 0.4], 'linewidth',0.5)

polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{3}), mean(ptcPiezoTuned{3}(:,1))]...
    + [sem(ptcPiezoTuned{3}), sem(ptcPiezoTuned{3}(:,1))], 'm-', 'linewidth',0.5),
polarplot(deg2rad([piezoAngles,0]), [mean(ptcPiezoTuned{3}), mean(ptcPiezoTuned{3}(:,1))]...
    - [sem(ptcPiezoTuned{3}), sem(ptcPiezoTuned{3}(:,1))], 'm-', 'linewidth',0.5)


legend({'Co-tuned','45','75-105', '135'})
title({'All piezo tuned neurons'})

ax = gca;
ax.ThetaZeroLocation = 'left';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = piezoAngles;






%%
%% Fig 5J
%%

piezoMeanAngle = cell2mat(cellfun(@(x,y) y(find(isfinite(x).*isfinite(y))), isTuned(:,2), isTuned(:,4), 'un', 0));
piezoMeanAmplitude = cell2mat(cellfun(@(x,y) y(find(isfinite(x).*isfinite(y))), isTuned(:,2), isTuned(:,5), 'un', 0));


piezoMeanAngle45 = cell2mat(cellfun(@(x,y) y(find((x==45).*isfinite(y))), isTuned(:,1), isTuned(:,4), 'un', 0));
piezoMeanAmplitude45 = cell2mat(cellfun(@(x,y) y(find((x==45).*isfinite(y))), isTuned(:,1), isTuned(:,5), 'un', 0));

piezoMeanAngleMid = cell2mat(cellfun(@(x,y) y(find((x>=75 & x<=105).*isfinite(y))), isTuned(:,1), isTuned(:,4), 'un', 0));
piezoMeanAmplitudeMid = cell2mat(cellfun(@(x,y) y(find((x>=75 & x<=105).*isfinite(y))), isTuned(:,1), isTuned(:,5), 'un', 0));

piezoMeanAngle135 = cell2mat(cellfun(@(x,y) y(find((x==135).*isfinite(y))), isTuned(:,1), isTuned(:,4), 'un', 0));
piezoMeanAmplitude135 = cell2mat(cellfun(@(x,y) y(find((x==135).*isfinite(y))), isTuned(:,1), isTuned(:,5), 'un', 0));

% Polar histogram of mean tuned angle

histBin = 30;
histRange = [-180:histBin:180];
distTemp = -180:histBin:180;
distPlot = [distTemp; distTemp];
distPlotRad = deg2rad(distPlot(:));
distPlotRad = distPlotRad(2:end);
distPiezoMeanTuned = histcounts(piezoMeanAngle, histRange);
distPiezoMeanTuned45 = histcounts(piezoMeanAngle45, histRange);
distPiezoMeanTunedMid = histcounts(piezoMeanAngleMid, histRange);
distPiezoMeanTuned135 = histcounts(piezoMeanAngle135, histRange);

distPiezoMeanTunedNorm = histcounts(piezoMeanAngle, histRange, 'norm', 'probability');
distPiezoMeanTunedNorm45 = histcounts(piezoMeanAngle45, histRange, 'norm', 'probability');
distPiezoMeanTunedNormMid = histcounts(piezoMeanAngleMid, histRange, 'norm', 'probability');
distPiezoMeanTunedNorm135 = histcounts(piezoMeanAngle135, histRange, 'norm', 'probability');

tempAll = [distPiezoMeanTunedNorm; distPiezoMeanTunedNorm];
tempAll = tempAll(:);
tempAll = [tempAll; tempAll(1)];

temp45 = [distPiezoMeanTunedNorm45; distPiezoMeanTunedNorm45];
temp45 = temp45(:);
temp45 = [temp45; temp45(1)];

tempMid = [distPiezoMeanTunedNormMid; distPiezoMeanTunedNormMid];
tempMid = tempMid(:);
tempMid = [tempMid; tempMid(1)];

temp135 = [distPiezoMeanTunedNorm135; distPiezoMeanTunedNorm135];
temp135 = temp135(:);
temp135 = [temp135; temp135(1)];

figure,
polarplot(distPlotRad, tempAll, 'k-', 'linewidth', 3), hold on
polarplot(distPlotRad-0.01, temp45, 'b-')
polarplot(distPlotRad+0.01, tempMid, 'g-')
polarplot(distPlotRad+0.02, temp135, 'r-')
legend({'All co-tuned', 'Awake 45', 'Awake 75-105', 'Awake 135'})
title('Mean passive tuned angle distribution')

ax = gca;
ax.ThetaZeroLocation = 'left';
ax.ThetaDir = 'clockwise';
% ax.ThetaTick = piezoAngles;




%%
%% 5K - Upward vs Downward deflection for 45 and 135 degree tuning
%%
% A little bit backwards. (270~300 i.e., -90 ~ -60; and 60~90) 

% stacked bar, but only considering co-tuned neurons
allObjectAngle = cell2mat(isTuned(:,1));
indObject45 = find(allObjectAngle == 45);
indObject135 = find(allObjectAngle == 135);

allPiezoAngle = cell2mat(isTuned(:,4));
% indPiezoNT = find(isnan(allPiezoAngle));
indPiezoDown = find(allPiezoAngle > -90 & allPiezoAngle < -60);
indPiezoUp = find(allPiezoAngle > 60 & allPiezoAngle < 90);
indPiezoTuned = find(isfinite(allPiezoAngle));
indPiezoOther = setdiff(indPiezoTuned, union(indPiezoUp, indPiezoDown));

props = zeros(3,2); % (:,1) 45, (:,1) 135; (1,:) up, (2,:) not-tuned, (3,:) down

props(1,1) = length(intersect(indObject45, indPiezoUp)) / length(intersect(indObject45, indPiezoTuned));
props(2,1) = length(intersect(indObject45, indPiezoOther)) / length(intersect(indObject45, indPiezoTuned));
props(3,1) = length(intersect(indObject45, indPiezoDown)) / length(intersect(indObject45, indPiezoTuned));

props(1,2) = length(intersect(indObject135, indPiezoUp)) / length(intersect(indObject135, indPiezoTuned));
props(2,2) = length(intersect(indObject135, indPiezoOther)) / length(intersect(indObject135, indPiezoTuned));
props(3,2) = length(intersect(indObject135, indPiezoDown)) / length(intersect(indObject135, indPiezoTuned));

numBS = 10000;
numShuffle = 10000;
ciDown45 = zeros(numBS,1);
ciUp45 = zeros(numBS,1);
ciDown135 = zeros(numBS,1);
ciUp135 = zeros(numBS,1);
shuffleDown45 = zeros(numShuffle,1);
shuffleUp45 = zeros(numShuffle,1); % just to check
shuffleDown135 = zeros(numShuffle,1); % just to check
shuffleUp135 = zeros(numShuffle,1);
for bsi = 1 : numBS
    bsInds = randi(length(allObjectAngle), length(allObjectAngle),1);
    bsObjectAngle = allObjectAngle(bsInds);
    bsIndObject45 = find(bsObjectAngle == 45);
    bsIndObject135 = find(bsObjectAngle == 135);
    
    bsPiezoAngle = allPiezoAngle(bsInds);
    bsIndPiezoDown = find(bsPiezoAngle > -90 & bsPiezoAngle < -60);
    bsIndPiezoUp = find(bsPiezoAngle > 60 & bsPiezoAngle < 90);
    bsIndPiezoTuned = find(isfinite(bsPiezoAngle));
    
    ciDown45(bsi) = length(intersect(bsIndObject45, bsIndPiezoDown)) / length(intersect(bsIndObject45, bsIndPiezoTuned));
    ciUp45(bsi) = length(intersect(bsIndObject45, bsIndPiezoUp)) / length(intersect(bsIndObject45, bsIndPiezoTuned));
    ciDown135(bsi) = length(intersect(bsIndObject135, bsIndPiezoDown)) / length(intersect(bsIndObject135, bsIndPiezoTuned));
    ciUp135(bsi) = length(intersect(bsIndObject135, bsIndPiezoUp)) / length(intersect(bsIndObject135, bsIndPiezoTuned));
end

for shi = 1 : numShuffle
    shPiezoAngle = allPiezoAngle(randperm(length(allPiezoAngle)));
    shIndPiezoDown = find(shPiezoAngle > -90 & shPiezoAngle < -60);
    shIndPiezoUp = find(shPiezoAngle > 60 & shPiezoAngle < 90);
    shIndPiezoTuned = find(isfinite(shPiezoAngle));
    shuffleDown45(shi) = length(intersect(shIndPiezoDown, indObject45)) / length(intersect(indObject45, shIndPiezoTuned));
    shuffleUp45(shi) = length(intersect(shIndPiezoUp, indObject45)) / length(intersect(indObject45, shIndPiezoTuned));
    shuffleDown135(shi) = length(intersect(shIndPiezoDown, indObject135)) / length(intersect(indObject135, shIndPiezoTuned));
    shuffleUp135(shi) = length(intersect(shIndPiezoUp, indObject135)) / length(intersect(indObject135, shIndPiezoTuned));
end

% draw
figure, hold on

bar(1,sum(props(1:3,1)), 'facecolor', 'r', 'edgecolor', 'r')
bar(1,sum(props(2:3,1)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', [0.6 0.6 0.6])
bar(1,props(3,1), 'facecolor', 'b', 'edgecolor', 'b')

plot([0.5 1.5], ones(1,2)*prctile(ciDown45, 97.5), 'c--')
plot([0.5 1.5], ones(1,2)*prctile(ciDown45, 2.5), 'c--')
plot([0.5 1.5], 1-ones(1,2)*prctile(ciUp45, 97.5), 'm--')
plot([0.5 1.5], 1-ones(1,2)*prctile(ciUp45, 2.5), 'm--')


plot([0.5 2.5], ones(1,2)*(prctile(shuffleDown45, 5)+prctile(shuffleDown135,5))/2, '--', 'color', [0.6 0.6 0.6])
plot([0.5 2.5], ones(1,2)*(mean(shuffleDown45)+mean(shuffleDown135))/2, 'k--')
plot([0.5 2.5], ones(1,2)*(prctile(shuffleDown45, 95)+prctile(shuffleDown135,95))/2, '--', 'color', [0.6 0.6 0.6])


bar(2,sum(props(1:3,2)), 'facecolor', 'r' , 'edgecolor', 'r')
bar(2,sum(props(2:3,2)), 'facecolor', [0.6 0.6 0.6], 'edgecolor', [0.6 0.6 0.6])
bar(2,props(3,2), 'facecolor', 'b', 'edgecolor', 'b')

plot([1.5 2.5], ones(1,2)*prctile(ciDown135, 97.5), 'c--')
plot([1.5 2.5], ones(1,2)*prctile(ciDown135, 2.5), 'c--')
plot([1.5 2.5], 1-ones(1,2)*prctile(ciUp135, 97.5), 'm--')
plot([1.5 2.5], 1-ones(1,2)*prctile(ciUp135, 2.5), 'm--')

plot([0.5 2.5], 1-ones(1,2)*(prctile(shuffleUp135, 5)+prctile(shuffleUp45, 5))/2, '--', 'color', [0.6 0.6 0.6])
plot([0.5 2.5], 1-ones(1,2)*(mean(shuffleUp135)+mean(shuffleUp45))/2, 'k--')
plot([0.5 2.5], 1-ones(1,2)*(prctile(shuffleUp135, 95)+prctile(shuffleUp45, 95))/2, '--', 'color', [0.6 0.6 0.6])

xlim([0 3]), xticks(1:2), xticklabels({'45', '135'})
yticks(0:0.2:1)

%% check if values from shuffling are the same across angle selection
mean(shuffleUp45)
prctile(shuffleUp45,95)

mean(shuffleUp135)
prctile(shuffleUp135,95)



mean(shuffleDown45)
prctile(shuffleDown45,95)

mean(shuffleDown135)
prctile(shuffleDown135,95)
















%%
%% S5C - Preferred passive deflection direction histogram
%%
awakeTunedAngle = cell2mat(cellfun(@(x,y) x(find(isfinite(x).*isfinite(y))), isTuned(:,1), isTuned(:,2), 'un', 0));
piezoTunedAngle = cell2mat(cellfun(@(x,y) y(find(isfinite(x).*isfinite(y))), isTuned(:,1), isTuned(:,2), 'un', 0));

angleRange = [piezoAngles,360];
polarDist = zeros(4,length(angleRange)-1);
polarDist(1,:) = histcounts(piezoTunedAngle(find(awakeTunedAngle == 45)), angleRange, 'norm', 'probability');
polarDist(2,:) = histcounts(piezoTunedAngle(find(awakeTunedAngle >= 75 & awakeTunedAngle <=105)), angleRange, 'norm', 'probability');
polarDist(3,:) = histcounts(piezoTunedAngle(find(awakeTunedAngle == 135)), angleRange, 'norm', 'probability');

polarDist(4,:) = histcounts(piezoTunedAngle(find(isfinite(piezoTunedAngle))), angleRange, 'norm', 'probability');

p = circ_wwtest( piezoTunedAngle(find(awakeTunedAngle == 45)), piezoTunedAngle(find(awakeTunedAngle == 135)), ones(length(find(awakeTunedAngle==45)),1),  ones(length(find(awakeTunedAngle==135)),1));


tempAngle = [angleRange;angleRange]-22.5;
polarRad = deg2rad(tempAngle(:));
polarRad = polarRad(2:end);

tempAll = [polarDist(4,:); polarDist(4,:)];
tempAll = tempAll(:);
tempAll = [tempAll; tempAll(1)];

temp45 = [polarDist(1,:); polarDist(1,:)];
temp45 = temp45(:);
temp45 = [temp45; temp45(1)];

tempMid = [polarDist(2,:); polarDist(2,:)];
tempMid = tempMid(:);
tempMid = [tempMid; tempMid(1)];

temp135 = [polarDist(3,:); polarDist(3,:)];
temp135 = temp135(:);
temp135 = [temp135; temp135(1)];

figure
polarplot(polarRad, tempAll, 'k-', 'linewidth', 2), hold on
polarplot(polarRad-0.01, temp45, 'b-'), hold on
polarplot(polarRad+0.01, tempMid, 'g-')
polarplot(polarRad+0.02, temp135, 'r-')

ax = gca;
ax.ThetaZeroLocation = 'left';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = piezoAngles;

legend({'All co-tuned', '45', '75-105','135'})
title(sprintf('p = %s', num2str(p,3)))










