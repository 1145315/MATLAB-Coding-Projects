ca= readcell('RData.xlsx');
data= ca(3:38,:);
[r,c]= size(data);
%% for logistic regression likelihood for malicious intent

% get intent columns
intent= data(:,34);
% binarize them
intentLog= strcmpi(intent, 'Malicious Intent');

% get embodiment col
body= data(:,1)
% binarize them
bodyLog= strcmpi(body, 'AI');

% get truthfulness col
truthfulness= data(:,2)
% binarize them
truthLog= strcmpi(truthfulness, 'Dishonest');

% get outcome col

% put last 3 in 3 x C matrix
varmat= [bodyLog truthLog];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Avg Trust scores%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trustVals= data(:,20:end-9);
trustVals= [trustVals(:,1:end-5) trustVals(:,1:end-3:end)];
[r,c]= size(trustVals);
avgTS= [];
score=0;
for x= 1:r;
    for y=1:c;
        prompt= trustVals{x,y};
        if y >= 1 & y <=5;
            switch prompt;
                case 'Strongly disagree'
                    score= score + 7;
                case 'Disagree'
                    score= score + 6;
                case 'Somewhat disagree'
                    score= score + 5;
                case 'Neither agree nor disagree'
                    score= score + 4;
                case 'Somewhat agree'
                    score= score + 3;
                case 'Agree'
                    score= score + 2;
                case 'Strongly agree'
                    score= score + 1;
            end 
        else
            switch prompt;
                case 'Strongly disagree'
                    score= score + 1;
                case 'Disagree'
                    score= score + 2;
                case 'Somewhat disagree'
                    score= score + 3;
                case 'Neither agree nor disagree'
                    score= score + 4;
                case 'Somewhat agree'
                    score= score + 5;
                case 'Agree'
                    score= score + 6;
                case 'Strongly agree'
                    score= score + 7;
            end
        end
    end
    avg= score ./ 12;
    avgTS= [avgTS avg;];
    score= 0;
end
avgTS= avgTS'    
GLM= fitglm(varmat,intentLog, 'linear','distribution','binomial','link', 'logit')
groups= [body truthfulness]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Arg= [body truthfulness num2cell(avgTS)]
r1 = Arg(:,1);
r2 = Arg(:,2);
R =cell2mat(Arg(:,3));
q1 = unique(r1);
q2 = unique(r2);
if isnumeric(q1)
 q1 = cellstr(num2str(q1));r1 = cellstr(num2str(r1));
end
if isnumeric(q2)
 q2 = cellstr(num2str(q2));r2 = cellstr(num2str(r2));
end
Z = zeros(size(R));
Yaligned1 = zeros(size(R));
Yaligned2 = zeros(size(R));
Yaligned3 = zeros(size(R));
mu = nanmean(R);
for i = 1:size(q1,1)
 ME1 = nanmean(R(strcmp(r1,q1{i}))) - mu;
 for j = 1:size(q2,1)
 Z(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) = R(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) - nanmean(R(strcmp(r1,q1{i}) & strcmp(r2,q2{j})));
 ME2 = nanmean(R(strcmp(r2,q2{j}))) - mu;
 ME12 = nanmean(R(strcmp(r1,q1{i}) & strcmp(r2,q2{j}))) - mu;
 MEI = ME12 - ME1 - ME2 + mu; 
 Yaligned1(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) = Z(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) + ME1;
 Yaligned2(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) = Z(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) + ME2;
 Yaligned3(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) = Z(strcmp(r1,q1{i}) & strcmp(r2,q2{j})) + MEI;
 end
end

R1 = tiedrank(Yaligned1);
R2 = tiedrank(Yaligned2);
R3 = tiedrank(Yaligned3);
Raligned = [R1,R2,R3]
[p,tbl,stats] = anovan(Raligned(:, 3), {body, truthfulness}, 'model', 'interaction', 'varnames', {'Embodiment', 'Truthfulness'})

% Post-hoc tests for Truthfulness and Embodiment
comparison_truthfulness = multcompare(stats, 'Dimension', 1);
comparison_embodiment = multcompare(stats, 'Dimension', 2);

% Group means and standard errors
[means_truthfulness, se_truthfulness] = grpstats(avgTS, truthfulness, {'mean', 'sem'});
[means_embodiment, se_embodiment] = grpstats(avgTS, body, {'mean', 'sem'});

% --- Plot 1: Effect of Truthfulness ---
figure;
bar(means_truthfulness, 'FaceColor', 'flat');
hold on;
errorbar(1:2, means_truthfulness, se_truthfulness, 'k.', 'LineWidth', 1.5);
xticks([1 2]);
xticklabels({'Honest', 'Dishonest'});
ylabel('Average Change in Trust');
title('Effect of Truthfulness on Trust');
set(gca, 'FontSize', 12);
text(1.5, max(means_truthfulness) + 0.1, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
hold off;

% --- Plot 2: Effect of Embodiment ---
figure;
bar(means_embodiment, 'FaceColor', 'flat');
hold on;
errorbar(1:length(means_embodiment), means_embodiment, se_embodiment, 'k.', 'LineWidth', 1.5);
xticks([1 2 3]);
xticklabels({'AI', 'Robot', 'Human'});
ylabel('Average Change in Trust');
title('Effect of Embodiment on Trust');
set(gca, 'FontSize', 12);
text(2, max(means_embodiment) + 0.1, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
hold off;

% --- Plot 3: Interaction Effect ---
interaction_means = grpstats(avgTS, {body, truthfulness}, 'mean');
interaction_se = grpstats(avgTS, {body, truthfulness}, 'sem');

figure;
group_labels = {'Honest-AI', 'Dishonest-AI', 'Honest-Robot', 'Dishonest-Robot', 'Honest-Human', 'Dishonest-Human'};
bar(interaction_means, 'FaceColor', 'flat');
hold on;
errorbar(1:length(interaction_means), interaction_means, interaction_se, 'k.', 'LineWidth', 1.5);
xticks(1:6);
xticklabels(group_labels);
ylabel('Average Change in Trust');
title('Interaction Effect of Truthfulness and Embodiment');
set(gca, 'FontSize', 10);
text(4, max(interaction_means) + 0.1, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
hold off;
