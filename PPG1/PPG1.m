%% Loading Data

clear;
clc;
close all;

% Initialization of the seed
rng('shuffle');
SyntDataBase = load ('data/SyntheticData.mat');

% Initializing a zeros matrix 
SyntData = zeros(64,105000);
j = 1;
% Extracting samples of all subjects and putting them in an array
for i=1:(numel(SyntDataBase.data))
    
    if SyntDataBase.data(i).group == "fm"
        tmp = SyntDataBase.data(i).ppg;
        
        % "fm" type needs to be padded
         if numel(tmp.v) == 104999
             tmp.v = [tmp.v 0];
         end
         
        SyntData(j,:) = tmp.v(1,:);
        j = j+1;
    end
end

% SyntData now represents the matrix of clean signals

%% Training Phase

%% Simulation movement PPG data: Adding White Gaussian Noise to clean data

num_signals = 50; % Number of signals used for the training phase

% Inizialization of the matrixs that will be the input of the SVM
BigSTDhr = zeros(1,30*num_signals);
BigSTDamp = zeros(1,30*num_signals);
BigSTDsd = zeros(1,30*num_signals);
BigSTDwav = zeros(1,30*num_signals);
BigLabels = zeros(1,30*num_signals);

% For each training signal:
for ggg=1:num_signals
    
% Flag to print or not the plots
    plot_lb = false;
    
%      if ggg == 5
%          plot_lb = true;
%      end
    
% 30 Rows, one for each segment, with the 4 parameters and labels needed for classification
    [STDhrTr,STDampTr,STDsdTr,STDwavTr,labels] = SVM_fun(SyntData(ggg,:),plot_lb);
   
% Saving the results in the big table for the SVM    
    BigSTDhr(1,1+(ggg-1)*30:(ggg*30)) = STDhrTr;
    BigSTDamp(1,1+(ggg-1)*30:(ggg*30)) = STDampTr;
    BigSTDsd(1,1+(ggg-1)*30:(ggg*30)) = STDsdTr;
    BigSTDwav(1,1+(ggg-1)*30:(ggg*30)) = STDwavTr;
    
% Plotting the parameters with the right format [210s]
    % STDsTot(STDhrTr,STDampTr,STDsdTr,STDwavTr,plot_lb);
    
    BigLabels(1,1+(ggg-1)*30:(ggg*30)) = labels;
end

%% Training SVM with segments of multiple signals
% Building the table needed by the SVM
tbl = table(BigSTDhr',BigSTDamp',BigSTDsd',BigSTDwav', ...
    'VariableNames', {'STDhr';'STDamp';'STDsd';'STDwav'});
% Columns of classifiers as input of the SVM
Y = BigLabels';

% Training a SVM with a linear kernel
Mdl = fitcsvm(tbl,Y,'KernelFunction', 'linear');

%% Test phase

num_test_sign = length(SyntData(:,1))-(num_signals+1); %The remains signals from the dataset
% Initializing the confusion matrix needed for the statistics
confusion_matrix = zeros(2,2); % [14 signals inside]

for c = num_signals+1:length(SyntData(:,1))

     plot_lb = false;
     
    if c == num_signals+2
       plot_lb = true;
    end
    
    % Computing all parameters for the SVM, labels here are needed to create
    % the confusion matrix
    [STDhrTest,STDampTest,STDsdTest,STDwavTest,labelTest] = SVM_fun(SyntData(c,:),plot_lb);

    %Creating the table for the input
    tbl_test = table(STDhrTest',STDampTest',STDsdTest',STDwavTest', ...
        'VariableNames', {'STDhr';'STDamp';'STDsd';'STDwav'});
    
    % Obtaining the predicted labels
    [Predict_label] = predict(Mdl,tbl_test);

    % Making the algorithm more robust w.r.t. a certain type of error. 
    for k = 2:length(Predict_label)-1
        if Predict_label(k) == 0
            %If the conseguente and antecedent segments are corrupted, for construction the one must be corrupted too.
           if Predict_label(k-1) == 1 && Predict_label(k+1) == 1 
               Predict_label(k) = 1;
           end
        end
    end
    
    % One way to save the labels for the 2-nd paper
    %      if c == num_signals+2
    %          save('Labels','Predict_label');
    %      end
    %    

    %Computing the Confusion matrix
    for k = 1:length(Predict_label)
        if Predict_label(k) == 0 && labelTest(k) == 0
            confusion_matrix(1,1) = confusion_matrix(1,1)+1;
        elseif Predict_label(k) == 1 && labelTest(k) == 0
            confusion_matrix(1,2) = confusion_matrix(1,2)+1;
        elseif Predict_label(k) == 0 && labelTest(k) == 1
            confusion_matrix(2,1) = confusion_matrix(2,1)+1;
        else
            confusion_matrix(2,2) = confusion_matrix(2,2)+1;
        end
    end
    
   % Optional print
   STDsTot(STDhrTest,STDampTest,STDsdTest,STDwavTest,Predict_label,plot_lb);

    % Optional print of the classified segments.
    %     STDhrCl= STDhrTest(1,Predict_label == 0)';
    %     STDhrCorr = STDhrTest(1,Predict_label == 1)';
    % 
    %     STDampCl= STDampTest(1,Predict_label == 0)';
    %     STDampCorr = STDampTest(1,Predict_label == 1)';
    % 
    %     STDsdCl= STDsdTest(1,Predict_label == 0)';
    %     STDsdCorr = STDsdTest(1,Predict_label == 1)';
    % 
    %     STDwavCl= STDwavTest(1,Predict_label == 0)';
    %     STDwavCorr = STDwavTest(1,Predict_label == 1)';
    % 
    %     figure();
    %     hold on;
    %     scatter(STDhrCl,STDampCl);
    %     scatter(STDhrCorr,STDampCorr);
    %     hold off;
    
end
%% Plots

% close all
x = linspace(min(tbl{:, 3}),max(tbl{:, 3}),10);
y = linspace(min(tbl{:, 4}),max(tbl{:, 4}),10);
z1 = mean(tbl{:, 1});
z2 = mean(tbl{:, 2});

[X1,X2] = meshgrid(x,y);

[score] = predict(Mdl,[z1*ones(size(X1(:))),z2*ones(size(X1(:))), X1(:),X2(:)]);
scoreGrid = reshape(score,size(X2,1),size(X1,1));

% Trying to plot the hyperplane computed by the SVM (Not as the paper shown)
figure
imagesc(x, y, scoreGrid, [0, 1])
colorbar
set(gca, 'Ydir', 'normal')
% title(sprintf('%.5f\n%.5f', z1, z2))
hold on
plot(tbl{Y==0, 3}, tbl{Y==0, 4}, 'cO', 'MarkerFaceColor', 'c') 
plot(tbl{Y==1, 3}, tbl{Y==1, 4}, 'rO', 'MarkerFaceColor', 'r') 
xlabel('STDsd')
ylabel('STDwav')

% Plotting statistical values of parameters
figure
c = predict(Mdl, tbl);
subplot(2,2,1)
boxplot(tbl{:, 1}, c)
title('STDhr');
xlabel('Computed labels');
ylabel('Parameter value');
subplot(2,2,2)
boxplot(tbl{:, 2}, c)
title('STDamp');
xlabel('Computed labels');
ylabel('Parameter value');
subplot(2,2,3)
boxplot(tbl{:, 3}, c)
title('STDsd');
xlabel('Computed labels');
ylabel('Parameter value');
subplot(2,2,4)
boxplot(tbl{:, 4}, c)
title('STDwav');
xlabel('Computed labels');
ylabel('Parameter value');
%% Evaluation of the training model 

tot_segm = sum(confusion_matrix,'all');
false_corrupted_percentage = (confusion_matrix(1,2)/tot_segm)*100;
false_clean_percentage = (confusion_matrix(2,1)/tot_segm)*100;
error_rate_percentage = false_corrupted_percentage + false_clean_percentage;
accuracy_percentage = 100-error_rate_percentage;
precision_percentage = (confusion_matrix(2,2)/(confusion_matrix(2,2)+confusion_matrix(2,1)))*100;