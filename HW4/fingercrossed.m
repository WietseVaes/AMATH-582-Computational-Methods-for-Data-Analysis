%% Clean workspace
clear all; close all; clc

%% Load training data

load('CP4_training_labels.mat')
load('CP4_training_images.mat')

train_ima = zeros(size(training_images,2)^2,size(training_images,3));
for index1 = 1:size(training_images,3)
    train_ima(:,index1) = reshape(training_images(:,:,index1),[],1);
end

Trainall_DWT = dc_wavelet(train_ima);

Ntrain = size(Trainall_DWT,2);

Training_DWT1 = Trainall_DWT(:,1:(Ntrain/2));
training_labels1 = training_labels(1:(Ntrain/2));
Training_DWT2 = Trainall_DWT(:,(Ntrain/2+1):end);
training_labels2 = training_labels(:,(Ntrain/2+1):end);

[U,S,V] = svd(Training_DWT1,'econ');

load('CP4_test_labels.mat')
load('CP4_test_images.mat')

test_ima = zeros(size(test_images,2)^2,size(test_images,3));
for index1 = 1:size(test_images,3)
    test_ima(:,index1) = reshape(test_images(:,:,index1),[],1);
end

Test_DWT = dc_wavelet(test_ima);

%% Combos
Allcombos = nchoosek(0:9,2);

nc = size(Allcombos,1);
nf = size(S,1);
sucRates = zeros(nc,1);
Thresholds = zeros(nc,1);
meannums = zeros(nc,2);
for index1 = 1:nc
    sucRateold = 0;
    for index2 = 10:nf
        [threshold, sucRate, meannum1, meannum2,w] = CP4cor(index2,U,S,V,Allcombos(index1,1),Allcombos(index1,2),Training_DWT2,test_labels,training_labels);
        if sucRate > sucRateold
            sucRateold = sucRate;
            sucRates(index1) = sucRate;
            Thresholds(index1) = threshold;
            meannums(index1,:) = [meannum1, meannum2];
            ws{index1} = w;
            Features(index1) = length(w);
            if sucRate >=.99
                break
            end
        end
    end
    index1/nc
end

%% Using thresholds
sortedsucRates = flip(sort(sucRates));
Resultcounter = zeros(size(Test_DWT,2),10);
for index1 = 1:nc
    clear feature w threshold Ufeat
    feature = Features(index1);
    w = ws{index1};
    threshold = Thresholds(index1);
    Ufeat = U(:,1:feature);


    testMat = Ufeat'*Test_DWT;
    w_projected_test = w'*testMat;
    Results = w_projected_test >= threshold;

    if meannums(index1,1)>meannums(index1,1)
        for index2 = 1:sum(Results)
        Resultcounter(Results,Allcombos(index1,1)+1) = Resultcounter(Results,Allcombos(index1,1)+1)+1;
        Resultcounter(~Results,Allcombos(index1,2)+1) = Resultcounter(~Results,Allcombos(index1,2)+1)+1;
        end
    else
        Resultcounter(Results,Allcombos(index1,2)+1) = Resultcounter(Results,Allcombos(index1,2)+1)+1;
        Resultcounter(~Results,Allcombos(index1,1)+1) = Resultcounter(~Results,Allcombos(index1,1)+1)+1;
    end
end
Resultlabels = zeros(size(Resultcounter,1),1);
posind = 1:10;
for index1 = 1:size(Resultcounter,1)
    maxind = posind(Resultcounter(index1,:)==max(Resultcounter(index1,:)));
    if length(maxind)>1
        Resultlabels(index1) = -1;
    else
        Resultlabels(index1) = maxind-1;
    end
end



function dcData = dc_wavelet(dcfile)
[m,n] = size(dcfile);
pxl = sqrt(m);
nw = m/4; % wavelet resolution cus downsampling
dcData = zeros(nw,n);

for k = 1:n
    X = im2double(reshape(dcfile(:,k),pxl,pxl));
    [~,cH,cV,~]=dwt2(X,'haar'); % only want horizontal and vertical
    cod_cH1 = rescale(abs(cH)); % horizontal rescaled
    cod_cV1 = rescale(abs(cV)); % vertical rescaled
    cod_edge = cod_cH1+cod_cV1; % edge detection
    dcData(:,k) = reshape(cod_edge,nw,1);
end
end