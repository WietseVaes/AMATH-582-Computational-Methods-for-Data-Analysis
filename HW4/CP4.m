%%% Clean workspace
clear all; close all; clc

%%% Load training data

load('CP4_training_labels.mat')
load('CP4_training_images.mat')


train_ima = zeros(size(training_images,2)^2,size(training_images,3));
for index1 = 1:size(training_images,3)
    train_ima(:,index1) = reshape(training_images(:,:,index1),[],1);
end


Training_DWT = dc_wavelet(train_ima);
%load('Training_DWT.mat')

[U,S,V] = svd(Training_DWT,'econ');

sig = diag(S);

for index1 = 1:length(sig)
    ener(index1) = sum(sig(1:index1).^2)/sum(sig.^2); % Shape: 6x1 double
end

%semilogy(1:20,ener(1:20));

A1 = 15; % The number of PCA modes we are going to project onto.

training_projected_temp = S*V';
training_projected = training_projected_temp(1:A1,:);
U = U(:,1:A1);

A2 = U;

zero_label = find(training_labels==0);
one_label = find(training_labels==1);

minlength = min(length(zero_label),length(one_label));
zero_label(minlength+1:end) = [];
one_label(minlength+1:end) = [];

zeros_training = training_projected(:,zero_label);
ones_training = training_projected(:,one_label);

mz = mean(zeros_training,2);
mo = mean(ones_training,2);

Varw = 0;
n=size(zeros_training,2);
for index1 = 1:n
    Varw = Varw + (zeros_training(:,index1)-mz)*(zeros_training(:,index1)-mz)';
end
n=size(ones_training,2);
for index1 = 1:n
    Varw = Varw + (ones_training(:,index1)-mo)*(ones_training(:,index1)-mo)';
end

Varb = (mz-mo)*(mz-mo)';
A3 = Varw;
A4 = Varb;

[V2, D] = eig(Varb,Varw); % linear disciminant analysis; i.e., generalized eval. prob.
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

A5 = w;

w_projected_zeros = w'*zeros_training;
w_projected_ones = w'*ones_training;

% Find the threshold value just like in Week7_LDA.m.  Save it as A6

sortzeros = sort(w_projected_zeros);
sortones = sort(w_projected_ones);
t1 = length(sortzeros); % start on the right
t2 = 1; % start on the left
while sortzeros(t1) > sortones(t2) 
    t1 = t1 - 1;
    t2 = t2 + 1;
end
% ^ go past each other
threshold = (sortzeros(t1) + sortones(t2))/2; % get the midpoint


A6 = threshold;

load('CP4_test_labels.mat')
load('CP4_test_images.mat')

test_ima = zeros(size(test_images,2)^2,size(test_images,3));
for index1 = 1:size(test_images,3)
    test_ima(:,index1) = reshape(test_images(:,:,index1),[],1);
end

Test_Ima_zo = test_ima(:,test_labels == 1 | test_labels == 0);
Labels_Test = test_labels(test_labels == 1 | test_labels == 0)';

Test_DWT= dc_wavelet(Test_Ima_zo);
%load('Test_DWT.mat')
testMat = U'*Test_DWT;
w_projected_test = w'*testMat;

A7 = (w_projected_test >= threshold)*1;

Results = w_projected_test >= threshold;
err = abs(Results - Labels_Test);
errNum = sum(err);
sucRate = 1 - errNum/length(Labels_Test)
save('vars','A1','A2','A3','A4','A5','A6','A7')
%%% Put any helper functions here
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
