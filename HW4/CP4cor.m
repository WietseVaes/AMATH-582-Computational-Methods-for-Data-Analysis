function [threshold, sucRate, meannum1, meannum2,w] = CP4cor(feature,U,S,V,num1,num2,Test_DWT,test_labels,training_labels)

training_projected_temp = S*V';
training_projected = training_projected_temp(1:feature,:);
U = U(:,1:feature);

zero_label = find(training_labels==num1);
one_label = find(training_labels==num2);

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

[V2, D] = eig(Varb,Varw); % linear disciminant analysis; i.e., generalized eval. prob.
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);

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

threshold = (sortzeros(t1) + sortones(t2))/2; 
meannum1 = mean(w_projected_zeros);
meannum2 = mean(w_projected_ones);

Test_Ima_zo = Test_DWT(:,test_labels == num2 | test_labels == num1);
Labels_Test = test_labels(test_labels == num2 | test_labels == num1)';

testMat = U'*Test_Ima_zo;
w_projected_test = w'*testMat;
Results = w_projected_test >= threshold;
if meannum1 > meannum2
    Testingggg = (Labels_Test-num2)/(num1-num2);
else
    Testingggg = (Labels_Test-num1)/(num2-num1);
end
err = abs(Results - Testingggg);
errNum = sum(err);
sucRate = 1 - errNum/length(Labels_Test);
plot(w_projected_test(Results&err==0),zeros(length(w_projected_test(Results&err==0)),1),'bo','Filled'); hold on;
plot(w_projected_test(~Results&err==0),ones(length(w_projected_test(~Results&err==0)),1),'ro','Filled')
plot(w_projected_test(Results&err==1),ones(length(w_projected_test(Results&err==1)),1),'r*')
plot(w_projected_test(~Results&err==1),zeros(length(w_projected_test(~Results&err==1)),1),'b*')
title('Projected test cases'); legend('Ones','Zeros','Incorrect guess')
ylim([-.5,1.5]);set(gca,'XTick',[], 'YTick', []);xline(threshold,'k:','LineWidth',2)

end

% Approxi = U(:,1:A1)*S(1:A1,1:A1)*V(:,1:A1)';
% indexes = [16,2];
% for index1 = 1:2
%     indexzo = indexes(index1);
%     figure(index1+1)
%     subplot(1,3,1)
%     imshow(training_images(:,:,indexzo))
%     title('Original')
%     subplot(1,3,2)
%     imshow(reshape(Training_DWT(:,indexzo),14,14))
%     title('Edges')
%     subplot(1,3,3)
%     imshow(reshape(Approxi(:,indexzo),14,14))
%     title("15-rank approximation of edges")
%     sgtitle(['Images through processing for ', num2str(index1-1)])
%     subplot(1,3,1)
%     set(gca,'YAxisLocation','right')
%     ylabel('\Downarrow', 'FontSize', 50)
%     subplot(1,3,2)
%     set(gca,'YAxisLocation','right')
%     ylabel('\Downarrow', 'FontSize', 50)
% 
% end

% for index1 = 1:length(sig)
%     ener(index1) = sum(sig(1:index1).^2)/sum(sig.^2);
% end
% figure(1);
% semilogy(1:120,ener(1:120));
% title('Ratio of energies of singular values.')
% xlabel('Rank n')
% ylabel('Energy ratio')
