clear all

%% Ideal case
load('Xt1_1.mat')
load('Xt2_1.mat')
load('Xt3_1.mat')

minsize = min(size(Xt1_1,2),min(size(Xt2_1,2),size(Xt3_1,2)));
Xt1_1(:,(minsize+1):end)=[];Xt2_1(:,(minsize+1):end)=[];Xt3_1(:,(minsize+1):end)=[];

for index1 = 1:size(Xt1_1,1)
    Xt1_1(index1,:) = Xt1_1(index1,:)-mean(Xt1_1(index1,:));
    Xt2_1(index1,:) = Xt2_1(index1,:)-mean(Xt2_1(index1,:));
    Xt3_1(index1,:) = Xt3_1(index1,:)-mean(Xt3_1(index1,:));
end

Xt1 = [Xt1_1;Xt2_1;Xt3_1];

[U,S,V] = svd(Xt1);

Y = transpose(U)*Xt1;
A1 = Y; 

sig = diag(S);

for index1 = 1:length(sig)
    A2(index1,1) = sum(sig(1:index1).^2)/sum(sig.^2); % Shape: 6x1 double
end
figure(2); hold on;
semilogy(1:6,A2,'*:', 'MarkerSize',15);
xlabel('# of singular values added')
ylabel('energy')
title('Log scaled energies')
set(gca,'xtick',[1,2,3,4,5,6])

nimpo = 2;


A3 = U(:,1:nimpo)*S(1:nimpo,1:nimpo)*(V(:,1:nimpo)'); % Shape: 6x226 double


%%% Plot the approximations from SVD for report (not for autograder)
% Plot each rank-n approximation (up to one more than A12) to observe the
% convergence

figure(1); hold on
for index2 = 1:6
    AA = U(:,1:index2)*S(1:index2,1:index2)*(V(:,1:index2)');
    plot(1:length(AA(1,:)),AA(1,:),'LineWidth',7-index2)
end
    xlabel('t')
    ylabel('Centered x-direction')
    title('n-rank approximations of the motion of mass on a spring in one direction')
    set(gca,'xtick',[])
    legend('n=1','n=2','n=3','n=4','n=5','n=6')
    xlim([1,length(AA(index2,:))])
    figure(4)
for index1 = 1:6
    subplot(2,3,index1)
    Tempinfo = U(:,1:index1)*S(1:index1,1:index1)*(V(:,1:index1)');
    %plot(1:length(Tempinfo(2,:)),Tempinfo(6,:))
    plot3(Tempinfo(6,:),Tempinfo(1,:),Tempinfo(4,:))
    xlabel('x-direction')
ylabel('y-direction')
zlabel('height')
title([num2str(index1),'-rank approximation'])

end
sgtitle('Motion of mass on spring')
figure(3)
subplot(1,3,1)
plot(1:length(Xt1_1(1,:)),Xt1_1(1,:))
    xlabel('t')
    ylabel('Centered x-direction')
    title('The data')
    set(gca,'xtick',[])
    xlim([1,length(Xt1_1(1,:))])
for index1 = [1,6]
    subplot(1,3,floor(index1/6)+2)
    AAA = U(:,index1)*S(index1,index1)*(V(:,index1)');
    plot(1:length(AAA(1,:)),AAA(1,:))
    xlabel('t')
    ylabel('Centered x-direction')
    title('Approximation using the highest singular value')
    set(gca,'xtick',[])
    xlim([1,length(AA(index2,:))])
end

figure(5)
subplot(1,2,1)
Tempinfo = U(:,1:2)*S(1:2,1:2)*(V(:,1:2)');
plot3(Tempinfo(6,:),Tempinfo(1,:),Tempinfo(4,:),'LineWidth',3); hold on
plot3(Tempinfo(6,2),Tempinfo(1,2),Tempinfo(4,2),'.g','MarkerSize',30)
plot3(Tempinfo(6,end),Tempinfo(1,end),Tempinfo(4,end),'.r','MarkerSize',30)
xlabel('x-direction')
ylabel('y-direction')
zlabel('height')
title('Motion of mass on a spring')
subplot(1,2,2)
Tempinfo = U(:,1:2)*S(1:2,1:2)*(V(:,1:2)');
plot3(Tempinfo(6,:),Tempinfo(1,:),Tempinfo(4,:),'LineWidth',3); hold on
plot3(Tempinfo(6,2),Tempinfo(1,2),Tempinfo(4,2),'.g','MarkerSize',30)
plot3(Tempinfo(6,end),Tempinfo(1,end),Tempinfo(4,end),'.r','MarkerSize',30)
xlabel('x-direction')
ylabel('y-direction')
zlabel('height')
title('Motion of mass on a spring')
clearvars -except A1 A2 A3
%% Test 2
load('Xt1_2.mat')
load('Xt2_2.mat')
load('Xt3_2.mat')

minsize = min(size(Xt1_2,2),min(size(Xt2_2,2),size(Xt3_2,2)));
Xt1_2(:,(minsize+1):end)=[];Xt2_2(:,(minsize+1):end)=[];Xt3_2(:,(minsize+1):end)=[];

for index1 = 1:size(Xt1_2,1)
    Xt1_2(index1,:) = Xt1_2(index1,:)-mean(Xt1_2(index1,:));
    Xt2_2(index1,:) = Xt2_2(index1,:)-mean(Xt2_2(index1,:));
    Xt3_2(index1,:) = Xt3_2(index1,:)-mean(Xt3_2(index1,:));
end

Xt2 = [Xt1_2;Xt2_2;Xt3_2];

[U,S,V] = svd(Xt2,'econ');

Y = transpose(U)*Xt2;
A4 = Y; 

sig = diag(S);
for index1 = 1:length(sig)
    A5(index1,1) = sum(sig(1:index1).^2)/sum(sig.^2); % Shape: 6x1 double
end

nimpo = 3;


A6 = U(:,1:nimpo)*S(1:nimpo,1:nimpo)*(V(:,1:nimpo)'); % Shape: 6x226 double

clearvars -except A1 A2 A3 A4 A5 A6
%% Test 3
load('Xt1_3.mat')
load('Xt2_3.mat')
load('Xt3_3.mat')

minsize = min(size(Xt1_3,2),min(size(Xt2_3,2),size(Xt3_3,2)));
Xt1_3(:,(minsize+1):end)=[];Xt2_3(:,(minsize+1):end)=[];Xt3_3(:,(minsize+1):end)=[];

for index1 = 1:size(Xt1_3,1)
    Xt1_3(index1,:) = Xt1_3(index1,:)-mean(Xt1_3(index1,:));
    Xt2_3(index1,:) = Xt2_3(index1,:)-mean(Xt2_3(index1,:));
    Xt3_3(index1,:) = Xt3_3(index1,:)-mean(Xt3_3(index1,:));
end

Xt3 = [Xt1_3;Xt2_3;Xt3_3];

[U,S,V] = svd(Xt3,'econ');

Y = transpose(U)*Xt3;
A7 = Y; 

sig = diag(S);
for index1 = 1:length(sig)
    A8(index1,1) = sum(sig(1:index1).^2)/sum(sig.^2); % Shape: 6x1 double
end

nimpo = 3;


A9 = U(:,1:nimpo)*S(1:nimpo,1:nimpo)*(V(:,1:nimpo)'); % Shape: 6x226 double

clearvars -except A1 A2 A3 A4 A5 A6 A7 A8 A9
%% Test 4
load('Xt1_4.mat')
load('Xt2_4.mat')
load('Xt3_4.mat')

minsize = min(size(Xt1_4,2),min(size(Xt2_4,2),size(Xt3_4,2)));
Xt1_4(:,(minsize+1):end)=[];Xt2_4(:,(minsize+1):end)=[];Xt3_4(:,(minsize+1):end)=[];

for index1 = 1:size(Xt1_4,1)
    Xt1_4(index1,:) = Xt1_4(index1,:)-mean(Xt1_4(index1,:));
    Xt2_4(index1,:) = Xt2_4(index1,:)-mean(Xt2_4(index1,:));
    Xt3_4(index1,:) = Xt3_4(index1,:)-mean(Xt3_4(index1,:));
end

Xt4 = [Xt1_4;Xt2_4;Xt3_4];

[U,S,V] = svd(Xt4,'econ');

Y = transpose(U)*Xt4;
A10 = Y; 

sig = diag(S);
for index1 = 1:length(sig)
    A11(index1,1) = sum(sig(1:index1).^2)/sum(sig.^2); % Shape: 6x1 double
end

nimpo = 4;


A12 = U(:,1:nimpo)*S(1:nimpo,1:nimpo)*(V(:,1:nimpo)'); % Shape: 6x226 double

Motions = {A6,A9,A12};
titles = {'noisy','off center','Rotation'};
    figure(8)
for index1 = 1:3
    subplot(1,3,index1)
    plot3(Motions{index1}(6,:),Motions{index1}(1,:),Motions{index1}(4,:),'LineWidth',3); hold on
    plot3(Motions{index1}(6,2),Motions{index1}(1,2),Motions{index1}(4,2),'.g','MarkerSize',30)
    plot3(Motions{index1}(6,end),Motions{index1}(1,end),Motions{index1}(4,end),'.r','MarkerSize',30)
    xlabel('x-direction')
ylabel('y-direction')
zlabel('height')
title(titles{index1})
end
sgtitle('Motion of mass on spring for other videos')
