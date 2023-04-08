% Clean workspace
 clear all; close all;

 % Preamble from CP1_sample.m.  The autograder will take care of any files
 % needed for this step.  Please only submit solution file.

 load('Kraken.mat') %load data matrix
 L = 10; % spatial domain
 n = 64; % Fourier modes
 x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x; %create 3D axis arrays with 64 points
 k = (2*pi/(2*L))*[0:(n/2 - 1), -n/2:-1]; %create frequency array and rescale them to be 2pi periodic 
 
 %Create 3D grids for both spatial domain and frequency domain 
 [X,Y,Z] = meshgrid(x,y,z);
 [kx,ky,kz] = meshgrid(k,k,k);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

 % sum all realizations in frequency space
 % after your for loop ends, save the sum as variable A1
 Unt = zeros(round(size(Kraken,1)^(1/3)),round(size(Kraken,1)^(1/3)),round(size(Kraken,1)^(1/3)));
 for index4 = 1:size(Kraken,2)
    ReshapedKrak(:,:,:,index4) = reshape(Kraken(:, index4), n, n, n);
    Ut(:, :, :,index4) = fftn(ReshapedKrak(:,:,:,index4)); % We need to reshape our data into a tensor, which represents a cube of Fourier modes in x-y-z space
    Unt(:, :, :) = Unt(:, :, :) + Ut(:, :, :,index4);
 end
 A1 = Unt;


 % Average the sum over the 49 realizations (i.e., A1/49) and save as A2

 A2 = A1/49;

 
 % find the peak frequencies in x, y, and z directions; i.e., find the
 % max in each direction of the averaged sum A2.
 % save these variables as A3, A4, and A5
locmax = find(max(max(max(abs(A2))))==abs(A2));

kxmax = ky(locmax); kymax = kx(locmax); kzmax = kz(locmax);
 A3 = kxmax;
 A4 = kymax;
 A5 = kzmax;
 
%create an appropriate Gaussian filter and save it as A6
Gaussfilter = exp(-((ky-A3).^2+(kx-A4).^2+(kz-A5).^2)/L);
A6 = Gaussfilter;

% Using the peak frequencies for the filtered signal, estimate the x, y, and z coordinates of the Kraken over time and save as A7, A8, A9
for index4 = 1:size(ReshapedKrak,4)
    Filteredfreq = Ut(:,:,:,index4).*Gaussfilter;
    FilteredKrak( :, :, :,index4) = ifftn(Filteredfreq);
    locmax = max(max(max(abs(FilteredKrak( :, :, :,index4)))))==abs(FilteredKrak( :, :, :,index4));
    xmaxreal(index4) = X(locmax); ymaxreal(index4) = Y(locmax); zmaxreal(index4) = Z(locmax);
end
A7 = ymaxreal;
A8 = xmaxreal;
A9 = zmaxreal;

%% plotting
% Plot the location in x-y-z space over time for your report (not for the autograder)
figure();
plot3(A7,A8,zmaxreal); hold on;
scatter3(A7(1),A8(1),zmaxreal(1),'green', 'filled')
scatter3(A7(end),A8(end),zmaxreal(end),'red', 'filled')
xlabel('x'); ylabel('y'); zlabel('z')
xlim([-10,10]);ylim([-10,10]);zlim([-10,10])
title('Path of the Kraken')
% Plot the projection onto the x-y plane for your reprot (not for the autograder)
figure();
plot(A7,A8); hold on;
scatter(A7(1),A8(1),'green', 'filled')
scatter(A7(end),A8(end),'red', 'filled')
plot(0*(-10:10),-10:10,'--k')
xlabel('x'); ylabel('y')
xlim([-10,10]);ylim([-10,10]);
title('Path of the Kraken')
% Save a table of the x-y-z coordinates for your report (not for the autograder)


