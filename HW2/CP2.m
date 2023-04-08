load 'CP2_SoundClip.mat'

Fs = 44100;
S = y'; 
w = length(y)/4;

S1 = S((1-1)*w+1:1*w); 
S2 = S((2-1)*w+1:2*w); 
S3 = S((3-1)*w+1:3*w); 
S4 = S((4-1)*w+1:4*w); 

L = length(S1)/Fs; 
n = length(S1);  
t = [0:1/Fs:L - 1/Fs]; 
tau = 0:0.1:L;
k = 2*pi*(1/L/2)*[0:n/2-1 -n/2:-1]; 
ks = fftshift(k);

[A1,i1] =  FilteredSgtcalc(t, k, tau, S1);   % Shape:  484560x110 double
[A2,i2] =  FilteredSgtcalc(t, k, tau, S2);   % Shape:  484560x110 double
[A3,i3] =  FilteredSgtcalc(t, k, tau, S3);   % Shape:  484560x110 double
[A4,i4] =  FilteredSgtcalc(t, k, tau, S4);   % Shape:  484560x110 double


% figure(1)
% subplot(121)
% pcolor(tau,ks,log(A1+1))
% shading interp
% set(gca,'Fontsize',16)
% colormap(hot)
% colorbar
% xlabel('time (t)'), ylabel('frequency (k)')
% title('Spectogram first segment.')
% subplot(122)
% pcolor(tau,ks,log(A1+1))
% shading interp
% set(gca,'ylim',[0 600],'Fontsize',16)
% xlabel('time (t)'), ylabel('frequency (k)')
% title('Spectogram first segment over frequencies 0-600.')
% figure(2)
% pcolor(tau,ks,log(A2+1))
% shading interp
% set(gca,'ylim',[0 600],'Fontsize',16)
% colormap(hot)
% colorbar
% xlabel('time (t)'), ylabel('frequency (k)')
% title('Spectogram second segment over frequencies 0-600.')
% figure(3)
% pcolor(tau,ks,log(A3+1))
% shading interp
% set(gca,'ylim',[0 600],'Fontsize',16)
% colormap(hot)
% colorbar
% xlabel('time (t)'), ylabel('frequency (k)')
% title('Spectogram third segment over frequencies 0-600.')
% figure(4)
% pcolor(tau,ks,log(A4+1))
% shading interp
% set(gca,'ylim',[0 600],'Fontsize',16)
% colormap(hot)
% colorbar
% xlabel('time (t)'), ylabel('frequency (k)')
% title('Spectogram fourth segment over frequencies 0-600')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the bassline
S = y'; 
L = length(S)/Fs; 
n = length(S); 
t = [0:1/Fs:L - 1/Fs]; 
k = (1/L)*[0:n/2-1 -n/2:-1];

% Take the Fourier transform of S, and in freq. space isolate all freqs. 
% (in absolute value) that you determine should be part of the baseline 
% according to spectrogram (or also just by listening); that is, all points
% in the transformed function not within the frequency range you determined
% should be set to zero (kind of like a Shannon filter, but simpler than
% what we did in lecture).
% You may have to do this part a few times with different thresholds to get
% it right.

Sf = fft(S);

spacing = 25000;
Filter1 = (k>= 0 & k<= 350);
Sff1 = Sf.*Filter1;

Filter2 = k<=1200 & k>= 350 ;
Sff2 = Sf.*Filter2;



% After thresholding the transformed function, take the inverse transform
% of the thresholded function and save it as A5.
fS1 = ifft(Sff1);

A5 =   fS1.'; %Shape:  1938240x1 double

fS2 = ifft(Sff2);

A6 =  fS2.';

%p9 = audioplayer(A6, Fs); playblocking(p9);

%Play sound (not for autograder)

%Plot the amplitude S over time (for the report, not for the autograder)

function [Sgt_spec, mindex] = FilteredSgtcalc(t, k, tau, S)

    Sgt_spec = zeros(length(k),length(tau));
    a = 400; 
    range = [1:1800]; 

    Fs = 44100;
    L = length(S)/Fs;
    for index1 = 1:length(tau)
        g = exp(-a*(t-tau(index1)).^2);
        Sg = g.*S;
        Sgt = fft(Sg);
        [mval, mindex] = max(abs(Sgt(range)));
        FilterSgt = exp(-((abs(k)-k(mindex)).^2)/L);
        fSgt = Sgt.*FilterSgt;
        Sgt_spec(:,index1) = fftshift(abs(fSgt));
    end
end
