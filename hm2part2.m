%% PART 2: Piano
clear all; close all; clc
tr_piano = 16; % record time in seconds
y=audioread('music1.wav'); 

p = y.';
pf = fft(p);  %% FFT of signal
pfs = fftshift(pf); %% FFT shifted of signal

L=tr_piano; n=701440; % set up time domain and fourier modes
t2=linspace(0,L,n+1); t=t2(1:n);  % setup time vector based on signal data, periodic
k=(2*pi/L)*[0:n/2-1 -n/2:-1];   % setup freq vector
ks=fftshift(k); % flip so ks has ordered modes with 0 at center

% plot signal and fft of signal
figure(1)
sgtitle('Mary had a little lamb:  Piano');
%subplot(2,1,1);
Fs = length(p)/tr_piano;
plot((1:length(p))/Fs,p); 
xlim([0 16]);
set(gca,'Fontsize',[9]) ;
xlabel('Time [sec]'); ylabel('Amplitude');
title ('Audio Signal');

subplot(2,1,2);
plot(ks,abs(pfs)/max(abs(pfs))); 
xlabel('freq (\omega)'); 
ylabel('fft'); 
title('Fourier Transform');

figure(2)
subplot(2,1,1);
plot(ks,abs(pfs)/max(abs(pfs))); 
xlabel('freq (\omega)'); 
ylabel('fft'); 
title('Fourier Transform close up near peaks');
xlim([0 5000]);

%find max value in three areas with peaks, this is the lowest note
[M1,I1]= max(pfs);
max1_wave = ks(1, I1);% first peak


%zero out anything above 2x this value (covers it's overtones in addition
%to any overtones of the three later peaks
ind = abs(ks(1,:)) >3000; ; %grab indices along second column where value greater than lowest multiple of peak1
pfs_clean = num2cell(pfs); % youre going to need to turn it into a cell array
pfs_clean(1,ind) = {0}; %set those indices to zero
pfs_clean = cell2mat(pfs_clean);


subplot(2,1,2);
plot(ks,abs(pfs_clean)/max(abs(pfs_clean))); 
xlabel('freq (\omega)'); 
ylabel('fft'); 
title('Fourier Transform close up near peaks, cleaned');
xlim([0 5000]);

% ifft and iffshift to pfs_clean data for gabor filtering
pf_clean = ifftshift(pfs_clean);
p_clean = ifft(pf_clean);



%Gaussian
figure(3)
sgtitle('Gaussian Filter  with a = 10 & b = 0.2')
pgt_spec=[]; 
tslide=[0:0.2: 16];
for ii=1:length(tslide)
    g=exp(-10*(t-tslide(ii)).^2); % Gabor 
    pg = g.*p_clean;
    pgt=fft(pg); 
    pgt_spec=[pgt_spec;  abs(fftshift(pgt))]; 
    %subplot(3,1,1), plot(t,p,'k', t,g, 'r') ;  
    %title('translation resolution= .2')
   % xlabel('time(t)', 'Fontsize', [10])
    % ylabel('Amplitude')
   % subplot(4,1,2), plot(t,pg,'k')
  %  subplot(4,1,3), plot(ks,abs(fftshift(pgt))/max(abs(pgt))) 
 
end
%
subplot(2,1,1)
pcolor(tslide,ks,pgt_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[10]) 
colormap(hot)
title('Spectrogram');
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 300])
ylim([-4000 4000])

xlabel('time(t)');
ylabel('frequency (\omega)')
%
subplot(2,1,2)
Hz = ks/(2*pi);
pcolor(tslide,Hz,pgt_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[10]) 
colormap(hot)
title('Spectrogram: close up (Hz)');
xlabel('time(t)');
ylabel('frequency (Hz)')
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 800])
ylim([100 400])


% Step Function
figure(4)
sgtitle('Shannon step filter w/ a = 0.2 and b = 0.2')
pst_spec=[]; 
tslide=[0:0.2: 16];
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 0.2/2; % Shannon filter
    ps = sh.*p_clean;
    pst=fft(ps); 
    pst_spec=[pst_spec;  abs(fftshift(pst))]; 
    %subplot(3,1,1), plot(t,p_clean,'k', t,sh, 'r') ;  
    %xlabel('time(t)', 'Fontsize', [10])
     %ylabel('Amplitude')
   % subplot(4,1,2), plot(t,pg,'k')
  %  subplot(4,1,3), plot(ks,abs(fftshift(pgt))/max(abs(pgt))) 
 
end

%
subplot(2,1,1)
pcolor(tslide,ks,pst_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[10]) 
colormap(hot)
title('Spectrogram');
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 800])
ylim([-3500 3500])
xlabel('time(t)');
ylabel('frequency (\omega)')
%
subplot(2,1,2)
Hz = ks/(2*pi);
pcolor(tslide,Hz,pst_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[10]) 
colormap(hot);
title('Spectrogram: close up (Hz)');
xlabel('time(t)');
ylabel('frequency (Hz)')
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([100 400])
ylim([200 550])


%% PART 2: Recorder %%
clear all;
tr_rec= 14; % record time in seconds
y=audioread('music2.wav'); 
r = y.';
rf = fft(r);  % FFT of signal
rfs = fftshift(rf); %% FFT shifted of signal

L=tr_rec; n=length(y) ; % set up time domain and fourier modes
t2=linspace(0,L,n+1); t=t2(1:n);  % setup time vector based on signal data, periodic
k=(2*pi/L)*[0:n/2-1 -n/2:-1];   % setup freq vector
ks=fftshift(k); % flip so ks has ordered modes with 0 at center

% plot signal and fft of signal
figure(5)
subplot(2,1,1)
sgtitle('Mary had a little lamb:  Recorder');
Fs = length(r)/tr_rec;
plot((1:length(r))/Fs,r); 
set(gca,'Fontsize',[10]) ;
xlabel('Time [sec]'); ylabel('Amplitude');
title ('Recorder:  Audio Signal', 'Fontsize', [12]);
%
subplot(2,1,2)
plot(ks,abs(rfs)/max(abs(rfs))); 
xlabel('freq (\omega)'); 
ylabel('fft'); 
title('Fourier Transform');

figure(6)
subplot(2,1,1);
plot(ks,abs(rfs)/max(abs(rfs))); 
set(gca,'Fontsize',[12]) ;
xlabel('freq (\omega)'); 
ylabel('fft'); 
title('Fourier Transform: Dominant Peaks (recorder)', 'Fontsize', [13]);

%zero out anything above 2x this value (covers it's overtones in addition
%to any overtones of the three later peaks
ind = abs(ks(1,:)) >10000; ; %grab indices along second column where value greater than lowest multiple of peak1
rfs_clean = num2cell(rfs); % youre going to need to turn it into a cell array
rfs_clean(1,ind) = {0}; %set those indices to zero
rfs_clean = cell2mat(rfs_clean);

subplot(2,1,2);
plot(ks,abs(rfs_clean)/max(abs(rfs_clean))); 
xlabel('freq (\omega)'); 
ylabel('fft'); 
title('Fourier Transform close up near peaks, cleaned');

% ifft and iffshift to pfs_clean data for gabor filtering
rf_clean = ifftshift(rfs_clean);
r_clean = ifft(rf_clean);

% Gaussian
figure(7)
sgtitle('Gaussian Filter  with a = 10 & translation = 0.2')
rgt_spec=[]; 
tslide=[0:0.2: 16];
for ii=1:length(tslide)
    g=exp(-12*(t-tslide(ii)).^2); % Gabor 
    rg = g.*r_clean;
    rgt=fft(rg); 
    rgt_spec=[rgt_spec;  abs(fftshift(rgt))]; 
    subplot(3,1,1), plot(t,r_clean,'k', t,g, 'r') ;  
    xlabel('time(t)', 'Fontsize', [10])
     ylabel('Amplitude')
   % subplot(4,1,2), plot(t,pg,'k')
  %  subplot(4,1,3), plot(ks,abs(fftshift(pgt))/max(abs(pgt))) 
 
end
%
subplot(3,1,2)
pcolor(tslide,ks,rgt_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[10]) 
colormap(hot)
title('Spectrogram');
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
ylim([-25000 25000])
caxis([0 300])
xlim([0 14])

subplot(3,1,3)
Hz = ks/(2*pi);
pcolor(tslide,Hz,rgt_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[10]) 
colormap(hot)
title('Spectrogram: close up (Hz)');
colorbar
ylabel('frequency (Hz)');
h = colorbar; ylabel(h, 'FFT of filtered signal')
ylim([400 1400])
caxis([0 400])
xlim([0 14])

% Step Function
figure(8)
sgtitle('Step filter (shannon)  w/ a = 0.2 & translation = 0.2', 'Fontsize', [14])
rst_spec=[]; 
tslide=[0:0.2: 16];
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 0.2/2; % Gabor 
    rs = sh.*r_clean;
    rst=fft(rs); 
    rst_spec=[rst_spec;  abs(fftshift(rst))]; 
end

%
subplot(2,1,1)
pcolor(tslide,ks,rst_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[12]) 
colormap(hot)
ylabel('frequency (\omega)');
title('Spectrogram');
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
ylim([-25000 25000])
caxis([0 200])
xlim([0 14])
%
subplot(2,1,2)
Hz = ks/(2*pi);
pcolor(tslide,Hz,rst_spec.'), 
shading interp
xlabel('time(t)')
set(gca,'Fontsize',[12]) 
colormap(hot);
ylabel('frequency (Hz)');
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram: close up (Hz)');
ylim([600 1400])
caxis([0 300])
xlim([0 14])

