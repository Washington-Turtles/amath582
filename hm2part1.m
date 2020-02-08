%% Part 1, Exploring Gabor window sizes %%
clear all; close all; clc
load handel
%%
v = y'/2; % new vector, transpose of y/2
vf = fft(v);  %% FFT of signal
vfs = fftshift(vf); %% FFT shifted of signal

L=9; n=73113; % set up time domain and fourier modes
t2=linspace(0,L,n+1); t= t2(1:n); % setup time vector based on signal data, periodic
k=(2*pi/L)*[0:n/2 -n/2:-1];  % setup freq vector for odd fourier mode (symmetrical)
ks=fftshift(k); % flip so ks has ordered modes with 0 at center

% plot signal and fft of signal
figure(1)
sgtitle('Signal of Interest v(t)');
subplot(2,1,1);
plot((1:length(v))/Fs,v); 
xlim([0 9]);
set(gca,'Fontsize',[10]) ;
xlabel('Time [sec]'); 
ylabel('Amplitude'); 
title('v(t)');
subplot(2,1,2);
plot(ks,abs(vfs)/max(abs(vfs))); 
xlabel('frequency (\omega)'); 
ylabel('fft(v)'); 
title('Fourier Transform - v(t)');


%%  Gabor Method with Gaussian fIlter
figure(2) 
sgtitle('Gaussian Filter with varying widths')
width=[20 2 0.2];  % using 3 different gabor widths

for j=1:3
    g=exp(-width(j)*(t-4.5).^2);  % gabor window is gaussian filter
    subplot(3,3,(j))
    plot(t, v,'k', 'Linewidth', [.25]), hold on 
    plot(t,g,'r','Linewidth',[1]) 
    set(gca,'Fontsize',[9]) 
    xlim([0 9]);
    ylabel('v(t), g(t)')
    xlabel('time (t)')
    title(['a = ', num2str(width(j))])
    vg = g.*v;
    vgf = fft(vg);
    vgfs_all(j , :) = abs(fftshift(vgf))/max(abs(vgf));  %% FFT of signal
    
    subplot(3,3, j+3), plot(t,vg,'k') 
    xlim([0 9]);
    set(gca,'Fontsize',[9]) 
    title('filtered signal (vg)')
    ylabel('v(t)g(t)'), xlabel('time (t)')
    
    
    subplot(3,3,(j)+6);
    title('Fourier Transform');
    plot(ks, vgfs_all(j, :), 'b', 'Linewidth', [1]);
    set(gca, 'Fontsize', [9]);
    ylabel('fft(vg)');
    xlabel('frequency (\omega)');
   
end
 
%  Sliding Gabor Window (Gaussian), wide window, coarse translation
figure(3)
sgtitle('Gaussian Filter  with a = 0.2 (Wide Width)')
vgt_spec=[]; 
tslide=[0:3: 9];
for ii=1:length(tslide)
    g=exp(-.2*(t-tslide(ii)).^2); % Gabor 
    vg = g.*v;
    vgt=fft(vg); 
    vgt_spec=[vgt_spec;  abs(fftshift(vgt))]; 
    title('Coarse Translation (b = 3)')
    subplot(2,2,1), plot(t,v,'k', t,g, 'r') ;  hold on;
    xlabel('time(t)', 'Fontsize', [10])
     ylabel('Amplitude')
end

subplot(2,2,3)
pcolor(tslide,ks,vgt_spec.'), 
shading interp
title('Spectrogram (coarse)')
xlabel('time(t)')
ylabel('frequency (\omega)');
set(gca,'Fontsize',[10]) 
colormap(hot)
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 200])


% Gaussian (wide window, fine translation)
vgt_spec=[]; 
tslide=[0:.5: 9];
for ii=1:length(tslide)
    g=exp(-.2*(t-tslide(ii)).^2); % Gabor 
    vg = g.*v;
    vgt=fft(vg); 
    vgt_spec=[vgt_spec;  abs(fftshift(vgt))]; 
    subplot(2,2,2), plot(t,v,'k', t,g, 'r') ;  hold on;
        title('Fine Translation (b = 0.5)')
     xlabel('time(t)', 'Fontsize', [10]);
     ylabel('Amplitude')
end

subplot(2,2,4)
pcolor(tslide,ks,vgt_spec.'), 
shading interp
 xlabel('time(t)', 'Fontsize', [10]);
ylabel('frequency (\omega)')
set(gca,'Fontsize',[10]) 
colormap(hot)
colorbar
title('Spectrogram (fine)');
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 200])

%  Sliding Gabor Window (Gaussian) - narrow window, coarse translation
figure(4)
sgtitle('Gaussian Filter  with a = 15 (Narrow Width)')
vgt_spec=[]; 
tslide=[0:1: 9];
for ii=1:length(tslide)
    g=exp(-15*(t-tslide(ii)).^2); % Gabor 
    vg = g.*v;
    vgt=fft(vg); 
    vgt_spec=[vgt_spec;  abs(fftshift(vgt))]; 
    subplot(2,2,1), plot(t,v,'k', t,g, 'r') ;  hold on;
      xlabel('time(t)', 'Fontsize', [10]);
     ylabel('Amplitude')
     title('Coarse Translation (b = 1)')
   
end

subplot(2,2,3)
pcolor(tslide,ks,vgt_spec.'), 
shading interp
xlabel('time(t)', 'Fontsize', [10]);
ylabel('frequency (\omega)')
colormap(hot)
colorbar
title('Spectrogram (coarse)');
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 200])

% Gaussian (narrow window, fine translation)
vgt_spec=[]; 
tslide=[0:.1: 9];
for ii=1:length(tslide)
    g=exp(-15*(t-tslide(ii)).^2); % Gabor 
    vg = g.*v;
    vgt=fft(vg); 
    vgt_spec=[vgt_spec;  abs(fftshift(vgt))]; 
    subplot(2,2,2), plot(t,v,'k', t,g, 'r') ;  hold on;
   title('Fine Translation (b = 0.1)')
    xlabel('time(t)', 'Fontsize', [10]);
     ylabel('Amplitude')
end

subplot(2,2,4)
pcolor(tslide,ks,vgt_spec.'), 
shading interp
xlabel('time(t)', 'Fontsize', [10]);
ylabel('frequencyuency (\omega)')
colormap(hot)
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram (fine)');
caxis([0 200])

%% best fit Gabot with Gaussian
figure(5)
vgt_spec=[]; 
tslide=[0:.3: 9];
for ii=1:length(tslide)
    g=exp(-9*(t-tslide(ii)).^2); % Gabor 
    vg = g.*v;
    vgt=fft(vg); 
    vgt_spec=[vgt_spec;  abs(fftshift(vgt))]; 
    subplot(2,1,1), plot(t,v,'k', t,g, 'r') ;  hold on;
   title('/sigma = 10 & translation = 0.3')
    xlabel('time(t)', 'Fontsize', [10]);
     ylabel('Amplitude')
end

subplot(2,1,2)
pcolor(tslide,ks,vgt_spec.'), 
shading interp
xlabel('time(t)', 'Fontsize', [10]);
ylabel('frequency (\omega)')
colormap(hot)
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram (fine)');
caxis([0 200])

%  Gabor Method with mexican hat filter
figure(6) 
sgtitle('Mexican Hat Filter of varying widths')
width=[1 .5 0.2];  % using 3 different gabor widths

for jj=1:3
    m=(1-((t-4.5)./width(jj)).^2).*exp(-(t-4.5).^2/(2.*width(jj).^2));  % gabor window is gaussian filter
    subplot(3,3,(jj))
    plot(t, v,'k'), hold on 
    plot(t,m,'r','Linewidth',[2]) 
    set(gca,'Fontsize',[9]) 
    xlim([0 9]);
    ylabel('v(t), m(t)')
    xlabel('time (t)')
    title( width(jj))
    vm = m.*v;
    vmf = fft(vm);
    vmfs_all(jj, :) = abs(fftshift(vmf))/max(abs(vmf));  %% FFT of signal
    
    subplot(3,3, jj+3), plot(t,vm,'k') 
    set(gca,'Fontsize',[9]) 
    xlim([0 9]);
    title('filtered signal (vm)')
    ylabel('v(t)m(t)'), xlabel('time (t)')
    
    
    subplot(3,3,(jj)+6);
    plot(ks, vmfs_all(jj, :), 'b', 'Linewidth', [1]);
    title('Fourier Transform (vm)')
    set(gca, 'Fontsize', [9]);
    ylabel('fft(vm)');
    xlabel('frequency (\omega)');
   
end


% gabor window (mexican hat) wide window coarse translation
figure(7)
sgtitle('Mexican Hat Filter  with a = 1.5 (wide width)')
vmt_spec=[]; 
tslide=0:3:9;
for ii=1:length(tslide)
    m=(1-((t-tslide(ii))./1.5).^2).*exp(-(t-tslide(ii)).^2/(2.*1.5.^2)); % Gabor 
    vm = m.*v;
    vmt=fft(vm); 
    vmt_spec=[vmt_spec;  abs(fftshift(vmt))]; 
    subplot(2,2,1), plot(t,v,'k', t,m, 'r'); hold on;
         ylabel('Amplitude');
     xlabel('time(t)')
    title('Coarse Translation (b = 3)')
end
%
subplot(2,2,3),
pcolor(tslide,ks,vmt_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
title('Spectrogram (coarse)');
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 200])
ylabel('frequency (\omega)')
xlabel('time(t)')

% mhat, wide window, fine translation
vmt_spec=[]; 
tslide=0:.5:9;
for ii=1:length(tslide)
    m=(1-((t-tslide(ii))./1.5).^2).*exp(-(t-tslide(ii)).^2/(2.*1.5.^2)); % Gabor 
    vm = m.*v;
    vmt=fft(vm); 
    vmt_spec=[vmt_spec;  abs(fftshift(vmt))]; 
    subplot(2,2,2), plot(t,v,'k', t,m, 'r'); hold on;
     title('Fine Translation (b = 1.5)')
     ylabel('Amplitude');
     xlabel('time(t)')
end

subplot(2,2,4),
pcolor(tslide,ks,vmt_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram (fine)');
caxis([0 200])
ylabel('frequency (\omega)')
xlabel('time(t)')

%% narrow window coarse translation
figure(8)
sgtitle('Mexican Hat Filter  with a = 0.2 (narrow width)')
vmt_spec=[]; 
tslide=0:1.5:9;
for ii=1:length(tslide)
    m=(1-((t-tslide(ii))./0.2).^2).*exp(-(t-tslide(ii)).^2/(2.*0.2.^2)); % Gabor 
    vm = m.*v;
    vmt=fft(vm); 
    vmt_spec=[vmt_spec;  abs(fftshift(vmt))]; 
    subplot(2,2,1), plot(t,v,'k', t,m, 'r'); hold on;
      title('Coarse Translation (b = 1.5)')
       ylabel('Amplitude');
     xlabel('time(t)')
end
subplot(2,2,3),
pcolor(tslide,ks,vmt_spec.'), 
title('Spectrogram (coarse)');
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
caxis([0 200])
h = colorbar; ylabel(h, 'FFT of filtered signal')
ylabel('frequency (\omega)')
xlabel('time(t)')

% mhat, wide window, fine translation
vmt_spec=[]; 
tslide=0:.3:9;
for ii=1:length(tslide)
    m=(1-((t-tslide(ii))./0.2).^2).*exp(-(t-tslide(ii)).^2/(2.*0.2.^2)); % Gabor 
    vm = m.*v;
    vmt=fft(vm); 
    vmt_spec=[vmt_spec;  abs(fftshift(vmt))]; 
    subplot(2,2,2), plot(t,v,'k', t,m, 'r'); hold on;
      title('Fine Translation (b = 0.3)')
       ylabel('Amplitude');
     xlabel('time(t)')
end
%
subplot(2,2,4),
pcolor(tslide,ks,vmt_spec.'), 
title('Spectrogram:  Mexican Hat Filter')
shading interp
set(gca,'Fontsize',[9]) 
title('Spectrogram (fine)');
colormap('hot')
colorbar
caxis([0 200])
h = colorbar; ylabel(h, 'FFT of filtered signal')
ylabel('frequency (\omega)')
xlabel('time(t)')


%%best fit
figure(9)
vmt_spec=[]; 
tslide=0:.3:9;
for ii=1:length(tslide)
    m=(1-((t-tslide(ii))./0.2).^2).*exp(-(t-tslide(ii)).^2/(2.*0.2.^2)); % Gabor 
    vm = m.*v;
    vmt=fft(vm); 
    vmt_spec=[vmt_spec;  abs(fftshift(vmt))]; 
    subplot(2,1,1), plot(t,v,'k', t,m, 'r'); hold on;
    title('Signal of Interest w// Mexican Hat filter, all translations')
  end
%
subplot(2,1,2),
pcolor(tslide,ks,vmt_spec.'), 
title('Spectrogram:  Mexican Hat Filter')
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
caxis([0 200])

%% gabor window (Step Function: Shannon) 
figure(10)
sgtitle('Shannon Step Filter with varying widths');
width=[3 1 0.2];  % using 3 different gabor widths

for jj=1:3
    sh = abs(t-4.5) <= width(jj)/2;  % gabor window with Shannon Step filter
    subplot(3,3,(jj))
    plot(t, v,'k'), hold on 
    plot(t,sh,'r','Linewidth',[1]) 
    set(gca,'Fontsize',[9]) 
    xlim([0 9]);
    ylabel('v(t), m(t)')
    xlabel('time (t)')
    title( width(jj))
    vs = sh.*v;
    vsf = fft(vs);
    vsfs_all(jj, :) = abs(fftshift(vsf))/max(abs(vsf));  %% FFT of signal
    
    subplot(3,3, jj+3), plot(t,vs,'k') 
    set(gca,'Fontsize',[9]) 
    xlim([0 9]);
    title('filtered signal (vs)')
    ylabel('v(t)s(t)'), xlabel('time (t)')
    
    
    subplot(3,3,(jj)+6);
    plot(ks, vsfs_all(jj, :), 'b', 'Linewidth', [1]);
    set(gca, 'Fontsize', [9]);
    ylabel('fft(vs)');
    xlabel('frequency (\omega)');
   
end
% wide window coarse translation
figure(11)
sgtitle('Shannon Step Filter with a = 3 (wide width)')
vst_spec=[]; 
tslide=0:1.5:9;
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 3/2; % shannon filter
    vs = sh.*v;
    vst=fft(vs); 
    vst_spec=[vst_spec;  abs(fftshift(vst))]; 
    subplot(2,2,1), plot(t,v,'k', t,sh, 'r'); hold on;
    title('Coarse Translation (b = 1.5)')
           ylabel('Amplitude');
     xlabel('time(t)')

end

subplot(2,2,3),
pcolor(tslide,ks,vst_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
caxis([0 200])
h = colorbar; ylabel(h, 'FFT of filtered signal')
ylabel('frequency (\omega)')
    title('Spectrogram (coarse)')
xlabel('time(t)')
% shannon, wide window, fine translation
vst_spec=[]; 
tslide=0:.5:9;
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 3/2; % Gabor 
    vs = sh.*v;
    vst=fft(vs); 
    vst_spec=[vst_spec;  abs(fftshift(vst))]; 
    subplot(2,2,2), plot(t,v,'k', t,sh, 'r'); hold on;
    title('Fine Translation (b = 0.5)')
     ylabel('Amplitude');
     xlabel('time(t)')
end
%
subplot(2,2,4),
pcolor(tslide,ks,vst_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
h = colorbar; ylabel(h, 'FFT of filtered signal')
caxis([0 200])
    title('Spectrogram (fine)')
ylabel('frequency (\omega)')
xlabel('time(t)')

% shannon, narrow window, coarse translation
figure(12)
sgtitle('Shannon Step Filter a = 0.5 (narrow width)')
vst_spec=[]; 
tslide=0:1.5:9;
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 0.5/2; % Gabor 
    vs = sh.*v;
    vst=fft(vs); 
    vst_spec=[vst_spec;  abs(fftshift(vst))]; 
    subplot(2,2,1), plot(t,v,'k', t,sh, 'r'); hold on;
    title('Coarse Translation (b = 1.5)')
     ylabel('Amplitude');
     xlabel('time(t)')
end

subplot(2,2,3),
pcolor(tslide,ks,vst_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
caxis([0 200])
ylabel('frequency (\omega)')
xlabel('time(t)')
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram (coarse)')
% shannon, narrow window, fine translation
vst_spec=[]; 
tslide=0:.1:9;
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 0.5/2; % Gabor 
    vs = sh.*v;
    vst=fft(vs); 
    vst_spec=[vst_spec;  abs(fftshift(vst))]; 
    subplot(2,2,2), plot(t,v,'k', t,sh, 'r'); hold on;
    title('Fine Translation (b = 0.1)')
      ylabel('Amplitude');
     xlabel('time(t)')
end
%
subplot(2,2,4),
pcolor(tslide,ks,vst_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
colorbar
caxis([0 200])
ylabel('frequency (\omega)')
xlabel('time(t)')
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram (fine)')

% shannon, narrow window, fine translation
figure(13)
vst_spec=[]; 
tslide=0:0.1:9;
for ii=1:length(tslide)
    sh = abs(t-tslide(ii)) <= 0.5/2; % Gabor 
    vs = sh.*v;
    vst=fft(vs); 
    vst_spec=[vst_spec;  abs(fftshift(vst))]; 
    subplot(2,1,1), plot(t,v,'k'); hold on;
    ylabel('Amplitude')
    xlabel('time(t)')
    title('original audio signal')
end

subplot(2,1,2);
pcolor(tslide,ks,vst_spec.'), 
shading interp
set(gca,'Fontsize',[9]) 
colormap('hot')
ylabel('frequency (\omega)')
xlabel('time(t)')
h = colorbar; ylabel(h, 'FFT of filtered signal')
title('Spectrogram:  Shannon Filter (a = 0.5, b = 0.1)')