clear all; close all; clc

load cam1_1.mat
load cam2_1.mat
load cam3_1.mat

xs = [1: 640];
ys =[1:480];

%% Cam 1
cam1 = vidFrames1_1; 
f1x =repelem((abs(xs-322) <= 32), 480, 1) ; %create filter along vert axis
f1y = repelem(abs(ys- 290) <= 90, 640, 1).'; %create filter along horiz axis
filter1 = f1x.*f1y; %final filter
frame_count = length(cam1(1,1,1,:));
ts = [1:frame_count];
% filter all frames (movement starts at 8 or 9)
for jj = 1: frame_count
    frame = double(rgb2gray(cam1(:,:,:,jj)));
    frame_f = filter1.*frame;
    frame_fn = frame_f/max(max(frame_f));
    cam1_fn(:,:,:,jj) = frame_fn;
end
% find max value in each frame (end with 1x226 plot)
c1 = zeros(2, frame_count).';
for kk = 1:frame_count
    frame = cam1_fn(:,:,:,kk);
    I = find(frame == 1);
    [y,x] = ind2sub(size(frame), I);
    y = min(y);
    x = mean(x);
   c1(kk, :) = [-x,-y]';
end

%% Cam 2
cam2= vidFrames2_1; 
f1x =repelem((abs(xs-305) <= 35), 480, 1) ; %create filter along vert axis
f1y = repelem(abs(ys- 200) <= 100, 640, 1).'; %create filter along horiz axis
filter1 = f1x.*f1y; %final filter
frame_count = length(cam2(1,1,1,:));
ts = [1:frame_count];
% filter all frames (movement starts at 8 or 9
for jj = 1: frame_count
    frame = double(rgb2gray(cam2(:,:,:,jj)));
    frame_f = filter1.*frame;
    frame_fn = frame_f/max(max(frame_f));
    cam2_fn(:,:,:,jj) = frame_fn;
end
% find max value in each frame (end with 1x226 plot)
c2 = zeros(2, frame_count).';
for kk = 1:frame_count
    frame = cam2_fn(:,:,:,kk);
    I = find(frame > .97);
    [y,x] = ind2sub(size(frame), I);
    y = min(y);
    x = median(x);
   c2(kk, :) = [-x,-y]';
end

%% Cam 3
cam3= vidFrames3_1; 
f1x =repelem((abs(xs-350) <= 80), 480, 1) ; %create filter along vert axis
f1y = repelem(abs(ys- 270) <= 25, 640, 1).'; %create filter along horiz axis
filter1 = f1x.*f1y; %final filter
frame_count = length(cam3(1,1,1,:));
ts = [1:frame_count];
% filter all frames (movement starts at 8 or 9
for jj = 1: frame_count
    frame = double(rgb2gray(cam3(:,:,:,jj)));
    frame_f = filter1.*frame;
    frame_fn = frame_f/max(max(frame_f));
    cam3_fn(:,:,:,jj) = frame_fn;
end
% find max value in each frame (end with 1x226 plot)
c3= zeros(2, frame_count).';
for kk = 1:frame_count
    frame = cam3_fn(:,:,:,kk);
    I = find(frame == 1);
    [y,x] = ind2sub(size(frame), I);
    y = min(y);
    x = min(x);
   c3(kk, :) = [-x,-y]';
end

%% truncate x,y from each cam to match (4.5 cycles)
x1 = c1(10:180, 1) ; y1 = c1(10:180, 2); %cam1  30-170
x2 = c2(20:190, 1) ; y2 = c2(20:190, 2); %cam2  40-180
x3 = c3(10:180, 1) ; y3 = c3(10:180, 2); %cam2 30-170
A = horzcat(x1, y1, x2, y2, x3, y3)';  % Matrix for SVD
%% calculate displacement
[m,n]=size(A); %compute data size
mn=mean(A,2); %compute mean for each row
A=A-repmat(mn, 1, n); % subtract mean
fs = 1:length(A); 
r1 = rank(A); %full rank due to noise
covA = cov(A); % covariance matrix
%
figure(1)
sgtitle('Ideal Case')
subplot(3,1,1); plot(fs', A(1,:), 'b', 'Linewidth', [2]);  hold on;
title('Camera 1')
plot(fs', A(2,:), 'r', 'Linewidth', [2])
legend('X direction','Y direction', 'Location', 'NorthEast') 
xlabel('Frame No.')
ylabel('displacement');
ylim([-200 300])
subplot(3,1,2); plot(fs', A(3,:), 'b', 'Linewidth', [2]);  hold on;
title('Camera 2')
plot(fs', A(4,:), 'r', 'Linewidth', [2])
legend('X direction','Y direction', 'Location', 'NorthEast') 
xlabel('Frame No.')
ylabel('displacement');
ylim([-200 300])
subplot(3,1,3); plot(fs', A(5,:), 'b', 'Linewidth', [2]);  hold on;
title('Camera 3')
plot(fs', A(6,:), 'r', 'Linewidth', [2])
legend('X direction','Y direction', 'Location', 'NorthEast') 
xlabel('Frame No.')
ylabel('displacement');
ylim([-200 300])

%%
[u,s,v]=svd(A');
sig = diag(s);
lambda = sig.^2;
Y = u*A';

energy1 = sig(1)/sum(sig);
energy2 = sum(sig(1:2))/sum(sig);
energy3 = sum(sig(1:3))/sum(sig);
energy4 = sum(sig(1:4))/sum(sig);
energy5 = sum(sig(1:5))/sum(sig);
energy6 = sum(sig(:))/sum(sig);
modes =1:6;
%%
figure(2)
subplot(3,2,1), plot(sig,'ko','Linewidth',[3])
sgtitle('Ideal Case', 'Fontsize', [15]);
title('Sigma (\sigma) values' ,'Fontsize', [10]);
set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6 7])
ylim([0 1300]); 
xlabel('mode')
ylabel('\sigma')

subplot(3,2,2), semilogy(sig,'ko','Linewidth',[3])
title('Sigma (\sigma) values - log y', 'Fontsize', [10]);
set(gca,'Fontsize',[13],'Ytick',[10^(0) 10^(1) 10^(2) 10^3 10^4 10^5 10^6 10^7],...
   'Xtick',[0 1 2 3 4 5 6]); 
ylim([0 2000]);
xlabel('mode')
ylabel('\sigma')
subplot(3,1,2) 
plot(fs,u(:,1),'k','Linewidth',[2])
ylim([-0.2 0.2])
ylabel('u');
xlabel('frame no.')
legend('Mode 1')
subplot(3,1,3) 
plot(modes,v(:,1),'k','Linewidth',[2])
ylabel('v');
xlabel('mode')
legend('Mode 1')
ylim([-0.1 1])



