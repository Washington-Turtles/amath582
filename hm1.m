clear all; close all; clc; 
load Testdata

L=15; % spatial domain 
n=64; % Fourier modes 
x2=linspace(-L,L,n+1); % spacial discretization
x=x2(1:n);  % dimension x
y=x;  % dimension y
z=x;  % dimension z
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];  % FFT frequency components
ks=fftshift(k);  % move 0 frequencies to center for visualization purposes

[X,Y,Z]=meshgrid(x,y,z);   % create spatial grid space
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);  % create frequency gridspace
Untp_sum = zeros(n,n,n);
Untp_all = zeros(length(Undata(:, 1,1,1)),n,n,n);

% Average data in frequency domain to remove white noise and identify central frequency 
for jj = 1:length(Undata(:, 1,1,1))
        Un = reshape(Undata(jj,:),n,n,n);
        Untp = fftshift(fftn(Un));
        Untp_sum = Untp_sum + Untp ;
        Untp_all(jj,:,:,:) = Untp; % save this for filter/inverse iter
        
end
Utpave=abs(Untp_sum)/jj;
Utpave_norm = Utpave/max(Utpave(:));

% viz average in 3D after normalization without filter
close all
isosurface(Kx,Ky,Kz,abs(Utpave_norm), .55) 
    axis([-10 10 -10 10 -10 10]), grid on, drawnow
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')

%% extract max value/indices
[M, I] = max(abs(Utpave(:)));
[kmy,kmx,kmz] = ind2sub(size(Utpave), I);
kx = ks(kmx); ky = ks(kmy); kz =  ks(kmz);

%  create 3D Gaussian filter given tau and CF
tau = 4; % standard dev for Gaussian filter (how narrow should the filter be)
filter = exp(-(tau.*(Kx - kx).^2 + tau.*(Ky - ky).^2 + tau.*(Kz - kz).^2));

%% Apply filter to average signal for viz of cleaned CF
close all
Utpf = filter.*Utpave;
isosurface(Kx,Ky,Kz,abs(Utpf), .75) 
    axis([-10 10 -10 10 -10 10]), grid on, drawnow
title('Center Frequency after filter')  
xlabel('Kx')
ylabel('Ky')
zlabel('Kz')
print(gcf, '-dpng', 'id_cf.png')

%% Generate cleaned Spatial Realizations
for mm = 1:jj
   Untp = squeeze(Untp_all(mm,:,:,:));
    Untpf = filter.*(Untp);
    Untf = ifftshift(Untpf);
    Unf = ifftn(Untf);     
    Unf_all(mm,:,:,:) = Unf;
end

close all
%%  Identify coordinates for each realization and create table figure
xs = zeros(jj, 1); ys = zeros(jj, 1); zs = zeros(jj, 1);
for pp = 1:jj
    real = squeeze(Unf_all(pp, :, :, :));
    [PM, PI] = max(abs(real(:)));
    [py,px,pz] = ind2sub(size(real), PI);
    xx = x(px); yy = y(py); zz =  z(pz);
    xs(pp, :) = x(px); ys(pp,:) = y(py); zs(pp,:) = z(pz); 

end

% plot Isosurface for each realization
close all
for vv = 1:jj
  Unf = squeeze(Unf_all(vv, :, :, :));
  Unf_norm = abs(Unf)/(max(abs(Unf(:))));
     
  isosurface(X,Y,Z,abs(Unf_norm), .9) 
  axis([-15 15 -15 15 -15 15]), grid on, drawnow   
  hold on
end
title('Marble Locations')
view([-40.5 5]);
xlabel('X')
ylabel('Y')
zlabel('Z')
print(gcf, '-dpng', 'marble_loc.png');


close all
plot3(xs, ys, zs, 'r', 'Linewidth', 3)
view([-40.5 5]);
xlim([-15 15]);
ylim([-15 15]);
zlim([-15 15]);
grid('on');
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Marble Path')
print(gcf, '-dpng', 'descent_path.png');








