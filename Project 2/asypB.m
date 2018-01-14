% ME C176/BioE C119, Fall 2015
% Asymmetrical Bending Analysis

%% PICTURE FOR NO IMPLANT
clear all;
A=imread('Femur.png');
A = rgb2gray(A);  % Convert the image to grayscale
A = double(A);    % Convert the uint8 (8-bit) to a double
                  %   This makes operations with A have double precision
BMD = 2*A/255;

Eyoung=zeros(86);
Eold=zeros(86);
                  
 for i=1:86         % CASE 1: create binary E matrix with E=17 or 13
     for j=1:86
         if BMD(i,j) > 1
            Eyoung(i,j)=17;
            Eold(i,j)=17;

         end
     end
 end


%% PICTURE FOR IMPLANT

B=imread('FemurA.png');
C = double(B);    % Convert the uint8 (8-bit) to a double
                  %   This makes operations with A have double precision

pixcolor=zeros(86);
n=0;

for i=1:86
    for j=1:86
    pixcolor=squeeze(C(i,j,:));
    if pixcolor(1,1) ~= pixcolor(2,1) || pixcolor(1,1) ~= pixcolor(3,1)
        n=n+1;
        implantZvalues(n)=i;
        implantYvalues(n)=j;
    end
    end
end

B = rgb2gray(B);  % Convert the image to grayscale
B = double(B);

BMD = 2*B/255;

Eyoung=zeros(86);
Eold=zeros(86);

 for i=1:86         % CASE 1: create binary E matrix with E=17 or 0
     for j=1:86
         if BMD(i,j) > 1
            Eyoung(i,j)=17;
            Eold(i,j)=13;
         end
     end
 end
 
 for i=1:n
     Eyoung(implantZvalues(i),implantYvalues(i))=200;
     Eold(implantZvalues(i),implantYvalues(i))=200;
 end

                  
%% Now calculate E and Go...

res=0.78e-3;	% resolution in meters per pixel
 
 My = 46; % moment in y [N-m]
 Mz = 28; % moment in z [N-m]
 stressXyoung = zeros (86);   % Create stress matrix young
 stressXold = zeros (86);   % Create stress matrix old


% Find BMD weighted centroid [You'll need to modify this for implants]
sumy=0;
sumz=0;
total_sum=0;

for i=1:86
	for j=1:86
            z=(j-1)*res+res/2;	% z position of pixel CENTER in m
            y=(i-1)*res+res/2;	% y position of pixel CENTER in m
            sumz=sumz+z*Eyoung(i,j);
            sumy=sumy+y*Eyoung(i,j);
            total_sum= total_sum + Eyoung(i,j);
            sumz=sumz+z*Eold(i,j);
            sumy=sumy+y*Eold(i,j);
            total_sum= total_sum + Eold(i,j);
	end
end


% Show image YOUNG
figure(1)
imagesc(Eyoung),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

% Centroid coordinates
y_hat=sumy/total_sum;	% in m
z_hat=sumz/total_sum;	% in m
fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

% Plot centroid axes on the image 
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])

% Show image OLD
figure(2)
imagesc(Eold),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

% Centroid coordinates
y_hat=sumy/total_sum;	% in m
z_hat=sumz/total_sum;	% in m
fprintf('\n y_hat = %3.2f mm'  ,y_hat*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat*1000);

% Plot centroid axes on the image 
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])

Iyystary=0;
Izzstary=0;
Iyzstary=0;
Iyystaro=0;
Izzstaro=0;
Iyzstaro=0;
Ipix = res^4/12;
maxstressyoung = 0;
maxstressold = 0;

for i=1:86
    for j=1:86
        z=(j-1)*res+res/2;	% z position of pixel CENTER in m
        y=(i-1)*res+res/2;	% y position of pixel CENTER in m
        Iyystary=(Eyoung(i,j)*(Ipix+(z-z_hat)^2*res^2))+Iyystary;
        Izzstary=(Eyoung(i,j)*(Ipix+(y-y_hat)^2*res^2))+Izzstary;
        Iyzstary=(Eyoung(i,j)*(0+(y-y_hat)*(z-z_hat)*res^2))+Iyzstary;
        Iyystaro=(Eold(i,j)*(Ipix+(z-z_hat)^2*res^2))+Iyystaro;
        Izzstaro=(Eold(i,j)*(Ipix+(y-y_hat)^2*res^2))+Izzstaro;
        Iyzstaro=(Eold(i,j)*(0+(y-y_hat)*(z-z_hat)*res^2))+Iyzstaro;

    end
end

for i=1:86
	for j=1:86
            z=(j-1)*res+res/2;	% z position of pixel CENTER in m
            y=(i-1)*res+res/2;	% y position of pixel CENTER in m
        if Eyoung(i,j)~=0
            stressXyoung(i,j) = Eyoung(i,j)*((My*Izzstary+Mz*Iyzstary)*(z-z_hat)-(My*Iyzstary+Iyystary*Mz)*(y-y_hat))/(Iyystary*Izzstary-Iyzstary^2); % stress in the X direction (tens+, comp-)
        end
        if Eold(i,j)~=0
            stressXold(i,j) = Eold(i,j)*((My*Izzstaro+Mz*Iyzstaro)*(z-z_hat)-(My*Iyzstaro+Iyystaro*Mz)*(y-y_hat))/(Iyystaro*Izzstaro-Iyzstaro^2); % stress in the X direction (tens+, comp-)
        end
	end
end

stressXnoimplantyoung = stressXyoung;
maxstressyoung = 0;
stressXnoimplantold = stressXold;
maxstressold = 0;

for i=1:n
        stressXnoimplantyoung(implantZvalues(i),implantYvalues(i))=0;
        stressXnoimplantold(implantZvalues(i),implantYvalues(i))=0;
end

for i=1:86
    for j=1:86
        if stressXnoimplantyoung(i,j) > maxstressyoung
        maxstressyoung = stressXnoimplantyoung(i,j);
        end
        if stressXnoimplantold(i,j) > maxstressold
        maxstressold = stressXnoimplantold(i,j);
        end
    end
end


% YOUNG
myoung=(My*Izzstary+Mz*Iyzstary)/(My*Iyzstary+Mz*Iyystary);
Xyoung=20:70;
Yyoung=myoung*(Xyoung-z_hat*1282)+y_hat*1282;

figure(3)
contourf(stressXyoung, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
plot(Xyoung,Yyoung);
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

figure(4)
contourf(stressXnoimplantyoung, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
plot(Xyoung,Yyoung);
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

%OLD 
mold=(My*Izzstaro+Mz*Iyzstaro)/(My*Iyzstaro+Mz*Iyystaro);
Xold=20:60;
Yold=mold*(Xold-z_hat*1282)+y_hat*1282;

figure(3)
contourf(stressXold, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
plot(Xold,Yold);
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

figure(4)
contourf(stressXnoimplantold, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
plot(Xold,Yold);
line([z_hat/res,z_hat/res],[1,86])
line([1,86],[y_hat/res,y_hat/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;

