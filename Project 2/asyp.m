% ME C176/BioE C119, Fall 2015
% Asymmetrical Bending Analysis

%% PICTURE IMPORT FOR NO IMPLANT

A=imread('Femur.png');
A = rgb2gray(A);  % Convert the image to grayscale
A = double(A);    % Convert the uint8 (8-bit) to a double
                  %   This makes operations with A have double precision
BMD = 2*A/255;

Eyoung=zeros(86);
Eold=zeros(86);
Ehet=zeros(86);

% Create E matrix for YOUNG, OLD, HETERO
 for i=1:86
     for j=1:86
         if BMD(i,j) > 1
            Eyoung(i,j)=17;
            Eold(i,j)=17;
         end
         Ehet(i,j)=6.72068*BMD(i,j)^1.33885;
     end
 end

 n=0;


%% PICTURE IMPORT FOR IMPLANT

B=imread('FemurA.png');
C = double(B);    % Convert the uint8 (8-bit) to a double
                  %   This makes operations with A have double precision

pixcolor=zeros(86);
n=0;

% Find where are RGB values are not equal, which would mean not grayscale.
% Then create arrays of Z,Y pixel values, representing the implant pixels.
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

% Convert image matrix to grayscale
B = rgb2gray(B);  % Convert the image to grayscale
B = double(B);

BMD = 2*B/255;

Eyoung=zeros(86);
Eold=zeros(86);
Ehet=zeros(86);

 for i=1:86         % Create E matrix for YOUNG, OLD, HETERO
     for j=1:86
         if BMD(i,j) > 1
            Eyoung(i,j)=17;
            Eold(i,j)=13;
         end
         Ehet(i,j)=6.72068*BMD(i,j)^1.33885;
     end
 end
 
 % Replace known pixel locations of implant with E=200
 for i=1:n
     Eyoung(implantZvalues(i),implantYvalues(i))=200;
     Eold(implantZvalues(i),implantYvalues(i))=200;
     Ehet(implantZvalues(i),implantYvalues(i))=200;
 end

                  
%% Now calculate...

res=0.78e-3;	% resolution in meters per pixel

 My = 46; % moment in y [N-m]
 Mz = 28; % moment in z [N-m]
 stressXyoung = zeros (86);   % Create stress matrix young
 stressXold = zeros (86);   % Create stress matrix old
 stressXhet = zeros (86);   % Create stress matrix hetero


% Find BMD weighted centroid
% YOUNG
sumyyoung=0;
sumzyoung=0;
total_sum_young=0;
% OLD
sumyold=0;
sumzold=0;
total_sum_old=0;
% HETERO
sumyhet=0;
sumzhet=0;
total_sum_het=0;

% Obtain zhat and yhat numerators and denominators
for i=1:86
	for j=1:86
            z=(j-1)*res+res/2;	% z position of pixel CENTER in m
            y=(i-1)*res+res/2;	% y position of pixel CENTER in m
            sumzyoung=sumzyoung+z*Eyoung(i,j);
            sumyyoung=sumyyoung+y*Eyoung(i,j);
            total_sum_young= total_sum_young + Eyoung(i,j);
            sumzold=sumzold+z*Eold(i,j);
            sumyold=sumyold+y*Eold(i,j);
            total_sum_old= total_sum_old + Eold(i,j);
            sumzhet=sumzhet+z*Ehet(i,j);
            sumyhet=sumyhet+y*Ehet(i,j);
            total_sum_het= total_sum_het + Ehet(i,j);
	end
end


% Show image YOUNG
figure(1)
imagesc(Eyoung),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;
% Centroid coordinates YOUNG
y_hat_young=sumyyoung/total_sum_young;	% in m
z_hat_young=sumzyoung/total_sum_young;	% in m
fprintf('\n y_hat = %3.2f mm'  ,y_hat_young*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat_young*1000);
% Plot centroid axes on the image YOUNG
line([z_hat_young/res,z_hat_young/res],[1,86])
line([1,86],[y_hat_young/res,y_hat_young/res])

% Show image OLD
figure(2)
imagesc(Eold),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;
% Centroid coordinates OLD
y_hat_old=sumyold/total_sum_old;	% in m
z_hat_old=sumzold/total_sum_old;	% in m
fprintf('\n y_hat = %3.2f mm'  ,y_hat_old*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat_old*1000);
% Plot centroid axes on the image OLD
line([z_hat_old/res,z_hat_old/res],[1,86])
line([1,86],[y_hat_old/res,y_hat_old/res])

% Show image HETERO
figure(3)
imagesc(Ehet),colormap(gray); hold on;
xlabel('z (pixels)'); ylabel('y (pixels)'); axis square;
% Centroid coordinates HETERO
y_hat_het=sumyhet/total_sum_het;	% in m
z_hat_het=sumzhet/total_sum_het;	% in m
fprintf('\n y_hat = %3.2f mm'  ,y_hat_het*1000);
fprintf('\n z_hat = %3.2f mm\n',z_hat_het*1000);
% Plot centroid axes on the image HETERO
line([z_hat_het/res,z_hat_het/res],[1,86])
line([1,86],[y_hat_het/res,y_hat_het/res])

% Create Istar terms
Iyystary=0;
Izzstary=0;
Iyzstary=0;
Iyystaro=0;
Izzstaro=0;
Iyzstaro=0;
Iyystarhet=0;
Izzstarhet=0;
Iyzstarhet=0;
Ipix = res^4/12;
maxstressyoung = 0;
maxstressold = 0;
maxstresshet = 0;

% Integrate Istar terms
for i=1:86
    for j=1:86
        z=(j-1)*res+res/2;	% z position of pixel CENTER in m
        y=(i-1)*res+res/2;	% y position of pixel CENTER in m
        Iyystary=(Eyoung(i,j)*(Ipix+(z-z_hat_young)^2*res^2))+Iyystary;
        Izzstary=(Eyoung(i,j)*(Ipix+(y-y_hat_young)^2*res^2))+Izzstary;
        Iyzstary=(Eyoung(i,j)*(0+(y-y_hat_young)*(z-z_hat_young)*res^2))+Iyzstary;
        Iyystaro=(Eold(i,j)*(Ipix+(z-z_hat_old)^2*res^2))+Iyystaro;
        Izzstaro=(Eold(i,j)*(Ipix+(y-y_hat_old)^2*res^2))+Izzstaro;
        Iyzstaro=(Eold(i,j)*(0+(y-y_hat_old)*(z-z_hat_old)*res^2))+Iyzstaro;
        Iyystarhet=(Ehet(i,j)*(Ipix+(z-z_hat_het)^2*res^2))+Iyystarhet;
        Izzstarhet=(Ehet(i,j)*(Ipix+(y-y_hat_het)^2*res^2))+Izzstarhet;
        Iyzstarhet=(Ehet(i,j)*(0+(y-y_hat_het)*(z-z_hat_old)*res^2))+Iyzstarhet;
    end
end

% Calculate Stresses
for i=1:86
	for j=1:86
            z=(j-1)*res+res/2;	% z position of pixel CENTER in m
            y=(i-1)*res+res/2;	% y position of pixel CENTER in m
        if Eyoung(i,j)~=0
            stressXyoung(i,j) = Eyoung(i,j)*((My*Izzstary+Mz*Iyzstary)*(z-z_hat_young)-(My*Iyzstary+Iyystary*Mz)*(y-y_hat_young))/(Iyystary*Izzstary-Iyzstary^2); % stress in the X direction (tens+, comp-)
        end
        if Eold(i,j)~=0
            stressXold(i,j) = Eold(i,j)*((My*Izzstaro+Mz*Iyzstaro)*(z-z_hat_old)-(My*Iyzstaro+Iyystaro*Mz)*(y-y_hat_old))/(Iyystaro*Izzstaro-Iyzstaro^2); % stress in the X direction (tens+, comp-)
        end
        if Ehet(i,j)~=0
            stressXhet(i,j) = Ehet(i,j)*((My*Izzstarhet+Mz*Iyzstarhet)*(z-z_hat_het)-(My*Iyzstarhet+Iyystarhet*Mz)*(y-y_hat_het))/(Iyystarhet*Izzstarhet-Iyzstarhet^2); % stress in the X direction (tens+, comp-)
        end

	end
end

stressXnoimplantyoung = stressXyoung;
maxstressyoung = 0;
stressXnoimplantold = stressXold;
maxstressold = 0;
stressXnoimplanthet = stressXhet;
maxstresshet = 0;

% Remove implant stresses in order to find max stress, only through bone.
for i=1:n
        stressXnoimplantyoung(implantZvalues(i),implantYvalues(i))=0;
        stressXnoimplantold(implantZvalues(i),implantYvalues(i))=0;
        stressXnoimplanthet(implantZvalues(i),implantYvalues(i))=0;
end

% ISOLATING MAX STRESSES
for i=1:86
    for j=1:86
        if stressXnoimplantyoung(i,j) > maxstressyoung
        maxstressyoung = stressXnoimplantyoung(i,j);
        locyoungz=j*res-res/2;
        locyoungy=i*res-res/2;
        pixyoungz=j;  % Pixel location used for normalized stress calc
        pixyoungy=i;  % Pixel location used for normalized stress calc
        end
        if stressXnoimplantold(i,j) > maxstressold
        maxstressold = stressXnoimplantold(i,j);
        locoldz=j*res-res/2;
        locoldy=i*res-res/2;
        pixoldz=j;  % Pixel location used for normalized stress calc
        pixoldy=i;  % Pixel location used for normalized stress calc
        end
        if stressXnoimplanthet(i,j) > maxstresshet
        maxstresshet = stressXnoimplanthet(i,j);
        lochetz=j*res-res/2;
        lochety=i*res-res/2;
        pixhetz=j;  % Pixel location used for normalized stress calc
        pixhety=i;  % Pixel location used for normalized stress calc
        end
    end
end


% YOUNG PLOTS
myoung=(My*Izzstary+Mz*Iyzstary)/(My*Iyzstary+Mz*Iyystary);
Xyoung=25:60;
Yyoung=myoung*(Xyoung-z_hat_young*1282)+y_hat_young*1282;

figure(4)
contourf(stressXyoung, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
caxis([-30*10^6,30*10^6]);
plot(Xyoung,Yyoung);
line([z_hat_young/res,z_hat_young/res],[1,86])
line([1,86],[y_hat_young/res,y_hat_young/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Stress in Young Bone w/ Implant'); axis square;

figure(5)
contourf(stressXnoimplantyoung, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
caxis([-30*10^6,30*10^6]);
plot(Xyoung,Yyoung);
line([z_hat_young/res,z_hat_young/res],[1,86])
line([1,86],[y_hat_young/res,y_hat_young/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Stress in Young Bone'); axis square;

% OLD PLOTS
mold=(My*Izzstaro+Mz*Iyzstaro)/(My*Iyzstaro+Mz*Iyystaro);
Xold=25:60;
Yold=mold*(Xold-z_hat_old*1282)+y_hat_old*1282;

figure(6)
contourf(stressXold, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
caxis([-30*10^6,30*10^6]);
plot(Xold,Yold);
line([z_hat_old/res,z_hat_old/res],[1,86])
line([1,86],[y_hat_old/res,y_hat_old/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Stress in Old Bone w/ Implant'); axis square;

figure(7)
contourf(stressXnoimplantold, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
caxis([-30*10^6,30*10^6]);
plot(Xold,Yold);
line([z_hat_old/res,z_hat_old/res],[1,86])
line([1,86],[y_hat_old/res,y_hat_old/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Stress in Old Bone'); axis square;

% HETERO PLOTS
mhet=(My*Izzstarhet+Mz*Iyzstarhet)/(My*Iyzstarhet+Mz*Iyystarhet);
Xhet=25:60;
Yhet=mhet*(Xold-z_hat_het*1282)+y_hat_het*1282;

figure(8)
contourf(stressXhet, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
caxis([-30*10^6,30*10^6]);
plot(Xhet,Yhet);
line([z_hat_het/res,z_hat_het/res],[1,86])
line([1,86],[y_hat_het/res,y_hat_het/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Stress in Hetero Bone w/ Implant'); axis square;

figure(9)
contourf(stressXnoimplanthet, 100, 'linecolor', 'none'); hold on;
set(gca, 'ydir','reverse');
caxis([-30*10^6,30*10^6]);
plot(Xhet,Yhet);
line([z_hat_het/res,z_hat_het/res],[1,86])
line([1,86],[y_hat_het/res,y_hat_het/res])
xlabel('z (pixels)'); ylabel('y (pixels)'); title('Stress in Hetero Bone'); axis square;

STRESSyoung = stressXyoung(30,61);  %Stress at pixel of interest
STRESSold = stressXold(30,61);  %Stress at pixel of interest
STRESShetero = stressXhet(29,58);  %Stress at pixel of interest