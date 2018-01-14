B=imread('FemurA.png');
% B = rgb2gray(B);  % Convert the image to grayscale
 B = double(B);    % Convert the uint8 (8-bit) to a double
                  %   This makes operations with A have double precision

imshow(B)

