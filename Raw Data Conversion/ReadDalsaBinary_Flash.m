%% 
% Summary: converts a .bin file into two image stacks.

% Inputs:       thefile - 

% Outputs:      CBV - cerebral blood volume images from green light. Used with 
%                       Arduino_Cam_to_LEDs_v3.ino to control the flashing
%                       of LEDs, blood volume should always be contained in
%                       the first frame.

%               oxy - oxygenation images from red light.

function [CBV,oxy]=ReadDalsaBinary_Flash(thefile, image_height, image_width)

pixels_per_frame=image_width*image_height;
%open the file , get file size , back to the begining
fid=fopen(thefile);
fseek(fid,0, 'eof');
thefile_size=ftell(fid);
fseek(fid,0, 'bof');

% identify the number of frames to read. Each frame has a previously
% defined width and height (as inputs), along with a grayscale "depth" of
% 2"
nframes_to_read=floor(thefile_size/(2*pixels_per_frame))
counter1 = 1;
counter2 = 1;
for n=1:nframes_to_read
    if n == 300;
        test = 1;
    end
    z=fread(fid, pixels_per_frame,'*int16','b');
    if rem(n,2) ~= 0;
        img=reshape(z(1:pixels_per_frame),image_height,image_width);
        CBV{counter1} = rot90(img',2);
        counter1 = counter1+1;
    else
        img=reshape(z(1:pixels_per_frame),image_height,image_width);
        oxy{counter2} = rot90(img',2);
        counter2 = counter2+1;
    end
end

fclose('all');