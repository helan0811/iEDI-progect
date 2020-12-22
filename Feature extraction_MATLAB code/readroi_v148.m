function [point]=readroi_v148(filename)
% in ImageJ v1.48 the format of ROI changes

in=fopen(filename,'rb','b');
% fetch border info
fseek(in,6,-1); 
roitype = fread(in,1,'int8');
fread(in,1,'int8');

high=fread(in,1,'int16');
left=fread(in,1,'int16');
low=fread(in,1,'int16');
right=fread(in,1,'int16');

if roitype == 7 % polygon roi
nPoint=fread(in,1,'int16');


fseek(in,23*2,0);
point = fread(in,[nPoint,2],'int16');

point(:,1)=point(:,1)+left;
point(:,2)=point(:,2)+high;


else if roitype == 0 % rect roi
    % to be changed
    point = zeros(4,2);
    %point = a=[left high;left low;right low;right high];
    point = [left high;left low;right low;right high];
    end
end
fclose(in);
