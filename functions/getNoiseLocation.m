function coords = getNoiseLocation(data, flag)

if nargin < 2
    flag = 1;
end

if flag
    figNoise = figure; 
    imshow(imadjust(im2double(uint16(data(:,:,50)))));
    
else
    figNoise = figure;
    imshow(imadjust(data(:,:,50)));
end


title('draw box around a region with noise only');
% k = waitforbuttonpress;
% point1 = get(gca,'CurrentPoint');    % button down detected
% finalRect = rbbox;                   % return figure units
% point2 = get(gca,'CurrentPoint');    % button up detected
% point1 = point1(1,1:2);              % extract x and y
% point2 = point2(1,1:2);
% p1 = min(point1,point2);             % calculate locations
% offset = abs(point1-point2);         % and dimensions
% x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
% y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
% hold on
% axis manual
% plot(x,y,'r','linewidth',2);
% hold off;
% 
% xmin = round(min(x));
% xmax = round(max(x));
% ymin = round(min(y));
% ymax = round(max(y));
% 
% coords = [xmin xmax ymin ymax];
coord = getrect; % [xmin ymin width height]
coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
close(figNoise);
