function [OptoROI,xvals,yvals,xcenter,ycenter]=DefineOptoROI(dalsaFrame)
Optofil=input('Input dalsa filename with optogenetic stimulus\n','s');
Optoimg=ReadDalsaBinary_Matrix(Optofil,256,256);
Flashimg=std(double(Optoimg),0,3);
figure(100);imagesc(Flashimg);axis image; axis off;
figure(99);imagesc(dalsaFrame);colormap gray; axis image; axis off;
title('Select center of optogenetic ROI');
[xcenter,ycenter]=ginput(1);
thecircle=drawcircle('Center',[round(xcenter,0),round(ycenter,0)],'Radius',15,'Color',[1 1 1], 'LineWidth', 0.25,'InteractionsAllowed','none');
ROIbounds=round(thecircle.Vertices,0);
OptoROI=roipoly(dalsaFrame,ROIbounds(:,1),ROIbounds(:,2));
xvals=ROIbounds(:,1);
yvals=ROIbounds(:,2);