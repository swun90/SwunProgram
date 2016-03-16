%makeMovie 
%   make a movie with 'fps' frames per second
%   i.e. the whole video would be N/fps seconds,
%   in which N is the number of all figures.
%   'quality' controls the video size, and it verifies from 0 to 100.
%   An example, fps=12, quality=50, 1000 figures and a mp4 video 
%   would be 24.8MB.
%
% change 'figPrefix' first!
% see also VideoWriter

figPrefix='02z'; % the several symbol of the figures' name.
fps=12;
quality=50;
myVideo=VideoWriter('myVideo','MPEG-4'); 
% options of second parameters: MPEG-4, Motion JPEG AVI, Motion JPEG 2000

Nfig=length(dir([figPrefix,'_*.png']));
myVideo.FrameRate=fps;
myVideo.Quality=quality;
open(myVideo);


for i=1:Nfig
    disp(i);
    p=imread(['02z_',sprintf('%4.4d',i),'.png']);
    imag=imresize(p,[902,1204]);
    imshow(imag);
    currFrame=getframe;
    writeVideo(myVideo,currFrame);
end
close(myVideo);