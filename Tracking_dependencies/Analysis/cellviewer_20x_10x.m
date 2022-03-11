function cellviewer_20x_10x(x,y,shots,pos,channel1,channel2,impth,app)
figure;
s(1)=subplot(1,3,1);
hold on;
s(2)=subplot(1,3,2);
hold on;
s(3)=subplot(1,3,3);
hold on;
subplot(s(1))
hax = dscatter(x,y);
handle=get(hax,'children');
set(handle,'ButtonDownFcn',{@showcell});



function showcell ( objectHandle , eventData)
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = coordinates(1,1:2);
subplot(s(1))
t=title([ ' x:' num2str(coordinates(1)) ', y:' num2str(coordinates(2))]);
set(t,'interpreter','none');

distance=sqrt((coordinates(1)-x).^2+(coordinates(2)-y).^2);
[~,index]=min(distance);
shot=strcat(num2str(shots(index,1)), '_', num2str(shots(index,2)),'_', num2str(shots(index,3)));
position=pos(index,:);
im=imread([impth shot '\' shot '_' channel1 app]);
subplot(s(2))
imshow(log2(single(im)),[]);
plot(position(1),position(2),'.r');
xlim([-100 100] + position(1)); ylim([-100 100] + position(2));
t=title([shot ' x:' num2str(position(1)) ', y:' num2str(position(2))]);
set(t,'interpreter','none');
set(gca,'YDir','Reverse')

im=imread([impth shot '\' shot '_' channel2 app]);
subplot(s(3))
imshow(log2(single(im)),[]);
plot(position(1),position(2),'.r');
xlim([-100 100] + position(1)); ylim([-100 100] + position(2));
t=title([shot ' x:' num2str(position(1)) ', y:' num2str(position(2))]);
set(t,'interpreter','none');
set(gca,'YDir','Reverse')

colormap(s(1),parula)
%imtool(im);

end

end

%{
figure, imshow(im,[]), hold on
scatter(pos(:,1),pos(:,2))
%}