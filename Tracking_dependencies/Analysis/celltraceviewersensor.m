function celltraceviewersensor(x,y,gem,cdk,xframes,frames_per_hour,IDs)
figure;
s(1)=subplot(1,3,1);
hold on;
s(2)=subplot(1,3,2);
%hold on;
s(3)=subplot(1,3,3);
%hold on;
subplot(s(1))
hax = dscatter(x,y);
xlabel('Hours since APCoff'); ylabel('Intensity');
handle=get(hax,'children');
set(handle,'ButtonDownFcn',{@showcell});
xmax=max(x);
ymax=max(y);


function showcell ( objectHandle , eventData)
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = coordinates(1,1:2);
distance=sqrt(((coordinates(1)-x)/xmax).^2+((coordinates(2)-y)/ymax).^2);
[~,index]=min(distance);
subplot(s(2))
plot(xframes/frames_per_hour,gem(index,:));
title(['Cell: ' num2str(IDs(index))]);
xlabel('Time (hrs)'); ylabel('Geminin');
xlim([min(xframes) max(xframes)]/frames_per_hour + 1);
ylim([0 max(gem(index,:))]);
subplot(s(3))
plot(xframes/frames_per_hour,cdk(index,:));
title(['Cell: ' num2str(IDs(index))]);
xlabel('Time (hrs)'); ylabel('Sensor');
xlim([min(xframes) max(xframes)]/frames_per_hour + 1);
ylim([0 max(cdk(index,:))]);




end

end