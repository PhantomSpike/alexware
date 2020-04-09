function imfuseQ(image1,varargin)

% imfuseQ(image1,image2,...) overlay image2 to 4 on image1. The sliders sets the transparency of

% each image. image1 is black and white, images 2 to 4 are red/green/blue

% respectively. Images must be of same size.



% Figure

hFig = gcf;

set(hFig,'DeleteFcn',@closeCallback);



hueNames = {'red' 'green' 'blue'};

hue{1} = cat(3, ones(size(image1)), zeros(size(image1)), zeros(size(image1)));

hue{2} = cat(3, zeros(size(image1)), ones(size(image1)), zeros(size(image1)));

hue{3} = cat(3, zeros(size(image1)), zeros(size(image1)), ones(size(image1)));



image1 = squeeze(image1 ./ max(max(image1)));

figure(hFig);

imshow(image1);

hold on



for i = 1:length(varargin),

    % Slider i

    sliderhandles(i) = figure('units','normalized','position',[0.65 0.6+0.1*(i-1) 0.25 0.05],'menubar','none','toolbar','none','name',sprintf('Transparency adjustment %s',hueNames{i}));

    hSlider(i) = uicontrol(sliderhandles(i),'Style','slider','units','normalized','position',[0 0 1 1],'value',0.5,'callback',@changeContrast);

    

    % Normalize image

    handles(i).image = squeeze(varargin{i}./ max(max(varargin{i})));

    

    % Overlay

    figure(hFig);

    handles(i).h = imshow(hue{i});

    set(handles(i).h, 'AlphaData', handles(i).image .* get(hSlider(i),'Value'));

    

    hSlider(i).UserData = handles(i);

end

hold off

hFig.UserData = sliderhandles;

end



function changeContrast(hObject, eventdata, handles)

set(hObject.UserData.h, 'AlphaData',hObject.UserData.image .* get(hObject,'Value'));

end



function closeCallback(hObject, eventdata, handles)

for i = 1:length(hObject.UserData),

    if isvalid(hObject.UserData(i)),

        close(hObject.UserData(i))

    end

end

end