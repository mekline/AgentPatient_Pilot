function JpgResize(image_type)
%Resizes a jpg image
    original_dir = pwd
    
    %PsychToolbox to get winHeight
    screenNum = max(Screen('Screens'));  %Highest screen number is most likely correct display
    SCREEN_ADJUST = 1.2;
    %windowInfo = PTBhelper('initialize',screenNum);
	%wPtr = windowInfo{1}; %pointer to window on screen that's being referenced
    wPtr = 10;
    %rect = windowInfo{2}; %dimensions of the window
    rect = [0 0 1280 800];
        winWidth = rect(3);
        winHeight = rect(4)*SCREEN_ADJUST;
    

    %Save file in image folder
    DESTINATION_IMAGE_DIR = fullfile(pwd, 'images', ['resized_jpg_' image_type]);
    SOURCE_IMAGE_DIR = fullfile(pwd, 'images', image_type);
    item_index = 1;
    
    cd(SOURCE_IMAGE_DIR)
    images = dir('*.jpg');
    for i=1:length(images)
        source_image = fullfile(images(i).name);
            %Sets up image and overlays text
            I = imread(source_image);
            %oldimsize = size(I)
            I = imresize(I, [winHeight, NaN]);
            %newimsize = size(I)

            %Sets up file to save; numbers indicate index at which stimulus
            %was presented
%           fileToSave = [condition '_' flip_word '_' itemNumber '.jpg'];
            fileToSave = fullfile(DESTINATION_IMAGE_DIR, images(i).name);

            %Displays text image
            %figure
            %imshow(I)

            %Saves text image
            imwrite(I,fileToSave,'jpg')

            %Increments counter
            item_index = item_index + 1;
    end
cd(original_dir) %brings you back to your original directory
    
end