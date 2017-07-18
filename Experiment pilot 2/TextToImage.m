function TextToImage(order, run)
%Takes in text and turns it into a jpeg file in images

    %Constants
    FONT_SIZE = 90;
    EMS_TO_PIXELS = FONT_SIZE*.59;
    
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
    
    %Reads in materials file
    ORDER_DIR = fullfile(pwd, 'reference_data');
    order_filename = ['AgentPatientStimuli_reference_' order num2str(run) '_data.csv'];
    order_filename = fullfile(ORDER_DIR, order_filename);
    all_materials = readtable(order_filename);
    sentences = all_materials.Sentence;
    conditions = all_materials.Condition;
    flips = all_materials.Flip;
    agentNames = all_materials.AgentName;
    patientNames = all_materials.PatientName;
    agentShapes = all_materials.AgentShape;
    patientShapes = all_materials.PatientShape;
    itemNumbers = all_materials.ItemNumber;
    numSentences = length(sentences);

    
    %Save file in image folder
    IMAGE_DIR = fullfile(pwd, 'images', 'resized_sentences');
    
    item_index = 1;
    
    for i=1:numSentences
        condition = conditions{i};
        if ~strcmp(char(condition),'NULL ')
            %Extract info for the item we're currently on
            text = sentences{i};
            flip = flips{i};
            agent = [agentNames{i} ' ' agentShapes{i}];
            patient = [patientNames{i} ' ' patientShapes{i}];
            itemNumber = itemNumbers{i};
            
            %Determine person to be highlighted
            switch condition
                case 'Agent'
                    person = agent;
                case 'Patient'
                    person = patient;
            end
            
            %Determine flip (for file naming purposes)
            switch flip
                case '0'
                    flip_word = 'active';
                case '1'
                    flip_word = 'passive';
            end

            %Spaces don't work out uniformly, so adjust for each name
            %Will want to switch for the one being highlighted
            %TWo adjustment options, one for 1st position in sentence and one for 2nd
            switch person
                case 'Melissa Oval'
                    %adjustment = 5;
                    adjustment1 = 2;
                case 'Lily Triangle'
                    %adjustment = 5;
                    adjustment1 = 0;
                case 'Kyle Square'
                    %adjustment = 4;
                    adjustment1 = 1;
                case 'Zach Star'
                    %adjustment = 4;
                    adjustment1 = 1;
            end
    
            
            %Determine whether person comes 1st or 2nd in sentence
            if xor(strcmp(char(condition), 'Agent'), strcmp(flip, num2str(1)))
                position = 1;
            else
                position = 2;
            end
            
            %Determine where highlight starts and ends, based on position
            switch position
                case 1
                    %Determine where to start the highlight
                    highlight_start = 275-adjustment1;

                case 2
                    highlight_start = 275 + EMS_TO_PIXELS*(length(text)-length(person));
            end
            
            %Make a string with as many spaces as there are characters in person
            highlight_box = blanks(length(person)-1);

            %Sets up image and overlays text
            I = imread('blank-white-rectangle.png');
            %I = imresize(I, [winHeight, NaN]);
            position_text = [250 1000];
            position_highlight = [highlight_start 1000];

            RGB = insertText(I,position_text,text,'FontSize',FONT_SIZE,'BoxOpacity',0,'Font','Courier');
%           RGB = insertText(RGB,position_highlight,highlight_box,'FontSize',FONT_SIZE,'BoxOpacity',.4,'Font','Courier');

            %Sets up file to save; numbers indicate index at which stimulus
            %was presented
%           fileToSave = [condition '_' flip_word '_' itemNumber '.jpg'];
            fileToSave = ['Base_' flip_word '_' itemNumber '.jpg'];
            fileToSave = fullfile(IMAGE_DIR, fileToSave);

            %Displays text image
            figure
            imshow(RGB)

            %Saves text image
%           imwrite(RGB,fileToSave,'jpg')
            if strcmp(condition, 'Agent')
                imwrite(RGB,fileToSave,'jpg')
            end

            %Increments counter
            item_index = item_index + 1;
        end
        
    end

    
end