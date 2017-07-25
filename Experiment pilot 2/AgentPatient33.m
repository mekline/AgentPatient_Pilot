%AgentPatientStimuli.m
% 
% DESCRIPTION
% Subjects view stimuli - either pictures or sentences ? about one shape
% performing an action on another (e.g., 'Melissa Oval is being bounced by
% Kyle Square'). Sometimes the agent is highlighted; other times the
% patient is highlighted.
% 
% Run PsychStartup; in command line, then call the function
% 
% Function call: AgentPatient33(subjID, image_type, order, run)
% 
% RUNTIME: should be 480 sec (8 min, 240 TRs) ##TOUPDATE
%          -runs 10 times faster with subjID 'debug'
%          -actually takes about 1475 sec due to lag
% 
% INPUTS:
%   -subjID: subject ID string
%   -image_type: 'sentences' for sentence stimuli, 'stills' for picture
%   stimuli
%   -order: A, B, C D, or E; determines which preset item ordering to use
%   -(Each participant should get at least XXX runs, each with a different
%   ordering)
% 
% OUTPUT:
%   -.csv file with information on what was presented when to whom for how
%   long (data/AgentPatient33_subjID_image_type_order_data.csv)
% 
% ADJUSTMENTS:
% To change font size, go to Set display options.
% If the aspect ratio is off, change SCREEN_ADJUST in Set experiment
% constants until it looks better.

function AgentPatient33(subjID, image_type, order)



    start_time = GetSecs;
    

    %% Make sure inputs are valid and raise an error otherwise
    %subjID is a string
    assert(ischar(subjID), 'subjID (first argument) must be a string');

    %image_type is 'sentences' or 'stills'
    assert(strcmp(image_type, 'sentences') | strcmp(image_type, 'stills'), sprintf('image_type (second argument) must be ''sentences'' or ''stills''\nExample of correct format: AgentPatientStimuli(''myID'',''sentences'',''A'',1)'));
    
    %order is a letter A-E
    assert(ismember(order, ['A','B','C','D','E']), sprintf('order (third argument) must be a letter A-E enclosed in quotes\nExample of correct format: AgentPatientStimuli(''myID'',''sentences'',''A'',1)'));
 
    %% Make sure we don't accidentally overwrite a data file
    %This is where the data file will go
    DATA_DIR = fullfile(pwd, 'data');
    %This is what we'll call the data file we're making
    fileToSave = ['AgentPatient33_' subjID '_' image_type '_' order '_data.csv'];
    %The file should be in the DATA_DIR folder
    fileToSave = fullfile(DATA_DIR, fileToSave);
    
    % Error message if data file already exists.
    if exist(fileToSave,'file') && ~strcmpi(subjID, 'debug')
        str = input('The data file already exists for this subject! Overwrite? (y/n)','s');
        if ~isequal(str,'y')
            error('The data file already exists for this subject!');
        end
	end
    
        
    %% Set experiment constants
    %Timing (in seconds)              
    FIX_DUR     = 0.0; %Length of trial-initial fixation
    ITI         = 0.0; %Inter-trial interval
    BLINK_DUR   = 0.1; %Length of one blink on or off
    TIME_TO_Q   = 4.0; %delay until question, for question trials
    TIME_Q      = 2.0; %time question visible, for question trials
    SCREEN_ADJUST = 1.2; %factor to adjust aspect ratio by
    
    %Based on image type (sentences or stills), decide what Flip means
    switch image_type
        case 'sentences'
            %all_materials.Flip is just 0 or 1; this is what it means in each case
            flip_word_0 = 'active';
            flip_word_1 = 'passive';
        case 'stills'
            flip_word_0 = 'orig';
            flip_word_1 = 'flipped';
    end
    
    %% Set display options
    %Font sizes
    sentFontSize = 40;      %question sentences
    fixFontSize = 40;       %fixation cross
    
    %% Set up orders
    ORDER_DIR = fullfile(pwd, 'orders'); %folder containing all orders
    
    order_filename = ['AgentPatient33_Order_' order '.csv']; %order for this run
    order_filename = fullfile(ORDER_DIR, order_filename); %full path for that order
    all_materials = readtable(order_filename); %order stored as a table
    
    numEvents = height(all_materials); %the total number of trials and fixations
    
    %% Make the experiment run faster if subjID is 'debug'
    if strcmpi(subjID, 'debug')
        scale = 0.02;
        FIX_DUR = FIX_DUR * scale;
        BLINK_DUR = BLINK_DUR * scale;
        ITI = ITI * scale;
        
        all_materials.IntendedOnset = all_materials.IntendedOnset * scale;
        all_materials.IntendedDuration = all_materials.IntendedDuration * scale;
    end
    
    %% Read in the stimuli materials
    %Save all_materials to a matfile
    MATERIALS_DIR = fullfile(pwd, 'data'); %where to put the saved materials
    mat_filename = ['AgentPatient33_' subjID '_' order '_materials.mat']; %materials to save
    mat_filename = fullfile(MATERIALS_DIR, mat_filename);
    save(mat_filename, 'all_materials'); %save all_materials as mat_filename
    
    %Load the all_materials table from the mat file
    %If the mat file doesn't exist, make sure the user entered the correct
    %inputs
    try
         load(mat_filename);
         
    catch errorInfo
         fprintf('%s%s\n\n', 'error message: ', errorInfo.message)
         
         error('\n%s\n\t%s\n\t%s\n', ...
               'Please make sure the following conditions are met:', ...
               '1) subjID is the same for run 1 and run 2', ...
               '2) run is 1 for the first run and 2 for the second');
    end
     
     %Set up the data that we want to save
     resultsHdr = {'Trial', 'SubjID',        'Run',       'Order',   'ItemNumber'   'IntendedOnset', ...
                   'ActualOnset', 'IntendedDuration', 'ActualDuration',      'Condition', 'Flip', 'FlipMeaning',  'Sentence', ...
                   'AgentName',     'AgentShape',          'PatientName'...
                   'PatientShape', 'BaseFilename', 'BlinkFilename', 'Response', 'Correct'};
 	
     %Set up results, the table that will hold this data
     results = cell(numEvents, length(resultsHdr));
     results = cell2table(results, 'VariableNames', resultsHdr);
    
    %Fill in the user input information
    results.SubjID(:) = {subjID};
    results.Order(:) = {order};
    
    %Fill in Condition anxd Flip in the results file, one by one
    for eventNum=1:numEvents
        results.Condition{eventNum} = char(all_materials.Condition(eventNum)); %Use base condition to get the right stim files!
        results.Flip{eventNum} = char(all_materials.Flip(eventNum));
        results.ItemNumber{eventNum} = char(all_materials.ItemNumber(eventNum));
        results.IntendedOnset{eventNum} = all_materials.IntendedOnset(eventNum);
        results.IntendedDuration{eventNum} = all_materials.IntendedDuration(eventNum);
    end

	%% Set up screen and keyboard for Psychtoolbox
    %Screen
    screenNum = max(Screen('Screens'));  %Highest screen number is most likely correct display
    windowInfo = PTBhelper('initialize',screenNum);
	wPtr = windowInfo{1}; %pointer to window on screen that's being referenced
    rect = windowInfo{2}; %dimensions of the window
        winWidth = rect(3);
        winHeight = rect(4)*SCREEN_ADJUST;
    oldEnableFlag = windowInfo{4};
    HideCursor;
    PTBhelper('stimImage',wPtr,'WHITE');
    PTBhelper('stimText',wPtr,'Loading images\n\n(Don''t start yet!)',30);
     
    %Keyboard
    keyboardInfo = [];
    keyboardInfo = PTBhelper('getKeyboardIndex');
    kbIdx = [keyboardInfo{1}];
    escapeKey = keyboardInfo{2};
    
    %% Set up cells containing image file data
    %Initialize the cells and IMAGE_DIR
    
    %images with highlight
    img_files = cell(1,72);
    img_stims = cell(1,72);
    
    %images with no highlight
    base_files = cell(1,72);
    base_stims = cell(1,72);
    
    IMAGE_DIR = fullfile(pwd, 'images', image_type); %the folder we're taking images from

    index = 1; %counter
    load_start_time = GetSecs;
    for eventNum=1:numEvents
        condition = all_materials.Condition(eventNum); %TOUPDATE: Is this coming from the right condition? 
        intendedDuration = all_materials.IntendedDuration(eventNum);
        flip = all_materials.Flip{eventNum}; 
        
        if ~strcmp(char(condition), 'Fixation') %if this trial isn't a fixation
            %Determine flip
            switch flip
                case '0'
                    flip_word = flip_word_0;
                case '1'
                    flip_word = flip_word_1;
            end

            %what is the name of the image we want
            img_name = [char(condition) '_' flip_word '_' char(all_materials.ItemNumber{eventNum}) '.jpg'];
            base_name = ['Base_' flip_word '_' char(all_materials.ItemNumber{eventNum}) '.jpg'];
            
            
            %put it into a cell with the other image files
            img_files{1,index} = fullfile(IMAGE_DIR, img_name);
            base_files{1,index} = fullfile(IMAGE_DIR, base_name);
                
            %read in the base and highlight file we want
            base = base_files{index};
            base = imread(base,'jpg');
            image = img_files{index};
            image = imread(image,'jpg');
            fclose('all');
            
            %make it a texture so PTBHelper will like it
            img_stims{index} = Screen('MakeTexture', wPtr, double(image));
            base_stims{index} = Screen('MakeTexture', wPtr, double(base));

        
        results.BlinkFilename{eventNum} = char(img_name);
        results.BaseFilename{eventNum} = char(base_name);
        results.FlipMeaning{eventNum} = char(flip_word);
        results.Trial{eventNum} = index;
        index = index+1; %increment counter, needed to keep images in the right spot on the image list
        
        PTBhelper('stimText',wPtr,['Loading images\n\n(Don''t start yet!)\n' num2str(index) '/120'],30);
 
        end
    end 
    load_end_time = GetSecs;
    load_time = load_end_time - load_start_time
  

    
    %% Present the experiment
	% Wait indefinitely until trigger
    PTBhelper('stimText',wPtr,'Waiting for trigger...',sentFontSize);
    PTBhelper('waitFor','TRIGGER',kbIdx,escapeKey);
    
    runOnset = GetSecs; %remains the same
    item_index = 1;
    
    %Present each event
    try
        for eventNum = 1:numEvents
            
            condition = all_materials.Condition(eventNum);
            intendedOnset = all_materials.IntendedOnset(eventNum);
            intendedDuration = all_materials.IntendedDuration(eventNum);
            intendedOffset = intendedOnset + intendedDuration;
            flip = all_materials.Flip{eventNum};
            optCond = all_materials.OptCondition(eventNum); %Conditions again, marked with whether it's a question. Ugly!!
            
            actualOnset = GetSecs-runOnset; %What time is it right now?
   
            %If it's fixation
            if strcmp(condition, 'Fixation')

                %Show fixation cross
                PTBhelper('stimText', wPtr, '+', fixFontSize);
                fixEndTime = runOnset + intendedOffset;
                PTBhelper('waitFor',fixEndTime,kbIdx,escapeKey);
                
                %Save data
                
                results.Sentence{eventNum} = 'NA';
                results.AgentName{eventNum} = 'NA';
                results.AgentShape{eventNum} = 'NA';
                results.PatientName{eventNum} = 'NA';
                results.PatientShape{eventNum} = 'NA';
                results.Response{eventNum} = 'NA';
                results.Correct{eventNum} = 'NA';
                
            %If there's a sentence to be presented (i.e., not NULL)   
            else
                if char(all_materials.Flip(eventNum)) == '0'
                    sentence = char(all_materials.ProgressiveActive(eventNum));
                elseif char(all_materials.Flip(eventNum)) == '1'
                    sentence = char(all_materials.ProgressivePassive(eventNum));
                end
                
                %A question trial?
                question = 0;
                if sum(strcmp(optCond, {'AgentQ','PatientQ'}))
                    question = 1;
                    theQ = char(all_materials.Question(eventNum)); %XXXXXSTART HERE!
                    theA = char(all_materials.Answer(eventNum));
                    
                    disp(theQ)
                    disp(theA)
                end
                
                eventEndTime = runOnset + intendedOffset;
                if question
                    sentEndTime = eventEndTime - TIME_TO_Q - TIME_Q;
                    questionTime = eventEndTime - TIME_TO_Q;
                else
                    sentEndTime = eventEndTime;
                end
                    
                %Blink sentence until !!ALMOST!! sentEndTime
                while GetSecs < (sentEndTime-3*BLINK_DUR) %Any fewer than this tends to lead to long bleedovers that have to be corrected/uneven trial lengths
                    PTBhelper('stimImage', wPtr, item_index, img_stims);
                    WaitSecs(BLINK_DUR);
                    PTBhelper('stimImage', wPtr, item_index, base_stims);
                    WaitSecs(BLINK_DUR);
                end
                
                %Play a question if applicable!
                if question
                    PTBhelper('stimText', wPtr, '+', fixFontSize);
                    PTBhelper('waitFor',questionTime,kbIdx,escapeKey);
                    %PTBhelper('stimText', wPtr, %%%QUESTION TEXT%%%%% , fixFontSize);
                    %AND RECORD INPUT!!
                end

                
                              
                %Save Sentence, agent info to results
                results.Sentence{eventNum} = sentence;
                results.AgentName{eventNum} = char(all_materials.AgentName(eventNum));
                results.AgentShape{eventNum} = char(all_materials.AgentShape(eventNum));
                results.PatientName{eventNum} = char(all_materials.PatientName(eventNum));
                results.PatientShape{eventNum} = char(all_materials.PatientShape(eventNum));
               
                
                %Update loop variables (only for active stim presentation
                %events!)
                item_index = item_index + 1;
                
                %Catch up to the end of the trial.
                PTBhelper('waitFor',eventEndTime,kbIdx,escapeKey);
                
            end
            
            %%Save actual onset and duration back to the results file!
            %%(v. useful for debugging)
            results.ActualOnset{eventNum} = actualOnset;
            actualDuration = GetSecs - (actualOnset + runOnset);
            results.ActualDuration{eventNum} = actualDuration;
            
        end

        ran_completely = true;
        
    catch errorInfo
        ran_completely = false;
        
        fprintf('%s%s\n\n', 'error message: ', errorInfo.message)
        for k=1:length(errorInfo.stack)
            disp(errorInfo.stack(k))
        end
    end
    
    %Save all data
	writetable(results, fileToSave);
    
     %Close the PTB screen
	Screen('CloseAll');
	ShowCursor;
    
    %Restore the old level.
    Screen('Preference','SuppressAllWarnings',oldEnableFlag);
    disp(['Run ' order ' finished; data for this run saved to ' fileToSave])
end

%% %% Debugging functions
function [wPtr, rect] = openDebugWindow(screenNum, rect)
    Screen('CloseAll');
    ShowCursor;
    clear Screen
    
    rect = rect / 2;
    rect(1) = 5;
    rect(2) = 5;

    java; %clear java cache
    KbName('UnifyKeyNames');
    warning('off','MATLAB:dispatcher:InexactMatch');
    AssertOpenGL;
    suppress_warnings = 1;
    Screen('Preference', 'SuppressAllWarnings', suppress_warnings);
    Screen('Preference', 'TextRenderer', 0);
    Screen('Preference', 'SkipSyncTests', 1);
    [wPtr,rect] = Screen('OpenWindow',screenNum,1,rect,[],[],[],[],[],kPsychGUIWindow,[]);
end

        

