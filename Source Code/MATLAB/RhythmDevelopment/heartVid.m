% Request the directory for saving the file
        savedir = uigetdir;
        % If the cancel button is selected cancel the function
        if savedir == 0
            return
        end
        % Request the desired name for the movie file
        filename = inputdlg('Enter Filename:');
        filename = char(filename);
        % Check to make sure a value was entered
        if isempty(filename)
            error = 'A filename must be entered! Function cancelled.';
            msgbox(error,'Incorrect Input','Error');
            return
        end
        filename = char(filename);
        % Create path to file
        movname = [savedir,'/',filename,'.avi'];
        % Start writing the video
        
        vidObj = VideoWriter(movname,'Motion JPEG AVI');
        open(vidObj);
        
        imgdir = uigetdir('Select Image Directory');
        cd(imgdir)
        
        fileList = dir;
            % grab files that are tiffs
            checkFiles = zeros(size(fileList,1),1);
            for n = 1:length(checkFiles)
                if length(fileList(n).name) > 4
                    checkFiles(n) = strcmp(fileList(n).name(end-3:end),'tiff');
                else
                    checkFiles(n) = 0;
                end
            end
            % grab indices of the files that are tiffs
            checkFiles = checkFiles.*(1:length(checkFiles))';
            checkFiles = unique(checkFiles);
            checkFiles = checkFiles(2:end);
            % remove directories from file list
            fileList = fileList(checkFiles);
            
            % identify period that separates the name and file type
            charCheck = zeros(length(fileList(1).name),1);
            for n = 1:length(charCheck)
                % char(46) is a period
                charCheck(n) = fileList(1).name(n) == char(46);
                if charCheck(n) == 1
                    middleInd = n;
                    break
                end
            end
            f = figure; 
            a = imread(fileList(1).name); 
            imshow(a);
            [a_cropped, rect] = imcrop(a);
            writeVideo(vidObj,a_cropped);
            writeVideo(vidObj,a_cropped);
            writeVideo(vidObj,a_cropped);
            writeVideo(vidObj,a_cropped);
            writeVideo(vidObj,a_cropped);
            for n = 2:size(fileList)
                a = imread(fileList(n).name);
                a_cropped = imcrop(a, rect);
                writeVideo(vidObj,a_cropped);
                writeVideo(vidObj,a_cropped);
                writeVideo(vidObj,a_cropped);
                writeVideo(vidObj,a_cropped);
                writeVideo(vidObj,a_cropped);
            end
            close(vidObj);
            