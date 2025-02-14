function [blueLocs, outputCellArray, bbox, bboxResize]= avgPixelMapFinal(filename1)
% gather all names into a cell array
files= [{filename1}];
for x=1:length(files)

    fname= files{x};
    img= imread(fname);
    img2= imread(fname);
    %size
    [r,c,l]= size(img);

    % create background for collage
    border= 20;
    collage= uint8(zeros(r + 2.*border, 3*c + 4*border,l));
    collage(:,:,1)= 255;
    collage(:,:,2)=255;
    collage(:,:,3)=255;

    %layers
    justR= img(:,:,1);
    justG= img(:,:,2);
    justB= img(:,:,3);

    % linearize
    redVec= justR(:);
    greenVec= justG(:);
    blueVec= justB(:);

    % averages
    redAvg= round(mean(redVec),2);
    bluAvg= round(mean(blueVec),2);
    greenAvg= round(mean(greenVec),2);

    % Image Name
    imageName= strtok(fname, '.');

    % Plt each Channel on a white background with borders
    allBlack = zeros(size(img,1,2),class(img));
    justRed = cat(3,justR,allBlack,allBlack);
    [r3,c3,l3]= size(justR);
    collage(border+1:r + border, border + 1: c + border, :)= justRed;
    justGreen = cat(3,allBlack,justG,allBlack);
    collage(border+1:r + border, 2*border + c + 1: 2*c + 2*border,:)= justGreen;
    justBlue = cat(3,allBlack,allBlack,justB);
    collage(border+1:r + border, 3*border + 2*c + 1: 3*c + 3*border,:)= justBlue;

    subplot(1,1,x)
    imshow(collage)
    title(sprintf('%s- [Red Avg: %d | Green Avg: %d | Blue Avg: %d]', ...
    imageName, redAvg, greenAvg, bluAvg))
end
[rowsJustB,colJustB,layersJustB]= size(justB);
imshow(justB)
mat= ones(rowsJustB,colJustB).*-1;
for x=1:rowsJustB
    for y=1:colJustB
        if x>2 && y>2 && x+2<=rowsJustB && y+2 <= colJustB
            tR= x-2;
            tC= y-2;
            bR= x+2;
            bC= y+2;
            mat2= justB(tR:bR,tC:bC);
            pixelAvg= mean(mat2(:));
            mat(x,y)= pixelAvg;
        else
            currentPixel = justB(x, y);
            pixels = currentPixel; % Start with the current pixel
           
            % Attempt to access each neighboring pixel with try-catch
            try pixelAbove = justB(x-1, y); catch, pixelAbove = []; end
            try pixel2Above = justB(x-2, y); catch, pixel2Above = []; end
            try pixelRight = justB(x, y+1); catch, pixelRight = []; end
            try pixel2Right = justB(x, y+2); catch, pixel2Right = []; end
            try pixelLeft = justB(x, y-1); catch, pixelLeft = []; end
            try pixel2Left = justB(x, y-2); catch, pixel2Left = []; end
            try pixelDown = justB(x+1, y); catch, pixelDown = []; end
            try pixel2Down = justB(x+2, y); catch, pixel2Down = []; end
            try pixelLeftDiag = justB(x-1, y-1); catch, pixelLeftDiag = []; end
            try pixelRightDiag = justB(x-1, y+1); catch, pixelRightDiag = []; end
            try pixelBotRightDiag = justB(x+1, y+1); catch, pixelBotRightDiag = []; end
            try pixelBotLDiag = justB(x+1, y-1); catch, pixelBotLDiag = []; end
    
            % Adding knight's move pixels
            try knightMove1 = justB(x-2, y+1); catch, knightMove1 = []; end % Up 2, Right 1
            try knightMove2 = justB(x-1, y+2); catch, knightMove2 = []; end % Up 1, Right 2
            try knightMove3 = justB(x+1, y+2); catch, knightMove3 = []; end % Down 1, Right 2
            try knightMove4 = justB(x+2, y+1); catch, knightMove4 = []; end % Down 2, Right 1
            try knightMove5 = justB(x+2, y-1); catch, knightMove5 = []; end % Down 2, Left 1
            try knightMove6 = justB(x+1, y-2); catch, knightMove6 = []; end % Down 1, Left 2
            try knightMove7 = justB(x-1, y-2); catch, knightMove7 = []; end % Up 1, Left 2
            try knightMove8 = justB(x-2, y-1); catch, knightMove8 = []; end % Up 2, Left 1
           
            % Collect all available pixels (excluding empty values)
            pixels = [pixels; pixelAbove; pixel2Above; pixelRight; pixel2Right; pixelLeft; pixel2Left; pixelDown; pixel2Down; pixelLeftDiag; pixelRightDiag; pixelBotRightDiag; pixelBotLDiag; knightMove1; knightMove2; knightMove3; knightMove4; knightMove5; knightMove6; knightMove7; knightMove8];
           
            % Calculate the average
            pixelAvg = mean(pixels);
            mat(x,y)= pixelAvg;
        end
    end
end
output= mat;
thresholdMax= output >= (floor(max(mat(:)))./8)
grayedImg = mat2gray(uint8(thresholdMax)*255);
newImg= imbinarize(grayedImg);
figure
imshow(grayedImg)

% use region props to get the coordinates of the boundings boxes (outputs
% structure)
    s2= regionprops(newImg, 'centroid');
    centBlue= cat(1,s2.Centroid);
    [aRows,aCols,l]= size(grayedImg);
    xs= max(min(floor(centBlue(:, 1)), aCols), 1);
    ys= max(min(floor(centBlue(:, 2)), aRows), 1);
    [cRows,cCols]= size(centBlue);
    pointTracker= {};
    nBlocs= [xs,ys]
    for h1= 1:cRows
        startPoint= floor(nBlocs(h1,:));
        xC= startPoint(2); %ROWS
        yC= startPoint(1); %COLS
        startX1= startPoint(2); % by dimensions of the image
        startY1= startPoint(1);
        startX2= startPoint(2);
        startY2= startPoint(1);

% then start counting in each direction (right) and collecting points
        va= startY1;
        found= false;
        ptCt1= 0;
        while va <= aCols && ~found
            if grayedImg(xC,va)>=0.6
                va = va + 1;
            else
                found= true;
                rightPoint= [va,xC];
                ptCt1= ptCt1 + 1;
            end
        end
        var= startY2;
        found2= false;
        ptCt2= 0;
        while var <= aCols & ~found2
            if grayedImg(xC,var)>=0.6
                var=var - 1;
            else
                found2= true;
                leftPoint= [var,xC];
                ptCt2= ptCt2 + 1;
            end
        end
% then start counting in each direction (top) and collecting points
        v= startX1;
        found3= false;
        ptCt3= 0
        while v <= aRows & ~found3
            if grayedImg(v,yC)>=0.6;
                v= v + 1;
            else
                found3= true;
                downPoint= [yC,v];
                ptCt3= ptCt3 + 1;
            end
        end
        %disp(storedPoint3);
       % plot(storedPoint3, 'b*')
    % then start counting in each direction (bottom) and collecting points
        vs= startX2;
        found4= false;
        ptCt4= 0
        while vs <= aRows & ~found4
            if grayedImg(vs,yC)>=0.6;
                vs=vs - 1;
            else
                found4= true;
                upPoint= [yC,vs];
                ptCt4= ptCt4+1;
            end
        end
        upRight = [rightPoint(1), upPoint(2)];
        botRight= [rightPoint(1), downPoint(2)];
        botLeft= [leftPoint(1), downPoint(2)];
        upLeft= [leftPoint(1), upPoint(2)];
        points= [upLeft;upRight;botRight;botLeft;upLeft];
        points2= [upLeft;upRight;botRight;botLeft];
        spx1s= points(:,1);
        spy1s= points(:,2);
        %subplot(4,aRows,h1)

        %subplot()
        hold on;
        plot(spx1s,spy1s, 'b-')
    end
    %%
    blueLocs = [ys, xs];
    [rowsBlueLocs, colsBlueLocs] = size(blueLocs);
    MagicRed = 30;
    MagicGreen = 30;
    allThree = {};
    redBlue = {};
    greenBlue = {};
    pureBlue= {};
    for h2= 1:cRows
        startPoint= floor(nBlocs(h2,:));
        xC= startPoint(2); %ROWS
        yC= startPoint(1); %COLS
        startX1= startPoint(2);
        startY1= startPoint(1);
        startX2= startPoint(2);
        startY2= startPoint(1);

% then start counting in each direction (right) and collecting points
        va= startY1;
        found= false;
        ptCt1= 0;
        % could maybe add to this conditional to include things from other
        % cell clqassifications. Ex if BW(xC,va,1) >= avgVec(2) from red
        % green
        while va <= aCols && ~found
            if grayedImg(xC,va)>=0.6
                va = va + 1;
            else
                found= true;
                rightPoint= [va,xC];
                ptCt1= ptCt1 + 1;
            end
        end
       % disp(storedPoint);
       % plot(storedPoint, 'b*')
% then start counting in each direction (left) and collecting points
        var= startY2;
        found2= false;
        ptCt2= 0;
        while var <= aCols & ~found2
            if grayedImg(xC,var)>=0.6
                var=var - 1;
            else
                found2= true;
                leftPoint= [var,xC];
                ptCt2= ptCt2 + 1;
            end
        end
% then start counting in each direction (top) and collecting points
        v= startX1;
        found3= false;
        ptCt3= 0;
        while v <= aRows & ~found3
            if grayedImg(v,yC)>=0.6;
                v= v + 1;
            else
                found3= true;
                downPoint= [yC,v];
                ptCt3= ptCt3 + 1;
            end
        end
        %disp(storedPoint3);
       % plot(storedPoint3, 'b*')
    % then start counting in each direction (bottom) and collecting points
        vs= startX2;
        found4= false;
        ptCt4= 0;
        while vs <= aRows & ~found4
            if grayedImg(vs,yC)>=0.6;
                vs=vs - 1;
            else
                found4= true;
                upPoint= [yC,vs];
                ptCt4= ptCt4+1;
            end
        end

        upRight = [rightPoint(1), upPoint(2)];
        botRight= [rightPoint(1), downPoint(2)];
        botLeft= [leftPoint(1), downPoint(2)];
        upLeft= [leftPoint(1), upPoint(2)];
        points= [upLeft;upRight;botRight;botLeft;upLeft];
        points2= [upLeft;upRight;botRight;botLeft];
        upPointMid= (upLeft(1)+botLeft(1))./2;
        distUp= ((upLeft(1)-upPointMid));
        distHorz= (830*(upLeft(2)-(upLeft(2)+upRight(2))./2));
        newTl= [upLeft(1)-distUp,upLeft(2)-distHorz];
        newTR= [upRight(1)-distUp,upRight(2)+distHorz];
        newBR= [botRight(1)+distUp,botRight(2)+distHorz];
        newBL= [botLeft(1)+distUp,botLeft(2)-distHorz];

        snippet= img2(newTl(2):newBL(2),newTl(1):newTR(1),:);
        snippetR= snippet(:,:,1)
        snippetG= snippet(:,:,2)
        snippetB= snippet(:,:,3)

        blackSquare= zeros(301,301);
        upDown= upLeft(1)-botLeft(1);
        leftRight= upLeft(2)-upRight(2);
        %backSquare(151-(upDown./2):151+(upDown./2),150-(leftRight./2):150+(leftRight./2))= snippet
        [rowSnip,colSnip,l]= size(snippet);
        rowAdded= 70-rowSnip;
        colAdded= 70-colSnip;
        newSnip= padarray(snippet,[rowAdded colAdded],uint8(0),'both');

        figure
        %% ask abt imlocalbrighten
        subplot(1,4,1)
        imshow(imlocalbrighten(newSnip,0.3));
    
        subplot(1,4,2)
        imshow(imlocalbrighten(newSnip(:,:,1),0.3));
    
        subplot(1,4,3)
        imshow(imlocalbrighten(newSnip(:,:,2),0.3));
    
        subplot(1,4,4)
        imshow(imlocalbrighten(newSnip(:,:,3),0.3));

        blueLocs = [ys, xs]
        MagicRed = 30;
        MagicGreen = 30;
        startPointBlueLocs = blueLocs(h2,:);
        startPointBlueLocsX = startPointBlueLocs(1);
        startPointBlueLocsY = startPointBlueLocs(2);

        if any(snippetR(:) >= MagicRed) && any(snippetG(:)>= MagicGreen) % 2. Corrected indexing
            sgtitle('This cell has been classified as meeting All Three Criteria')
            allThree = [allThree {blueLocs(h2,:)}]
        elseif any(snippetR(:) >= MagicRed)
            sgtitle('This cell has been classified as meeting the Red Blue Criteria')
            redBlue = [redBlue {blueLocs(h2,:)}]
        elseif any(snippetG(:)>= MagicGreen)
            sgtitle('This cell has been classified as meeting the Green Blue Criteria')
            greenBlue = [greenBlue {blueLocs(h2,:)}]
        else
            sgtitle('This cell has been classified as meeting the Pure Blue criteria')
            pureBlue = [pureBlue {blueLocs(h2,:)}]
        end
    end
    allThree= [{'All Three'};allThree']
    redBlue=  [{'Red Blue'};redBlue']
    greenBlue= [{'Green Blue'};greenBlue']
    pureBlue= [{'Pure Blue'}; pureBlue']
        %{
    %%
     while indexX<=rowsBlueLocs % 1. Corrected indexing
        startPointBlueLocs = blueLocs(indexX,:);
        startPointBlueLocsX = startPointBlueLocs(1);
        startPointBlueLocsY = startPointBlueLocs(2);
        if img(startPointBlueLocsX, startPointBlueLocsY, 1) >= MagicRed && img(startPointBlueLocsX, startPointBlueLocsY, 2) >= MagicGreen % 2. Corrected indexing
            allThree = [allThree {[startPointBlueLocsX,startPointBlueLocsY]}]; % 3. Corrected appending
            indexX= indexX + 1;
        elseif img(startPointBlueLocsX,startPointBlueLocsY,1) >= MagicRed
            redBlue= [redBlue {[startPointBlueLocsX,startPointBlueLocsY]}]; % 3. Corrected appending
            indexX= indexX + 1;
        elseif img(startPointBlueLocsX, startPointBlueLocsY, 2) >= MagicGreen
            greenBlue= [greenBlue {[startPointBlueLocsX,startPointBlueLocsY]}]; % 3. Corrected appending
            indexX= indexX + 1;
        else
            pureBlue= [pureBlue {[startPointBlueLocsX,startPointBlueLocsY]}];
        end
        indexX= indexX+1;
    end
%}
    %% Take a snapshot of the surrounding pixels
    % take distance from centroid to tope side of the box

end