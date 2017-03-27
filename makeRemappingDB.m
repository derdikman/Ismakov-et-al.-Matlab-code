function makeRemappingDB(excel_file,start_num)
%
%
% makeRemappingDB
%
%
% read an excel file containing the file and position data, and generate a
% data-base containing this information.
%
% excel_file - name of file containing unit data
% start_num - starting from the <start_num> line in the excel_sheet
% (deafult = 1)
%
% 2nd version, which keeps the data structure simple
%
% Excel contains the following columns:
% cell x trial id: rat cell date id 


% rat       - number of rat -- first 3 letters
% date      - run date -- letters 4 to 7
% A1        - suffix of first set of files (open field)
% A2        - suffix of 4th set of files (open field)
% exp_type - info on the type of experiment
% food      - info on food administration
% tetrode   - tetrode number (i.e. t6)
% cell      - cell number (i.e. t6c1)
% arena     - name of arena (room10old, room9 etc...)
% histology 
% ... additional optional fields ...
%
% Created by Dori Derdikman, CBM, Feb 2008
% edited for Kate Jefferey remapping data by Rebekkah, Apr 2016

dbstop if error 

if ~exist('start_num')
    start_num = 1;    
end

% various parameters
excel_file= 'Sum.xlsx'; 


parms.data_dir = '\\192.114.21.198\Dori_Data\data\rebekkah\All rats';  % root of data directory
parms.dir_save_data= '\\192.114.21.198\Dori_Data\data\rebekkah\All rats\all cells with no pos data fixing';
% parms.dim = 160; % dimensions of box + 10   % ???
parms.excel_file = excel_file;  % excel file containing the data-base data

cd(parms.data_dir);

% read excel sheet

excel_sheet = read_excel(excel_file);

% create db structure from excel sheet

db_all = build_db_from_excel(excel_sheet,parms);

% read spike data and pos data into db structure

ncells = length(db_all);

% hProg = progbar('start'); Progress = 0;

for cell_num = start_num: ncells %261:283 %97:ncells %start_num: ncells

    rec = db_all(cell_num);
    
    % read pos data
    rec = HP_read_pos_data(rec,parms);
            
    % read spike data
    rec = HP_read_spike_data(rec,parms);

    
%     db_all(cell_num) = rec;
%     save('HP_db','db_all');
    
    
   % ind = find(excel_file == '.');
   % excel_prefix = excel_file(1:ind-1);
        
    key = [rec.rat '_' rec.date '_t'  ...
           num2str(rec.tetrode) 'c' num2str(rec.cell) '_tr' rec.trial_num];
    db = rec;
    if length(rec.rat) == 3
        file_name = key;
    elseif length(rec.rat) >3    
        file_name= key(2:end);
    end

    cd(parms.dir_save_data)
    save(file_name,'db','parms','db_all');
    cd(parms.data_dir)
    disp(['===> save cell in file ' file_name ' <===']);
    
%    Progress = cell_num/ncells; progbar( 'update', hProg, Progress );
 
    pack;
end               

% Progress = 1; progbar( 'close', hProg, Progress );

disp('finished');



disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function db = build_db_from_excel(excel_sheet,parms)
%
% copy fields from excel sheet to DB
%
for i = 1:length(excel_sheet)
    
   % db(i) = excel_sheet(i);
    db(i).rat = excel_sheet(i).Rat;
    db(i).date = excel_sheet(i).Date;
    db(i).trial_num = excel_sheet(i).Trial_no;
    db(i).envt_type = excel_sheet(i).Envt_type;
    db(i).tetrode = excel_sheet(i).Tetrode;
    db(i).cell = excel_sheet(i).Cluster;
    db(i).arena = excel_sheet(i).TriaL1_type;
    
    %db(i).tetrode = str2num(db(i).tetrode(2)); % this is done because of the strange way it is kept in the excel
    %db(i).cell = str2num(db(i).cell(4)); % assume < 10 cells (which is ok for my data).
    
    %
    % determine files, cut files, sessions, and cut_file in structure(AB)
    % later pos_data and spike_data will also be put in this structure.
    %

    % open rat folder
    % open date folder
    % open appropriate tetrode and cluster folder
    
    rat_num = sprintf('r%s', excel_sheet(i).Rat);
    date_num = sprintf('%s', excel_sheet(i).Date);
    
    if isempty(date_num)
        cut_file= '';
    else
        cut_file= sprintf('%s%s%s%s_%s.cut',excel_sheet(i).Rat, date_num(4:5), ...
        date_num(1:2),excel_sheet(i).Trial_no,excel_sheet(i).Tetrode); 
    end 
    
    if length(rat_num) == 4
    db(i).cut_file = cut_file;
    elseif length(rat_num) >4
         db(i).cut_file = cut_file(2:end);
    end
    
    disp('')
%     for AB = {'A1' 'A2'}
%         AB = AB{:};
%         db(i).(AB).sessions = excel_sheet(i).(AB); % A1 in DB is a structure etc.
%         db(i).(AB).cut_file = [db(i).date db(i).(AB).sessions(1:2) ...
%             '_' num2str(db(i).tetrode) '.cut'];
%         if length(db(i).(AB).sessions) == 5
%             db(i).(AB).files{1} = [db(i).date db(i).(AB).sessions(1:2)];
%             db(i).(AB).files{2} = [db(i).date db(i).(AB).sessions(4:5)]; % assuming sessions is of the
%             % form 03,04
%         elseif length(db(i).(AB).sessions) == 2 % whitlock
%             db(i).(AB).files{1} = [db(i).date db(i).(AB).sessions(1:2)];
%         end
%     end

end % enumeration on excel sheet


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function excel_sheet = read_excel(excel_file)
%
% read excel sheet into structure
%

% read excel sheet

[num, txt, raw] = xlsread(excel_file);

fields = deblank(raw(2,:)); % excel file starts at line 2
raw = raw(3:end,:); % file really ends at 654 (rest are crossed out for some reason?)

% % get rid of empty rows % does not work because of a matlab bug
% 
% full_rows = find(~isnan([raw{:,1}])); % find non-empty rows (in which first cell is a NaN after reading)
% raw = raw(full_rows,:);

% get rid of empty columns

keep_columns = find(cellfun('isclass',fields,'char'));
fields = fields(keep_columns);
raw = raw(:,keep_columns);
       
% turn all numeric values into text
num_vals = cellfun('isclass',raw,'double');
num_inds = find(num_vals);

for ind = num_inds'
        if isnan(raw{ind})
            raw{ind} = '';
        else
            raw{ind} = num2str(raw{ind});
        end
end

% turn blank spaces into underscores for fields
fields= strrep(fields, ' ', '_');

excel_sheet = cell2struct(deblank(raw),fields,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function db = HP_read_pos_data(db,parms)
%
% loop on data-base records and read relevant pos-data,
% using different routines for the open-field and hairpin cases
%

%
% get open-field position data
%
   
    % read open field data
    
    dim_name= db.envt_type;

    if strfind(dim_name, '60')         
        parms.dim= 70;
    elseif strfind(dim_name, '100')
        parms.dim= 110;
    elseif strfind(dim_name, '120')
        parms.dim= 130;
    end
    
    db.pos_data = get_pos_data(db,parms);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos_data = get_pos_data(rec,parms,rotation)

if isempty(rec.cut_file)
    pos_data = struct();
    return;
end

if ~exist('rotation')
    rotation = 0;
end
%
% get position data for the open-field case
%

dim = parms.dim;
data_dir = parms.data_dir;

rat = rec.rat;

time_offset = [];

combPosX1 = []; combPosY1 = []; combPosT = []; combPosCT = []; 
combPosX2 = []; combPosY2 = [];
offset = 0;

    posFile = [data_dir '\r' rat '\' rec.date '\' rec.cut_file(1:end-6) '.pos'];
   
    disp(['reading from ' posFile]);
    
    x1 = []; y1 = []; t = []; ct = [];
    x2 = []; y2 = []; t = []; ct = [];
    % read the position data (assume single red led)

    [x1,y1,t] = GR_getPos(posFile,'red LED');
    [x2,y2,t] = GR_getPos(posFile,'green LED');

    % Copy of the postion timestamps before illegal times are removed

    ct = t;

   %%%%%%%%%%% REMOVED TO SEE WHAT HAPPENS
  %  [x1,y1,t] = GR_deletePathOutliers(x1,y1,t,dim);

    % Treshold for how far a rat can move (100cm/s), in one sample (sampFreq =
    % 50 Hz)
%     treshold = 100/50;
% 
%     [x1,y1,t] = GR_remBadTrack(x1,y1,t,treshold,dim);
%     [x1,y1] = GR_interporPos(x1,y1);
%     [x1,y1] = GR_posMeanFilter(x1,y1);

    % Remove NaN from the position samples if present
    
    x1(1) = NaN;   % remove the first sample to disconnect between combined files
    indNaN = isnan(x1);
    x1(indNaN) = [];    y1(indNaN) = [];    t(indNaN) = [];
    
    x2(1) = NaN;   % remove the first sample to disconnect between combined files

    x2(indNaN) = [];
    y2(indNaN) = [];

    indNaN=[];
    indNaN =isnan(x2);
    
    x2(indNaN) = [];
    y2(indNaN) = [];
    x1(indNaN) = [];
    y1(indNaN) = [];
    t(indNaN) = [];
    
    % rotate path so that it fits with hairpin data rotation
    
     rAngle = rotation*2*pi/360;
      [x1,y1] = rotatePath(x1,y1,rAngle);
     [x2,y2] = rotatePath(x2,y2,rAngle);
    
    % append new pos data to combined pos data
    
    time_offset = [time_offset; offset];
    combPosX1 = [combPosX1; x1];  % the first point is discarded to disconnect between the files
    combPosY1 = [combPosY1; y1];
    combPosT = [combPosT; t+offset];
    combPosCT = [combPosCT; ct + offset];

    % Time offset
    offset = combPosCT(end) + 0.02;
    

% set mapAxis

centre = GR_centreBox(combPosX1,combPosY1,dim);
combPosX1 = combPosX1 - centre(1);
combPosY1 = combPosY1 - centre(2);

obsSLength = max(max(combPosX1)-min(combPosX1),max(combPosY1)-min(combPosY1));
% We use bin size 3cm x 3cm
binWidth = 3;
%sprintf('%s%1.1f%s','Binwidth: ',binWidth,' cm.');
bins = ceil(obsSLength/binWidth);
sLength = binWidth * bins;
mapAxis = (-sLength/2+binWidth/2):binWidth:(sLength/2-binWidth/2);

mapAxis = GR_setMapAxis(combPosX1,combPosY1,mapAxis,binWidth,dim);

 combPosX2 = [combPosX2; x2];  % the first point is discarded to disconnect between the files
    combPosY2 = [combPosY2; y2];
%     combPosT = [combPosT; t+offset];
%     combPosCT = [combPosCT; ct + offset];

    % Time offset
%     offset = combPosCT(end) + 0.02;
    

% set mapAxis

centre = GR_centreBox(combPosX2,combPosY2,dim);
combPosX2 = combPosX2 - centre(1);
combPosY2 = combPosY2 - centre(2);

obsSLength = max(max(combPosX2)-min(combPosX2),max(combPosY2)-min(combPosY2));
% We use bin size 3cm x 3cm
binWidth = 3;
%sprintf('%s%1.1f%s','Binwidth: ',binWidth,' cm.');
bins = ceil(obsSLength/binWidth);
sLength = binWidth * bins;
mapAxis = (-sLength/2+binWidth/2):binWidth:(sLength/2-binWidth/2);

mapAxis = GR_setMapAxis(combPosX2,combPosY2,mapAxis,binWidth,dim);



pos_data.x1 = single(combPosX1);
pos_data.y1 = single(combPosY1);
pos_data.x2 = single(combPosX2);
pos_data.y2 = single(combPosY2);
pos_data.t = single(combPosT);
pos_data.ct = single(combPosCT);
pos_data.dim = parms.dim;
pos_data.mapAxis = single(mapAxis);
pos_data.time_offset = time_offset;
pos_data.trial_type = 'open_field';
pos_data.axis = single(mapAxis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated the position angles according to the angle tAngle [radians]
function [newX,newY] = rotatePath(x,y,tAngle)

newX = x * cos(tAngle) - y * sin(tAngle); % Tilted x-coord
newY = x * sin(tAngle) + y * cos(tAngle); % Tilted y-coord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapAxis = GR_setMapAxis(posx,posy,mapAxis,binWidth,dim)
%           Function for setting the axis for drawing

% Check for asymmetri in the path. If so correct acount for it in
% mapAxis
minX = min(posx);
maxX = max(posx);
minY = min(posy);
maxY = max(posy);
if minX < mapAxis(1)
    nXtra = ceil(abs(minX-mapAxis(1))/binWidth);
else
    nXtra = 0;
end
if maxX > mapAxis(end)
    pXtra = ceil(abs(maxX-mapAxis(end))/binWidth);
else
    pXtra = 0;
end
if nXtra
    for nn =1:nXtra
        tmp = mapAxis(1) - binWidth;
        mapAxis = [tmp; mapAxis'];
        mapAxis = mapAxis';
    end
end
if pXtra
    tmp = mapAxis(end) + binWidth;
    mapAxis = [mapAxis'; tmp];
    mapAxis = mapAxis';
end

if minY < mapAxis(1)
    nXtra = ceil(abs(minX-mapAxis(1))/binWidth);
else
    nXtra = 0;
end
if maxY > mapAxis(end)
    pXtra = ceil(abs(maxX-mapAxis(end))/binWidth);
else
    pXtra = 0;
end
if nXtra
    for nn =1:nXtra
        tmp = mapAxis(1) - binWidth;
        mapAxis = [tmp; mapAxis'];
        mapAxis = mapAxis';
    end
end
if pXtra
    tmp = mapAxis(end) + binWidth;
    mapAxis = [mapAxis'; tmp];
    mapAxis = mapAxis';
end

% Put on 3 extra cm on each side.
mapAxis = [mapAxis(1)-1.5;mapAxis'];
mapAxis = [mapAxis; mapAxis(end)+1.5];
mapAxis = mapAxis';
mapAxis = [mapAxis(1)-1.5;mapAxis'];
mapAxis = [mapAxis; mapAxis(end)+1.5];
mapAxis = mapAxis';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function arm_vec = GR_getArmVec(x,xBorders)
% %
% % for every point in x, determine to which arm it belongs
% %
% 
% arm_vec = zeros(size(x));
% 
% narms = length(xBorders)+1;
% 
% % leftmost arm (1st arm)
% 
% cArm = 1;
% arm_vec( x < xBorders(cArm) ) = cArm;
% 
% % arms 2-9
% 
% for cArm = 2:narms-1
%     arm_vec( x > xBorders(cArm-1) & x < xBorders(cArm) ) = cArm;
% end
% 
% % 10th arm (rightmost arm)
% 
% cArm = narms;
% arm_vec( x > xBorders(cArm-1) ) = cArm;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,y1,t] = GR_deletePathOutliers(x1,y1,t,dim)
% %
% % get rid of really strange values in path (use box dimension)
% %
% good_inds = x1 > nanmedian(x1)-dim & ...
%             x1 < nanmedian(x1)+dim & ...
%             y1 > nanmedian(y1)-dim & ...
%             y1 < nanmedian(y1)+dim;
%         
% x1 = x1(good_inds);
% y1 = y1(good_inds);
% t = t(good_inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the centre of the box
function centre = GR_centreBox(posx,posy,dim)

% posx = posx(abs(posx) < 2*dim);  % restrict the possible values of x and y according to the box dimension
% posy = posy(abs(posy) < 2*dim);

% Find border values for path and box
maxX = max(posx);
minX = min(posx);
maxY = max(posy);
minY = min(posy);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
centre = GR_findCentre(NE,NW,SW,SE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the centre of the box from the corner coordinates
function centre = GR_findCentre(NE,NW,SW,SE)

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posx,posy,post] = GR_getPos(posfile,colour)
%  
%   [posx,posy,post] = getpos(posfile,colour,arena)
%
%   Copyright (C) 2004 Sturla Molden
%   Centre for the Biology of Memory
%   NTNU
%   Modified by Raymond Skjerpeng 2004

[tracker,trackerparam] = HP_ImportVideoTracker(posfile);
if (trackerparam.num_colours ~= 4)
    error('getpos requires 4 colours in video tracker file.');
end    
post = zeros(trackerparam.num_pos_samples,1);


N = size(tracker(1).xcoord,2);
if N == 2 % A two point tracking has been done 
    temp = zeros(trackerparam.num_pos_samples,4);
else % Normal tracking
    temp = zeros(trackerparam.num_pos_samples,8);
end
for ii = 1:trackerparam.num_pos_samples
    post(ii) = tracker(ii).timestamp;
    temp(ii,:) = [tracker(ii).xcoord tracker(ii).ycoord];
end

[didFix,fixedPost] = GR_fixTimestamps(post);
if didFix
    post = fixedPost;
    disp('Continue to read data');
end

if N == 4
    switch colour
        case {'red LED'}
            posx = temp(:,1) + trackerparam.window_min_x;
            posy = temp(:,5) + trackerparam.window_min_y;
        case {'green LED'}
            posx = temp(:,2) + trackerparam.window_min_x;
            posy = temp(:,6) + trackerparam.window_min_y;
        case {'blue LED'}
            posx = temp(:,3) + trackerparam.window_min_x;
            posy = temp(:,7) + trackerparam.window_min_y;
        case {'black on white'}
            posx = temp(:,4) + trackerparam.window_min_x;
            posy = temp(:,8) + trackerparam.window_min_y;
        otherwise
            error(sprintf('unknown colour "%s"',colour));
    end    
end

if N == 2
    M = length(temp(:,1));
    n1 = sum(isnan(temp(:,1)));
    n2 = sum(isnan(temp(:,2)));
    if n1 > n2
        posx = temp(:,2);
        posy = temp(:,4);
    else
        posx = temp(:,1);
        posy = temp(:,3);
    end
end

numPos = length(posx);
numPost = length(post);
if numPos ~= numPost
    posx = posx(1:numPost);
    posy = posy(1:numPost);
end


% if N == 2
%     switch colour
%         case {'red LED'}
%             posx = temp(:,1) + trackerparam.window_min_x;
%             posy = temp(:,3) + trackerparam.window_min_y;
%         case {'green LED'}
%             posx = temp(:,2) + trackerparam.window_min_x;
%             posy = temp(:,4) + trackerparam.window_min_y;
%         
%         otherwise
%             error(sprintf('unknown colour "%s"',colour));
%     end    
% end

% index = find ( (posx==0) & (posy==511) ); % this is internal to CBM's equipment.
% posx(index) = NaN;                        % 
% posy(index) = NaN;                        % 


% if (nargin > 2)
%     [posx, posy] = HP_arenaConfig(posx,posy,arena);
% end
post = post - post(1);

% Remove NaN in position if these appear at the end of the file
% [posx,posy,post] = GR_removeNaN_1(posx,posy,post);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes NaN from the position samples if these occure at the end of the
% position file
function [posx,posy,post] = GR_removeNaN_1(posx,posy,post)

while 1
    if isnan(posx(end)) | isnan(posy(end))
        posx(end) = [];
        posy(end) = [];
        post(end) = [];
    else
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removes position "jumps", i.e position samples that imply that the rat is
% moving quicker than physical possible.
function [x,y,t] = GR_remBadTrack(x,y,t,treshold,dim)

% Indexes to position samples that are to be removed
remInd = [];
% remInd = find(abs(x)> 2* dim); % remove x which are much out of box
% remInd = [remInd; find(abs(y) > 2* dim)];

diffX = diff(x);
diffY = diff(y);
diffR = sqrt(diffX.^2 + diffY.^2);
ind = find(diffR > treshold);

if ind(end) == length(x)
    offset = 2;
else
    offset = 1;
end

for ii = 1:length(ind)-offset
    if ind(ii+1) == ind(ii)+1
        % A single sample position jump, tracker jumps out one sample and
        % then jumps back to path on the next sample. Remove bad sample.
        remInd = [remInd; ind(ii)+1];
        ii = ii+1;
        continue
    else
        % Not a single jump. 2 possibilities:
        % 1. Tracker jumps out, and stay out at the same place for several
        % samples and then jumps back.
        % 2. Tracker just has a small jump before path continues as normal,
        % unknown reason for this. In latter case the samples are left
        % untouched.
        idx = find(x(ind(ii)+1:ind(ii+1)+1)==x(ind(ii)+1));
        if length(idx) == length(x(ind(ii)+1:ind(ii+1)+1));
            remInd = [remInd; (ind(ii)+1:ind(ii+1)+1)'];
        end
    end
end
% Remove the samples
x(remInd) = [];
y(remInd) = [];
t(remInd) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If 1-5 positon samples are missing, this function will interpolate
% the missing samples. If more than 5 samples are missing in a row they are
% left as NaN.
function [x1,y1] = GR_interporPos(x1,y1)

% Find the indexes to the missing samples, for the red tracking diode
ind = isnan(x1);
ind = find(ind==1);
N = length(ind);
if N == 0
    return
end
% Set start and stop points for the loop
if ind(N) >= length(x1)-5
    endLoop = N-5;
else
    endLoop = N;
end
if ind(1) == 1
    startLoop = 2;
else
    startLoop = 1;
end

for ii = startLoop:endLoop
    if length(find(ind==ind(ii)+1)) == 0 & length(find(ind==ind(ii)-1))==0
        % Only one missing sample in a row
        x1(ind(ii)) = (x1(ind(ii)-1)+x1(ind(ii)+1))/2;
        y1(ind(ii)) = (y1(ind(ii)-1)+y1(ind(ii)+1))/2;
    else
        if length(find(ind==ind(ii)+2)) == 0 & length(find(ind==ind(ii)-1))==0
            % 2 missing samples in a row
            xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+2));
            yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+2));
            x1(ind(ii)) = x1(ind(ii)-1) + 1/3*xDist;
            y1(ind(ii)) = y1(ind(ii)-1) + 1/3*yDist;
            x1(ind(ii)+1) = x1(ind(ii)-1) + 2/3*xDist;
            y1(ind(ii)+1) = y1(ind(ii)-1) + 2/3*yDist;
            ii = ii+1;
        else
            if length(find(ind==ind(ii)+3)) == 0 & length(find(ind==ind(ii)-1))==0
                % 3 missing samples in a row
                xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+3));
                yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+3));
                x1(ind(ii)) = x1(ind(ii)-1) + 1/4*xDist;
                y1(ind(ii)) = y1(ind(ii)-1) + 1/4*yDist;
                x1(ind(ii)+1) = x1(ind(ii)-1) + 1/2*xDist;
                y1(ind(ii)+1) = y1(ind(ii)-1) + 1/2*yDist;
                x1(ind(ii)+2) = x1(ind(ii)-1) + 3/4*xDist;
                y1(ind(ii)+2) = y1(ind(ii)-1) + 3/4*yDist;
                ii = ii+2;
            else
                if length(find(ind==ind(ii)+4)) == 0 & length(find(ind==ind(ii)-1))==0
                    % 4 missing samples in a row
                    xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+4));
                    yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+4));
                    x1(ind(ii)) = x1(ind(ii)-1) + 1/5*xDist;
                    y1(ind(ii)) = y1(ind(ii)-1) + 1/5*yDist;
                    x1(ind(ii)+1) = x1(ind(ii)-1) + 2/5*xDist;
                    y1(ind(ii)+1) = y1(ind(ii)-1) + 2/5*yDist;
                    x1(ind(ii)+2) = x1(ind(ii)-1) + 3/5*xDist;
                    y1(ind(ii)+2) = y1(ind(ii)-1) + 3/5*yDist;
                    x1(ind(ii)+3) = x1(ind(ii)-1) + 4/5*xDist;
                    y1(ind(ii)+3) = y1(ind(ii)-1) + 4/5*yDist;
                    ii = ii+3;
                else
                    if length(find(ind==ind(ii)+5)) == 0 & length(find(ind==ind(ii)-1))==0
                        % 5 missing samples in a row
                        xDist = abs(x1(ind(ii)-1)-x1(ind(ii)+5));
                        yDist = abs(y1(ind(ii)-1)-y1(ind(ii)+5));
                        x1(ind(ii)) = x1(ind(ii)-1) + 1/6*xDist;
                        y1(ind(ii)) = y1(ind(ii)-1) + 1/6*yDist;
                        x1(ind(ii)+1) = x1(ind(ii)-1) + 2/6*xDist;
                        y1(ind(ii)+1) = y1(ind(ii)-1) + 2/6*yDist;
                        x1(ind(ii)+2) = x1(ind(ii)-1) + 3/6*xDist;
                        y1(ind(ii)+2) = y1(ind(ii)-1) + 3/6*yDist;
                        x1(ind(ii)+3) = x1(ind(ii)-1) + 4/6*xDist;
                        y1(ind(ii)+3) = y1(ind(ii)-1) + 4/6*yDist;
                        x1(ind(ii)+4) = x1(ind(ii)-1) + 5/6*xDist;
                        y1(ind(ii)+4) = y1(ind(ii)-1) + 5/6*yDist;
                        ii = ii+4;
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posx,posy] = GR_posMeanFilter(posx,posy)

% Smooth samples with a mean filter over 15 samples
for cc = 8:length(posx)-7
    posx(cc) = nanmean(posx(cc-7:cc+7));   
    posy(cc) = nanmean(posy(cc-7:cc+7));
end

%__________________________________________________________________________
%
%           Function for fixing the position timestamps
%__________________________________________________________________________

function [didFix,fixedPost] = GR_fixTimestamps(post)

% First time stamp in file
first = post(1);
% Number of timestamps
N = length(post);
uniqePost = unique(post);

if length(uniqePost)~=N
    disp('Position timestamps are corrected');
    didFix = 1;
    numZeros = 0;
    % Find the number of zeros at the end of the file
    while 1
        if post(end-numZeros)==0
            numZeros = numZeros + 1;
        else
            break;
        end
    end
    
    last = first + (N-1-numZeros) *0.02;
    fixedPost = first:0.02:last;
    fixedPost = fixedPost';
else
    didFix = 0;
    fixedPost = [];
end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function db = HP_read_spike_data(db,parms)
%
% for every cell and every cut_file read spike data
%
% (AB) stands for A1,B1,B2,A2

datafile = [parms.data_dir '\r' db.rat '\' db.date '\' db.cut_file(1:end-6) '.' db.tetrode];

combTS = GR_getspikes(datafile);

                   
    ts = GR_getspikes(datafile);
    
    combTS = [combTS; ts+db.pos_data.time_offset];
    

db.pos_data.ts = combTS;

% get cut file

cut_file = db.cut_file;
disp(['cut file: ' cut_file]);

db.pos_data.cut = GR_getcut([parms.data_dir '\r' db.rat '\' db.date '\' cut_file]);
cut_vals = db.pos_data.cut;
% get rid of the 0 values in cut file (we never look at cell 0)

non_zero_inds = find(cut_vals ~= 0);
db.pos_data.cut = ...
    db.pos_data.cut(non_zero_inds);
db.pos_data.ts = ...
    db.pos_data.ts(non_zero_inds);

cell_form= str2num(db.cell);    % string to number 
cut_inds = find(db.pos_data.cut == cell_form);
db.spike_data.cut = db.pos_data.cut(cut_inds);
db.spike_data.ts = db.pos_data.ts(cut_inds);

if isempty(db.spike_data.ts)
    db.spike_data.x = [];
    db.spike_data.y = [];
    if strcmp(db.pos_data.trial_type,'hairpin')
        db.spike_data.arm = [];
        db.spike_data.aim = [];
        db.spike_data.linear = [];
    end
    return;
end
        
% Calculate the spike positions

db.spike_data.x = GR_interp(db.spike_data.ts,db.pos_data.x1, ...
                                 db.pos_data.t,db.pos_data.ct);
db.spike_data.y = GR_interp(db.spike_data.ts,db.pos_data.y1, ...
                                 db.pos_data.t,db.pos_data.ct);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the position of the spikes if they are in the legal time zone
function spkx = GR_interp(ts,posx,post,cPost)

if isempty(posx) 
    spkx = []; 
    return
end

N = length(ts);
post_ind = interp1(ts,1:N,post,'linear','extrap');
cPost_ind = interp1(ts,1:N,cPost,'linear','extrap');
bad_inds = unique(round(setdiff(cPost_ind,post_ind)));
bad_inds = bad_inds(bad_inds > 0);  % bad inds must be positive
bad_inds = bad_inds(bad_inds <= N);
[tmp,unique_inds,rev_inds] = unique(post);

warning off MATLAB:interp1:NaNinY

spkx = interp1(post(unique_inds),posx(unique_inds),ts,'linear','extrap');

warning on MATLAB:interp1:NaNinY

spkx(bad_inds) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clust = GR_getcut(cutfile);
clust = [];
fid = fopen(cutfile, 'rt');
if fid == -1
    return
end

while ~feof(fid)
    string = fgetl(fid);
    if (length(string))
    if (string(1) == 'E') 
        break;
    end
    end
end
while ~feof(fid)
  string = fgetl(fid);
  if length(string)
     content = sscanf(string,'%u')';
     clust = [clust content];
  end
end
fclose(fid);
clust = clust';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ts,ch1,ch2,ch3,ch4] = GR_getspikes(filename)
%
%   [ts,ch1,ch2,ch3,ch4] = getspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

[spikes,spikeparam] = importspikes(filename);
ts = [spikes.timestamp1]';
nspk = spikeparam.num_spikes;
spikelen = spikeparam.samples_per_spike;
ch1 = reshape([spikes.waveform1],spikelen,nspk)';
ch2 = reshape([spikes.waveform2],spikelen,nspk)';
ch3 = reshape([spikes.waveform3],spikelen,nspk)';
ch4 = reshape([spikes.waveform4],spikelen,nspk)';


function [spikes,spikeparam] = importspikes(filename)
%
%   [spikes,spikeparam] = importspikes(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

fid = fopen(filename,'r','ieee-be');
if (fid < 0)
   error(sprintf('Could not open %s\n',filename)); 
end    

% read all bytes, look for 'data_start'
fseek(fid,0,-1);
sresult = 0;
[bytebuffer, bytecount] = fread(fid,inf,'uint8');
for ii = 10:length(bytebuffer)
    if strcmp( char(bytebuffer((ii-9):ii))', 'data_start' )
        sresult = 1;
        headeroffset = ii;
        break;
    end
end
if (~sresult)
    fclose(fid);
    error(sprintf('%s does not have a data_start marker', filename));
end

% count header lines
fseek(fid,0,-1);
headerlines = 0;
while(~feof(fid))
    txt = fgetl(fid);
    tmp = min(length(txt),10);
    if (length(txt))
        if (strcmp(txt(1:tmp),'data_start'))
            break;
        else
            headerlines = headerlines + 1;
        end
    else
        headerlines = headerlines + 1;
    end   
end    

% find timebase
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^timebase.*')))
        timebase = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 96 kHz');   
    timebase = 96000;    
end

% find duration
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^duration.*')))
        duration = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Duration not reported, defaulting to last time stamp');   
    duration = inf;    
end

% find number of spikes
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_spikes.*')))
        num_spikes = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of spikes not reported, using all that can be found');   
    num_spikes = inf;    
end

% find bytes per sample
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_sample.*')))
        bytes_per_sample = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per sample not reported, defaulting to 1');   
    bytes_per_sample = 1;    
end

% find bytes per timestamp
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_timestamp.*')))
        bytes_per_timestamp = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per timestamp not reported, defaulting to 4');   
    bytes_per_timestamp = 4;    
end

% find samples per spike
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^samples_per_spike.*')))
        samples_per_spike = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Samples per spike not reported, defaulting to 50');   
    samples_per_spike = 50;    
end

% check spike format
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^spike_format.*')))
        if (length(regexp(txt,'^spike_format t,ch1,t,ch2,t,ch3,t,ch4')))
            sresult = 1;
            break;
        else
           fclose(fid);
           error(sprintf('Unknown spike format, cannot read spikes from %s',filename));   
        end
    end
end    
if (~sresult)
    fclose(fid);
    error(sprintf('No spike format reported, cannot read spikes from %s.\nAre you sure this is a spike file?',filename));   
end

% close the file
fclose(fid);

% count the number of spikes in the file
spikelen = 4 * (bytes_per_sample * samples_per_spike + bytes_per_timestamp);
num_spikes_in_file = floor((bytecount - headeroffset)/spikelen);
if (isfinite(num_spikes))
    if (num_spikes_in_file > num_spikes)
        warning(sprintf('%d spikes reported in header, but %s seems to contain %d spikes.',num_spikes,filename,num_spikes_in_file));
    elseif (num_spikes_in_file < num_spikes)
        warning(sprintf('%d spikes reported in header, but %s can contain have %d spikes.',num_spikes,filename,num_spikes_in_file));
        num_spikes = num_spikes_in_file;    
    end
else
    num_spikes = num_spikes_in_file;
end
    
% allocate memory for return values

spikestruct = struct('timestamp1',0,'waveform1',zeros(samples_per_spike,1), ...
                     'timestamp2',0,'waveform2',zeros(samples_per_spike,1), ...
                     'timestamp3',0,'waveform3',zeros(samples_per_spike,1), ...
                     'timestamp4',0,'waveform4',zeros(samples_per_spike,1));

spikes = repmat(spikestruct,num_spikes,1);
                        
% out the spikes into the struct, one by one

big_endian_vector =  (256.^((bytes_per_timestamp-1):-1:0))';
little_endian_matrix = repmat(256.^(0:(bytes_per_sample-1))',1,samples_per_spike);

for ii = 1:num_spikes
   % sort the bytes for this spike
   spikeoffset = headeroffset + (ii-1)*spikelen;
   t1_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w1_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w1_bytes( w1_bytes > 127 ) = w1_bytes( w1_bytes > 127 ) - 256;
   w1_bytes = reshape(w1_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t2_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w2_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w2_bytes( w2_bytes > 127 ) = w2_bytes( w2_bytes > 127 ) - 256;
   w2_bytes = reshape(w2_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t3_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w3_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w3_bytes( w3_bytes > 127 ) = w3_bytes( w3_bytes > 127 ) - 256;
   w3_bytes = reshape(w3_bytes,bytes_per_sample,samples_per_spike);
   spikeoffset = spikeoffset + bytes_per_sample*samples_per_spike;
   t4_bytes = bytebuffer((spikeoffset+1):(spikeoffset+bytes_per_timestamp));
   spikeoffset = spikeoffset + bytes_per_timestamp;
   w4_bytes = bytebuffer((spikeoffset+1):(spikeoffset+(bytes_per_sample*samples_per_spike)));
   w4_bytes( w4_bytes > 127 ) = w4_bytes( w4_bytes > 127 ) - 256;
   w4_bytes = reshape(w4_bytes,bytes_per_sample,samples_per_spike);
   % interpret the bytes for this spike
   spikes(ii).timestamp1 = sum(t1_bytes .* big_endian_vector) / timebase; % time stamps are big endian
   spikes(ii).timestamp2 = sum(t2_bytes .* big_endian_vector) / timebase;
   spikes(ii).timestamp3 = sum(t3_bytes .* big_endian_vector) / timebase;
   spikes(ii).timestamp4 = sum(t4_bytes .* big_endian_vector) / timebase;
   spikes(ii).waveform1 =  sum(w1_bytes .* little_endian_matrix, 1); % signals are little-endian
   spikes(ii).waveform2 =  sum(w2_bytes .* little_endian_matrix, 1);
   spikes(ii).waveform3 =  sum(w3_bytes .* little_endian_matrix, 1);
   spikes(ii).waveform4 =  sum(w4_bytes .* little_endian_matrix, 1);
end
if (~isfinite(duration))
    duration = ceil(spikes(end).timestamp1);
end
spikeparam = struct('timebase',timebase,'bytes_per_sample',bytes_per_sample,'samples_per_spike',samples_per_spike, ...
                    'bytes_per_timestamp',bytes_per_timestamp,'duration',duration,'num_spikes',num_spikes);

