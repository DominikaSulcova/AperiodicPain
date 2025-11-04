%% AperiodicPain - individual data pre-processing
% author:   Dominika Sulcova, MSH Hamburg
% created:  2025   
% 
% This script runs the EEG pre-processing pipeline for a single subject 
% taking part in the AperiodicPain (AP) study.
% 
% Data: 
%   1) resting state EEG (RS) 
%           --> 1.5 mins eyes open + 1.5 mins eyes closed
%           --> recorded only from S007 (38 subjects in total)
%   2) laser-evoked potentials (LEP)
%           --> fixed Yap-laser parameters: 3ms, 5mm, 1,75J  
%           --> left and right hand stimulated, repeated in ABBA format
%           --> in total 4 EEG blocks of 30 stimuli
%   3) somatosensory-evoked potentials (SEP)
%           --> median nerve stimulation
%           --> 0.5ms, intensity above motor threshold (visible thumb twitch)
%           --> same structure as LEPs
%   4) visual-evoked potentials (VEP)
%           --> pattern-onset VEPs - 1° checkerboard stimuli, 200ms flash
%           --> 2 blocks of 30 stimuli
% 
% Processing steps:
% 
% Output:
%   1) AP_info      --> structure containing all information about
%                       the subject, experimental session, data pre-processing 
%                       including individual processing parameters...
%                   --> one row per subject
%   2) AP_data      --> structure containing final processed data  
%                   --> single-trial evoked potentials
%                   --> single-trial PSD of pre-stimulus EEG
%                   --> single-trial PSD of source activities
%                   --> one row per subject
%   3) AP_measures  --> structure containing extracted measures
%                   --> amplitudes and latencies of late LEP/SEP/VEP components (single-trial)
%                   --> single-trial estimates of pre-stimulus aperiodic exponent
%                   --> pain intensity ratings
%                   --> one row per subject  
% 
%% parameters
clear all; clc

% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
folder.input = uigetdir(pwd, 'Coose the input folder');             % raw data
folder.data = uigetdir(pwd, 'Choose the data folder');              % processed data
folder.output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, loutput file, exports 
cd(folder.data)

% output
study = 'AP';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);

% ask for subject_idx if necessary
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% load the info structure
if exist(output_file) == 2
    output_vars = who('-file', sprintf('%s', output_file));
    if ismember('AP_info', output_vars)
        load(output_file, 'AP_info')
    else
        AP_info = struct;
        save(output_file, 'AP_info','-append')
    end
else
    AP_info = struct;
    save(output_file, 'AP_info')
end
clear output_vars

%% 1) fill in subject & session information 
% ----- section input -----
params.laser.intensity = 1.75;
params.laser.pulse = 3;
params.laser.diameter = 5;
params.electrical.duration = 200;
params.electrical.target = 'median nerve';
params.visual.stimulus = 'checkerboard 1°';
params.visual.duration = 200;
% -------------------------
fprintf('section 1: fill in subject & session info\n\n')

% get subject & session info
prompt = {'date:', 'subject:', 'age (years):', 'male:', 'handedness score:', ...
    'height (cm):', 'weight (kg):', 'right arm (cm):', 'left arm (cm)', ...
    'starting modality:', 'starting side:', ...
    'ES intensity (mA) - right:', 'ES intensity (mA) - left:'};
dlgtitle = 'Subject & session';
dims = [1 50];
definput = {date, 'S200', '18', '1', '50', ...
    '170', '70', '60', '60', ...
    'laser/electrical', 'right/left', ...
    '8.0', '8.0'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% verify subject index
if subject_idx ~= str2num(session_info{2}(end-1:end))
    fprintf('WARNING: manually-encoded subject index (%d) does not match the index extracted from subject information (%d)!\n', subject_idx, str2num(session_info{2}(end-1:end)))
    while true
        subject_idx = input('Please enter the correct subject index: ');
        if isnumeric(subject_idx) && subject_idx > 0 && mod(subject_idx, 1) == 0
            break;
        else
            fprintf('⚠Please enter a valid, positive whole number.\n');
        end
    end
end

% fill in the metadata structure
AP_info(subject_idx).date = session_info{1};
AP_info(subject_idx).ID = session_info{2};
AP_info(subject_idx).age = str2num(session_info{3});
AP_info(subject_idx).male = str2num(session_info{4});
AP_info(subject_idx).handedness = str2num(session_info{5});
AP_info(subject_idx).body.height = str2num(session_info{6});
AP_info(subject_idx).body.weight = str2num(session_info{7});
AP_info(subject_idx).body.arm_right = str2num(session_info{8});
AP_info(subject_idx).body.arm_left = str2num(session_info{9});
AP_info(subject_idx).session.site = 'hand';
AP_info(subject_idx).session.starting_modality = session_info{10};
AP_info(subject_idx).session.starting_side = session_info{11};
AP_info(subject_idx).stimulation.laser.intensity = params.laser.intensity;
AP_info(subject_idx).stimulation.laser.pulse = params.laser.pulse;
AP_info(subject_idx).stimulation.laser.diameter = params.laser.diameter;
AP_info(subject_idx).stimulation.electrical.intensity.right = str2num(session_info{12});
AP_info(subject_idx).stimulation.electrical.intensity.left = str2num(session_info{13});
AP_info(subject_idx).stimulation.electrical.duration = params.electrical.duration;
AP_info(subject_idx).stimulation.electrical.target = params.electrical.target;
AP_info(subject_idx).stimulation.visual.stimulus = params.visual.stimulus;
AP_info(subject_idx).stimulation.visual.duration = params.visual.duration;
clear session_info

% get pre-LEP temperature
prompt = {'right hand, block 1:','right hand, block 2:', 'left hand, block 1:','left hand, block 2:'};
dlgtitle = 'Pre-LEP temperature';
dims = [1 35];
definput = {'35.5', '35.5', '35.5', '35.5'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% fill in the metadata structure
AP_info(subject_idx).temperature.right_block1 = session_info{1};
AP_info(subject_idx).temperature.right_block2 = session_info{2};
AP_info(subject_idx).temperature.left_block1 = session_info{3};
AP_info(subject_idx).temperature.left_block2 = session_info{4};
clear session_info 

% save to the output file
save(output_file, 'AP_info', '-append');
fprintf('section 1 finished.\n\n')

%% 2) import continuous EEG data, pre-process and save for letswave
% ----- section input -----
params.data = {'RS', 'LEP', 'SEP', 'VEP'};
params.folder = 'EEG';
params.data_n = [4, 4, 4, 2];
params.crop_margin = 5;
params.downsample = 5;
params.suffix = {'crop' 'ds' 'dc'};
% -------------------------
fprintf('section 2: import & pre-process continuous EEG data\n\n')

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% cycle though datasets and import in letswave format
fprintf('loading:\n')
for a = 1:length(params.data)
    % identify the appropriate files
    file2import = dir(sprintf('%s\\%s\\%s\\*%s*.vhdr', folder.input, AP_info(subject_idx).ID, params.folder, params.data{a}));

    % remove average files if necessary
    file2rmv = logical([]);
    for b = 1:length(file2import)
        if contains(file2import(b).name, 'avg') 
            file2rmv(b) = true;
        else
            file2rmv(b) = false;
        end
    end
    file2import(file2rmv) = [];
    
    % check that the number of files matches
    if size(file2import, 1) ~= params.data_n(a)
        fprintf('WARNING: incorrect number (%d) of %s recordings was found!\n', size(file2import, 1), params.data{a})
    end

    % identify datasets
    datanames = []; blocks = [];
    for b = 1:length(file2import)
        % prepare the dataname
        dataname = replace(file2import(b).name, ' ', '');
        underscores = strfind(dataname, '_');
        dataname = dataname(underscores(1) + 1 : end - 5);
        dataname = replace(dataname, '_', ' ');
        if contains(dataname, 'LEP')
            dataname = replace(dataname, 'LEP', 'laser');
        elseif contains(dataname, 'SEP')
            dataname = replace(dataname, 'SEP', 'electric');
        elseif contains(dataname, 'VEP')
            dataname = replace(dataname, 'VEP', 'visual');
        end
        datanames{b} = dataname;

        % identify the block
        block = regexp(file2import(b).name, 'b(\d+)', 'tokens');
        block = str2num(block{1}{1});
        blocks(b) = block;        
    end

    % identify repetition block
    if a > 1 && a < 4
        blocks_sorted = sort(blocks);
        blocks_rep.b1 = blocks_sorted([1, 2]);
        blocks_rep.b2 = blocks_sorted([3, 4]);
    elseif a == 4
        blocks_sorted = sort(blocks);
        blocks_rep.b1 = blocks_sorted(1);
        blocks_rep.b2 = blocks_sorted(2);
    end

    % cycle through datasets
    for c = 1:length(file2import)
        % identify the filename
        filename = sprintf('%s\\%s', file2import(c).folder, file2import(c).name);

        % deal with blocks
        breaks = strfind(datanames{c}, ' ');
        if contains(dataname, 'RS')
            datanames{c}(breaks(end) + 1 : end) = [];
            if blocks(c) < 4
                datanames{c}(end + 1 : end + 3) = 'pre';
            else
                datanames{c}(end + 1 : end + 4) = 'post';
            end
        else    
            datanames{c}(breaks(end) + 1 : end) = [];
            if ismember(blocks(c), blocks_rep.b1)                
                datanames{c}(end + 1 : end + 2) = 'b1';
            elseif ismember(blocks(c), blocks_rep.b2)
                datanames{c}(end + 1 : end + 2) = 'b2';
            end
        end

        % provide update
        fprintf('%s ...\n', datanames{c})
        
        % encode 
        AP_info(subject_idx).dataset(blocks(c) - 1).block = blocks(c);
        AP_info(subject_idx).dataset(blocks(c) - 1).filename = file2import(c).name(1:end-5);
        AP_info(subject_idx).dataset(blocks(c) - 1).data = datanames{c}(6:end);
    
        % import the dataset
        dataset.raw(blocks(c) - 1).name = datanames{c}(6:end);
        [dataset.raw(blocks(c) - 1).header, dataset.raw(blocks(c) - 1).data, ~] = RLW_import_VHDR(filename);
    
        % rename in the header
        dataset.raw(blocks(c) - 1).header.name = datanames{c};
    end
end
fprintf('done.\n%d datasets imported.\n\n', length(dataset.raw))

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% pre-process continuous data and save for letswave
fprintf('pre-processing:\n')
for d = 1:length(dataset.raw)
    % provide update
    fprintf('%d. dataset: %s\n', d, AP_info(subject_idx).dataset(d).data)

    % select data
    lwdata.header = dataset.raw(d).header;
    lwdata.data = dataset.raw(d).data; 

    % assign electrode coordinates
    fprintf('assigning electrode coordinates...\n')
    option = struct('filepath', sprintf('%s\\letswave 7\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', folder.toolbox), ...
        'suffix', '', 'is_save', 0);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(1).process = 'electrode coordinates assigned';
        AP_info(subject_idx).preprocessing(1).params.layout = 'standard 10-20-cap81';
        AP_info(subject_idx).preprocessing(1).suffix = [];
        AP_info(subject_idx).preprocessing(1).date = sprintf('%s', date);
    end

    % select events to keep
    fprintf('checking events...\n')
    if ~isempty(lwdata.header.events)
        events = unique({lwdata.header.events.code}); 
    else
        error('ERROR: no events/triggers found in the recording!')
    end
    [selection, ok] = listdlg( ...
        'PromptString', 'Select the events to KEEP:', ...
        'SelectionMode', 'multiple', ...
        'ListString', events, ...
        'InitialValue', 1:numel(events));      
    if ok
        events = events(selection);
        fprintf('keeping following event codes: ');
        for e = 1:length(events)
            fprintf('%s ', events{e})            
        end
        fprintf('\n')
    else
        fprintf('you cancelled the selection!\n--> keeping all available event codes: ');
        for e = 1:length(events)
            fprintf('%s ', events{e})            
        end
        fprintf('\n')
    end

    % remove all other events, count and encode
    event_idx = logical([]);
    event_cnt = zeros(1, length(events));
    for e = 1:length(lwdata.header.events)
        if ismember(lwdata.header.events(e).code, events)
            event_idx(e) = true;
            for f = 1:length(events)
                if strcmp(lwdata.header.events(e).code, events{f})
                    event_cnt(f) = event_cnt(f) + 1;
                end
            end
        else
            event_idx(e) = false;
        end
    end
    lwdata.header.events = lwdata.header.events(event_idx);
    if d == 1
        AP_info(subject_idx).preprocessing(2).process = 'EEG elvents identified';
        AP_info(subject_idx).preprocessing(2).suffix = [];
        AP_info(subject_idx).preprocessing(2).date = sprintf('%s', date);
    end
    AP_info(subject_idx).preprocessing(2).params(d).data = AP_info(subject_idx).dataset(d).data;
    AP_info(subject_idx).preprocessing(2).params(d).eventcodes = events;
    AP_info(subject_idx).preprocessing(2).params(d).eventcounts = event_cnt;
    for e = 1:length(events)
        fprintf('--> %d event(s) labeled ''%s'' identified.\n', event_cnt(e), events{e})
    end

    % crop around events
    fprintf('cropping around events...\n')
    if lwdata.header.events(1).latency - params.crop_margin > 0
        params.crop(1) = lwdata.header.events(1).latency - params.crop_margin;
    else
        params.crop(1) = 0;
    end
    if length(events) == 1 & event_cnt(1) == 1
        params.crop(2) = floor(lwdata.header.datasize(6) * lwdata.header.xstep); 
    else
        params.crop(2) = lwdata.header.events(end).latency + params.crop_margin;
    end    
    option = struct('xcrop_chk', 1, 'xstart', params.crop(1), 'xend', params.crop(2), ...
        'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(3).process = 'continuous data cropped';
        AP_info(subject_idx).preprocessing(3).suffix = params.suffix{1};
        AP_info(subject_idx).preprocessing(3).date = sprintf('%s', date);
    end
    AP_info(subject_idx).preprocessing(3).params(d).data = AP_info(subject_idx).dataset(d).data;
    AP_info(subject_idx).preprocessing(3).params(d).start = params.crop(1);
    AP_info(subject_idx).preprocessing(3).params(d).end = params.crop(2);

    % downsample 
    fprintf('downsampling...\n')
    option = struct('x_dsratio', params.downsample, 'suffix', params.suffix{2}, 'is_save', 0);
    lwdata = FLW_downsample.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(4).process = sprintf('downsampled');
        AP_info(subject_idx).preprocessing(4).params.ratio = params.downsample;
        AP_info(subject_idx).preprocessing(4).params.fs_original = 1/lwdata.header.xstep * params.downsample;
        AP_info(subject_idx).preprocessing(4).params.fs_final = 1/lwdata.header.xstep;
        AP_info(subject_idx).preprocessing(4).suffix = params.suffix{2};
        AP_info(subject_idx).preprocessing(4).date = sprintf('%s', date);
    end

    % remove DC + linear detrend, save for letswave
    fprintf('removing DC and applying linear detrend...\nsaving for letswave...\n')
    option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(5).process = sprintf('DC + linear detrend on continuous data');
        AP_info(subject_idx).preprocessing(5).suffix = params.suffix{3};
        AP_info(subject_idx).preprocessing(5).date = sprintf('%s', date);
    end

    % update dataset
    dataset.raw(d).header = lwdata.header;
    dataset.raw(d).data = lwdata.data; 
    if d == length(dataset.raw)
        fprintf('done.\n\n')
    else
        fprintf('\n')
    end
end

% save info structure and move on
save(output_file, 'AP_info', '-append');
clear a b c d e f file2import file2rmv filename dataname datanames block blocks blocks_rep blocks_sorted breaks underscores ...
    lwdata option events selection ok event_idx event_cnt
fprintf('section 2 finished.\n\n')

%% 3) pre-process stimulation data
% ----- section input -----
params.prefix = 'dc ds crop';
params.data = {'laser', 'electric', 'visual'};
params.triggers = {'L  1', 'E  1', 'V  1'};
params.bandpass = [0.1 80];
params.notch = 50;
params.epoch.long = [-1.5, 1.5];
params.interpolate_n = 4;
params.suffix = {'bandpass' 'notch' 'reref' 'ep' 'dc' 'ar'};
% ------------------------- 
fprintf('section 3: pre-process stimulation data\n\n')

% re-load dataset if needed, remove RS data
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, AP_info(subject_idx).ID));
    dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'raw');
    clear data2load
end
fprintf('removing RS data...\n')
data_idx = logical([]);
for a = 1:length(dataset.raw)
    if contains(dataset.raw(a).name, 'RS') 
        data_idx(a) = true;
    else
        data_idx(a) = false;
    end
end
dataset.raw(data_idx) = [];
fprintf('done.\n\n')

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% apply frequency filters
%!make sure the notch filter code is fixed
for b = 1:length(dataset.raw)
    % provide update
    fprintf('%d. dataset: %s\n', b, dataset.raw(b).name)

    % select data
    lwdata.header = dataset.raw(b).header;
    lwdata.data = dataset.raw(b).data; 

    % bandpass filter
    fprintf('applying Butterworth bandpass filter...\n')
    option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2),'low_cutoff', params.bandpass(1),...
        'filter_order', 4, 'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
    if b == 1
        AP_info(subject_idx).preprocessing(6).process = 'bandpass filtered';
        AP_info(subject_idx).preprocessing(6).params.filter = 'Butterworth';
        AP_info(subject_idx).preprocessing(6).params.order = 4;
        AP_info(subject_idx).preprocessing(6).params.limits = [params.bandpass(1), params.bandpass(2)];
        AP_info(subject_idx).preprocessing(6).suffix = params.suffix{1};
        AP_info(subject_idx).preprocessing(6).date = sprintf('%s', date);
    end

    % notch filter
    fprintf('applying FFT notch filter...\n')
    option = struct('filter_type', 'notch', 'notch_fre', params.notch, 'notch_width', 2, 'slope_width', 2,...
        'harmonic_num', 2,'suffix', params.suffix{2},'is_save', 1);
    lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);
    if b == 1
        AP_info(subject_idx).preprocessing(7).process = 'notch filtered';
        AP_info(subject_idx).preprocessing(7).params.filter = 'FFT';
        AP_info(subject_idx).preprocessing(7).params.width = 2;
        AP_info(subject_idx).preprocessing(7).params.slope = 2;
        AP_info(subject_idx).preprocessing(7).suffix = params.suffix{2};
        AP_info(subject_idx).preprocessing(7).date = sprintf('%s', date);
    end

    % update dataset
    dataset.preprocessed(b).name = dataset.raw(b).name;
    dataset.preprocessed(b).header = lwdata.header;
    dataset.preprocessed(b).data = lwdata.data; 
    if b == length(dataset.raw)
        fprintf('done.\n\n')
    else
        fprintf('\n')
    end
end

% visual check + adjust triggers if necessary
fprintf('please verify visually the quality of data and adjust triggers if necessary.\n')
addpath(genpath([folder.toolbox '\letswave 7']));
letswave

% ask for channels to interpolate
while true
    answer = input('do you want to interpolate any channel? (y/n): ', 's');
    if any(strcmpi(answer, {'y','n'}))
        break;
    else
        disp('please enter "y" or "n"!');
    end
end

% update dataset
fprintf('updating dataset...\n')
params.prefix = regexp(dataset.preprocessed(1).header.name, ['(.*)\s+' AP_info(subject_idx).ID], 'tokens');
params.prefix = params.prefix{1}{1};
dataset_old = dataset;
data2load = dir(sprintf('%s*%s*', params.prefix, AP_info(subject_idx).ID));
dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'preprocessed');
dataset_new = dataset;
[dataset_old.preprocessed] = dataset_new.preprocessed;
dataset = dataset_old;
clear data2load dataset_old dataset_new

% interpolate bad channels if necessary
params.channels = {dataset.preprocessed(1).header.chanlocs.labels};
params.chanlocs = dataset.preprocessed(1).header.chanlocs;
switch answer
    case 'n'
        % provide update
        fprintf('no channels were interpolated.\n');

        % encode 
        AP_info(subject_idx).preprocessing(8).process = sprintf('no channels interpolated');
        AP_info(subject_idx).preprocessing(8).date = sprintf('%s', date);
    case 'y'
        % select channels to interpolate
        [selection, ok] = listdlg( ...
        'PromptString', 'Select channels:', ...
            'SelectionMode', 'multiple', ...
            'ListString', params.channels, ...
            'InitialValue', 1:numel(params.channels));      
        if ok
            % provide update
            channels2interp = params.channels(selection);
            fprintf('interpolating following channel(s): ');
            for c = 1:length(channels2interp)
                fprintf('%s ', channels2interp{c})            
            end
            fprintf('\n')

            % interpolate identified channels in all datsets
            for c = 1:length(channels2interp)    
                % calculate distances with other electrodes
                chan_n = selection(c);
                chan_dist = -ones(length(params.channels), 1);
                for d = setdiff(1:length(params.channels), chan_n)
                    if params.chanlocs(d).topo_enabled == 1
                        chan_dist(d) = sqrt((params.chanlocs(d).X - params.chanlocs(chan_n).X)^2 + ...
                            (params.chanlocs(d).Y - params.chanlocs(chan_n).Y)^2 + ...
                            (params.chanlocs(d).Z - params.chanlocs(chan_n).Z)^2);
                    end
                end
                chan_dist((chan_dist==-1)) = max(chan_dist);
                [~, chan_idx] = sort(chan_dist);
    
                % identify neighbouring channels
                chan_idx = chan_idx(1:params.interpolate_n);
                channels2use = params.channels(chan_idx);
    
                % cycle through all datasets
                for d = 1:length(dataset.preprocessed)
                    % select data
                    lwdata.header = dataset.preprocessed(d).header;
                    lwdata.data = dataset.preprocessed(d).data;
        
                    % interpolate using the neighboring electrodes
                    option = struct('channel_to_interpolate', channels2interp{c}, 'channels_for_interpolation_list', {channels2use}, ...
                        'suffix', '', 'is_save', 0);
                    lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);
        
                    % update dataset
                    dataset.preprocessed(d).header = lwdata.header;
                    dataset.preprocessed(d).data = lwdata.data;  
                end
                
                % encode
                if c == 1
                    AP_info(subject_idx).preprocessing(8).process = sprintf('bad channel(s) interpolated');
                    AP_info(subject_idx).preprocessing(8).date = sprintf('%s', date);
                end
                AP_info(subject_idx).preprocessing(8).params.bad{c} = channels2interp{c};
                AP_info(subject_idx).preprocessing(8).params.channels_used{c} = strjoin(channels2use, ' ');  
            end
        else
            % encode 
            fprintf('you cancelled the selection!\n--> no channels will be interpolated.\n');
            AP_info(subject_idx).preprocessing(8).process = sprintf('no channels interpolated');
            AP_info(subject_idx).preprocessing(8).date = sprintf('%s', date);
        end
end
fprintf('\n')

% re-reference and segment into long epochs
for d = 1:length(dataset.preprocessed)
    % provide update
    fprintf('%d. dataset: %s\n', d, dataset.preprocessed(d).name)

    % select data
    lwdata.header = dataset.preprocessed(d).header;
    lwdata.data = dataset.preprocessed(d).data; 

    % re-reference to common average
    fprintf('re-referencing to common average...\n')
    option = struct('reference_list', {params.channels}, 'apply_list', {params.channels}, ...
        'suffix', params.suffix{3}, 'is_save', 0);
    lwdata = FLW_rereference.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(9).process = sprintf('re-referenced to common average');
        AP_info(subject_idx).preprocessing(9).suffix = params.suffix{3};
        AP_info(subject_idx).preprocessing(9).date = sprintf('%s', date);
    end

    % select the trigger
    for e = 1:length(params.data)
        if contains(lwdata.header.name, params.data{e})
            trigger = params.triggers{e};
        end
    end

    % epoch
    fprintf('segmenting into long epochs...\n')
    option = struct('event_labels', {trigger}, 'x_start', params.epoch.long(1), 'x_end', params.epoch.long(2), ...
        'x_duration', params.epoch.long(2) - params.epoch.long(1), 'suffix', params.suffix{4}, 'is_save', 0);
    lwdata = FLW_segmentation.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(10).process = sprintf('segmented to long epochs');
        AP_info(subject_idx).preprocessing(10).params.trigger = trigger;
        AP_info(subject_idx).preprocessing(10).params.limits = params.epoch.long;
        AP_info(subject_idx).preprocessing(10).suffix = params.suffix{4};
        AP_info(subject_idx).preprocessing(10).date = sprintf('%s', date);
    end

    % remove DC + linear detrend, save for letswave
    fprintf('removing DC and applying linear detrend...\nsaving for letswave...\n')
    option = struct('linear_detrend', 1, 'suffix', params.suffix{5}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        AP_info(subject_idx).preprocessing(11).process = sprintf('DC + linear detrend on epoched data');
        AP_info(subject_idx).preprocessing(11).suffix = params.suffix{5};
        AP_info(subject_idx).preprocessing(11).date = sprintf('%s', date);
    end

    % update dataset
    dataset.preprocessed(d).header = lwdata.header;
    dataset.preprocessed(d).data = lwdata.data; 
    if d == length(dataset.preprocessed)
        fprintf('done.\n\n')
    else
        fprintf('\n')
    end
end

% remove bad epochs
fprintf('please check data visually and discard bad epochs.\n')
addpath(genpath([folder.toolbox '\letswave 6']));
letswave
for f = 1:length(dataset.preprocessed)
    filenames{f} = sprintf('%s %s.mat', params.suffix{6}, dataset.preprocessed(f).header.name);
end
wait4files(filenames);
fprintf('\n')

% load dataset with bad epochs removed
fprintf('updating dataset...\n')
params.prefix = regexp(dataset.preprocessed(1).header.name, ['(.*)\s+' AP_info(subject_idx).ID], 'tokens');
params.prefix = params.prefix{1}{1};
dataset_old = dataset;
data2load = dir(sprintf('%s*%s*', params.suffix{6}, AP_info(subject_idx).ID));
dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'filtered');
dataset_new = dataset;
[dataset_old.filtered] = dataset_new.filtered;
dataset = dataset_old;
clear data2load dataset_old dataset_new

% encode bad epochs
fprintf('encoding bad epochs...\n')
AP_info(subject_idx).preprocessing(12).process = sprintf('bad epochs discarded');
AP_info(subject_idx).preprocessing(12).params.GUI = 'letswave';
AP_info(subject_idx).preprocessing(12).suffix = params.suffix{6};
AP_info(subject_idx).preprocessing(12).date = sprintf('%s', date);
for a = 1:length(dataset.filtered)
    % subset header
    header = dataset.filtered(a).header;

    % extract discarded expochs
    if ~isempty(header.history(end).configuration)
        if ~isempty(header.history(end).configuration.parameters.rejected_epochs)
            discarded = header.history(end).configuration.parameters.rejected_epochs;
        else
            discarded = [];
        end
    end

    % encode 
    AP_info(subject_idx).preprocessing(12).params.discarded(a).dataset = dataset.filtered(a).name;
    AP_info(subject_idx).preprocessing(12).params.discarded(a).trials = discarded;
    AP_info(subject_idx).preprocessing(12).params.kept(a).dataset = dataset.filtered(a).name;
    AP_info(subject_idx).preprocessing(12).params.kept(a).trials = header.datasize(1);
end
fprintf('done.\n\n')

% save info structure and move on
save(output_file, 'AP_info', '-append');
clear a b c d e f answer chan_dist chan_idx chan_n channels2interp channels2use data_idx ...
    discarded filenames header lwdata ok option selection trigger
fprintf('section 3 finished.\n\n')

%% 4) ICA - remove ocular & muscle artifacts
% ----- section input -----
params.prefix = 'ar dc ep reref notch bandpass dc ds crop';
params.ICA_comp = 30;
params.suffix = {'ica_all' 'icfilt'};
% ------------------------- 
fprintf('section 4: remove ocular & muscle artifacts using ICA\n\n')

% re-load dataset if needed
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, AP_info(subject_idx).ID));
    dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'filtered');
    clear data2load
    fprintf('done.\n\n')
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% compute ICA matrix 
fprintf('computing ICA matrix...\n')
lwdataset = dataset.filtered;  
addpath(genpath([folder.toolbox '\letswave 7']));
option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', params.ICA_comp, 'suffix', params.suffix{1}, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
fprintf('done.\n\n')

% extract ICA parameters
fprintf('extracting ICA parameters...\n')
matrix.mix = lwdataset(1).header.history(end).option.mix_matrix;
matrix.unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
params.chanlocs = lwdataset(1).header.chanlocs;
for i = 1:size(matrix.mix, 2)
    params.ICA_labels{i} = ['IC',num2str(i)];
end
params.ICA_SR = 1/lwdataset(1).header.xstep;

% encode
AP_info(subject_idx).preprocessing(13).process = 'ICA matrix computed';
AP_info(subject_idx).preprocessing(13).params.method = 'runica';
AP_info(subject_idx).preprocessing(13).params.components = params.ICA_comp;
AP_info(subject_idx).preprocessing(13).params.chanlocs = params.chanlocs;
AP_info(subject_idx).preprocessing(13).params.labels = params.ICA_labels;
AP_info(subject_idx).preprocessing(13).params.SR = params.ICA_SR;
AP_info(subject_idx).preprocessing(13).params.matrix = matrix;
AP_info(subject_idx).preprocessing(13).suffix = params.suffix{1};
AP_info(subject_idx).preprocessing(13).date = sprintf('%s', date);

% update dataset and adjust for letswave 6
fprintf('updating dataset...\n')
for a = 1:length(lwdataset)
    % update filtred dataset
    dataset.filtered(a).header = lwdataset(a).header;
    dataset.filtered(a).data = lwdataset(a).data;

    % adjust header for letswave 6
    dataset.filtered(a).header.history(end).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
    dataset.filtered(a).header.history(end).configuration.parameters = dataset.filtered(a).header.history(end).option;  
    [dataset.filtered(a).header.history(end).configuration.parameters.ICA_um] = dataset.filtered(a).header.history(end).configuration.parameters.unmix_matrix; 
    [dataset.filtered(a).header.history(end).configuration.parameters.ICA_mm] = dataset.filtered(a).header.history(end).configuration.parameters.mix_matrix; 
    dataset.filtered(a).header.history(end).configuration.parameters = rmfield(dataset.filtered(a).header.history(end).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
    header = dataset.filtered(a).header;
    save(sprintf('%s.lw6', dataset.filtered(a).header.name), 'header');
end

% unmix data
fprintf('unmixing data...\n')
for b = 1:length(dataset.filtered)
    for e = 1:size(dataset.filtered(b).data, 1)
        dataset.unmixed(b).name = dataset.filtered(b).name;
        dataset.unmixed(b).header = dataset.filtered(b).header;
        dataset.unmixed(b).data(e, :, 1, 1, 1, :) = matrix.unmix * squeeze(dataset.filtered(b).data(e, :, 1, 1, 1, :));        
    end
end

% calculate PSD across all datasets and trials 
fprintf('estimating spectral content of ICs...\n')
psd = [];
for c = 1:params.ICA_comp
    for d = 1:length(dataset.unmixed)
        for e = 1:size(dataset.unmixed(d).data, 1)
            [psd(c, d, e, :), freq] = pwelch(squeeze(dataset.unmixed(d).data(e, c, 1, 1, 1, :)), ...
                [], [], [], AP_info(subject_idx).preprocessing(13).params.SR);  
        end
    end
end
AP_info(subject_idx).preprocessing(13).params.PSD = squeeze(mean(psd, [2, 3]));

% plot IC spectral content (in 2 figures to keep the size readable)
addpath(genpath([folder.toolbox '\letswave 6']));
for a = 1:2
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    components2plot = (a-1)*(params.ICA_comp/2) + [1:ceil(params.ICA_comp/2)];
    for f = 1:ceil(params.ICA_comp/2)
        % plot the topography
        matrix = AP_info(subject_idx).preprocessing(13).params.matrix.mix;
        subplot(ceil(length(components2plot)/3), 6, (f-1)*2 + 1);
        topoplot(double(matrix(:, components2plot(f))'), params.chanlocs, 'maplimits', [-4 4], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color',[1 1 1]);
        title(params.ICA_labels{components2plot(f)})
    
        % plot the psd
        subplot(ceil(length(components2plot)/3), 6, (f-1)*2 + 2);
        plot(freq(1:max(find(freq <= 45))), log10(AP_info(subject_idx).preprocessing(13).params.PSD(components2plot(f), 1:max(find(freq <= 45)))));
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
    end
    saveas(gcf, sprintf('%s\\figures\\ICA_%s_%d.png', folder.output, AP_info(subject_idx).ID, a));
end

% remove artifactual ICs
fprintf('please perform ICA manually in letswave 6.\n')
letswave
for f = 1:length(dataset.preprocessed)
    filenames{f} = sprintf('%s %s.mat', params.suffix{2}, dataset.filtered(f).header.name);
end
wait4files(filenames);
fprintf('\n')

% load dataset with artifactual ICs removed
fprintf('updating dataset...\n')
params.prefix = regexp(dataset.filtered(1).header.name, ['(.*)\s+' AP_info(subject_idx).ID], 'tokens');
params.prefix = params.prefix{1}{1};
dataset_old = dataset;
data2load = dir(sprintf('%s*%s*', params.suffix{2}, AP_info(subject_idx).ID));
dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'filtered');
dataset_new = dataset;
[dataset_old.filtered] = dataset_new.filtered;
dataset = dataset_old;
clear data2load dataset_old dataset_new

% ask for the input
prompt = {'blinks:', 'horizontal:', 'muscles:', 'electrodes:'};
dlgtitle = 'ICs removed';  
dims = [1 60];
definput = {'', '', '', ''};
input = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput 

% encode ICA 
AP_info(subject_idx).preprocessing(14).process = 'artifactual ICs discarded';
AP_info(subject_idx).preprocessing(14).suffix = params.suffix{2};
AP_info(subject_idx).preprocessing(14).date = sprintf('%s', date);
AP_info(subject_idx).preprocessing(14).params.kept = params.ICA_comp - length([str2num(input{1}), str2num(input{2}), str2num(input{3}), str2num(input{4})]);
AP_info(subject_idx).preprocessing(14).params.removed.blinks = str2num(input{1});
AP_info(subject_idx).preprocessing(14).params.removed.horizontal = str2num(input{2});
AP_info(subject_idx).preprocessing(14).params.removed.muscles = str2num(input{3});
AP_info(subject_idx).preprocessing(14).params.removed.electrodes = str2num(input{4});

% save info structure and move on
save(output_file, 'AP_info', '-append');
clear a b c d e f i components2plot filenames freq input lwdataset matrix option psd 
fprintf('section 4 finished.\n\n')

%% 5) split stimulation data into conditions
% ----- section input -----
params.prefix = 'icfilt ica_all ar dc ep reref notch bandpass dc ds crop';
params.data = {'laser', 'electric', 'visual'};
params.triggers = {'L  1', 'E  1', 'V  1'};
params.ready = 'G  1';
params.epoch.ERP = [-0.3, 1];
params.epoch.prestim = [-1, 0];
params.suffix = {'ERP' 'prestim'};
% ------------------------- 
fprintf('section 5: split stimulation data into conditions\n\n')

% re-load dataset if needed
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, AP_info(subject_idx).ID));
    dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'filtered');
    clear data2load
    fprintf('done.\n\n')
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% cycle through data
for a = 1:length(params.data)
    % select all the appropriate datasets
    idx_data = logical([]);
    for b = 1:length(dataset.filtered)
        if contains(dataset.filtered(b).name, params.data{a})
            idx_data(b) = true;
        else
            idx_data(b) = false;
        end
    end
    data2process = dataset.filtered(idx_data);

    % select the trigger
    trigger = params.triggers{a};

    % cycle through datasets
    for b = 1:length(data2process)
        % provide update
        fprintf('dataset: %s\n', data2process(b).name)

        % cycle through epoch types
        for c = 1:length(params.suffix)   
            % subset the data
            lwdata.header = data2process(b).header;
            lwdata.data = data2process(b).data;

            % determine epoch limits
            epoch = params.epoch.(params.suffix{c});
            
            % epoch
            fprintf('saving the %s epoch...\n', params.suffix{c})
            option = struct('event_labels', {trigger}, 'x_start', epoch(1), 'x_end', epoch(2), ...
                'x_duration', epoch(2) - epoch(1), 'suffix', params.suffix{c}, 'is_save', 1);
            lwdata = FLW_segmentation.get_lwdata(lwdata, option);
            if b == 1 
                AP_info(subject_idx).preprocessing(15).process = sprintf('segmented to short epochs');
                AP_info(subject_idx).preprocessing(15).params.trigger{a} = trigger;
                AP_info(subject_idx).preprocessing(15).params.epochs(c).data = params.suffix{c};
                AP_info(subject_idx).preprocessing(15).params.epochs(c).limits = epoch;
                AP_info(subject_idx).preprocessing(15).suffix{c} = params.suffix{c};
                AP_info(subject_idx).preprocessing(15).date = sprintf('%s', date);
            end

            % save to the dataset
            dataset.(params.suffix{c}).(params.data{a})(b).name = data2process(b).name;
            dataset.(params.suffix{c}).(params.data{a})(b).header = lwdata.header;
            dataset.(params.suffix{c}).(params.data{a})(b).data = lwdata.data;
        end
        fprintf('\n')
    end
end
fprintf('done.\n\n')

% save info structure and move on
save(output_file, 'AP_info', '-append');
clear a b c d data2process epoch header idx_data option lwdata trigger
fprintf('section 5 finished.\n\n')

%% to add:
% single-trial ERP pre-processing:
%   - LEP pre-processing per component --> N1, N2, P2
%   - SEP N2 and P2 estimation
%   - VEP N2 and P2 estimation
% baseline correcton after final ICA for all conditions!
% single-subject average plotting

%% in the next script AperiodicPain_analysis:
% extraction of ERP amplitudes and latencies
% single-trial prestim pre-processing at sensor level ==> aperiodic
% exponent extraction
% single-trial source activity estimation ==> processing at source level

%% 6) LEPs: pre-process for single-trial analysis
% ----- section input -----
params.prefix = 'ERP icfilt ica_all ar dc ep reref notch bandpass dc ds crop';
params.suffix = {'reref_AFz' 'ica_N1' 'icfilt' 'bl'};
params.ref = 'AFz';
params.baseline = [-0.25 0];
% -------------------------
fprintf('section 6: LEPs - pre-process for single-trial analysis\n\n')

% re-load dataset if needed
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, AP_info(subject_idx).ID));
    dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'ERP');
    clear data2load
    fprintf('done.\n\n')
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% subset LEP datasets
for d = 1:length(dataset.ERP)
    if contains(dataset.ERP(d).name, 'laser')
        data_idx(d) = true;
    else
        data_idx(d) = false;
    end
end
subset = dataset.ERP(data_idx);
lwdataset = subset;

% re-reference to chosen frontal central electrode
fprintf('re-referencing to %s...\n', params.ref)
params.chanlocs = lwdataset(1).header.chanlocs;
for d = 1:length(lwdataset)
    % choose dataset
    lwdata.data = lwdataset(d).data;
    lwdata.header = lwdataset(d).header;

    % re-reference
    option = struct('reference_list', {{params.ref}}, 'apply_list', {{params.chanlocs.labels}},...
    'suffix', params.suffix{1}, 'is_save', 1);
    lwdata = FLW_rereference.get_lwdata(lwdata, option); 

    % update subset
    lwdataset(d).data = lwdata.data;
    lwdataset(d).header = lwdata.header;
end
AP_info(subject_idx).preprocessing(16).process = 'LEPs: re-reference to create dataset for N1 analysis';
AP_info(subject_idx).preprocessing(16).params.ref = params.ref;
AP_info(subject_idx).preprocessing(16).suffix = params.suffix{1};
AP_info(subject_idx).preprocessing(16).date = sprintf('%s', date);

% compute ICA matrix 
fprintf('computing ICA matrix (%s components) ...\n', AP_info(subject_idx).preprocessing(14).params.kept)
option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', AP_info(subject_idx).preprocessing(14).params.kept, 'suffix', params.suffix{2}, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
fprintf('done.\n\n')

% extract ICA parameters
fprintf('extracting ICA parameters...\n')
matrix.mix = lwdataset(1).header.history(end).option.mix_matrix;
matrix.unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
params.chanlocs = lwdataset(1).header.chanlocs;
for i = 1:size(matrix.mix, 2)
    params.ICA_labels{i} = ['IC',num2str(i)];
end
params.ICA_SR = 1/lwdataset(1).header.xstep;

% encode
AP_info(subject_idx).preprocessing(17).process = 'LEPs: second ICA matrix computed';
AP_info(subject_idx).preprocessing(17).params.method = 'runica';
AP_info(subject_idx).preprocessing(17).params.components = AP_info(subject_idx).preprocessing(14).params.kept;
AP_info(subject_idx).preprocessing(17).params.chanlocs = params.chanlocs;
AP_info(subject_idx).preprocessing(17).params.labels = params.ICA_labels;
AP_info(subject_idx).preprocessing(17).params.SR = params.ICA_SR;
AP_info(subject_idx).preprocessing(17).params.matrix = matrix;
AP_info(subject_idx).preprocessing(17).suffix = params.suffix{2};
AP_info(subject_idx).preprocessing(17).date = sprintf('%s', date);

% update dataset and adjust for letswave 6
fprintf('updating dataset...\n')
for a = 1:length(lwdataset)
    % update filtred dataset
    dataset.LEP_N1(a).name = subset(a).name;;
    dataset.LEP_N1(a).header = lwdataset(a).header;
    dataset.LEP_N1(a).data = lwdataset(a).data;

    % adjust header for letswave 6
    dataset.LEP_N1(a).header.history(end).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
    dataset.LEP_N1(a).header.history(end).configuration.parameters = dataset.LEP_N1(a).header.history(end).option;  
    [dataset.LEP_N1(a).header.history(end).configuration.parameters.ICA_um] = dataset.LEP_N1(a).header.history(end).configuration.parameters.unmix_matrix; 
    [dataset.LEP_N1(a).header.history(end).configuration.parameters.ICA_mm] = dataset.LEP_N1(a).header.history(end).configuration.parameters.mix_matrix; 
    dataset.LEP_N1(a).header.history(end).configuration.parameters = rmfield(dataset.LEP_N1(a).header.history(end).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
    header = dataset.LEP_N1(a).header;
    save(sprintf('%s.lw6', dataset.LEP_N1(a).header.name), 'header');
end

% unmix data
fprintf('unmixing data...\n')
for b = 1:length(dataset.LEP_N1)
    for e = 1:size(dataset.LEP_N1(b).data, 1)
        dataset.unmixed_N1(b).name = dataset.LEP_N1(b).name;
        dataset.unmixed_N1(b).header = dataset.LEP_N1(b).header;
        dataset.unmixed_N1(b).data(e, :, 1, 1, 1, :) = matrix.unmix * squeeze(dataset.LEP_N1(b).data(e, :, 1, 1, 1, :));        
    end
end
visual.data = [];
for c = 1:length(dataset.unmixed_N1)
    for e = 1:size(dataset.unmixed_N1(c).data, 1)
        visual.data(end + 1, :, :) = squeeze(dataset.unmixed_N1(c).data(e, :, 1, 1, 1, :));        
    end
end
visual.data = permute(visual.data, [2, 1, 3]);

% plot component topographies and timecourse
addpath(genpath([folder.toolbox '\letswave 6']));
visual.x = (dataset.unmixed_N1(1).header.xstart : dataset.unmixed_N1(1).header.xstep : dataset.unmixed_N1(1).header.xstart + dataset.unmixed_N1(1).header.datasize(6)*dataset.unmixed_N1(1).header.xstep - dataset.unmixed_N1(1).header.xstep)*1000; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for f = 1:size(visual.data, 1)
    % plot the topography
    subplot(ceil(size(visual.data, 1)/2), 6, (f-1)*3 + 1);
    topoplot(double(matrix.mix(:, f)'), params.chanlocs, 'maplimits', [-3 3], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
    set(gca,'color',[1 1 1]);
    title(params.ICA_labels{f})

    % plot average timecourse 
    subplot(ceil(size(visual.data, 1)/2), 6, (f-1)*3 + 2);
    plot(visual.x, squeeze(mean(visual.data(f, :, :), 2)), 'Color', 'black', 'LineWidth', 2)
    xlim([-50 500])
    xlabel('time (ms)');
    ylabel('trial');

    % plot timecourse per trial 
    subplot(ceil(size(visual.data, 1)/2), 6, (f-1)*3 + 3);
    imagesc(visual.x, 1:size(visual.data, 2), squeeze(visual.data(f, :, :)));
    xlabel('time (ms)');
    ylabel('trial');
end
saveas(gcf, sprintf('%s\\figures\\ICA_N1_%s.png', folder.output, AP_info(subject_idx).ID));

% remove N2P2 IC(s)
fprintf('please perform ICA manually in letswave.\n')
letswave
for f = 1:length(dataset.LEP_N1)
    filenames{f} = sprintf('%s %s.mat', params.suffix{3}, dataset.LEP_N1(f).header.name);
end
wait4files(filenames);
fprintf('\n')

% ask for the input
prompt = {'N2P2 component(s):'};
dlgtitle = 'N2P2 removed';  
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput 

% encode ICA 
AP_info(subject_idx).preprocessing(18).process = 'LEPs: IC(s)s containing N2P2 component removed';
AP_info(subject_idx).preprocessing(18).params.kept = AP_info(subject_idx).preprocessing(17).params.components - length(str2num(input{1}));
AP_info(subject_idx).preprocessing(18).params.removed = str2num(input{1});
AP_info(subject_idx).preprocessing(18).suffix = params.suffix{3};
AP_info(subject_idx).preprocessing(18).date = sprintf('%s', date);

% load filtered N1 dataset
fprintf('updating dataset...\n')
params.prefix = regexp(dataset.LEP_N1(1).header.name, ['(.*)\s+' AP_info(subject_idx).ID], 'tokens');
params.prefix = params.prefix{1}{1};
dataset_old = dataset;
data2load = dir(sprintf('%s %s*%s*', params.suffix{3}, params.prefix, AP_info(subject_idx).ID));
dataset = reload_dataset(AP_info(subject_idx).ID, data2load, 'LEP_N1');
dataset_new = dataset;
[dataset_old.LEP_N1] = dataset_new.LEP_N1;
dataset = dataset_old;
fprintf('\n')
clear data2load dataset_old dataset_new

% baseline correct both LEP datasets
fprintf('orrecting for baseline...\n')
addpath(genpath([folder.toolbox '\letswave 7']));
for d = 1:length(subset)
    % N2P2 data
    lwdata.data = subset(d).data;
    lwdata.header = subset(d).header;
    option = struct('operation', 'substract', 'xstart', params.baseline(1), 'xend', params.baseline(2), ...
        'suffix', params.suffix{4},'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata,option);

    % update the dataset
    dataset.LEP_N2P2(d).name = subset(d).name; 
    dataset.LEP_N2P2(d).data = lwdata.data; 
    dataset.LEP_N2P2(d).header = lwdata.header; 

    % N1 data
    lwdata.data = dataset.LEP_N1(d).data;
    lwdata.header = dataset.LEP_N1(d).header;
    option = struct('operation', 'substract', 'xstart', params.baseline(1), 'xend', params.baseline(2), ...
        'suffix', params.suffix{4},'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata,option);

    % update the dataset
    dataset.LEP_N1(d).data = lwdata.data; 
    dataset.LEP_N1(d).header = lwdata.header; 
end
AP_info(subject_idx).preprocessing(19).process = 'LEPs: N1 and N2P2 datasets baseline corrected';
AP_info(subject_idx).preprocessing(19).params.method = 'subtraction';
AP_info(subject_idx).preprocessing(19).params.baseline = params.baseline;
AP_info(subject_idx).preprocessing(19).suffix = params.suffix{4};
AP_info(subject_idx).preprocessing(19).date = sprintf('%s', date);
fprintf('\n')

% select EOIs for CWT
fprintf('selecting EOIs...\n')
for d = 1:length(dataset.LEP_N1)
    % identify EOI
    % split data per components
end

data.N1.cond1 = []; data.N2P2.cond1 = []; 
data.N1.cond2 = []; data.N2P2.cond2 = []; 
if length(file2process) == params.n_files
    % prepare data
    for f = 1:length(file2process)
        % load the dataset
        option = struct('filename', sprintf('%s\\%slw6', file2process(f).folder, file2process(f).name(1:end-3)));
        lwdata = FLW_load.get_lwdata(option);

        % bandpass filter 
        option = struct('filter_type', 'bandpass', 'high_cutoff', params.bandpass(2), 'low_cutoff', params.bandpass(1),'filter_order', 4, 'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata,option);

        % split data according to conditions    
        if contains(file2process(f).name, params.prefix{2})          % N1
            % identify EOI
            if contains(file2process(f).name, 'right')
                eoi = find(strcmp({lwdata.header.chanlocs.labels}, params.EOI{1}));
                eoi_all{f} = params.EOI{1};
            else
                eoi = find(strcmp({lwdata.header.chanlocs.labels}, params.EOI{2}));
                eoi_all{f} = params.EOI{2};
            end
            
            % subset the data and save to a datase                  
            if contains(file2process(f).name, params.conditions{1})
                for d = 1:lwdata.header.datasize(1)
                    data.N1.cond1(end+1, :) = squeeze(lwdata.data(d, eoi, 1, 1, 1, :));
                end
            else
                for d = 1:lwdata.header.datasize(1)
                    data.N1.cond2(end+1, :) = squeeze(lwdata.data(d, eoi, 1, 1, 1, :));
                end
            end
        else                                                        % N2 and P2
            % identify EOI
            eoi = find(strcmp({lwdata.header.chanlocs.labels}, params.EOI{3}));
            eoi_all{f} = params.EOI{3};

            % subset the data and save to a datase                  
            if contains(file2process(f).name, params.conditions{1})
                for d = 1:lwdata.header.datasize(1)
                    data.N2P2.cond1(end+1, :) = squeeze(lwdata.data(d, eoi, 1, 1, 1, :));
                end
            else
                for d = 1:lwdata.header.datasize(1)
                    data.N2P2.cond2(end+1, :) = squeeze(lwdata.data(d, eoi, 1, 1, 1, :));
                end
            end
        end
    end
else
    fprintf('Incorrect number of datasets found in the directory: %d\n', length(file2process));
    return;
end

% perform CWT filtering
fprintf('filtering using CWT masking ...\n')
addpath(genpath([folder.toolbox '\STEP']));
params.freq = params.bandpass(1):1:params.bandpass(2);
params.fs = 1/lwdata.header.xstep;
params.timepoints = lwdata.header.xstart + (0:lwdata.header.datasize(6) + lwdata.header.xstart) * lwdata.header.xstep;
fig1 = figure('Name', 'CWT filters', 'Position', [100, 100, 750, 550]);
for p = 1:length(params.dataset)
    for c = 1:length(params.conditions)
        % define inputs
        statement = sprintf('data_in = data.%s.cond%d(:, :);', params.dataset{p}, c);
        eval(statement)
        data_in = permute(data_in, [2,1]);         

        % compute the mask
        P_mask(p, c, :, :) = model_generation(data_in, params.freq, params.fs, params.timepoints', params.mask_threshold);
        P_mask(p, c, :, 1:find(params.timepoints == 0)) = 0;

        % filter using the mask
        data_out = tf_filtering(data_in, params.freq, params.fs, squeeze(P_mask(p, c, :, :)));
        statement = sprintf('data_filtered.%s.cond%d(:, :) = data_out'';', params.dataset{p}, c);
        eval(statement)

        % plot grand average filtered vs. unfiltered data
        visual = struct;
        t = tinv(0.975, size(data_in, 2) - 1); 
        visual.data(1, :) = mean(data_in, 2)'; visual.data(2, :) = mean(data_out, 2)';
        visual.sem(1, :) = visual.data(1, :) + std(data_in', 0, 1) / sqrt(size(data_in, 2)); 
        visual.sem(2, :) = visual.data(2, :) + std(data_out', 0, 1) / sqrt(size(data_out, 2)); 
        visual.CI_upper(1, :) = visual.data(1, :) + t * visual.sem(1, :); visual.CI_lower(1, :) = visual.data(1, :) - t * visual.sem(1, :); 
        visual.CI_upper(2, :) = visual.data(2, :) + t * visual.sem(2, :); visual.CI_lower(2, :) = visual.data(2, :) - t * visual.sem(2, :); 
        subplot(length(params.dataset), length(params.conditions), (p-1)*length(params.conditions) + c)
        hold on
        plot_ERP(visual.data, visual.CI_upper, visual.CI_lower, params.timepoints, 'labels', {'unfiltered' 'CWT filtered'}, ...
            'colours', params.colours, 'alpha', params.alpha, 'shading', 'off', 'legend_loc', 'southwest');
        title(sprintf('%s: %s', params.dataset{p}, params.conditions{c}), 'fontsize', 16, 'fontweight', 'bold')
        if (p-1)*length(params.dataset) + c ~= length(params.dataset) * length(params.conditions)
            legend('off');
        end
    end
end
AP_info(subject_idx).preprocessing(20).process = 'LEPs: target electrode selected and CWT filtered';
AP_info(subject_idx).preprocessing(20).params.EOIs = unique(eoi_all);
AP_info(subject_idx).preprocessing(20).params.P_mask = P_mask;
AP_info(subject_idx).preprocessing(20).date = sprintf('%s', date);

% save the output figure
saveas(fig1, sprintf('%s\\figures\\CWT_%s.png', folder.output, AP_info.single_subject(subject_idx).ID))

% save info structure and move on
save(output_file, 'AP_info', '-append');
clear a b c d e f i data_idx filenames header input lwdata lwdataset matrix option subset visual
fprintf('section 6 finished.\n\n')

%% functions
function dataset = reload_dataset(ID, data2load, fieldname)
% =========================================================================
% reloads pre-processed EEG data of a single subject for following 
% processing steps 
% input:    - ID = subject identifier
%           - data2load = list of datasets to load
%           - fieldname
% =========================================================================  
    % initiate output
    dataset = struct;
    
    % sort datasets
    header_idx = logical([]);
    data_idx = logical([]);
    for d = 1:length(data2load)
        if contains(data2load(d).name, 'lw6') 
            header_idx(d) = true;
            data_idx(d) = false;
        elseif contains(data2load(d).name, 'mat') 
            header_idx(d) = false;
            data_idx(d) = true;
        end
    end
    headers = data2load(header_idx);
    datas = data2load(data_idx);
    
    % load 
    if length(datas) == length(headers) 
        for d = 1:length(datas)
            % header
            load(sprintf('%s\\%s', headers(d).folder, headers(d).name), '-mat')
            dataname = regexp(header.name, [ID '\s*(.*)'], 'tokens');
            dataname = dataname{1}{1};
            dataset.(fieldname)(d).name = dataname;
            dataset.(fieldname)(d).header = header; 
    
            % data
            load(sprintf('%s\\%s', datas(d).folder, datas(d).name))
            dataset.(fieldname)(d).data = data;
        end
    else
        error('ERROR: Wrong number of available datasets to load! Check manually.')
    end
end
function wait4files(filenames)
% =========================================================================
% waitForFiles pauses the script until all required files appear in the
% working working directory
% --> file names are specified in a cell array
% =========================================================================    
    % loop to wait for files
    while true
        allFilesExist = true;
        
        % check for each file
        for i = 1:length(filenames)
            if isempty(dir(filenames{i}))
                allFilesExist = false;
                break;
            end
        end
        
        % if all files exist, break the loop
        if allFilesExist
            break;
        end
        
        % pause for 2s
        pause(2);
    end
end
function select_EEG(data_name, time_vector, data_visual, s, duration, trigpos, varargin)
    % check for colours
    if ~isempty(varargin)
        a = find(strcmpi(varargin, 'colour'));
        if ~isempty(a)
            col_default = false;   
            colour = varargin{a + 1};
        else
            col_default = true;   
        end     
    else
        col_default = true;
    end

    % create the figure
    fig = figure('Name', data_name, 'NumberTitle', 'off');
    hold on;    

    % plot the EEG data
    n_plots = size(data_visual, 1);
    for c = 1:n_plots
        if col_default
            plot(time_vector, data_visual(c, :), 'LineWidth', 1.2); 
        else
            plot(time_vector, data_visual(c, :), 'Color', colour(c, :), 'LineWidth', 1.2); 
        end
    end

    % configure axes
    xlabel('time (s)');
    ylabel('amplitude (µV)');
    xlim([time_vector(1), time_vector(end)]);
    y_lim = ylim;

    % plot the trigger
    line([trigpos, trigpos], y_lim, 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--')

    % add draggable window for selection
    rect = rectangle('Position', [0, y_lim(1), duration, y_lim(2) - y_lim(1)], ...
        'EdgeColor', 'r', 'LineWidth', 4);
    draggable_rect = make_draggable(rect, time_vector, s);

    % add 'Proceed' button
    proceed_btn = uicontrol('Style', 'pushbutton', 'String', 'Proceed', ...
        'Position', [470, 5, 70, 45], 'Callback', @proceed_callback);

    % callback function for the 'Proceed' button
    drawnow;
    uiwait(fig);
    function proceed_callback(~, ~)
        uiresume(fig);
        close(fig);
    end

    hold off;
end
function draggable_rect = make_draggable(rect, x, s)
    draggable_rect = struct('Rectangle', rect, 'Position', rect.Position);
    set(rect, 'ButtonDownFcn', @startDragFcn);
    function startDragFcn(~, ~)
        set(gcf, 'WindowButtonMotionFcn', @draggingFcn);
        set(gcf, 'WindowButtonUpFcn', @endDragFcn);
    end
    function draggingFcn(~, ~)
        pt = get(gca, 'CurrentPoint');
        new_pos = [pt(1,1), rect.Position(2), rect.Position(3), rect.Position(4)];
        draggable_rect.Position = new_pos;
        set(rect, 'Position', new_pos);
    end
    function endDragFcn(~, ~)
        set(gcf, 'WindowButtonMotionFcn', '');
        set(gcf, 'WindowButtonUpFcn', '');
        
        % output final position
        pos = draggable_rect.Position;
        start_time = round(pos(1), 3);
        end_time = round(start_time + pos(3), 3);
        limits = [start_time, end_time];
        assignin('base', 'limits', limits);
    end
end
function [peak_x, peak_y] = gfp_plot(x, y, xstep, labeled, varargin)
    % check whether to plot labels (default)
    if ~isempty(varargin)
        a = find(strcmpi(varargin, 'max_peaks'));
        if ~isempty(a)
            max_peaks = varargin{a + 1};
        end
    end
    
    % launch the figure  
    plot(x, y)
    yl = get(gca, 'ylim');
    cla
    
    % plot interpolated part
    hold on
    xlim([x(1), x(end)])
    
    % plot data, mark stimulus
    plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
    line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    
    % find peaks 
    [pks, locs] = findpeaks(y, 'MinPeakDistance', 10, 'MinPeakProminence', 0.05);
    for a = 1:length(locs)
        if x(1) + locs(a)*xstep <= 100
            idx(a) = false;
        elseif x(1) + locs(a)*xstep > 600
            idx(a) = false;        
        else
            idx(a) = true;
        end
    end
    pks = pks(idx); locs = locs(idx);
    if length(pks) > max_peaks
        pks = pks(1:max_peaks); 
        locs = locs(1:max_peaks);
    end
    
    % calculate peak coordinations
    for a = 1:length(locs)
        peak_x(a) = x(1) + locs(a)*xstep;
        peak_y(a) = pks(a);
    end
    peak_y = double(peak_y);
    
    % plot peaks
    for a = 1:length(locs)
        plot(peak_x(a), peak_y(a), ...
            'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
        line([peak_x(a), peak_x(a)], [yl(1), peak_y(a)], 'LineStyle', ':', 'Color', [1, 0, 0], 'LineWidth', 1.5)
        
        % label with latency (ms)
        if strcmp(labeled, 'on') 
            text(peak_x(a) - 180, peak_y(a), sprintf('%1.0f ms', peak_x(a)), 'Color', [1 0 0], 'FontSize', 14)
        end
    end
    
    % add parameters
    set(gca, 'fontsize', 14)
    ylim(yl)
    xlabel('time (s)')
    ylabel('potential (\muV)')
end
function topo_plot(header, data, x_pos, map_lims)
    varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};
    
    % fetch data to display
    x_visual = ceil((x_pos/1000 - header.xstart)/header.xstep);
    vector = data(1, :, 1, 1, 1, x_visual);
    
    %fetch chanlocs
    chanlocs = header.chanlocs;
    
    %parse data and chanlocs 
    i=1;
    for chanpos=1:size(chanlocs,2);
        vector2(i)=double(vector(chanpos));
        chanlocs2(i)=chanlocs(chanpos);
        i=i+1;
    end;
    
    topoplot(vector2,chanlocs2,varargin{:});
    set(gcf,'color',[1 1 1]);
end
function plot_ERP(visual_data, visual_CI_upper, visual_CI_lower, x, varargin)
    % default colour palette
    colour_palette = [
        0, 0, 1;       % Blue
        1, 0, 0;       % Red
        0, 1, 0;       % Green
        1, 1, 0;       % Yellow
        1, 0, 1;       % Magenta
        0, 1, 1;       % Cyan
        0.5, 0.5, 0.5; % Gray
        1, 0.5, 0;     % Orange
        0, 0.5, 0.5;   % Teal
        0.5, 0, 0.5    % Purple
    ];

    % check for varargins
    if ~isempty(varargin)
        % y limits
        a = find(strcmpi(varargin, 'ylim'));
        if ~isempty(a)
            y_limits = varargin{a + 1};
        else
            y_limits = [0,0];
        end

        % labels
        b = find(strcmpi(varargin, 'labels'));
        if ~isempty(b)
            labels = varargin{b + 1};
        else
            for c = 1:size(visual_data, 1)
                labels{c} = sprintf('condition %d', c);
            end
        end

        % colours
        d = find(strcmpi(varargin, 'colours'));
        if ~isempty(d)
            col = varargin{d + 1};
        else
            for c = 1:size(visual_data, 1)
                col(c, :) = colour_palette(c, :);
            end
        end

        % alpha
        e = find(strcmpi(varargin, 'alpha'));
        if ~isempty(e)
            alpha = varargin{e + 1};
        else
            alpha = 0.2;
        end

        % shading - default on
        f = find(strcmpi(varargin, 'shading'));
        if ~isempty(f) && strcmp(varargin{f + 1}, 'off')
            shading = false;
        else
            shading = true;
        end

        % legend location
        g = find(strcmpi(varargin, 'legend_loc'));
        if ~isempty(g) 
            legend_loc = varargin{g + 1};
        else
            legend_loc = 'southeast';
        end        
    end

    % loop through datasets to plot
    for t = 1:size(visual_data, 1) 
        P(t) = plot(x, visual_data(t, :), 'Color', col(t, :), 'LineWidth', 2.5);
        hold on
        if shading
            F(t) = fill([x fliplr(x)],[visual_CI_upper(t, :) fliplr(visual_CI_lower(t, :))], ...
            col(t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
            hold on
        end
    end

    % check y limits
    if y_limits(1) == 0 && y_limits(2) == 0
        y_limits = ylim;
    end

    % plot stimulus
    line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 3, 'LineStyle', '--')

    % % legend
    % legend(P, labels, 'Location', legend_loc, 'fontsize', 9)
    % legend('boxoff');

    % axes
    box off;
    ax = gca;
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    % ax.TickDir = 'out'; 
    ax.XColor = [0.5020    0.5020    0.5020]; 
    ax.YColor = [0.5020    0.5020    0.5020]; 

    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 13)
    ylim(y_limits)
    xlim([x(1), x(end)])  
    set(gca, 'Layer', 'Top')
    set(gca, 'YDir', 'reverse');
end
function plot_box(data_visual, datatype, col, labels)
    % determine x limits
    xl = [0.25, size(data_visual, 2) + 0.5];

    % % add zero line
    % line(xl, [0, 0], 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)

    switch datatype
        case 'amplitude' 
            % boxplot
            for t = 1:size(data_visual, 2)
                P(t) = boxchart(t * ones(size(data_visual, 1), 1), data_visual(:, t), ...
                    'BoxFaceColor', col(t, :), 'BoxWidth', 0.3, 'WhiskerLineColor', [0 0 0], 'MarkerColor', [0 0 0]);
                hold on
            end 
    
            % add legend
            lgd = legend(P, labels, 'Location', 'southwest');
            lgd.FontSize = 14;
            legend('boxoff')
        
            % y label
            ylabel('\Delta amplitude (\muV \pm SEM)')
    
            % other parameters
            xlim(xl)
            set(gca, 'XTick', []);
            set(gca, 'FontSize', 14) 
            set(gca, 'layer', 'top');
        case 'latency'
            % boxplot
            for t = 1:size(data_visual, 2)
                P(t) = boxchart(t * ones(size(data_visual, 1), 1), data_visual(:, t), 'Orientation', 'horizontal', ...
                    'BoxFaceColor', col(t, :), 'BoxWidth', 0.3, 'WhiskerLineColor', [0 0 0], 'MarkerColor', [0 0 0]);
                hold on
            end 
    
            % add legend
            lgd = legend(P, labels, 'Location', 'southwest');
            lgd.FontSize = 14;
            legend('boxoff')
        
            % y label
            xlabel('\Delta latency (ms \pm SEM')
    
            % other parameters
            set(gca, 'YTick', []);
            set(gca, 'FontSize', 14) 
            set(gca, 'layer', 'top');
    end
end