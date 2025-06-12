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
%   1) fill in subject & session information 
%           --> allows user to manually encode information about 
%               the subject, session, stimulation parameters...
%           --> all this and following information is encoded in AP_info 
%   2) import EEG data
%           --> searches in the input directory to identify all datasets
%           --> imports to MATLAB as variable 'dataset'
%           --> encodes all datasets' names to the metadata structure 
%   3) pre-processing block 1 
%           --> assigns electrode coordinates
%           --> removes DC + detrends the whole recording
%           --> bandpass filters (Butterworth, 4th order) - default
%               [0.1,80]Hz (modifiable)
%           --> notch filters at 50Hz
%           --> downsamples by a ratio - default 5 (modifiable)
%           --> updates 'dataset'
%           --> encodes processing steps to 'NLEP_info'
%   4) visual check 1
%           --> opens letswave - check for bad channels before
%               re-referencing to common average
%   5) pre-process LEPs and SEPs
%           --> re-references to common average
%           --> checks the number of events and re-labels them 'LEP'/'SEP'
%           --> segments relative to the event - default [-300 1000]ms (modifiable)
%           --> removes DC + detrends each epoch
%           --> encodes processing steps to 'NLEP_info'
%   6) pre-process RS-EEG 
%           --> RS-EEG: 
%                   - re-references to common average
%                   - selects 1 min of clean recording
%                   - chunks into non-overlapping epochs - default 1s
%                   - removes DC + detrends each epoch
%           --> pre-stimulus RS-EEG:
%                   - re-references to common average
%                   - selects 1s epochs from 'relaxed' pre-stimulus EEG -
%                       default [-2 -1]s preceding 'get ready' trigger
%                   - selects 1s epochs from 'ready' pre-stimulus EEG - 
%                       default 1s before the stimulus  
%                   - removes DC + detrends each epoch
%   7) visual check 2 
%           --> remove bad epochs in letswave
%   8) encode deleted ERP epochs, discard associated RS-EEG epochs
%   9) compute ICA - ERPs and RS-EEG together
%           --> ICA is run manually in letswave - recommended parameters:
%               merged all datasets, restrict to 25 ICs
%           --> plots component topographies and spectral content
%   10) encode info and finish pre-processing
%           --> user input: interpolated channels, removed epochs, removed
%           --> corrects to baseline - subtracts the average [-250 0]ms
%           --> averages across trials (saves before and after)
%   11) LEPs: identify N2P2 component and subtract it for N1 analysis
%           --> re-references to chosen frontal electrode - default AFz
%           --> second ICA is run manually in letswave - IC number
%           restricted to components leftover after the first ICA
%           --> plots component topographies and spectral content
%           --> waits for user to select and remove N2P2 component
%           --> extracta average N1 before and after filtering, plots
%           --> encodes to the output structure 
%   12) single-trial LEP analysis: time domain
%           --> selects target electrodes based on LEP component
%           --> performs CWT filtering - default mask threshold 0.85
%           --> identifies and saves average peak measures (cleaned)
%           --> prepares template data based on average
%           --> computes regressors based on the template, plots them
%           --> performs the regression on single trial data
%           --> extracts single-trial peak measures
%   13) average LEP analysis: time domain
%           --> extracts peak measures from data averaged across all trials 
%               for each subject, condition, and block 
%   14) RS-EEG analysis: spectral decomposition, sensor space
%           --> compute PSD for each trial and electrode, 
%               extract aperiodic component parameters 
%           --> compute PSD for averaged data at each electrode, 
%               extract aperiodic component parameters 
%           --> computes averages for target and control regions 
%               (pre-stimulus data only)
%   15) exports to R for statistics and visualization
%           --> long-format .csv table
%   16) ERP group-average visualization
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
addpath(genpath([folder.toolbox '\letswave 6']));
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
% extraction of ERP amplitudes and latencies
% single-trial prestim pre-processing at sensor level
% single-trial source activity estimation --> processing at source level


%% 11) LEPs: identify N2P2 component and subtract it for N1 analysis
% ----- section input -----
param.prefix = 'icfilt ica_all dc ep reref ds notch bandpass dc';
param.suffix = {'reref_AFz' 'ica_N1' 'icfilt' 'bl'};
param.n_files = 4;
param.ref = 'AFz';
param.EOI = {'C5' 'C6'};
param.TOI = [0.1 0.22];
param.baseline = [-0.25 0];
% -------------------------

% ask for subject number
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% update the info structure
load(output_file, 'AP_info');

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% identify LEP datasets
cd(folder.input)
file2process = dir(sprintf('%s*%s LEP*.mat', param.prefix, AP_info.single_subject(subject_idx).ID));
cd(folder.data)

% re-reference to chosen frontal central electrode
fprintf('re-referencing to %s...\n', param.ref)
for d = 1:length(file2process)
    % load dataset
    filepath = sprintf('%s\\%s.lw6', file2process(d).folder, file2process(d).name(1:end-4));
    option = struct('filename', filepath);
    lwdata = FLW_load.get_lwdata(option);

    % re-reference
    option = struct('reference_list', {{param.ref}}, 'apply_list', {{lwdata.header.chanlocs(1:length(lwdata.header.chanlocs)).labels}},...
    'suffix', param.suffix{1}, 'is_save', 1);
    lwdata = FLW_rereference.get_lwdata(lwdata, option);            
end

% identify outcome filenames
if length(file2process) == param.n_files
    for i = 1:length(file2process)
        filenames{i} = sprintf('%s %s',[param.suffix{2} ' ' param.suffix{1}], file2process(i).name);
    end
else
    fprintf('Incorrect number of datasets for ICA found in the directory: %d\n', length(file2process));
    return;
end

% open letswave and manually run ICA
fprintf('You can run the ICA now - for subject %d, compute %d components!\n', subject_idx, AP_info.single_subject(subject_idx).preprocessing.ICA.ICs_kept)
addpath(genpath([folder.toolbox '\letswave 6']));
letswave

% wait until all files are processed
wait4files(filenames);

% extract ICA matrices 
fprintf('extracting ICA matrices...\n')
load(sprintf('%slw6', filenames{1}(1:end-3)), '-mat');
ICA.matrix = header.history(end).configuration.parameters.ICA_mm;
ICA.unmix = header.history(end).configuration.parameters.ICA_um;
ICA.chanlocs = header.chanlocs;
for i = 1:size(ICA.unmix, 1)
    ICA.labels{i} = ['IC',num2str(i)];
end
ICA.fs = 1/header.xstep;

% unmix data
ICA.data = [];
fprintf('unmixing the data: dataset ')
for d = 1:length(filenames)
    fprintf('%d ...', d)
    load(filenames{d});
    for e = 1:size(data, 1)
        ICA.data(end + 1, :, :) = ICA.unmix * squeeze(data(e, :, 1, 1, 1, :));        
    end
end
ICA.data = permute(ICA.data, [2, 1, 3]);
fprintf('\n') 

% determine plotting parameters
x = (header.xstart : header.xstep : header.xstart + header.datasize(6)*header.xstep - header.xstep)*1000; 
x_start = round((param.TOI(1) - header.xstart)/header.xstep);
x_end = round((param.TOI(2) - header.xstart)/header.xstep);

% plot component topographies and spectral content
figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for f = 1:size(ICA.unmix, 1)
    % plot the topography
    subplot(ceil(size(ICA.unmix, 1)/2), 4, (f-1)*2 + 1);
    topoplot(ICA.matrix(:, f), ICA.chanlocs, 'maplimits', [-3 3], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
    set(gca,'color',[1 1 1]);
    title(ICA.labels{f})

    % plot the timecourse 
    subplot(ceil(size(ICA.unmix, 1)/2), 4, (f-1)*2 + 2);
    imagesc(x, 1:size(ICA.data, 2), squeeze(ICA.data(f, :, :)));
    xlabel('time (ms)');
    ylabel('trial');
end
saveas(gcf, sprintf('%s\\figures\\ICA_N1_%s.png', folder.output, AP_info.single_subject(subject_idx).ID));

% identify outcome filenames and wait for them to appear
for i = 1:length(filenames)
    filenames_filtered{i} = sprintf('%s %s',param.suffix{3}, filenames{i});
end
wait4files(filenames_filtered);
close(gcf)

% baseline correct all LEP data
addpath(genpath([folder.toolbox '\letswave 7']));
for d = 1:length(filenames_filtered)
    % raw data
    option = struct('filename', sprintf('%slw6', filenames{d}(18:end-3)));
    lwdata = FLW_load.get_lwdata(option);
    option = struct('operation', 'substract', 'xstart', param.baseline(1), 'xend', param.baseline(2), ...
        'suffix', param.suffix{4},'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata,option);

    % raw re-referenced data
    option = struct('filename', sprintf('%slw6', filenames{d}(8:end-3)));
    lwdata = FLW_load.get_lwdata(option);
    option = struct('operation', 'substract', 'xstart', param.baseline(1), 'xend', param.baseline(2), ...
        'suffix', param.suffix{4},'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata,option);

    % filtered data
    option = struct('filename', sprintf('%slw6', filenames_filtered{d}(1:end-3)));
    lwdata = FLW_load.get_lwdata(option);
    option = struct('operation', 'substract', 'xstart', param.baseline(1), 'xend', param.baseline(2), ...
        'suffix', param.suffix{4},'is_save', 1);
    lwdata = FLW_baseline.get_lwdata(lwdata,option);
end

% extract average N1 before and after filtering
for d = 1:length(filenames_filtered)
    % identify EOI
    if contains(filenames_filtered{d}, 'right')
        eoi(d) = find(strcmp({ICA.chanlocs.labels}, param.EOI{1}));
    else
        eoi(d) = find(strcmp({ICA.chanlocs.labels}, param.EOI{2}));
    end

    % extract N1 from unfiltered data 
    load(sprintf('%s %s', param.suffix{4}, filenames{d}(8:end)))
    data_in(d, :) = squeeze(mean(data(:, eoi(d), 1, 1, 1, :), 1));
    [peak_raw.amplitude(d), peak_raw.latency(d)] = min(data_in(d, x_start:x_end));
    peak_raw.latency(d) = peak_raw.latency(d) + x_start;
    data_topo_raw(d, :) = squeeze(mean(data(:, :, 1, 1, 1, peak_raw.latency(d) - 5: peak_raw.latency(d) + 5), [1, 6]));

    % extract N1 from filtered data 
    load(sprintf('%s %s', param.suffix{4}, filenames_filtered{d}(1:end)))
    data_filt(d, :) = squeeze(mean(data(:, eoi(d), 1, 1, 1, :), 1));
    [peak_filt.amplitude(d), peak_filt.latency(d)] = min(data_filt(d, x_start:x_end));
    peak_filt.latency(d) = peak_filt.latency(d) + x_start;
    data_topo_filt(d, :) = squeeze(mean(data(:, :, 1, 1, 1, peak_filt.latency(d) - 5: peak_filt.latency(d) + 5), [1, 6]));
end

% plot N1
addpath(genpath([folder.toolbox '\letswave 6']));
figure('Position', [300, 30, 600, 750]);
for d = 1:length(filenames_filtered)
    % plot the timecourse
    subplot(4, 4, (d-1)*4 + [1, 2])
    hold on
    plot(x, data_in(d, :), 'Color', [0.6510    0.6510    0.6510], 'LineWidth', 1.5);
    plot(x, data_filt(d, :), 'Color', [0.8314    0.1647    0.1647], 'LineWidth', 1.5);
    ylim = get(gca, 'ylim');
    plot([0, 0], ylim,  'Color', 'black', 'LineStyle', '--', 'LineWidth', 2)
    set(gca, 'YDir', 'reverse');
    xlim([-100 650])
    title(sprintf('%s', filenames{d}([end-16:end-4])))

    % plot N1 topography in unfiltered signal
    subplot(4, 4, (d-1)*4 + 3)
    hold on
    topoplot(double(data_topo_raw(d, :)'), ICA.chanlocs, 'maplimits', [-5 5], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
    set(gca,'color',[1 1 1]);
    text(0, -0.7, sprintf('%d ms', peak_raw.latency(d) + header.xstart * 1000), 'HorizontalAlignment', 'center');
    if d == 1
        title('unfiltered');
    end

    % plot N1 topography in unfiltered signal
    subplot(4, 4, (d-1)*4 + 4)
    hold on
    topoplot(double(data_topo_filt(d, :)'), ICA.chanlocs, 'maplimits', [-5 5], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
    set(gca,'color',[1 1 1]);
    text(0, -0.7, sprintf('%d ms', peak_filt.latency(d) + header.xstart * 1000), 'HorizontalAlignment', 'center');
    if d == 1
        title('after ICA');
    end
end
saveas(gcf, sprintf('%s\\figures\\ICA_N1_filtered_%s.png', folder.output, AP_info.single_subject(subject_idx).ID));

% % update output dataset
% for d = 1:length(filenames)
%     if contains(filenames{d}, 'left')
%         LEP_average(subject_idx).conditions{d} = filenames{d}(end - 15 : end - 4); 
%     else
%         LEP_average(subject_idx).conditions{d} = filenames{d}(end - 16 : end - 4); 
%     end
% end
% LEP_average(subject_idx).N1.raw = peak_raw;
% LEP_average(subject_idx).N1.filtered = peak_filt;
% save(output_file, 'LEP_average', '-append');

% update NLEP info
prompt = {'discarded ICs'};
dlgtitle = 'N1 ICA';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
AP_info.single_subject(subject_idx).preprocessing.LEP(1).process = '1 - data re-referenced for N1 analysis'; 
AP_info.single_subject(subject_idx).preprocessing.LEP(1).params.ref = param.ref;
AP_info.single_subject(subject_idx).preprocessing.LEP(1).date = sprintf('%s', date);
AP_info.single_subject(subject_idx).preprocessing.LEP(2).process = '2 - N2P2 component filtered out with ICA'; 
AP_info.single_subject(subject_idx).preprocessing.LEP(2).params.n_components = AP_info.single_subject(subject_idx).preprocessing.ICA.ICs_kept;
AP_info.single_subject(subject_idx).preprocessing.LEP(2).params.matrix = ICA.matrix;
AP_info.single_subject(subject_idx).preprocessing.LEP(2).params.unmix = ICA.unmix;
AP_info.single_subject(subject_idx).preprocessing.LEP(2).params.removed = input;
AP_info.single_subject(subject_idx).preprocessing.LEP(2).date = sprintf('%s', date);
AP_info.single_subject(subject_idx).preprocessing.LEP(3).process = '3 - all LEP datasets baseline corrected'; 
AP_info.single_subject(subject_idx).preprocessing.LEP(3).params.baseline = param.baseline;
AP_info.single_subject(subject_idx).preprocessing.LEP(3).date = sprintf('%s', date);
save(output_file, 'AP_info', '-append');

clear param file2process d i filenames filepath option lwdata data header ICA x x_start x_end e f filenames_filtered ...
    data_in peak_raw data_topo_raw data_filt peak_filt data_topo_filt ylim prompt dlgtitle dims definput input eoi

%% 12) single-trial LEP analysis: time domain
% ----- section input -----
param.prefix = {'bl icfilt ica_all dc ep reref ds notch bandpass dc' 'bl icfilt ica_N1 reref_AFz'};  
param.suffix = {'butt'};
param.n_files = 8;
param.bandpass = [1 30];
param.peak = {'N1' 'N2' 'P2'};
param.dataset = {'N1' 'N2P2'};
param.EOI = {'C5' 'C6' 'Cz'};
param.mask_threshold = 0.85;
param.alpha = 0.2;
param.inverse_default = 300;
param.colours = [0.5020    0.5020    0.5020; 0.8314    0.1647    0.1647];
param.TOI_span = [0.14, 0.16, 0.20];
param.TOI_estimate = [0, 0.300; 0, 0.450; 0, 0.600];
param.gaussian_size = 20;
param.gaussian_sigma = 5;
% ------------------------- 

% ask for subject number
if ~exist('subject_idx')
    prompt = {'subject number:'};
    dlgtitle = 'subject';
    dims = [1 40];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    subject_idx = str2num(input{1,1});
end
clear prompt dlgtitle dims definput input

% load the output structures
%load(output_file, 'NLEP_info', 'NLEP_measures', 'NLEP_data');

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% identify LEP datasets to process
folder.input = uigetdir(pwd, 'Choose the input folder');            
file2process = [dir(sprintf('%s\\%s_%s\\%s*%s LEP*.mat', folder.input, study, AP_info.single_subject(subject_idx).ID, param.prefix{1}, AP_info.single_subject(subject_idx).ID)); ...
    dir(sprintf('%s\\%s_%s\\%s*%s LEP*.mat', folder.input, study, AP_info.single_subject(subject_idx).ID, param.prefix{2}, AP_info.single_subject(subject_idx).ID))];

% determine conditions
for c = 1:length(AP_info.single_subject(subject_idx).condition)
    param.conditions{c} = replace(AP_info.single_subject(subject_idx).condition{c}, '_', ' ');
end
AP_info.single_subject(subject_idx).LEP(1).process = '1 - preprocessed data selected';
AP_info.single_subject(subject_idx).LEP(1).params.conditions = param.conditions;
AP_info.single_subject(subject_idx).LEP(1).params.datasets.N1 = {file2process(5:8).name};
AP_info.single_subject(subject_idx).LEP(1).params.datasets.N2P2 = {file2process(1:4).name};
AP_info.single_subject(subject_idx).LEP(1).date = sprintf('%s', date);
NLEP_measures(subject_idx).LEP_avg.conditions = param.conditions;

% prepare data
fprintf('preparing LEP data ...\n')
data.N1.cond1 = []; data.N2P2.cond1 = []; 
data.N1.cond2 = []; data.N2P2.cond2 = []; 
if length(file2process) == param.n_files
    % prepare data
    for f = 1:length(file2process)
        % load the dataset
        option = struct('filename', sprintf('%s\\%slw6', file2process(f).folder, file2process(f).name(1:end-3)));
        lwdata = FLW_load.get_lwdata(option);

        % bandpass filter 
        option = struct('filter_type', 'bandpass', 'high_cutoff', param.bandpass(2), 'low_cutoff', param.bandpass(1),'filter_order', 4, 'suffix', param.suffix{1}, 'is_save', 0);
        lwdata = FLW_butterworth_filter.get_lwdata(lwdata,option);

        % split data according to conditions    
        if contains(file2process(f).name, param.prefix{2})          % N1
            % identify EOI
            if contains(file2process(f).name, 'right')
                eoi = find(strcmp({lwdata.header.chanlocs.labels}, param.EOI{1}));
                eoi_all{f} = param.EOI{1};
            else
                eoi = find(strcmp({lwdata.header.chanlocs.labels}, param.EOI{2}));
                eoi_all{f} = param.EOI{2};
            end
            
            % subset the data and save to a datase                  
            if contains(file2process(f).name, param.conditions{1})
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
            eoi = find(strcmp({lwdata.header.chanlocs.labels}, param.EOI{3}));
            eoi_all{f} = param.EOI{3};

            % subset the data and save to a datase                  
            if contains(file2process(f).name, param.conditions{1})
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
param.freq = param.bandpass(1):1:param.bandpass(2);
param.fs = 1/lwdata.header.xstep;
param.timepoints = lwdata.header.xstart + (0:lwdata.header.datasize(6) + lwdata.header.xstart) * lwdata.header.xstep;
fig1 = figure('Name', 'CWT filters', 'Position', [100, 100, 750, 550]);
for p = 1:length(param.dataset)
    for c = 1:length(param.conditions)
        % define inputs
        statement = sprintf('data_in = data.%s.cond%d(:, :);', param.dataset{p}, c);
        eval(statement)
        data_in = permute(data_in, [2,1]);         

        % compute the mask
        P_mask(p, c, :, :) = model_generation(data_in, param.freq, param.fs, param.timepoints', param.mask_threshold);
        P_mask(p, c, :, 1:find(param.timepoints == 0)) = 0;

        % filter using the mask
        data_out = tf_filtering(data_in, param.freq, param.fs, squeeze(P_mask(p, c, :, :)));
        statement = sprintf('data_filtered.%s.cond%d(:, :) = data_out'';', param.dataset{p}, c);
        eval(statement)

        % plot grand average filtered vs. unfiltered data
        visual = struct;
        t = tinv(0.975, size(data_in, 2) - 1); 
        visual.data(1, :) = mean(data_in, 2)'; visual.data(2, :) = mean(data_out, 2)';
        visual.sem(1, :) = visual.data(1, :) + std(data_in', 0, 1) / sqrt(size(data_in, 2)); 
        visual.sem(2, :) = visual.data(2, :) + std(data_out', 0, 1) / sqrt(size(data_out, 2)); 
        visual.CI_upper(1, :) = visual.data(1, :) + t * visual.sem(1, :); visual.CI_lower(1, :) = visual.data(1, :) - t * visual.sem(1, :); 
        visual.CI_upper(2, :) = visual.data(2, :) + t * visual.sem(2, :); visual.CI_lower(2, :) = visual.data(2, :) - t * visual.sem(2, :); 
        subplot(length(param.dataset), length(param.conditions), (p-1)*length(param.conditions) + c)
        hold on
        plot_ERP(visual.data, visual.CI_upper, visual.CI_lower, param.timepoints, 'labels', {'unfiltered' 'CWT filtered'}, ...
            'colours', param.colours, 'alpha', param.alpha, 'shading', 'off', 'legend_loc', 'southwest');
        title(sprintf('%s: %s', param.dataset{p}, param.conditions{c}), 'fontsize', 16, 'fontweight', 'bold')
        if (p-1)*length(param.dataset) + c ~= length(param.dataset) * length(param.conditions)
            legend('off');
        end
    end
end
AP_info.single_subject(subject_idx).LEP(2).process = '2 - target electrode selected and CWT filtered';
AP_info.single_subject(subject_idx).LEP(2).params.EOIs = unique(eoi_all);
AP_info.single_subject(subject_idx).LEP(2).params.P_mask = P_mask;
AP_info.single_subject(subject_idx).LEP(2).date = sprintf('%s', date);

% save the output figure
saveas(fig1, sprintf('%s\\figures\\CWT_%s.png', folder.output, AP_info.single_subject(subject_idx).ID))

% extract single-trial peak values
fprintf('extracting peak measures ...\n')
fig2 = figure('Name', 'PCA regressors', 'Position', [50, 30, 550, 750]);
for p = 1:length(param.peak)
    for c = 1:length(param.conditions)
        % extract single-trial data and average for template, determine
        % peak polarity
        switch p
            case 1
                % subset data
                statement = sprintf('data_ST = data_filtered.N1.cond%d;', c);
                eval(statement)
                data_avg = mean(data_ST, 1);

                % identify polarity
                negative = true;   
            case 2
                % subset data
                statement = sprintf('data_ST = data_filtered.N2P2.cond%d;', c);
                eval(statement)
                data_avg = mean(data_ST, 1);

                % identify polarity
                negative = true;  
            case 3
                % subset data
                statement = sprintf('data_ST = data_filtered.N2P2.cond%d;', c);
                eval(statement)
                data_avg = mean(data_ST, 1);

                % identify polarity
                negative = false; 
        end

        %% identify and save average peak measures 
        % select a narrower window around the average peak 
        select_EEG(sprintf('%s: %s', param.peak{p}, param.conditions{c}), param.timepoints, data_avg, ...
            1/lwdata.header.xstep, param.TOI_span(p), 0);

        % compute window limits and extract the average peak values 
        x_peak_start = round((limits(1) - lwdata.header.xstart)/lwdata.header.xstep);
        x_peak_end = round((limits(2) - lwdata.header.xstart)/lwdata.header.xstep);
        if negative
            [y_peak, x] = min(data_avg(x_peak_start:x_peak_end));
        else
            [y_peak, x] = max(data_avg(x_peak_start:x_peak_end));
        end
        peak_sample = x_peak_start + x;

        % encode the average peak values and parameters to output
        % structures
        switch p
            case 1
                NLEP_measures(subject_idx).LEP_avg.N1.CWT_filtered.latency(c) = peak_sample*lwdata.header.xstep + lwdata.header.xstart;
                NLEP_measures(subject_idx).LEP_avg.N1.CWT_filtered.amplitude(c) = y_peak;
                AP_info.single_subject(subject_idx).LEP(3).params.TOI_peak.N1(c, :) = [x_peak_start, x_peak_end];
                AP_info.single_subject(subject_idx).LEP(3).params.peak_avg.N1(c) = peak_sample*lwdata.header.xstep + lwdata.header.xstart;
            case 2 
                NLEP_measures(subject_idx).LEP_avg.N2.CWT_filtered.latency(c) = peak_sample*lwdata.header.xstep + lwdata.header.xstart;
                NLEP_measures(subject_idx).LEP_avg.N2.CWT_filtered.amplitude(c) = y_peak;
                AP_info.single_subject(subject_idx).LEP(3).params.TOI_peak.N2(c, :) = [x_peak_start, x_peak_end];
                AP_info.single_subject(subject_idx).LEP(3).params.peak_avg.N2(c) = peak_sample*lwdata.header.xstep + lwdata.header.xstart;
            case 3
                NLEP_measures(subject_idx).LEP_avg.P2.CWT_filtered.latency(c) = peak_sample*lwdata.header.xstep + lwdata.header.xstart;
                NLEP_measures(subject_idx).LEP_avg.P2.CWT_filtered.amplitude(c) = y_peak;
                AP_info.single_subject(subject_idx).LEP(3).params.TOI_peak.P2(c, :) = [x_peak_start, x_peak_end];
                AP_info.single_subject(subject_idx).LEP(3).params.peak_avg.P2(c) = peak_sample*lwdata.header.xstep + lwdata.header.xstart;
        end       
        
        %% prepare template data
        % make sure the peak window is in the main analysed window 
        x_start = round((param.TOI_estimate(p, 1) - lwdata.header.xstart)/lwdata.header.xstep);
        x_end = round((param.TOI_estimate(p, 2) - lwdata.header.xstart)/lwdata.header.xstep);
        if x_start > x_peak_start
            x_start = x_peak_start;
        end
        if x_end < x_peak_end
            x_end = x_peak_end;
        end

        % crop the data 
        x_peak = peak_sample - x_start;
        data_template = data_avg(x_start:x_end); 

        % identify inversion point for N2 and P2
        if p ~= 1
            idx_flip = logical([]);
            for i = 1:size(data_template, 2)-1                
                if data_template(i) < 0 && data_template(i+1) > 0
                    idx_flip(i) = true;
                else
                    idx_flip(i) = false;
                end
            end 
            idx_flip(i+1) = false;
            inverse_point = find(idx_flip);
            for i = 1:length(inverse_point)
                if p == 2 && inverse_point(i) < x_peak 
                    idx_flip(inverse_point(i)) = false;
                elseif p == 3 && inverse_point(i) < NLEP_measures(subject_idx).LEP_avg.N2.CWT_filtered.latency(c)/lwdata.header.xstep
                    idx_flip(inverse_point(i)) = false;
                end
            end
            inverse_point = min(find(idx_flip));
            if isempty(inverse_point)
                inverse_point = param.inverse_default;
            end
        end
                       
        % select N2 or P2 component, if necessary, and smoothen using a
        % Gaussian filter
        x_gaussian = linspace(-param.gaussian_size / 2, param.gaussian_size / 2, param.gaussian_size);       
        gaussian_filter = exp(-x_gaussian .^ 2 / (2 * param.gaussian_sigma ^ 2));       % sigma = std
        gaussian_filter = gaussian_filter / sum(gaussian_filter);                       % normalize the filter
        if p == 2                                                                       % suppress data relative to the inversion, depending on peak 
            idx_flip(inverse_point:end) = 1;
            data_template(idx_flip) = 0;
            data_template = conv(data_template, gaussian_filter, 'same');
        elseif p == 3
            idx_flip(1:inverse_point) = 1;
            data_template(idx_flip) = 0;
            data_template = conv(data_template, gaussian_filter, 'same');
        end

        % encode to NLEP info
        switch p
            case 1
                AP_info.single_subject(subject_idx).LEP(3).params.TOI_window.N1(c, :) = [x_start, x_end];
            case 2 
                AP_info.single_subject(subject_idx).LEP(3).params.TOI_window.N2(c, :) = [x_start, x_end];
                AP_info.single_subject(subject_idx).LEP(3).params.inversion_point(c) = inverse_point*lwdata.header.xstep;
                if c == 1
                    AP_info.single_subject(subject_idx).LEP(3).params.gaussian.width = param.gaussian_size;
                    AP_info.single_subject(subject_idx).LEP(3).params.gaussian.sigma = param.gaussian_sigma;
                end
            case 3
                AP_info.single_subject(subject_idx).LEP(3).params.TOI_window.P2(c, :) = [x_start, x_end];
        end

        %% compute regressors based on the template
        % prepare inputs
        t = param.timepoints(x_start + 1 : x_end + 1);  % time vector
        fs = 1/lwdata.header.xstep;                     % sampling rate
        number_trial = round(fs*0.05);                  % the default latency jitter limit
        if fs > 500
            idx_reg = 3;                                 % if out of memory, this number should be increased
        elseif Fs <= 500
            idx_reg = 1; 
        end
        k = 1:idx_reg/(number_trial-1):2;                % the default compression limit
        template = zeros(length(k), length(data_template));

        % prepare the template
        T = []; template = [];
        for j = 1:length(k)
            T(j,:) = t/k(j);
            ur = resample(double(data_template), round(fs/k(j)), fs);
            peak_point = floor(x_peak * round(fs/k(j))/fs);
            if x_peak - peak_point <= 0
                template(j, 1:length(ur)) = ur; 
            else
                template(j, x_peak - peak_point + 1 : x_peak - peak_point + length(ur)) = ur;  
            end
        end

        % obtain regressors using PCA
        data_start = []; data_PCA = [];
        for j = 1:length(k)
            data_start = zeros(length(k), length(t) + length(k));
            for i = 1:length(k)
                data_start(i, idx_reg*i : length(t)+idx_reg*i-1) = template(j,:);    
            end
            data_start = data_start(:,round(length(k)/2)*idx_reg+1:length(t)+round(length(k)/2)*idx_reg);
            data_PCA = [data_PCA; data_start];
        end
        [W,Y,val] = pca(data_PCA);
        for r = 1:3  
            regressors{p, c}(:, r) = real(Y(r, :));
        end
        PCA(p, c) = sum(val(1:3))/sum(val(1:end));

        % plot regressors
        figure(fig2)
        subplot(length(param.peak), length(param.conditions), (p-1)*length(param.conditions) + c)
        hold on
        h1 = plot(t, squeeze(regressors{p, c}(:, 2))','g','linewidth',0.8);
        h2 = plot(t, squeeze(regressors{p, c}(:, 3))','b','linewidth',0.8);
        h3 = plot(t, squeeze(regressors{p, c}(:, 1))','r','linewidth',2.5);
        title(sprintf('%s: %s', param.peak{p}, param.conditions{c}))
        if p == 3
            xlabel('time (s)')
        end
        if c == 1
            ylabel('amplitude (\muV)')
        end
        if (p-1)*length(param.conditions) + c == length(param.peak)*length(param.conditions)
            lgd = legend([h3, h1, h2], {'waveform' 'temporal jitter' 'morphology'});
            lgd.Orientation = 'horizontal'; 
            lgd.Position = [0.5, 0.99, 0.1, 0.1];
            lgd.Position(1) = 0.5 - lgd.Position(3)/2;
            title(lgd, sprintf('subject &d - MLRd regressors'));
            lgd.Box = 'off';
        end
        hold off

        %% perform the regression on single trial data
        % prepare inputs
        data_ST = data_ST(:, x_start:x_end)';
        reg = regressors{p, c}(:, :);
        reg(:, end+1) = ones;
        for i = 1:size(data_ST,2)                                                         
            coefficients{p, c}(i, :) = (regress(data_ST(:, i), reg))';         
        end
        n_func = size(reg, 2);
        n_trials = size(data_ST, 2);

        % multiply design matrix y by parameters b for each epoch
        fmult = [];
        for i = 1:n_func          
            for j = 1:n_trials             
                fmult(:, i, j) = coefficients{p, c}(j, i).*reg(:,i);
            end
        end
        fits = squeeze(sum(fmult,2)); % sum up the contributions

        %% extract single-trial peak measures
        % identify peak window
        peak_duration = round(param.fs*param.TOI_span(p));
        peak_limits(1) = max([1 round(x_peak - peak_duration/2)]);
        peak_limits(2) = min([size(fits, 1) - 1, round(x_peak + peak_duration/2)]);

        % cycle through trials
        for i = 1:size(fits, 2)
            % identify minima and maxima within the peak window
            [tmax tmin vmax vmin] = extreme_point(fits(:,i)');
            Pmax = [tmax(intersect(find(tmax >= peak_limits(1)), find(tmax <= peak_limits(2)))); vmax(intersect(find(tmax >= peak_limits(1)),find(tmax <= peak_limits(2))))];
            Pmin = [tmin(intersect(find(tmin >= peak_limits(1)), find(tmin <= peak_limits(2)))); vmin(intersect(find(tmin >= peak_limits(1)),find(tmin <= peak_limits(2))))];
            
            % get the fit and correlation with real data
            fit{p, c}(i,:) = sum(fmult(:, 1:3, i), 2);
            corr_i = corrcoef(mean(data_ST, 2), fit{p, c}(i,:));
            fit{p, c}(i) = corr_i(1,2);

            % extract peak values
            if negative
                if corr_i(1,2) > 0
                    if length(Pmin) > 0
                        y_peak_i = min(Pmin(2,:));
                        x_peak_i = t(find(fits(:,i) == y_peak_i));
                        [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                    elseif length(Pmin) == 0
                        if length(Pmax) == 0
                            y_peak_i = NaN; x_peak_i = NaN; FWHM = NaN; slope = NaN; energy = NaN;
                        else
                            y_peak_i = max(Pmax(2,:));
                            x_peak_i = t(find(fits(:,i) == y_peak_i));
                            [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                        end
                    end 
                elseif corr_i(1,2) < 0
                    if length(Pmax) > 0
                        y_peak_i = max(Pmax(2,:));
                        x_peak_i = t(find(fits(:,i) == y_peak_i));
                        [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                    elseif length(Pmax) == 0
                        if length(Pmin) == 0
                            y_peak_i = NaN; x_peak_i = NaN; FWHM = NaN; slope = NaN; energy = NaN;
                        else
                            y_peak_i = min(Pmin(2,:));
                            x_peak_i = t(find(fits(:,i) == y_peak_i));
                            [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                        end
                    end
                end            
            else
                if corr_i(1,2) > 0
                    if length(Pmax) > 0
                        y_peak_i = max(Pmax(2,:));
                        x_peak_i = t(find(fits(:,i) == y_peak_i));
                        [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                    elseif length(Pmax)==0
                        if length(Pmin)==0
                            y_peak_i = NaN; x_peak_i = NaN; FWHM = NaN; slope = NaN; energy = NaN;
                        else
                            y_peak_i = min(Pmin(2,:));
                            x_peak_i = t(find(fits(:,i) == y_peak_i));
                            [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                        end
                    end   
                elseif corr_i(1,2) < 0
                    if length(Pmin) > 0
                        y_peak_i = min(Pmin(2,:));
                        x_peak_i = t(find(fits(:,i) == y_peak_i));
                        [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                    elseif length(Pmin) == 0
                        if length(Pmax) == 0
                            y_peak_i = NaN; x_peak_i = NaN; FWHM = NaN; slope = NaN; energy = NaN;
                        else
                            y_peak_i = max(Pmax(2,:));
                            x_peak_i = t(find(fits(:,i) == y_peak_i));
                            [FWHM, slope, energy] = width_FWHM(t, param.fs, fits(:,i), y_peak_i, x_peak_i);
                        end
                    end
                end
            end
            
            % encode to output variable
            LEP_ST.amplitude(p, c, i) = y_peak_i;  
            LEP_ST.latency(p, c, i) = x_peak_i;
            LEP_ST.FWHM(p, c, i) = FWHM;
            LEP_ST.slope(p, c, i) = slope;
            LEP_ST.energy(p, c, i) = energy; 
        end
    end
end

% save figure with regressors
saveas(fig2, sprintf('%s\\figures\\regressors_%s.png', folder.output, AP_info.single_subject(subject_idx).ID))

% endcode and save NLEP info
AP_info.single_subject(subject_idx).LEP(3).process = '3 - single-trial peak values extracted';
AP_info.single_subject(subject_idx).LEP(3).params.conditions = param.conditions;
AP_info.single_subject(subject_idx).LEP(3).params.method = 'MLR with dispersion term';
AP_info.single_subject(subject_idx).LEP(3).params.regressors = regressors;
AP_info.single_subject(subject_idx).LEP(3).params.coefficients = coefficients;
AP_info.single_subject(subject_idx).LEP(3).params.fits = fit;
AP_info.single_subject(subject_idx).LEP(3).date = sprintf('%s', date);
save(output_file, 'AP_info', '-append');

% endcode and save NLEP_measures
NLEP_measures(subject_idx).LEP_ST.conditions = param.conditions;
NLEP_measures(subject_idx).LEP_ST.amplitude = LEP_ST.amplitude;
NLEP_measures(subject_idx).LEP_ST.latency = LEP_ST.latency;
NLEP_measures(subject_idx).LEP_ST.FWHM = LEP_ST.FWHM;
NLEP_measures(subject_idx).LEP_ST.slope = LEP_ST.slope;
NLEP_measures(subject_idx).LEP_ST.energy = LEP_ST.energy;
save(output_file, 'NLEP_measures', '-append');

% save data to output structure
NLEP_data.LEP(subject_idx).ID = AP_info.single_subject(subject_idx).ID;
NLEP_data.LEP(subject_idx).conditions = param.conditions;
NLEP_data.LEP(subject_idx).unfiltered.N1 = data.N1;
NLEP_data.LEP(subject_idx).unfiltered.N2P2 = data.N2P2;
NLEP_data.LEP(subject_idx).CWT_filtered.N1 = data_filtered.N1;
NLEP_data.LEP(subject_idx).CWT_filtered.N2P2 = data_filtered.N2P2;
save(output_file, 'NLEP_data', '-append');

clear c d f i r j k p t x y file2process option lwdata eoi statement data_in data_out visual x_peak_start x_peak_end x_start x_end negative ...
    peak_sample data_template fs number_trial idx_reg T ur peak_point data_PCA data_start reg n_func n_trials fmult fits ...
    peak_duration peak_limits limits corr_i x_gaussian gaussian_filter data_avg idx_flip tmax tmin vmax vmin data_ST energy  eoi_all ...
    fig1 fig2 FWHM P_mask PCA Pmax Pmin slope template val W Y x_peak x_peak_i y_peak y_peak_i param coefficients regressors fit ...
    LEP_ST data data_filtered inverse_point lgd h1 h2 h3
fprintf('done.\n')

% ask if the subject is done
answer = questdlg(sprintf('Do you want to continue processing the data of subject %d?', subject_idx), 'Continue with the subject?', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
    case 'NO'
    	clear subject_idx answer
end

%% 13) average LEP analysis: time domain
% ----- section input -----
param.prefix_N1 = 'bl reref_Afz';   
param.EOI_N1 = {'C5' 'C6'};
param.n_files = 4;
param.peak = {'N1' 'N2' 'P2'}; 
param.processed = {'unfiltered' 'CWT_filtered'};
param.dataset = {{'N1' 'N2P2' 'N1_raw'}, {'N1' 'N2P2'}}; 
param.TOI_span = [0.14, 0.16, 0.20];
% ------------------------- 

% load output structures
if ~exist("folder")
    folder.output = uigetdir(pwd, 'Choose the OneDrive folder'); 
end
cd(folder.input)
load(output_file, 'AP_info', 'NLEP_measures', 'NLEP_data');

% add raw N1 data to the dataset, append block numbers
fprintf('adding raw N1 data to the output structure: subject')
for subject_idx = 1:length(AP_info.single_subject)
    % provide update
    fprintf(' %d ...', subject_idx)

    % identify datasets and conditions 
    file2process = dir(sprintf('%s\\%s_%s\\%s*.mat', folder.input, study, AP_info.single_subject(subject_idx).ID, param.prefix_N1));
    conditions = NLEP_data.LEP(subject_idx).conditions;
    if ~length(file2process) == param.n_files
        fprintf('ERROR: incorrect number of datasets (%d) found in the directory!\n', length(file2process));
        return;
    end

    % cycle through datasets
    NLEP_data.LEP(subject_idx).blocks{1} = []; NLEP_data.LEP(subject_idx).blocks{2} = [];
    NLEP_data.LEP(subject_idx).unfiltered.N1_raw.cond1 = []; NLEP_data.LEP(subject_idx).unfiltered.N1_raw.cond2 = [];
    load(sprintf('%s\\%slw6', file2process(1).folder, file2process(1).name(1:end-3)), '-mat')
    for d = 1:length(file2process)
        % load the dataset
        load(sprintf('%s\\%s', file2process(d).folder, file2process(d).name));

        % identify block number
        if contains(file2process(d).name, conditions{1})
            NLEP_data.LEP(subject_idx).blocks{1}(end+1 : end+size(data,1)) = str2num(file2process(d).name(end-4));
        else
            NLEP_data.LEP(subject_idx).blocks{2}(end+1 : end+size(data,1)) = str2num(file2process(d).name(end-4));
        end

        % identify EOI
        if contains(file2process(d).name, 'right')
            eoi(d) = find(strcmp({header.chanlocs.labels}, param.EOI_N1{1}));
        else
            eoi(d) = find(strcmp({header.chanlocs.labels}, param.EOI_N1{2}));
        end

        % append N1 from raw data 
        if contains(file2process(d).name, conditions{1})
            NLEP_data.LEP(subject_idx).unfiltered.N1_raw.cond1(end+1 : end+size(data, 1), :) = squeeze(data(:, eoi(d), 1, 1, 1, :)); 
        else
            NLEP_data.LEP(subject_idx).unfiltered.N1_raw.cond2(end+1 : end+size(data, 1), :) = squeeze(data(:, eoi(d), 1, 1, 1, :)); 
        end
    end
end
save(output_file, 'NLEP_data', '-append');
fprintf('\n')
fprintf('done.\n')

% calculate average values 
fprintf('calculating mean values: subject')
for subject_idx = 1:length(AP_info.single_subject)
    % provide update
    fprintf(' %d ...', subject_idx)

    % cycle through all datasets and conditions
    for c = 1:length(param.processed)
        for d = 1:length(param.dataset{c})
            for e = 1:length(NLEP_data.LEP(subject_idx).conditions)
                % subset data
                statement = sprintf('data = NLEP_data.LEP(subject_idx).%s.%s.cond%d;', param.processed{c}, param.dataset{c}{d}, e);
                eval(statement)

                % calculate average values per condition 
                average_cond(e).mean = mean(data, 1);
                average_cond(e).sd = std(data, 0, 1);
                average_cond(e).SEM = average_cond(e).sd / sqrt(size(data, 1));
                t = tinv(0.975, size(data, 1) - 1);                            
                average_cond(e).CI_upper = average_cond(e).mean + t * average_cond(e).SEM;
                average_cond(e).CI_lower = average_cond(e).mean - t * average_cond(e).SEM;

                % calculate average values per block
                for b = 1:2
                    data_block = data(NLEP_data.LEP(subject_idx).blocks{e} == b, :);
                    average_block((e-1)*2+b).mean = mean(data_block, 1);
                    average_block((e-1)*2+b).sd = std(data_block, 0, 1);
                    average_block((e-1)*2+b).SEM = average_block((e-1)*2+b).sd / sqrt(size(data_block, 1));
                    t = tinv(0.975, size(data_block, 1) - 1);                            
                    average_block((e-1)*2+b).CI_upper = average_block((e-1)*2+b).mean + t * average_block((e-1)*2+b).SEM;
                    average_block((e-1)*2+b).CI_lower = average_block((e-1)*2+b).mean - t * average_block((e-1)*2+b).SEM;
                end
            end

            % append to the data structure
            statement = sprintf('NLEP_data.LEP(subject_idx).%s.%s.average_cond = average_cond;', param.processed{c}, param.dataset{c}{d});
            eval(statement)
            statement = sprintf('NLEP_data.LEP(subject_idx).%s.%s.average_block = average_block;', param.processed{c}, param.dataset{c}{d});
            eval(statement)
        end
    end
end
save(output_file, 'NLEP_data', '-append');
fprintf('\n')
fprintf('done.\n')

% extract average peak values for all peaks
figure_counter = 1; 
for subject_idx = 7%1:length(NLEP_info.single_subject)
    % launch the output figure
    fig = figure(figure_counter);
    set(fig, 'Position', [50, 30, 750, 750])

    % extract peak values for every dataset and condition
    for p = 1:length(param.peak)
        for c = 1:length(NLEP_data.LEP(subject_idx).conditions)
            % subset data
            data = []; colour = [];
            if p == 1
                % raw data
                for b = 1:2
                    data(end+1, :) = NLEP_data.LEP(subject_idx).unfiltered.N1_raw.average_block((c-1)*2 + b).mean;
                    colour(end+1, :) = [0.6510    0.6510    0.6510];
                end

                % ICA filtered data
                for b = 1:2
                    data(end+1, :) = NLEP_data.LEP(subject_idx).unfiltered.N1.average_block((c-1)*2 + b).mean;
                    colour(end+1, :) = [0    0.4471    0.7412];
                end

                % CWT filtered data
                for b = 1:2
                    data(end+1, :) = NLEP_data.LEP(subject_idx).CWT_filtered.N1.average_block((c-1)*2 + b).mean;
                    colour(end+1, :) = [0 0 0];
                end
            else
                % raw data
                for b = 1:2
                    data(end+1, :) = NLEP_data.LEP(subject_idx).unfiltered.N2P2.average_block((c-1)*2 + b).mean;
                    colour(end+1, :) = [0    0.4471    0.7412];
                end

                % CWT filtered data
                for b = 1:2
                    data(end+1, :) = NLEP_data.LEP(subject_idx).CWT_filtered.N2P2.average_block((c-1)*2 + b).mean;
                    colour(end+1, :) = [0 0 0];
                end
            end
            
            % manually choose TOI
            data_name = sprintf('%s: %s - %s', NLEP_data.LEP(subject_idx).ID, param.peak{p}, NLEP_data.LEP(subject_idx).conditions{c});
            time_vector = header.xstart + (1:header.datasize(6))*header.xstep; 
            select_EEG(data_name, time_vector, data, 1/header.xstep, param.TOI_span(p), 0, 'colour', colour)

            % extract peak measures from all datasets 
            x_start = round((limits(1) - header.xstart)/header.xstep);
            x_end = round((limits(2) - header.xstart)/header.xstep);
            for d = 1:size(data, 1)
                if contains(param.peak{p}, 'N')
                    [y_peak(d), x_sample] = min(data(d, x_start:x_end));
                else
                    [y_peak(d), x_sample] = max(data(d, x_start:x_end));
                end
                x_peak(d) = (x_start + x_sample)*header.xstep + header.xstart;
            end
            
            % fill in NLEP_measures
            if p == 1
                % raw data
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.raw.amplitude(c, [1:2]) = y_peak([1:2]);', param.peak{p});
                eval(statement)
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.raw.latency(c, [1:2]) = x_peak([1:2]);', param.peak{p});
                eval(statement)

                % ICA filtered data
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.ICA_filtered.amplitude(c, [1:2]) = y_peak([3:4]);', param.peak{p});
                eval(statement)
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.ICA_filtered.latency(c, [1:2]) = x_peak([3:4]);', param.peak{p});
                eval(statement)

                % CWT filtered data
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.CWT_filtered.amplitude(c, [1:2]) = y_peak([5:6]);', param.peak{p});
                eval(statement)
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.CWT_filtered.latency(c, [1:2]) = x_peak([5:6]);', param.peak{p});
                eval(statement)                    
            else
                % raw data
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.raw.amplitude(c, [1:2]) = y_peak([1:2]);', param.peak{p});
                eval(statement)
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.raw.latency(c, [1:2]) = x_peak([1:2]);', param.peak{p});
                eval(statement)

                % CWT filtered data
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.CWT_filtered.amplitude(c, [1:2]) = y_peak([3:4]);', param.peak{p});
                eval(statement)
                statement = sprintf('NLEP_measures(subject_idx).LEP_avg.%s.CWT_filtered.latency(c, [1:2]) = x_peak([3:4]);', param.peak{p});
                eval(statement) 
            end

            % update NLEP_info
            statement = sprintf('NLEP_info.single_subject(subject_idx).LEP(4).params.TOI_%s(c, :) = limits;', param.peak{p});
            eval(statement)

            % update output figure
            figure(figure_counter);
            subplot(3, 2, (p-1)*2 + c)
            hold on
            for d = 1:size(data, 1)
                plot(time_vector, data(d, :), 'Color', colour(d, :), 'LineWidth', 1.2)
            end
            xlim([time_vector(1) time_vector(end)])
            y_lims = get(gca, 'ylim');
            ylim(y_lims)
            plot([0, 0], y_lims, 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2.2)
            r = rectangle('Position', [limits(1) y_lims(1) limits(2) - limits(1) y_lims(2) - y_lims(1)], ...
                'FaceColor', [0.9804    0.9804    0.5961], 'EdgeColor', 'none');
            uistack(r, 'bottom');
            set(gca, 'Layer', 'top');
            title(data_name(7:end))
            hold off
        end
    end

    % update and save output structures 
    AP_info.single_subject(subject_idx).LEP(4).process = '4 - peak measures extracted from average data';
    AP_info.single_subject(subject_idx).LEP(4).date = sprintf('%s', date);
    save(output_file, 'AP_info', 'NLEP_measures', '-append');

    % save figure, update counter
    saveas(fig, sprintf('%s\\figures\\LEP_peaks_avg_%s.png', folder.output, AP_info.single_subject(subject_idx).ID))
    figure_counter = figure_counter + 1;
end

clear param subject_idx file2process b c d e p t data header eoi average_cond data_block average_block conditions ...
    fig figure_counter statement data_name time_vector limits x_start x_end x_peak y_peak x_sample y_lims r colour

%% 14) RS-EEG analysis: spectral decomposition, sensor space
% ----- section input -----
param.prefix_data = {'icfilt ica_all chunked' 'icfilt ica_all RS'};
param.foi_limits = [5, 80];
param.eoi_visual = 'F3'; 
param.eoi_target = {'AF3' 'AFz' 'AF4' 'F3' 'F1' 'F2' 'F4'};
param.eoi_ctrl = {'PO3' 'POz' 'PO4' 'P3' 'P1' 'P2' 'P4'};
% ------------------------- 
% set directories and load info structure
if ~exist("folder")
    folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
    folder.input = uigetdir(pwd, 'Coose the input folder');             % pre-processed data
    folder.data = uigetdir(pwd, 'Choose the data folder');              % processed data
    folder.output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> figures, logfiles, output .mat file
    study = 'NLEP';
    output_file = sprintf('%s\\%s_output.mat', folder.output, study);
end
cd(folder.data)
load(output_file, 'AP_info', 'NLEP_measures', 'NLEP_data', 'slopes');

% load default header
load(sprintf('%s\\dataset_default.lw6', folder.output), '-mat');

% add FieldTrip to the top of search path
addpath(genpath([folder.toolbox '\FieldTrip']));
ft_defaults;

% extract PSD and measures
figure_counter = 1; 
for subject_idx = 1:length(AP_info.single_subject)  
    % identify RS-EEG data 
    files2process = dir(sprintf('%s\\%s_%s\\%s*.mat', folder.input, study, AP_info.single_subject(subject_idx).ID, param.prefix_data{1}));
    files2process = [files2process; dir(sprintf('%s\\%s_%s\\%s*.mat', folder.input, study, AP_info.single_subject(subject_idx).ID, param.prefix_data{2}))];

    % provide update 
    fprintf('subject %d - %d datasets found in the input directory\n', subject_idx, length(files2process))
    
    % launch the output figure
    fig = figure(figure_counter);
    set(fig, 'Position', [50, 30, 750, 900])

    % cycle through RS-EEG datasets 
    for f = 1:length(files2process)
        % provide update 
        fprintf('dataset %d/%d:\n', f, length(files2process))

        % load the data and header
        load(sprintf('%s\\%s', files2process(f).folder, files2process(f).name))
        load(sprintf('%s\\%s.lw6', files2process(f).folder, files2process(f).name(1:end-4)), '-mat')
        data = squeeze(data);

        % extract dataset name
        name_idx = strfind(files2process(f).name, 'RS');
        if name_idx > length(files2process(f).name)/2
            name = files2process(f).name(name_idx:end-4);
        else            
            name = sprintf('%s - %s', files2process(f).name(strfind(files2process(f).name, 'LEP'):end-4), files2process(f).name(name_idx+3:name_idx+7));
        end

        % compute power spectra for each trial 
        fprintf('computing single-trial power spectra ...\n')
        PSD = struct; PSD_st = struct;
        for e = 1:size(data,1)
            % create a FieldTrip data structure
            cfg = [];
            cfg.trial = {squeeze(data(e, :, :))};       
            cfg.time = {0 : header.xstep : (size(data, 3)-1)*header.xstep};      
            cfg.fsample = 1/header.xstep;
            cfg.label = {header.chanlocs.labels};  
            cfg.sampleinfo = [1 1000];
            data_trial = ft_datatype_raw(cfg);
            ft_checkdata(data_trial);

            % extract original spectra
            cfg = [];
            cfg.output = 'pow';
            cfg.foilim = param.foi_limits;  
            cfg.pad = 'nextpow2'; 
            cfg.method = 'irasa';    
            cfg.output = 'original';
            PSD(e).original = ft_freqanalysis(cfg, data_trial);

            % extract fractal specra
            cfg.output = 'fractal';
            PSD(e).fractal = ft_freqanalysis(cfg, data_trial);

            % compute oscillatory spectra
            cfg = [];
            cfg.parameter = 'powspctrm';
            cfg.operation     = 'x2-x1';
            PSD(e).oscillatory = ft_math(cfg, PSD(e).fractal, PSD(e).original);
        end
        for e = 1:length(PSD)
            PSD_st.original(e, :, :) = PSD(e).original.powspctrm(:, :);
            PSD_st.fractal(e, :, :) = PSD(e).fractal.powspctrm(:, :);
            PSD_st.oscillatory(e, :, :) = PSD(e).oscillatory.powspctrm(:, :);
        end
        PSD_freq = PSD(e).original.freq;
        
        % extract slope and offset of the aperiodic component 
        % consider only signal above 20 Hz
        for a = 1:size(PSD_st.fractal, 1)
            for b = 1:size(PSD_st.fractal, 2)
                PSD_data = squeeze(PSD_st.fractal(a, b, PSD_freq > 20));
                p = polyfit(log10(PSD_freq(PSD_freq > 20)), log10(PSD_data), 1);
                PSD_measures.slope(a, b) = p(1); 
                PSD_measures.offset(a, b) = 10^p(2); 
            end
        end
        
        % average PSD - prepare avg data
        fprintf('computing average power spectra ...\n')
        PSD = struct; PSD_avg = struct;
        cfg = [];
        cfg.trial = {squeeze(mean(data, 1))};       
        cfg.time = {0 : header.xstep : (size(data, 3)-1)*header.xstep};      
        cfg.fsample = 1/header.xstep;
        cfg.label = {header.chanlocs.labels};  
        cfg.sampleinfo = [1 1000];
        data_avg = ft_datatype_raw(cfg);
        ft_checkdata(data_avg);

        % average PSD - extract original spectra
        cfg = [];
        cfg.output = 'pow';
        cfg.foilim = param.foi_limits;  
        cfg.pad = 'nextpow2'; 
        cfg.method = 'irasa';    
        cfg.output = 'original';
        PSD.original = ft_freqanalysis(cfg, data_avg);
        PSD_avg.original = PSD.original.powspctrm;

        % average PSD - extract fractal specra
        cfg.output = 'fractal';
        PSD.fractal = ft_freqanalysis(cfg, data_avg);
        PSD_avg.fractal = PSD.fractal.powspctrm;

        % average PSD - compute oscillatory spectra
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.operation     = 'x2-x1';
        PSD.oscillatory = ft_math(cfg, PSD.fractal, PSD.original);
        PSD_avg.oscillatory = PSD.oscillatory.powspctrm;
        
        % compute average slope and offset
        for b = 1:size(PSD_avg.fractal, 1)
            PSD_data = squeeze(PSD_avg.fractal(b, PSD_freq > 20));
            p = polyfit(log10(PSD_freq(PSD_freq > 20)), log10(PSD_data), 1);
            PSD_measures.slope_avg(b) = p(1); 
            PSD_measures.offset_avg(b) = 10^p(2); 
        end

        % append data to NLEP_data and save
        fprintf('saving data ...\n')
        NLEP_data.RSEEG(subject_idx).dataset{f} = name;
        if f == 1
            NLEP_data.RSEEG(subject_idx).freq = PSD_freq;
        end
        NLEP_data.RSEEG(subject_idx).PSD_avg(f).original = PSD_avg.original;
        NLEP_data.RSEEG(subject_idx).PSD_avg(f).fractal = PSD_avg.fractal;
        NLEP_data.RSEEG(subject_idx).PSD_avg(f).oscillatory = PSD_avg.oscillatory;
        NLEP_data.RSEEG(subject_idx).PSD_st(f).original = PSD_st.original;
        NLEP_data.RSEEG(subject_idx).PSD_st(f).fractal = PSD_st.fractal;
        NLEP_data.RSEEG(subject_idx).PSD_st(f).oscillatory = PSD_st.oscillatory;
        save(output_file, 'NLEP_data', '-append');

        % append extracted measures to NLEP_measures and save
        NLEP_measures(subject_idx).RS_avg(f).dataset = name;
        NLEP_measures(subject_idx).RS_avg(f).slope = PSD_measures.slope_avg;
        NLEP_measures(subject_idx).RS_avg(f).offset = PSD_measures.offset_avg;
        NLEP_measures(subject_idx).RS_ST(f).dataset = name;
        NLEP_measures(subject_idx).RS_ST(f).slope = PSD_measures.slope;
        NLEP_measures(subject_idx).RS_ST(f).offset = PSD_measures.offset;
        save(output_file, 'NLEP_measures', '-append');

        % identify eoi for plotting
        eoi = find(strcmp({header.chanlocs.labels}, param.eoi_visual));

        % update output figure
        figure(figure_counter);
        subplot(ceil(length(files2process)/2), 2, f)
        hold on
        for e = 1:size(PSD_st.fractal, 1)
            h(e) = plot(log10(PSD_freq), log10(squeeze(PSD_st.fractal(e, eoi, :))), 'Color', [0.7451    0.9098    0.6000], 'LineWidth', 1);
        end
        h(end + 1) = plot(log10(PSD_freq), log10(squeeze(PSD_avg.fractal(eoi, :))), 'Color', [0.8314    0.2824    0.2824], 'LineWidth', 2);
        h(end + 1) = plot(log10(PSD_freq), log10(squeeze(mean(PSD_st.fractal(:, eoi, :), 1))), 'Color', [0.3373    0.6196    0.0863], 'LineWidth', 2);
        g(1) = plot(log10(PSD_freq(PSD_freq > 20)), log10(PSD_measures.offset_avg(eoi) * (PSD_freq(PSD_freq > 20).^PSD_measures.slope_avg(eoi))), ...
            'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', ':');
        g(2) = plot(log10(PSD_freq(PSD_freq > 20)), log10(mean(PSD_measures.offset(:, eoi)) * (PSD_freq(PSD_freq > 20).^mean(PSD_measures.slope(:, eoi)))), ...
            'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', ':');
        set(gca, 'Layer', 'top');
        title(name)   
        if f == length(files2process)
            lgd = legend([h(1), h(end), h(end-1)], {'single-trial data' 'mean single-trial' 'average data'});
            lgd.Location = 'northoutside';
            lgd.Orientation = 'horizontal'; 
            lgd.Position(1) = 0.5 - lgd.Position(3)/2;
            lgd.Position(2) = 0.95;
            title(lgd, sprintf('subject %d - aperiodic component', subject_idx));
            lgd.Box = 'off';
        end
        hold off        
    end

    % save output figure
    saveas(fig, sprintf('%s\\figures\\RS_aperiodic_%s.png', folder.output, AP_info.single_subject(subject_idx).ID))
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
save(output_file, 'NLEP_data', 'NLEP_measures', '-append');

% select relax and ready datasets
for subject_idx = 1:length(AP_info.single_subject)
    for c = 1:length(NLEP_measures(subject_idx).LEP_ST.conditions)
        % extract condition
        condition = NLEP_measures(subject_idx).LEP_ST.conditions{c};

        % select ready datasets
        ready = [];
        for d = 1:length(NLEP_measures(subject_idx).RS_ST)
            if contains(NLEP_measures(subject_idx).RS_ST(d).dataset, condition) && contains(NLEP_measures(subject_idx).RS_ST(d).dataset, 'ready')
                idx(d) = true;
            else
                idx(d) = false;
            end
        end
        if sum(idx) == 2
            data_ready = NLEP_measures(subject_idx).RS_ST(idx);
            for e = 1:2
                n_trials = sum(NLEP_data.LEP(subject_idx).blocks{c} == e);
                ready(end + 1 : end+n_trials, :) = data_ready(e).slope(1:n_trials, :); 
            end
        else
            fprintf('Subject %d - wrong number of ''ready'' datasets found: %d\n', subject_idx, sum(idx));
            ready = [];
        end

        % select relax datasets
        relax = []; 
        for d = 1:length(NLEP_measures(subject_idx).RS_ST)
            if contains(NLEP_measures(subject_idx).RS_ST(d).dataset, condition) && contains(NLEP_measures(subject_idx).RS_ST(d).dataset, 'relax')
                idx(d) = true;
            else
                idx(d) = false;
            end
        end
        if sum(idx) == 2
            data_relax = NLEP_measures(subject_idx).RS_ST(idx);
            for e = 1:2
                n_trials = sum(NLEP_data.LEP(subject_idx).blocks{c} == e);
                relax(end + 1 : end+n_trials, :) = data_relax(e).slope(1:n_trials, :); 
            end
        else
            fprintf('Subject %d - wrong number of ''relax'' datasets found: %d\n', subject_idx, sum(idx));
            relax = [];
        end

        % fill output variable
        slopes(subject_idx).ready{c} = ready;
        slopes(subject_idx).relax{c} = relax;
    end
end

% compute averages for target and control regions
if length(param.eoi_target) == length(param.eoi_ctrl)
    for e = 1:length(param.eoi_target)
        eoi_target(e) = find(strcmp({header.chanlocs.labels}, param.eoi_target{e}));
        eoi_ctrl(e) = find(strcmp({header.chanlocs.labels}, param.eoi_ctrl{e}));
    end
else
    fprintf('Numbers of target and control electrodes do not match!')
end
for subject_idx = 1:length(slopes)
    for c = 1:length(NLEP_data.LEP(subject_idx).conditions)
        for d = 1:size(slopes(subject_idx).ready{c}, 1)
            slopes(subject_idx).ready_target{c}(d) = mean(slopes(subject_idx).ready{c}(d, eoi_target));
            slopes(subject_idx).ready_ctrl{c}(d) = mean(slopes(subject_idx).ready{c}(d, eoi_ctrl));
            if ~isempty(slopes(subject_idx).relax{c})
                slopes(subject_idx).relax_target{c}(d) = mean(slopes(subject_idx).relax{c}(d, eoi_target));
                slopes(subject_idx).relax_ctrl{c}(d) = mean(slopes(subject_idx).relax{c}(d, eoi_ctrl));
            end
        end
    end
end

% save slopes and clean
save(output_file, 'slopes', '-append');
clear param figure_counter next_subject files2process name_idx a b c d e f g h p PSD PSD_st PSD_avg PSD_data PSD_freq answer ...
    cfg data data_avg data_ft data_trial eoi fig freq header lgd name PSD_measures ready relax subject_idx idx n_trials idx ...
    eoi_ctrl eoi_target data_ready data_relax condition slopes

%% 15) export to R for statistics and visualization
% ----- section input -----
param.LEP.dataset = {'N1' 'N2P2'};
param.LEP.processing = {'unfiltered' 'CWT_filtered'};
param.LEP.peak = {'N1' 'N2' 'P2'}; 
param.LEP.EOI = {'C5' 'C6' 'Cz'};
param.LEP.ref = {'AFz' 'avg'};
param.topo_prefix = {'bl icfilt ica_N1 reref_AFz' 'bl icfilt ica_all'};
param.RS.condition = {'relaxed' 'ready'};
param.RS.area = {'target' 'control'};
% ------------------------- 

% re-set directories and load output structures
folder.input = uigetdir(pwd, 'Choose the shared data folder'); 
folder.output = uigetdir(pwd, 'Choose the OneDrive folder'); 
study = 'NLEP';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
cd(folder.output)
% load(output_file, 'NLEP_info', 'NLEP_measures', 'NLEP_data');
load(output_file, 'AP_info', 'NLEP_measures', 'NLEP_data', 'slopes')

% check if the output tables already exists
output_variables = who('-file', output_file);
if ismember('NLEP_table4stats', output_variables)
    load(output_file, 'NLEP_table4stats');    
else
    NLEP_table4stats = table;
end
if ismember('NLEP_table4data', output_variables)
    load(output_file, 'NLEP_table4data');    
else
    NLEP_table4data = table;
end
if ismember('NLEP_table4topo', output_variables)
    load(output_file, 'NLEP_table4topo');    
else
    NLEP_table4topo = table;
end

% load default header
load('dataset_default.lw6', '-mat')

% put measured variables in a long-format output table
fprintf('exporting LEP peak measures to a table ...\n')
row_counter = height(NLEP_table4stats) + 1;
for a = 1:length(AP_info.single_subject)
    % determine contrast = subgroup  
    if contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'hand')
        contrast = 1;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'foot') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 2;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 3;
    end

    % condition = stimulated area
    for b = 1:length(NLEP_measures(a).LEP_avg.conditions)       
        % verify that data and measures match in condition sequence
        if ~strcmp(NLEP_data.LEP(a).conditions{b}, NLEP_measures(a).LEP_avg.conditions{b})
            fprintf('Condition sequences do not match - check output structures!\n')
            return
        end
        % time = block 1 / block 2
        for c = 1:2  
            % LEP peak
            for d = 1:length(param.LEP.peak)
                % single trials
                for e = find(NLEP_data.LEP(a).blocks{b} == c)
                    % subject 
                    NLEP_table4stats.subject(row_counter) = a;
                    NLEP_table4stats.ID{row_counter} = AP_info.single_subject(a).ID;
                    NLEP_table4stats.contrast(row_counter) = contrast;
                    NLEP_table4stats.age(row_counter) = AP_info.single_subject(a).age;
                    NLEP_table4stats.male(row_counter) = AP_info.single_subject(a).male;
                    NLEP_table4stats.handedness(row_counter) = AP_info.single_subject(a).handedness; 
                    if isempty(AP_info.single_subject(a).body) 
                        NLEP_table4stats.height(row_counter) = NaN;
                        NLEP_table4stats.weight(row_counter) = NaN;
                        NLEP_table4stats.limb_length(row_counter) = NaN;
                    else
                        NLEP_table4stats.height(row_counter) = AP_info.single_subject(a).body.height;
                        NLEP_table4stats.weight(row_counter) = AP_info.single_subject(a).body.weight;
                        if contains(NLEP_measures(a).LEP_avg.conditions{b}, 'hand')
                            statement = sprintf('NLEP_table4stats.limb_length(row_counter) = NLEP_info.single_subject(a).body.arm_%s;', NLEP_data.LEP(a).conditions{b}(6:end));
                            eval(statement)
                        else
                            statement = sprintf('NLEP_table4stats.limb_length(row_counter) = NLEP_info.single_subject(a).body.leg_%s;', NLEP_data.LEP(a).conditions{b}(6:end));
                            eval(statement)
                        end
                    end
                    if ismember('LEP', fieldnames(AP_info.single_subject(a).temperature))
                        statement = sprintf('NLEP_table4stats.temperature(row_counter) = NLEP_info.single_subject(a).temperature.LEP.%s_b%d;', replace(NLEP_measures(a).LEP_avg.conditions{b}, ' ', '_'), c);
                        eval(statement)
                    else
                        NLEP_table4stats.temperature(row_counter) = NaN;
                    end
    
                    % stimulation 
                    NLEP_table4stats.area{row_counter} = NLEP_data.LEP(a).conditions{b}(1:4);
                    NLEP_table4stats.side{row_counter} = NLEP_data.LEP(a).conditions{b}(6:end);
                    NLEP_table4stats.time(row_counter) = c;
    
                    % average LEP measures
                    NLEP_table4stats.peak{row_counter} = param.LEP.peak{d};
                    if d == 1
                        if contains(NLEP_data.LEP(a).conditions{b}, 'right')
                            NLEP_table4stats.EOI{row_counter} = param.LEP.EOI{1};
                        else
                            NLEP_table4stats.EOI{row_counter} = param.LEP.EOI{2};
                        end
                        NLEP_table4stats.reference{row_counter} = param.LEP.ref{1};
                    else
                        NLEP_table4stats.EOI{row_counter} = param.LEP.EOI{3};
                        NLEP_table4stats.reference{row_counter} = param.LEP.ref{2};
                    end
                    statement = sprintf('NLEP_table4stats.amplitude_avg_raw(row_counter) = NLEP_measures(a).LEP_avg.%s.raw.amplitude(b, c);', param.LEP.peak{d});
                    eval(statement)
                    statement = sprintf('NLEP_table4stats.latency_avg_raw(row_counter) = NLEP_measures(a).LEP_avg.%s.raw.latency(b, c);', param.LEP.peak{d});
                    eval(statement)
                    if d == 1
                        NLEP_table4stats.amplitude_avg_ICA(row_counter) = NLEP_measures(a).LEP_avg.N1.ICA_filtered.amplitude(b, c);
                        NLEP_table4stats.latency_avg_ICA(row_counter) = NLEP_measures(a).LEP_avg.N1.ICA_filtered.latency(b, c);
                    else
                        NLEP_table4stats.amplitude_avg_ICA(row_counter) = NaN;
                        NLEP_table4stats.latency_avg_ICA(row_counter) = NaN;
                    end
                    statement = sprintf('NLEP_table4stats.amplitude_avg_CWT(row_counter) = NLEP_measures(a).LEP_avg.%s.CWT_filtered.amplitude(b, c);', param.LEP.peak{d});
                    eval(statement)
                    statement = sprintf('NLEP_table4stats.latency_avg_CWT(row_counter) = NLEP_measures(a).LEP_avg.%s.CWT_filtered.latency(b, c);', param.LEP.peak{d});
                    eval(statement)

                    % single trial LEP measures
                    NLEP_table4stats.trial(row_counter) = e;
                    NLEP_table4stats.amplitude_ST(row_counter) = NLEP_measures(a).LEP_ST.amplitude(d, c, e);
                    NLEP_table4stats.latency_ST(row_counter) = NLEP_measures(a).LEP_ST.latency(d, c, e);

                    % single-trial RS-EEg measures
                    NLEP_table4stats.ready_target(row_counter) = slopes(a).ready_target{b}(e);
                    NLEP_table4stats.ready_ctrl(row_counter) = slopes(a).ready_ctrl{b}(e);
                    if ~isempty(slopes(a).relax_target)
                        NLEP_table4stats.relax_target(row_counter) = slopes(a).relax_target{b}(e);
                        NLEP_table4stats.relax_ctrl(row_counter) = slopes(a).relax_ctrl{b}(e);
                    else
                        NLEP_table4stats.relax_target(row_counter) = NaN;
                        NLEP_table4stats.relax_ctrl(row_counter) = NaN;
                    end

                    % evoked pain
                    dataset_pain = find(strcmp(sprintf('%s b%d', NLEP_data.LEP(a).conditions{b}, c), NLEP_measures(a).conditions));
                    trials_b1 = max(find(NLEP_data.LEP(a).blocks{b} == 1));
                    if e > trials_b1
                        trial_pain = e - trials_b1;
                    else
                        trial_pain = e;
                    end
                    NLEP_table4stats.pain(row_counter) = NLEP_measures(a).pain(dataset_pain, trial_pain);

                    % update row counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
fprintf('done.\n')
save(output_file, 'NLEP_table4stats', '-append');

%  add
fprintf('exporting LEP peak measures to a table ...\n')
row_counter = height(NLEP_table4stats) + 1;
for a = 1:length(AP_info.single_subject)
    % determine contrast = subgroup  
    if contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'hand')
        contrast = 1;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'foot') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 2;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 3;
    end

    % condition = stimulated area
    for b = 1:length(NLEP_measures(a).LEP_avg.conditions)       
        % verify that data and measures match in condition sequence
        if ~strcmp(NLEP_data.LEP(a).conditions{b}, NLEP_measures(a).LEP_avg.conditions{b})
            fprintf('Condition sequences do not match - check output structures!\n')
            return
        end
        % time = block 1 / block 2
        for c = 1:2  
            % LEP peak
            for d = 1:length(param.LEP.peak)
                % single trials
                for e = find(NLEP_data.LEP(a).blocks{b} == c)
                    % subject 
                    NLEP_table4stats.subject(row_counter) = a;
                    NLEP_table4stats.ID{row_counter} = AP_info.single_subject(a).ID;
                    NLEP_table4stats.contrast(row_counter) = contrast;
                    NLEP_table4stats.age(row_counter) = AP_info.single_subject(a).age;
                    NLEP_table4stats.male(row_counter) = AP_info.single_subject(a).male;
                    NLEP_table4stats.handedness(row_counter) = AP_info.single_subject(a).handedness; 
                    if isempty(AP_info.single_subject(a).body) 
                        NLEP_table4stats.height(row_counter) = NaN;
                        NLEP_table4stats.weight(row_counter) = NaN;
                        NLEP_table4stats.limb_length(row_counter) = NaN;
                    else
                        NLEP_table4stats.height(row_counter) = AP_info.single_subject(a).body.height;
                        NLEP_table4stats.weight(row_counter) = AP_info.single_subject(a).body.weight;
                        if contains(NLEP_measures(a).LEP_avg.conditions{b}, 'hand')
                            statement = sprintf('NLEP_table4stats.limb_length(row_counter) = NLEP_info.single_subject(a).body.arm_%s;', NLEP_data.LEP(a).conditions{b}(6:end));
                            eval(statement)
                        else
                            statement = sprintf('NLEP_table4stats.limb_length(row_counter) = NLEP_info.single_subject(a).body.leg_%s;', NLEP_data.LEP(a).conditions{b}(6:end));
                            eval(statement)
                        end
                    end
                    if ismember('LEP', fieldnames(AP_info.single_subject(a).temperature))
                        statement = sprintf('NLEP_table4stats.temperature(row_counter) = NLEP_info.single_subject(a).temperature.LEP.%s_b%d;', replace(NLEP_measures(a).LEP_avg.conditions{b}, ' ', '_'), c);
                        eval(statement)
                    else
                        NLEP_table4stats.temperature(row_counter) = NaN;
                    end
    
                    % stimulation 
                    NLEP_table4stats.area{row_counter} = NLEP_data.LEP(a).conditions{b}(1:4);
                    NLEP_table4stats.side{row_counter} = NLEP_data.LEP(a).conditions{b}(6:end);
                    NLEP_table4stats.time(row_counter) = c;
    
                    % average LEP measures
                    NLEP_table4stats.peak{row_counter} = param.LEP.peak{d};
                    if d == 1
                        if contains(NLEP_data.LEP(a).conditions{b}, 'right')
                            NLEP_table4stats.EOI{row_counter} = param.LEP.EOI{1};
                        else
                            NLEP_table4stats.EOI{row_counter} = param.LEP.EOI{2};
                        end
                        NLEP_table4stats.reference{row_counter} = param.LEP.ref{1};
                    else
                        NLEP_table4stats.EOI{row_counter} = param.LEP.EOI{3};
                        NLEP_table4stats.reference{row_counter} = param.LEP.ref{2};
                    end
                    statement = sprintf('NLEP_table4stats.amplitude_avg_raw(row_counter) = NLEP_measures(a).LEP_avg.%s.raw.amplitude(b, c);', param.LEP.peak{d});
                    eval(statement)
                    statement = sprintf('NLEP_table4stats.latency_avg_raw(row_counter) = NLEP_measures(a).LEP_avg.%s.raw.latency(b, c);', param.LEP.peak{d});
                    eval(statement)
                    if d == 1
                        NLEP_table4stats.amplitude_avg_ICA(row_counter) = NLEP_measures(a).LEP_avg.N1.ICA_filtered.amplitude(b, c);
                        NLEP_table4stats.latency_avg_ICA(row_counter) = NLEP_measures(a).LEP_avg.N1.ICA_filtered.latency(b, c);
                    else
                        NLEP_table4stats.amplitude_avg_ICA(row_counter) = NaN;
                        NLEP_table4stats.latency_avg_ICA(row_counter) = NaN;
                    end
                    statement = sprintf('NLEP_table4stats.amplitude_avg_CWT(row_counter) = NLEP_measures(a).LEP_avg.%s.CWT_filtered.amplitude(b, c);', param.LEP.peak{d});
                    eval(statement)
                    statement = sprintf('NLEP_table4stats.latency_avg_CWT(row_counter) = NLEP_measures(a).LEP_avg.%s.CWT_filtered.latency(b, c);', param.LEP.peak{d});
                    eval(statement)

                    % single trial LEP measures
                    NLEP_table4stats.trial(row_counter) = e;
                    NLEP_table4stats.amplitude_ST(row_counter) = NLEP_measures(a).LEP_ST.amplitude(d, c, e);
                    NLEP_table4stats.latency_ST(row_counter) = NLEP_measures(a).LEP_ST.latency(d, c, e);

                    % single-trial RS-EEg measures
                    NLEP_table4stats.ready_target(row_counter) = slopes(a).ready_target{b}(e);
                    NLEP_table4stats.ready_ctrl(row_counter) = slopes(a).ready_ctrl{b}(e);
                    if ~isempty(slopes(a).relax_target)
                        NLEP_table4stats.relax_target(row_counter) = slopes(a).relax_target{b}(e);
                        NLEP_table4stats.relax_ctrl(row_counter) = slopes(a).relax_ctrl{b}(e);
                    else
                        NLEP_table4stats.relax_target(row_counter) = NaN;
                        NLEP_table4stats.relax_ctrl(row_counter) = NaN;
                    end

                    % evoked pain
                    dataset_pain = find(strcmp(sprintf('%s b%d', NLEP_data.LEP(a).conditions{b}, c), NLEP_measures(a).conditions));
                    trials_b1 = max(find(NLEP_data.LEP(a).blocks{b} == 1));
                    if e > trials_b1
                        trial_pain = e - trials_b1;
                    else
                        trial_pain = e;
                    end
                    NLEP_table4stats.pain(row_counter) = NLEP_measures(a).pain(dataset_pain, trial_pain);

                    % update row counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
fprintf('done.\n')
save(output_file, 'NLEP_table4stats', '-append');

% export table to .csv for R 
writetable(NLEP_table4stats, 'NLEP_LEP_measures.csv');

% transform average EEG signals into a table
row_counter = height(NLEP_table4data) + 1;
fprintf('exporting average LEP signals: ')
for a = 1:length(AP_info.single_subject)
    % update
    fprintf('%d .. ', a)

    % determine contrast = subgroup  
    if contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'hand')
        contrast = 1;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'foot') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 2;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 3;
    end

    % condition = stimulated area
    for b = 1:length(NLEP_measures(a).LEP_avg.conditions)       
        % verify that data and measures match in condition sequence
        if ~strcmp(NLEP_data.LEP(a).conditions{b}, NLEP_measures(a).LEP_avg.conditions{b})
            fprintf('Condition sequences do not match - check output structures!\n')
            return
        end
        % time = block 1 / block 2
        for c = 1:2  
            % dataset 
            for d = 1:length(param.dataset)
                % determine electrode and reference
                if d == 1
                    if contains(NLEP_data.LEP(a).conditions{b}, 'right')
                        eoi = param.EOI{1};
                    else
                        eoi = param.EOI{2};
                    end
                    ref = param.ref{1};
                else
                    eoi = param.EOI{3};
                    ref = param.ref{2};
                end
                % processing
                for e = 1:length(param.processing)
                    % subject info
                    NLEP_table4data.subject(row_counter) = a;
                    NLEP_table4data.ID{row_counter} = AP_info.single_subject(a).ID;
                    NLEP_table4data.contrast(row_counter) = contrast;

                    % dataset info
                    NLEP_table4data.response{row_counter} = 'LEP';
                    NLEP_table4data.area{row_counter} = NLEP_data.LEP(a).conditions{b}(1:4);
                    NLEP_table4data.side{row_counter} = NLEP_data.LEP(a).conditions{b}(6:end);
                    NLEP_table4data.time(row_counter) = c;
                    NLEP_table4data.dataset{row_counter} = param.dataset{d};
                    NLEP_table4data.processing{row_counter} = param.processing{e};

                    % average sampled voltage
                    for f = 1:size(NLEP_data.LEP(a).unfiltered.N1.cond1, 2)
                        statement = sprintf('NLEP_table4data.v%d(row_counter) = NLEP_data.LEP(a).%s.%s.average_block((b-1)*2 + c).mean(f);', f, param.processing{e}, param.dataset{d});   
                        eval(statement)
                    end

                    % update row counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
fprintf('\n')
fprintf('done.\n')
save(output_file, 'NLEP_table4data', '-append');

% export for R
writetable(NLEP_table4data, 'NLEP_LEP_data.csv');

% preapre a dataset with peak topoplots
param.processing{1} = 'raw';
row_counter = height(NLEP_table4topo) + 1;
fprintf('exporting average LEP peak topographies: ')
for a = 1:length(AP_info.single_subject)
    % update
    fprintf('%d .. ', a)

    % determine contrast = subgroup  
    if contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'hand')
        contrast = 1;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'foot') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 2;
    elseif contains(NLEP_data.LEP(a).conditions{1}, 'hand') && contains(NLEP_data.LEP(a).conditions{2}, 'foot')
        contrast = 3;
    end

    % condition = stimulated area
    for b = 1:length(NLEP_measures(a).LEP_avg.conditions)       
        % verify that data and measures match in condition sequence
        if ~strcmp(NLEP_data.LEP(a).conditions{b}, NLEP_measures(a).LEP_avg.conditions{b})
            fprintf('Condition sequences do not match - check output structures!\n')
            return
        end
        % time = block 1 / block 2
        for c = 1:2  
            % peak 
            for d = 1:length(param.peak)
                % determine electrode and reference
                if d == 1
                    if contains(NLEP_data.LEP(a).conditions{b}, 'right')
                        eoi = param.EOI{1};
                    else
                        eoi = param.EOI{2};
                    end
                    ref = param.ref{1};
                else
                    eoi = param.EOI{3};
                    ref = param.ref{2};
                end
                % processing
                for e = 1:length(param.processing)
                    % subject info
                    NLEP_table4topo.subject(row_counter) = a;
                    NLEP_table4topo.ID{row_counter} = AP_info.single_subject(a).ID;
                    NLEP_table4topo.contrast(row_counter) = contrast;

                    % dataset info
                    NLEP_table4topo.response{row_counter} = 'LEP';
                    NLEP_table4topo.area{row_counter} = NLEP_data.LEP(a).conditions{b}(1:4);
                    NLEP_table4topo.side{row_counter} = NLEP_data.LEP(a).conditions{b}(6:end);
                    NLEP_table4topo.time(row_counter) = c;
                    NLEP_table4topo.peak{row_counter} = param.peak{d};
                    NLEP_table4topo.processing{row_counter} = param.processing{e};
                    statement = sprintf('peak_latency = NLEP_measures(a).LEP_avg.%s.%s.latency((b-1)*2 + c);', param.peak{d}, param.processing{e});
                    eval(statement)
                    NLEP_table4topo.latency{row_counter} = peak_latency;

                    % load the data
                    if  d == 1
                        file2load = dir(sprintf('%s\\%s_%s\\%s*LEP*%s b%d*', folder.input, study, NLEP_measures(a).ID, param.prefix_topo{1}, NLEP_data.LEP(a).conditions{b}, c));
                    else
                        file2load = dir(sprintf('%s\\%s_%s\\%s*LEP*%s b%d*', folder.input, study, NLEP_measures(a).ID, param.prefix_topo{2}, NLEP_data.LEP(a).conditions{b}, c));
                    end
                    if length(file2load) == 2
                        load(sprintf('%s\\%s', file2load(1).folder, file2load(1).name), '-mat')
                        load(sprintf('%s\\%s', file2load(2).folder, file2load(2).name))
                    else
                        fprintf('Wrong number of files were found in the input directory: &d', length(file2load))
                        return
                    end
                    
                    % extract data at peak timepoint
                    peak_sample = round((peak_latency - header.xstart)/header.xstep);
                    peak_data = mean(data(:, :, 1, 1, 1, peak_sample), 1);

                    % extract channel labels
                    labels = {header.chanlocs.labels};

                    % append to the table
                    for f = 1:size(peak_data, 2)
                        statement = sprintf('NLEP_table4topo.%s(row_counter) = peak_data(f);', labels{f});   
                        eval(statement)
                    end

                    % update row counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
fprintf('\n')
fprintf('done.\n')
save(output_file, 'NLEP_table4topo', '-append');

% export for R
writetable(NLEP_table4topo, 'NLEP_LEP_topo.csv');

clear param output_variables a b c d e f row_counter statement dataset_pain trials_b1 trial_pain ref eoi ...
    header data peak_latency peak_sample peak_data labels file2load contrast

%% 16) ERP group-average visualization
% to fix, bugging!
% ----- section input -----
param.response = {'LEP' 'SEP'};
param.contrast = {'hand_hand' 'foot_foot' 'hand_foot'};
param.EOI = 'Cz';
param.window.LEP = [-0.2, 0.6];
param.window.SEP = [-0.1, 0.4];
param.ylim.LEP = [-15 15];
param.ylim.SEP = [-5 5];
param.alpha = 0.2;
param.peaks = {'N1' 'N2' 'P2'};
figure_counter = 1;
% ------------------------- 

% update output variables
load(output_file, 'AP_info', 'NLEP_data', 'NLEP_measures', 'data_area');

% load default header
load('dataset_default.lw6', '-mat')

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% --- plot mean N1 and N2P2 signal from the hand area ---
% load hand data
data_N1 = []; data_N2P2 = [];
for subject_idx = 1:length(AP_info.single_subject)
    for c = 1:length(NLEP_data.LEP(subject_idx).conditions)
        if contains(NLEP_data.LEP(subject_idx).conditions{c}, 'hand')
            data_N1(end+1, :) = NLEP_data.LEP(subject_idx).unfiltered.N1.average_cond(c).mean;  
            data_N2P2(end+1, :) = NLEP_data.LEP(subject_idx).unfiltered.N2P2.average_cond(c).mean; 
        end
    end
end

% plot unfiltered N1 signal
avg_N1.mean = mean(data_N1, 1);
avg_N1.sd = std(data_N1, 0, 1);
avg_N1.sem = avg_N1.sd / sqrt(size(data_N1, 1));
t = tinv(0.975, size(data_N1, 1) - 1);                            
avg_N1.CI_upper = avg_N1.mean + t * avg_N1.sem;
avg_N1.CI_lower = avg_N1.mean - t * avg_N1.sem;
x = (1 : length(avg_N1.mean))*header.xstep + header.xstart;
fig = figure(figure_counter);
plot_ERP(avg_N1.mean, avg_N1.CI_upper, avg_N1.CI_lower, x, 'shading', 'on', 'colours', [0    0.4471    0.7412])
saveas(fig, sprintf('%s\\figures\\LEP_hand_N1.svg', folder.output))
figure_counter = figure_counter + 1;

% plot unfiltered N2P2 signal
avg_N2P2.mean = mean(data_N2P2, 1);
avg_N2P2.sd = std(data_N2P2, 0, 1);
avg_N2P2.sem = avg_N2P2.sd / sqrt(size(data_N2P2, 1));
t = tinv(0.975, size(data_N2P2, 1) - 1);                            
avg_N2P2.CI_upper = avg_N2P2.mean + t * avg_N2P2.sem;
avg_N2P2.CI_lower = avg_N2P2.mean - t * avg_N2P2.sem;
x = (1 : length(avg_N2P2.mean))*header.xstep + header.xstart;
fig = figure(figure_counter);
plot_ERP(avg_N2P2.mean, avg_N2P2.CI_upper, avg_N2P2.CI_lower, x, 'shading', 'on', 'colours', [0.4667    0.6745    0.1882])
saveas(fig, sprintf('%s\\figures\\LEP_hand_N2P2.svg', folder.output))
figure_counter = figure_counter + 1;

% --- plot peak topographies from the hand area ---
% extract electrode labels
for i = 1:length(header.chanlocs)
    labels{i} = header.chanlocs(i).labels;
end   
labels_flipped = labels;
for i = 1:length(labels)
    electrode_n = str2num(labels{i}(end));
    if isempty(electrode_n)
    else
        if mod(electrode_n,2) == 1              % odd number --> left hemisphere                    
            label_new = labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n + 1)];
            a = find(strcmpi(labels, label_new));
            if isempty(a)
            else
                labels_flipped{i} = label_new;
            end
        else                                    % even number --> right hemisphere 
            label_new = labels{i}(1:end-1);
            label_new = [label_new num2str(electrode_n - 1)];
            a = find(strcmpi(labels, label_new));
            if isempty(a)
            else
                labels_flipped{i} = label_new;
            end
        end
    end
end
labels_dict = cat(1, labels, labels_flipped)';

% identify electrode positions for flipping
electrode_idx = 1:length(labels);
electrode_idx_flipped = electrode_idx;
for c = 1:length(labels)
    idx1 = find(strcmpi(labels ,labels_dict{c,1}));
    idx2 = find(strcmpi(labels,labels_dict{c, 2}));
    if isempty(idx1)
    else
        if isempty(idx2)
        else
            electrode_idx_flipped(idx1) = electrode_idx(idx2);
        end
    end
end

% plot N1 topography
[~, peak_N1] = min(avg_N1.mean);
topo_N1 = data_area(3).mean(:, peak_N1);
fig = figure(figure_counter);
topoplot(topo_N1',header.chanlocs, 'maplimits', [-4 4], 'shading', 'interp', 'whitebk', 'on');
colorbar;
set(gcf,'color',[1 1 1]);
saveas(fig, sprintf('%s\\figures\\LEP_topo_hand_N1.svg', folder.output))
figure_counter = figure_counter + 1;
        
% plot N2 topography
[~, peak_N2] = min(avg_N2P2.mean);
topo_N2 = data_area(1).mean(:, peak_N2);
fig = figure(figure_counter);
topoplot(topo_N2',header.chanlocs, 'maplimits', [-10 10], 'shading', 'interp', 'whitebk', 'on');
colorbar;
set(gcf,'color',[1 1 1]);
saveas(fig, sprintf('%s\\figures\\LEP_topo_hand_N2.svg', folder.output))
figure_counter = figure_counter + 1;

% plot P2 topography
[~, peak_P2] = max(avg_N2P2.mean);
topo_P2 = data_area(1).mean(:, peak_P2);
fig = figure(figure_counter);
topoplot(topo_P2',header.chanlocs, 'maplimits', [-10 10], 'shading', 'interp', 'whitebk', 'on');
colorbar;
set(gcf,'color',[1 1 1]);
saveas(fig, sprintf('%s\\figures\\LEP_topo_hand_P2.svg', folder.output))
figure_counter = figure_counter + 1;

% split according to conditions 
data_area(1).condition = 'LEP N2P2 hand'; data_area(1).data = [];
data_area(2).condition = 'LEP N2P2 foot'; data_area(2).data = [];
data_area(3).condition = 'LEP N1 hand'; data_area(3).data = [];
data_area(4).condition = 'LEP N1 foot'; data_area(4).data = [];
data_area(5).condition = 'SEP hand'; data_area(5).data = [];
data_area(6).condition = 'SEP foot'; data_area(6).data = [];
for d = 1:length(dataset)
    % load the data
    data = dataset(d).data;

    % flip electrodes if stimulated on the left
    if strcmp(dataset(d).side, 'left')
        data = data(electrode_idx_flipped, :);
    end

    % encode
    if d == 1
        AP_info.group_analysis.preliminary(4).process = sprintf('4 - GFP peak localization');
        for a = 1:length(data_area)
            AP_info.group_analysis.preliminary(4).params.datasets{a} = data_area(a).condition;
        end
        AP_info.group_analysis.preliminary(4).params.labels = labels;
        AP_info.group_analysis.preliminary(4).params.flip_idx = electrode_idx_flipped;
        AP_info.group_analysis.preliminary(4).date = sprintf('%s', date);
    end

    % append to the right datset
    if strcmp(dataset(d).response, 'LEP') && strcmp(dataset(d).peak, 'N2P2')
        switch dataset(d).area
            case 'hand'
                data_area(1).data(end+1, :, :) = data;
            case 'foot'
                data_area(2).data(end+1, :, :) = data;
        end
    elseif strcmp(dataset(d).response, 'LEP') && strcmp(dataset(d).peak, 'N1')
        switch dataset(d).area
            case 'hand'
                data_area(3).data(end+1, :, :) = data;
            case 'foot'
                data_area(4).data(end+1, :, :) = data;
        end
    elseif strcmp(dataset(d).response, 'SEP') 
        switch dataset(d).area
            case 'hand'
                data_area(5).data(end+1, :, :) = data;
            case 'foot'
                data_area(6).data(end+1, :, :) = data;
        end
    end
end

% calculate mean, variance measures and gfp
for a = 1:length(data_area)
    data_area(a).mean = squeeze(mean(data_area(a).data, 1));
    data_area(a).SD = squeeze(std(data_area(a).data, 0, 1));
    data_area(a).SEM = data_area(a).SD / sqrt(size(data_area(a).data, 1));
    t = tinv(0.975, sqrt(size(data_area(a).data, 1)) - 1);                            
    data_area(a).CI_upper = data_area(a).mean + t * data_area(a).SEM;
    data_area(a).CI_lower = data_area(a).mean - t * data_area(a).SEM;
    data_area(a).GFP = std(data_area(a).mean, 1);
end

% plot GFP, identify peaks
x = -dataset(1).header.xstart*1000:1:699;
for a = 1:length(data_area)
    % launch the figure
    fig = figure;
    hold on
    
    % extract peak latencies
    h_axis(1) = subplot(3, param.gfp.max_peaks, [1 : 2*param.gfp.max_peaks]);    
    [data_area(a).gfp_peaks.latency, data_area(a).gfp_peaks.amplitude] = gfp_plot(x, data_area(a).GFP, 1, param.gfp.labeled, 'max_peaks', param.gfp.max_peaks);
    title(sprintf('grand average GFP: %s', data_area(a).condition), 'fontsize', 16, 'fontweight', 'bold')

    % choose data for topoplots 
    for e = 1:size(data_area(a).mean, 1)
        for i = 1:size(data_area(a).mean, 2)
            data_topoplot(1, e, 1, 1, 1, i) = data_area(a).mean(e, i);
        end
    end

    % add topoplots
    for k = 1:length(data_area(a).gfp_peaks.latency)
        % plot the topoplot
        h_axis(1 + k) = subplot(3, param.gfp.max_peaks, 2*param.gfp.max_peaks + k);
        topo_plot(dataset(1).header, data_topoplot, data_area(a).gfp_peaks.latency(k), [-3, 3])
    
        % shift down
        pos = get(h_axis(1 + k), 'Position');
        pos(2) = pos(2) - 0.05;
        set(h_axis(1 + k), 'Position', pos);
    
        % add timing
        text(-0.3, -0.8, sprintf('%1.0f ms', data_area(a).gfp_peaks.latency(k)), 'Color', [1 0 0], 'FontSize', 14)
    end
    hold off

    % save figure
    saveas(fig, sprintf('%s\\figures\\GFP %s.png', folder.output, data_area(a).condition))
end

% identify EOI
eoi = find(strcmp(AP_info.group_analysis.preliminary(4).params.labels, param.EOI));

% plot average ERPs
for a = 1:length(param.response)
    % determine axes properties    
    switch param.response{a}
        case 'LEP'
            x_start = round((param.window.LEP(1) - dataset(1).header.xstart)/dataset(1).header.xstep);
            x_end = round((param.window.LEP(2) - dataset(1).header.xstart)/dataset(1).header.xstep);
            x = (param.window.LEP(1) + (0:x_end - x_start) * dataset(1).header.xstep)*1000;
            y_limits = param.ylim.LEP;
        case 'SEP'
            x_start = round((param.window.SEP(1) - dataset(1).header.xstart)/dataset(1).header.xstep);
            x_end = round((param.window.SEP(2) - dataset(1).header.xstart)/dataset(1).header.xstep);
            x = (param.window.SEP(1) + (0:x_end - x_start) * dataset(1).header.xstep)*1000;
            y_limits = param.ylim.SEP;
    end   

    % loop through contrasts
    for b = 1:length(param.contrast)
        % get the data
        data = [];
        switch param.contrast{b}
            case 'hand_hand'
                % define plotting parameters
                conditions = {'right hand' 'left hand'};
                colours = [0.9216    0.4157    0.5059; 0.3373    0.6196    0.0863];

                % select the data
                statement = sprintf('data(1,:,:,:) = data_contrast.%s.%s.right;', param.response{a}, param.contrast{b});
                eval(statement)
                statement = sprintf('data(2,:,:,:) = data_contrast.%s.%s.left;', param.response{a}, param.contrast{b});
                eval(statement)

            case 'foot_foot'
                % define plotting parameters
                conditions = {'right foot' 'left foot'};
                colours = [0.9216    0.4157    0.5059; 0.3373    0.6196    0.0863];

                % select the data
                statement = sprintf('data(1,:,:,:) = data_contrast.%s.%s.right;', param.response{a}, param.contrast{b});
                eval(statement)
                statement = sprintf('data(2,:,:,:) = data_contrast.%s.%s.left;', param.response{a}, param.contrast{b});
                eval(statement)

            case 'hand_foot'
                % define plotting parameters
                conditions = {'hand' 'foot'};
                colours = [0.8314    0.2824    0.2824; 0.1843    0.6627    0.8706];

                % select the data
                statement = sprintf('data(1,:,:,:) = data_contrast.%s.%s.hand;', param.response{a}, param.contrast{b});
                eval(statement)
                statement = sprintf('data(2,:,:,:) = data_contrast.%s.%s.foot;', param.response{a}, param.contrast{b});
                eval(statement)
        end
        
        % calculate mean and CI
        t = tinv(0.975, size(data, 2) - 1); 
        visual_data = []; visual_sem = []; visual_CI_upper = []; visual_CI_lower = [];
        for c = 1:size(data, 1)
            visual_data(c,:) = mean(squeeze(data(c, :, eoi, x_start:x_end)), 1);
            visual_sem = std(squeeze(data(c, :, eoi, x_start:x_end)), 0, 1) / sqrt(size(data, 2));
            visual_CI_upper(c,:) = visual_data(c,:) + t * visual_sem;
            visual_CI_lower(c,:) = visual_data(c,:) - t * visual_sem;
        end

        % plot and save
        fig = figure('Position', [100, 100, 700, 550]);
        plot_ERP(visual_data, visual_CI_upper, visual_CI_lower, x, 'ylim', y_limits, 'labels', conditions, 'colours', colours, 'alpha', param.alpha);
        title(sprintf('grand average %s: %s vs. %s', param.response{a}, conditions{1}, conditions{2}), 'fontsize', 16, 'fontweight', 'bold')
        saveas(fig, sprintf('%s\\figures\\%s %s.png', folder.output, param.response{a}, param.contrast{b}))
    end
end

% add N1 peak values to output variable
values_contrast.N1.raw.hand_hand.right.amplitude = []; values_contrast.N1.raw.hand_hand.left.amplitude = [];
values_contrast.N1.raw.hand_hand.right.latency = []; values_contrast.N1.raw.hand_hand.left.latency = [];
values_contrast.N1.raw.foot_foot.right.amplitude = []; values_contrast.N1.raw.foot_foot.left.amplitude = [];
values_contrast.N1.raw.foot_foot.right.latency = []; values_contrast.N1.raw.foot_foot.left.latency = [];
values_contrast.N1.raw.hand_foot.hand.amplitude = []; values_contrast.N1.raw.hand_foot.foot.amplitude = [];
values_contrast.N1.raw.hand_foot.hand.latency = []; values_contrast.N1.raw.hand_foot.foot.latency = [];
values_contrast.N1.ICA_filtered.hand_hand.right.amplitude = []; values_contrast.N1.ICA_filtered.hand_hand.left.amplitude = [];
values_contrast.N1.ICA_filtered.hand_hand.right.latency = []; values_contrast.N1.ICA_filtered.hand_hand.left.latency = [];
values_contrast.N1.ICA_filtered.foot_foot.right.amplitude = []; values_contrast.N1.ICA_filtered.foot_foot.left.amplitude = [];
values_contrast.N1.ICA_filtered.foot_foot.right.latency = []; values_contrast.N1.ICA_filtered.foot_foot.left.latency = [];
values_contrast.N1.ICA_filtered.hand_foot.hand.amplitude = []; values_contrast.N1.ICA_filtered.hand_foot.foot.amplitude = [];
values_contrast.N1.ICA_filtered.hand_foot.hand.latency = []; values_contrast.N1.ICA_filtered.hand_foot.foot.latency = [];
for s = 1:length(AP_info.single_subject)
    % identify contrast and condition
    for c = 1:length(LEP_average(s).conditions)
        if contains(LEP_average(s).conditions{c}, 'hand')
            idx_contrast(c) = true;            
        else
            idx_contrast(c) = false;
        end
    end
    switch sum(idx_contrast)
        case 4
            contrast = 'hand_hand';
            condition = {'right' 'left'};
        case 2
            contrast = 'hand_foot';
            condition = {'hand' 'foot'};
        case 0
            contrast = 'foot_foot';
            condition = {'right' 'left'};
    end

    % join N1 measures to the output variable
    for c = 1:length(LEP_average(s).conditions)
        % split according to condition
        if contains(LEP_average(s).conditions{c}, condition{1})
            % raw amplitude
            statement = sprintf('values_contrast.N1.raw.%s.%s.amplitude(end+1) = LEP_average(s).N1.raw.amplitude(c);', contrast, condition{1});
            eval(statement)

            % raw latency
            statement = sprintf('values_contrast.N1.raw.%s.%s.latency(end+1) = LEP_average(s).N1.raw.latency(c);', contrast, condition{1});
            eval(statement)

            % ICA_filtered amplitude
            statement = sprintf('values_contrast.N1.ICA_filtered.%s.%s.amplitude(end+1) = LEP_average(s).N1.filtered.amplitude(c);', contrast, condition{1});
            eval(statement)

            % ICA_filtered latency
            statement = sprintf('values_contrast.N1.ICA_filtered.%s.%s.latency(end+1) = LEP_average(s).N1.filtered.latency(c);', contrast, condition{1});
            eval(statement)
        elseif contains(LEP_average(s).conditions{c}, condition{2})
            % raw amplitude
            statement = sprintf('values_contrast.N1.raw.%s.%s.amplitude(end+1) = LEP_average(s).N1.raw.amplitude(c);', contrast, condition{2});
            eval(statement)

            % raw latency
            statement = sprintf('values_contrast.N1.raw.%s.%s.latency(end+1) = LEP_average(s).N1.raw.latency(c);', contrast, condition{2});
            eval(statement)

            % ICA_filtered amplitude
            statement = sprintf('values_contrast.N1.ICA_filtered.%s.%s.amplitude(end+1) = LEP_average(s).N1.filtered.amplitude(c);', contrast, condition{2});
            eval(statement)

            % ICA_filtered latency
            statement = sprintf('values_contrast.N1.ICA_filtered.%s.%s.latency(end+1) = LEP_average(s).N1.filtered.latency(c);', contrast, condition{2});
            eval(statement)
        end
    end
end
save(output_file, 'values_contrast', '-append');
param.N1_cond = {'raw' 'ICA_filtered'};

% plot LEP peak values - bugging!
for c = 1:length(param.contrast) 
    % plotting parameters  
    switch param.contrast{c}
        case 'hand_hand'
            colours = [0.9216    0.4157    0.5059; 0.3373    0.6196    0.0863]; 
            labels = {'right hand' 'left hand'};
            n_subjects = length(values_contrast.N2.hand_hand.right.amplitude);
        case 'foot_foot'
            colours = [0.9216    0.4157    0.5059; 0.3373    0.6196    0.0863]; 
            labels = {'right foot' 'left foot'};
            n_subjects = length(values_contrast.N2.foot_foot.right.amplitude);
        case 'hand_foot'  
            colours = [0.8314    0.2824    0.2824; 0.1843    0.6627    0.8706];
            labels = {'hand' 'foot'};
            n_subjects = length(values_contrast.N2.hand_foot.hand.amplitude);
    end
    statement = sprintf('LEP_stats.%s.n_subjects = n_subjects;', param.contrast{c});
    eval(statement)
    clear data_amplitude data_latency

    for p = 1:length(param.peaks)
        if p == 1
            for n = 1:length(param.N1_cond)
                % select the data
                statement = sprintf('data_amplitude(1,:) = values_contrast.%s.%s.%s.right.amplitude;', param.peaks{p}, param.N1_cond{n}, param.contrast{c});
                eval(statement)
                statement = sprintf('data_amplitude(2,:) = values_contrast.%s.%s.%s.right.amplitude;', param.peaks{p}, param.N1_cond{n}, param.contrast{c});
                eval(statement)
                statement = sprintf('data_latency(1,:) = values_contrast.%s.%s.%s.left.latency;', param.peaks{p}, param.N1_cond{n}, param.contrast{c});
                eval(statement)
                statement = sprintf('data_latency(2,:) = values_contrast.%s.%s.%s.left.latency;', param.peaks{p}, param.N1_cond{n}, param.contrast{c});
                eval(statement)

                % append mean values
                for a = 1:2
                % mean
                statement = sprintf('LEP_stats.%s.%s.%s.amplitude.mean(a) = mean(data_amplitude(a,:));', param.contrast{c}, param.peaks{p}, param.N1_cond{n});
                eval(statement)
                statement = sprintf('LEP_stats.%s.%s.%s.latency.mean(a) = mean(data_latency(a,:));', param.contrast{c}, param.peaks{p}, param.N1_cond{n});
                eval(statement)

                % SD
                statement = sprintf('LEP_stats.%s.%s.%s.amplitude.SD(a) = std(data_amplitude(a,:));', param.contrast{c}, param.peaks{p}, param.N1_cond{n});
                eval(statement)
                statement = sprintf('LEP_stats.%s.%s.%s.latency.SD(a) = std(data_latency(a,:));', param.contrast{c}, param.peaks{p}, param.N1_cond{n});
                eval(statement)

                % SEM
                statement = sprintf('LEP_stats.%s.%s.%s.amplitude.SEM(a) = LEP_stats.%s.%s.%s.amplitude.SD(a) / sqrt(n_subjects);', param.contrast{c}, param.peaks{p}, param.N1_cond{n}, param.contrast{c}, param.peaks{p}, param.N1_cond{n});
                eval(statement)
                statement = sprintf('LEP_stats.%s.%s.%s.latency.SEM(a) = LEP_stats.%s.%s.%s.latency.SD(a) / sqrt(n_subjects);', param.contrast{c}, param.peaks{p}, param.N1_cond{n}, param.contrast{c}, param.peaks{p}, param.N1_cond{n});
                eval(statement)
            end
            
                % plot and save amplitude
                fig = figure('Position', [100, 100, 700, 550]);
                plot_box(data_amplitude', 'amplitude', colours, labels)
                title(sprintf('%s average amplitude - %: %s vs. %s', param.peaks{p}, param.N1_cond{n}, labels{1}, labels{2}), 'fontsize', 16, 'fontweight', 'bold')
                saveas(fig, sprintf('%s\\figures\\LEP_avg_unfiltered_%s_%s_%s_amplitude.png', folder.output, param.contrast{c}, param.peaks{p}, param.N1_cond{n}))
    
                % plot and save latency
                fig = figure('Position', [100, 100, 700, 550]);
                plot_box(data_latency', 'latency', colours, labels)
                title(sprintf('%s average latency - %s: %s vs. %s', param.peaks{p}, param.N1_cond{n}, labels{1}, labels{2}), 'fontsize', 16, 'fontweight', 'bold')
                saveas(fig, sprintf('%s\\figures\\LEP_avg_unfiltered_%s_%_%s_latency.png', folder.output, param.contrast{c}, param.peaks{p}, param.N1_cond{n}))
            end
        else
            data_amplitude = []; data_latency = [];
            % select the data
            statement = sprintf('data_amplitude(1,:) = values_contrast.%s.%s.right.amplitude;', param.peaks{p}, param.contrast{c});
            eval(statement)
            statement = sprintf('data_amplitude(2,:) = values_contrast.%s.%s.right.amplitude;', param.peaks{p}, param.contrast{c});
            eval(statement)
            statement = sprintf('data_latency(1,:) = values_contrast.%s.%s.left.latency;', param.peaks{p}, param.contrast{c});
            eval(statement)
            statement = sprintf('data_latency(2,:) = values_contrast.%s.%s.left.latency;', param.peaks{p}, param.contrast{c});
            eval(statement)

            % append mean values
            for a = 1:2
                % mean
                statement = sprintf('LEP_stats.%s.%s.amplitude.mean(a) = mean(data_amplitude(a,:));', param.contrast{c}, param.peaks{p});
                eval(statement)
                statement = sprintf('LEP_stats.%s.%s.latency.mean(a) = mean(data_latency(a,:));', param.contrast{c}, param.peaks{p});
                eval(statement)

                % SD
                statement = sprintf('LEP_stats.%s.%s.amplitude.SD(a) = std(data_amplitude(a,:));', param.contrast{c}, param.peaks{p});
                eval(statement)
                statement = sprintf('LEP_stats.%s.%s.latency.SD(a) = std(data_latency(a,:));', param.contrast{c}, param.peaks{p});
                eval(statement)

                % SEM
                statement = sprintf('LEP_stats.%s.%s.amplitude.SEM(a) = LEP_stats.%s.%s.amplitude.SD(a) / sqrt(n_subjects);', param.contrast{c}, param.peaks{p}, param.contrast{c}, param.peaks{p});
                eval(statement)
                statement = sprintf('LEP_stats.%s.%s.latency.SEM(a) = LEP_stats.%s.%s.latency.SD(a) / sqrt(n_subjects);', param.contrast{c}, param.peaks{p}, param.contrast{c}, param.peaks{p});
                eval(statement)
            end

            % plot and save amplitude
            fig = figure('Position', [100, 100, 700, 550]);
            plot_box(data_amplitude', 'amplitude', colours, labels)
            title(sprintf('%s average amplitude: %s vs. %s', param.peaks{p}, labels{1}, labels{2}), 'fontsize', 16, 'fontweight', 'bold')
            saveas(fig, sprintf('%s\\figures\\LEP_avg_unfiltered_%s_%s_amplitude.png', folder.output, param.contrast{c},  param.peaks{p}))

            % plot and save latency
            fig = figure('Position', [100, 100, 700, 550]);
            plot_box(data_latency', 'latency', colours, labels)
            title(sprintf('%s average latency: %s vs. %s', param.peaks{p}, labels{1}, labels{2}), 'fontsize', 16, 'fontweight', 'bold')
            saveas(fig, sprintf('%s\\figures\\LEP_avg_unfiltered_%s_%_latency.png', folder.output, param.contrast{c}, param.peaks{p}))
        end                  
    end
end

% get some averages
mean([AP_info.single_subject.age])
std([AP_info.single_subject.age])
sum([AP_info.single_subject.male])
mean(LEP_stats.hand_hand.N1.ICA_filtered.amplitude.mean)
mean(LEP_stats.hand_hand.N1.ICA_filtered.amplitude.SD)
mean(LEP_stats.hand_hand.N1.ICA_filtered.latency.mean)
mean(LEP_stats.hand_hand.N1.ICA_filtered.latency.SD)
mean(LEP_stats.hand_hand.N2.amplitude.mean)
mean(LEP_stats.hand_hand.N2.amplitude.SD)
mean(LEP_stats.hand_hand.N2.latency.mean)
mean(LEP_stats.hand_hand.N2.latency.SD)

clear param eoi x_start x_end x a b c t p conditios statement data visual_data visual_CI_upper visual_CI_lower...
    visual_sem colours idx_contrast contrast condition n_subjects n labels fig data_amplitude data_latency s 

%% final data check 
% check for extra/missing events --> generate .txt report
name = sprintf('%s\\NLEP_event_analysis.txt', folder.output);
fileID = fopen(name, 'a');
data = NLEP_data;
for a = 1:length(data.LEP)
    if ~isempty(data.LEP(a).ID)
        fprintf(fileID, sprintf('subject %d (%s) - LEP\r\n', a, data.LEP(a).ID));
        % unfiltered LEP data
        for b = fieldnames(data.LEP(a).unfiltered)'
            for c = {'cond1' 'cond2'}
                if size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1) > 60
                    fprintf(fileID, sprintf('     --> unfiltered %s %s: - too many events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1)));
                elseif size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1) < 60
                    fprintf(fileID, sprintf('     --> unfiltered %s %s: - missing events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).unfiltered.(b{1}).(c{1}), 1)));
                end
            end
        end
        % CWT-filtered LEP data
        for b = fieldnames(data.LEP(a).CWT_filtered)'
            for c = {'cond1' 'cond2'}
                if size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1) > 60
                    fprintf(fileID, sprintf('     --> CWT_filtered %s %s: - too many events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1)));
                elseif size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1) < 60
                    fprintf(fileID, sprintf('     --> CWT_filtered %s %s: - missing events - %d\r\n', ...
                        b{1}, c{1}, size(data.LEP(a).CWT_filtered.(b{1}).(c{1}), 1)));
                end
            end
        end
        fprintf(fileID, sprintf('subject %d (%s) - RSEEG\r\n', a, data.LEP(a).ID));
        % PSD data
        for c = 1:length(data.RSEEG(a).PSD_st)
            if contains(data.RSEEG(a).dataset{c}, 'RS')
                if size(data.RSEEG(a).PSD_st(c).original, 1) > 59
                    fprintf(fileID, sprintf('     --> %s: - too many events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                elseif size(data.RSEEG(a).PSD_st(c).original, 1) < 59
                    fprintf(fileID, sprintf('     --> %s: - missing events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                end
            elseif contains(data.RSEEG(a).dataset{c}, 'LEP')
                if size(data.RSEEG(a).PSD_st(c).original, 1) > 30
                    fprintf(fileID, sprintf('     --> %s: - too many events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                elseif size(data.RSEEG(a).PSD_st(c).original, 1) < 30
                    fprintf(fileID, sprintf('     --> %s: - missing events - %d\r\n', ...
                        data.RSEEG(a).dataset{c}, size(data.RSEEG(a).PSD_st(c).original, 1)));
                end
            end
        end
        fprintf(fileID, '\r\n');
    end
end
fclose(fileID);

% remove specific events from specific datasets
a = 18; 
dataset = 'cond1';
% dataset = {'LEP hand left b2 - ready'};
trial = 31;
% if RSEEG
for d = 1:length(dataset)
    c = find(strcmp(NLEP_data_1to35.RSEEG(a).dataset, dataset{d}));
    for b = fieldnames(NLEP_data_1to35.RSEEG(a).PSD_st)'
        NLEP_data_1to35.RSEEG(a).PSD_st(c).(b{1})(trial, :, :) = [];
    end
end
% if LEP
for b = {'unfiltered' 'CWT_filtered'}   
    for c = fieldnames(NLEP_data_1to35.LEP(a).(b{1}))'
        NLEP_data_1to35.LEP(a).(b{1}).(c{1}).(dataset)(trial, :) = [];
    end
end
NLEP_data_1to35.LEP(a).blocks{str2double(regexp(dataset, '\d+', 'match', 'once'))}(trial) = [];
save(input_file, 'AP_info', 'NLEP_data_1to35', '-append')

% remove specific measures
s = 17; 
c = 1;
trial = 32;
for a = fieldnames(NLEP_measures(s).LEP_ST)'
    if ~strcmp(a{1}, 'conditions') 
        for b = 1:3
            NLEP_measures(s).LEP_ST.(a{1})(b, c, trial) = 0;
            idx = NLEP_measures(s).LEP_ST.(a{1})(b, c, :) ~= 0;
            data = NLEP_measures(s).LEP_ST.(a{1})(b, c, idx);
            data(1, 1, end+1:60) = 0;
            NLEP_measures(s).LEP_ST.(a{1})(b, c, :) = data;
        end
    end
end
save(input_file, 'NLEP_measures', '-append')

% check data in letswave
addpath(genpath([folder.toolbox '\letswave 7']));
letswave

%% helper lines
% update NLEP_info
NLEP_info_Domi = AP_info;
load(output_file, 'AP_info');
for i = [33,34,35]
    NLEP_info_Domi.single_subject(i).preprocessing = AP_info.single_subject(i).preprocessing;
end
AP_info = NLEP_info_Domi;
save('NLEP_info_backup.mat', 'AP_info');
save(output_file, 'AP_info', '-append');
clear i NLEP_info_Domi

% load again the 'dataset' variable
file2import = dir(sprintf('ds*S017*.mat'));
seq = 1; %[8, 7, 6, 5, 3, 4, 1, 2, 11, 12, 9, 10];
for d = 1:length(seq)
    load(file2import(seq(d)).name)
    load(sprintf('%s.lw6', file2import(seq(d)).name(1:end-4)), '-mat')
    % header.name = replace(header.name, '_', ' ');
    % header.name(1:5) = [];
    dataset(d).data = data;
    dataset(d).header = header;
end
clear file2import seq d header data header

% save again some data from 'dataset' variable
for d=3:length(dataset)
    LW_save(dataset(d).header.name, '', dataset(d).header, dataset(d).data)
end

% get rid of events with no-zero latency
for d = 1:length(dataset)
    % identify non-zero epochs
    for e = 1:length(dataset(d).header.events)
        if dataset(d).header.events(e).latency == 0
            idx(e) = false;
        else
            idx(e) = true;
        end
    end
    % remove these events 
    dataset(d).header.events(idx) = [];
    clear idx
    % save header to letswave
    header = dataset(d).header;
    save(sprintf('%s.lw6', dataset(d).header.name), 'header')
end
clear d e header

% re-organize output structure
for s = 1:length(AP_info.single_subject)-1
    % channel interpolation
    interp = AP_info.single_subject(s).preprocessing.ERP(4).process; 
    AP_info.single_subject(s).preprocessing.ERP(4).process = interp(1:25);
    AP_info.single_subject(s).preprocessing.ERP(4).params.interpolated = interp(28:end);
    AP_info.single_subject(s).preprocessing.ERP(4).params.sources = interp(28:end);

    % removed epochs
    stim = {'LEP' 'SEP'};
    block = {'b1' 'b2'};
    epochs = AP_info.single_subject(s).preprocessing.ERP(5).process;
    for t = 1:2
        for c = 1:2
            cond = sprintf('%s', replace(AP_info.single_subject(s).condition{c}, '_', ' '));
            for b = 1:2
                AP_info.single_subject(s).preprocessing.ERP(5).params(end+1).dataset = sprintf('%s %s %s', stim{t}, cond, block{b});        
            end
        end
    end
    AP_info.single_subject(s).preprocessing.ERP(5).params(1).discarded = epochs(23:end);

    % ICA
    AP_info.single_subject(s).preprocessing.ERP(6).process = '6 - ICA';
    AP_info.single_subject(s).preprocessing.ERP(6).params = AP_info.single_subject(s).preprocessing.ICA;
    AP_info.single_subject(s).preprocessing = rmfield(AP_info.single_subject(s).preprocessing, 'ICA');

    % order the structure
    AP_info.single_subject(s).preprocessing.ERP = orderfields(AP_info.single_subject(s).preprocessing.ERP, {'process' 'params' 'date'}); 
    save(output_file, 'AP_info', '-append');
end
clear s t c b interp stim block cond epochs 
AP_info.single_subject = orderfields(AP_info.single_subject, {'ID' 'date' 'age' 'male' 'handedness' 'body' 'QST' 'condition' 'temperature' 'ES' 'dataset' 'preprocessing' 'LEP'}); 
NLEP_measures = orderfields(NLEP_measures, {'ID' 'conditions' 'BDI' 'PCS' 'SES' 'RT' 'pain'}); 

% quick means
mean([AP_info.single_subject.age])
epochs = 0;
for a = 1:length(AP_info.single_subject)
    for b = 1:8
        epochs = epochs + numel(AP_info.single_subject(a).preprocessing.ERP(5).params(b).discarded);
    end
end
epochs/(length(AP_info.single_subject) * 8)
clear a b epochs

% next subject loop
next_subject = true;
while next_subject  
    % ask for subject number
    if ~exist('subject_idx')
        prompt = {'subject number:'};
        dlgtitle = 'subject';
        dims = [1 40];
        definput = {''};
        input = inputdlg(prompt,dlgtitle,dims,definput);
        subject_idx = str2num(input{1,1});
    end
    clear prompt dlgtitle dims definput input

    % ask for continuation
    answer = questdlg('Do you want to continue with next subject?', 'Continue?', 'YES', 'NO', 'YES'); 
    switch answer
        case 'YES'
            subject_idx = subject_idx + 1;
        case 'NO'
    	    next_subject = false;
            clear subject_idx answer
    end
end

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