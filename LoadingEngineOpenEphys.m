function [out1,out2] = LoadingEngineOpenEphys(varargin)


% Loading engine requirements
% A loading engine must be callable as a Matlab function. This means it must be either a .m
% function-file or a compiled mex function-file. It must take as input one or three inputs and
% provide one or two outputs. MClust-3.0 should work with any loading engine that supplies this
% functionality.
% INPUTS
% ∑ fn = file name string
% ∑ records_to_get = a range of values
% ∑ record_units = a flag taking one of 5 cases (1,2,3,4 or 5)
% 1. implies that records_to_get is a timestamp list.
% 2. implies that records_to_get is a record number list
% 3. implies that records_to_get is range of timestamps (a vector with 2 elements: a
% start and an end timestamp)
% 4. implies that records_to_get is a range of records (a vector with 2 elements: a
% start and an end record number)
% 5. asks to return the count of spikes (records_to_get should be [] in this case)
% In addition, if only fn is passed in then the entire file should be read.
% OUTPUT
% ∑ [t, wv]
% ∑ t = n x 1: timestamps of each spike in file
% ∑ wv = n x 4 x 32 waveforms
% EXAMPLES
% ∑ [t,wv] = myLoadingEngine('myfile.dat', 1:10, 2) should return the time and waveforms
% for the first 10 spikes in the file.
% ∑ t = myLoadingEngine('myfile.dat') should return all the timestamps from the file.
% ∑ n = myLoadingEngine('myfile.dat', [], 5) should return the number of spikes in the file.


fn = varargin{1};

fid = fopen(fn);
% constants
NUM_HEADER_BYTES = 1024;
hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
eval(char(hdr'));
num_channels = header.num_channels;
num_samples = 40; % **NOT CURRENTLY WRITTEN TO HEADER**
sampleRate = header.sampleRate;

global MClust_AverageWaveform_ylim
MClust_AverageWaveform_ylim = [-500 500];

% uint8 eventType
% int64 timestamp (to align with timestamps from the continuous records)
% uint16 electrodeID
% uint16 number of channels (N)
% uint16 number of samples per spike (M)
% N*M uint16 samples (individual channels are contiguous)
% N uint16 gains (actually gain*1000, to increase resolution)
% N uint16 thresholds used for spike extraction
% One uint16 recording number (version 0.2 and higher)

event_type = fread(fid, 1, 'uint8','l'); % always equal to 4; ignore 1


OneRecordSize = 17+2*num_channels*(num_samples+2);

if nargin == 1
    t = fread(fid,inf,'int64',9+num_channels*(num_samples*2+2+2),'l');
    t = 10000*t/sampleRate;
    fseek(fid, 1024+15, 'bof')
    data = fread(fid,inf,[num2str(num_channels*num_samples) '*uint16'],...
        17+4*num_channels, 'l');
    data = (data-32768)/5;
    fclose(fid); 
    wv = reshape(data,[],num_samples,num_channels);
    wv = permute(wv,[1 3 2]);
    out1 = t;
    out2 = wv;
    
elseif nargin == 3
    records_to_get = varargin{3};
    record_units = varargin{2};
    switch records_to_get
    case 1 % timestamp list
        wv = zeros(length(record_units),num_channels,num_samples);
        t = fread(fid,inf,'int64',9+num_channels*(num_samples*2+2+2),'l');
        t = 10000*t/sampleRate;
        for i = 1 : length(record_units)
            offset = find(t>=record_units(i),1,'first');
            fseek(fid, 1024+(offset-1)*OneRecordSize+15, 'bof');
            data = fread(fid,num_channels*num_samples,'*uint16','l');
            data = (data-32768)/5;
            wv(i,:,:) = reshape(data,num_samples,num_channels)';
        end
        fclose(fid);
        out1= record_units;
        out2 = wv;
                
    case 2  % record number list        
        t = zeros(record_units,1);
        wv = zeros(length(record_units),num_channels,num_samples);
        for i = 1 : length(record_units)
            fseek(fid, 1024+(record_units(i)-1)*OneRecordSize+1, 'bof');
            t(i) = fread(fid,1,'*int64','l');
            fseek(fid, 1024+(record_units(i)-1)*OneRecordSize+15, 'bof');
            data = fread(fid,num_channels*num_samples,'*uint16','l');
            data = (data-32768)/5;
            wv(i,:,:) = reshape(data,num_samples,num_channels)';
        end
        t = 10000*t/sampleRate;
        fclose(fid);
        out1= t;
        out2 = wv;

        
    case 3  % range of time stamp
        t = fread(fid,inf,'int64',9+num_channels*(num_samples*2+2+2),'l');
        t = 10000*t/sampleRate;
        offset = find(t>=record_units(1),1,'first');
        l = find(t<=record_units(2),1,'last')-offset+1;
        t = t(offset:offset+l-1);
        fseek(fid,1024+OneRecordSize*(offset-1)+15,'bof');
        data = fread(fid,l,[num2str(num_channels*num_samples) '*uint16'],...
            17+4*num_channels, 'l');
        data = (data-32768)/5;
        fclose(fid); 
        wv = reshape(data,[],num_samples,num_channels);
        wv = permute(wv,[1 3 2]);
        out1 = t;
        out2 = wv;
        
    case 4  % range of record number
        fseek(fid,1024+OneRecordSize*(record_units(1)-1)+1,'bof');
        l = record_units(2) - record_units(1)+1;
        t = fread(fid,l,'int64',9+num_channels*(num_samples*2+2+2),'l');
        fseek(fid,1024+OneRecordSize*(record_units(1)-1)+15,'bof');
        data = fread(fid,l*num_channels*num_samples,[num2str(num_channels*num_samples) '*uint16'],...
            17+4*num_channels, 'l');
        data = (data-32768)/5;
        fclose(fid); 
        wv = reshape(data,[],num_samples,num_channels);
        wv = permute(wv,[1 3 2]);
        out1 = t;
        out2 = wv;
    case 5  % the count of spikes
        fid = fopen(fn);
        fseek(fid, 0, 'eof');
        filelength = ftell(fid);
        fclose(fid);
        out1 = (filelength-1024)/OneRecordSize;
        out2 = [];
    end
end
