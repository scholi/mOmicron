classdef Omicron
    %   Detailed explanation goes here
    
    properties
        Path % Path of the folder containing all the mtrx files
        fid=-1; % File ID
        IDs={}; % IDs found
        Params % List of the parameters
        Images % List of the images found
        nImages = 0;
        EOF = 0;
    end
    
    methods
        function self = Omicron(Path)
            self.Path = Path;
            Files = dir(self.Path);
            for i = 1:length(Files)
                FileName = Files(i).name;
                if length(FileName)>10
                    Ending = FileName(end-9:end);
                    if strcmp(Ending,'_0001.mtrx')
                        self.fid = fopen(strcat(self.Path,filesep,FileName),'rb');
                        fi=dir(strcat(self.Path,filesep,FileName));
                        FileLength = fi.bytes;
                    end
                end
            end
            if self.fid==-1
                error('The file ..._0001.mtrx was nor found');
            end
            if ~strcmp(fread(self.fid,8,'*char')','ONTMATRX')
                error(strcat('Unknown header! Wrong Matrix file format for file ',FileName));
            end
            fread(self.fid,4,'*char'); % Should be 0101
            h = waitbar(0,'Parsing data...');
            while ~self.EOF
                self=self.read_block();
                waitbar(ftell(self.fid)/FileLength);
            end
            close(h);
            fclose(self.fid);
        end  
        function [xvalues,ydata] = getSTS(self, ID, num)
            i = self.IDs{ID}{num};
            img = self.Images{i};
            f = fopen(strcat(self.Path,filesep,img.filename),'rb');
            if ~strcmp(fread(f,8,'*char')','ONTMATRX')
                error('Invalid STS file format');
            end
            if ~strcmp(fread(f,4,'*char')','0101')
                error('Invalid STS version');
            end
            fread(f,24);
            ss = fread(f,15,'uint32');
            if ~strcmp(fread(f,4,'*char')','ATAD')
                error('Data should be here, but are not. Please debug script.');
            end
            fread(f,4);
            ydata = fread(f, ss(8), 'int32')';
            % Reconstruct the x-axis. Take the start and end volatege (v1,v2) with the correct number of points and pad it to the data length. Padding is in 'reflect' mode in the case of Forward/backward scans.
            v1 = img.params.Spectroscopy.Device_1_Start.value; % Get the start voltage used for the scan
            v2 = img.params.Spectroscopy.Device_1_End.value; % Get the end voltage for the scan
            N = img.params.Spectroscopy.Device_1_Points.value;
            xvalues = linspace(v1,v2,N);
            if length(xvalues)<ss(7)
                xvalues=[xvalues,fliplr(xvalues)];
            end

            if length(ydata)<length(xvalues)
                ydata = [ydata, NaN(1, length(xvalues)-length(ydata))];
            end
            
        end
        function self = read_block(self, sub)
            indent = fread(self.fid,4,'*char')'; % 4bytes forming the header. Those are capital letters between A-Z
            if length(indent)<4
                self.EOF=1; % Enf-Of-File (EOF) reached
                return
            end
            bs = fread(self.fid, 1,'uint32'); % Block size
            if nargin<2
                bs = bs + 8;
            end
            p = ftell(self.fid);
            if strcmp(indent,'DOMP') % Block storing the parameters changed during an experiment
                fread(self.fid,12); % Skip 12 bytes
                inst = self.read_string();
                prop = self.read_string();
                unit = self.read_string();
                fread(self.fid,4);
                value=self.read_value();
                self.Params.(inst).(prop).value = value;
                self.Params.(inst).(prop).unit = unit;
%             elseif strcmp(indent,'CORP') % Processor of scanning window. Useless in this script for the moment
%                 fread(self.fid,12);
%                 self.read_string();
%                 self.read_string();
            elseif strcmp(indent,'FERB') % A file was stored
                fread(self.fid, 12);
                FileName = self.read_string();
                self.nImages = self.nImages + 1;
                self.Images{self.nImages}.filename = FileName;
                self.Images{self.nImages}.params = self.Params(:); % Store the actual parameters uded to record the file
                [token,~] = regexp(FileName, '^(.*?)--([0-9]*)_([0-9]*)\.([^_]+)_mtrx$','tokens','match');
                ID = str2num(token{1}{2});
                Num = str2num(token{1}{3});
                self.IDs{ID}{Num}=self.nImages;
%             elseif strcmp(indent,'SPXE')
%                 fread(self.fid,12);
%                 self.read_block(1); % read sub-block
%                 self.read_block(1);
%                 self.read_block(1);
%             elseif strcmp(indent,'LNEG')
%                 self.read_string();
%                 self.read_string();
%                 self.read_string();
%             elseif strcmp(indent, 'TSNI')
%                 anz = fread(self.fid,1,'uint32');
%                 rr = {};
%                 for ai = 1:anz
%                     rr{ai}.a = self.read_string();
%                     rr{ai}.b = self.read_string();
%                     rr{ai}.c = self.read_string();
%                     count = fread(self.fid,1,'uint32');
%                     rr{ai}.content={};
%                     for i = 1:count
%                         rr{ai}.content{i}.x = self.read_string();
%                         rr{ai}.content{i}.y = self.read_string();
%                     end
%                 end
%             elseif strcmp(indent,'SXNC')
%                 count = fread(self.fid,1,'uint32');
%                 rr={};
%                 for i=1:count
%                     rr{i}.a = self.read_string();
%                     rr{i}.b = self.read_string();
%                     rr{i}.i = i;
%                     k = fread(self.fid,1,'uint32');
%                     for j = 1:k
%                         rr{i}.kk{j}.x = self.read_string();
%                         rr{i}.kk{j}.y = self.read_string();
%                     end
%                 end
            elseif strcmp(indent,'APEE')
                fread(self.fid,12);
                num = fread(self.fid,1,'uint32');
                for i = 1:num    
                    inst = self.read_string();
                    grp = fread(self.fid,1,'uint32');
                    for j = 1:grp
                        prop = self.read_string();
                        self.Params.(inst).(prop).unit =self.read_string();
                        fread(self.fid,4);
                        self.Params.(inst).(prop).value = self.read_value();
                    end
                end
            end
            fseek(self.fid,p,-1); % Go back to the beginning of the block
            fseek(self.fid,bs,0); % Go to the next block
        end
        function out = read_string(self)
            N = fread(self.fid,1,'uint32');
            if N==0
                out='';
                return
            end
            bytes = fread(self.fid,2*N)';
            out = native2unicode(bytes,'UTF-16LE');
        end
        function out = read_value(self)
            t = fread(self.fid,4,'*char')';
            if strcmp(t,'BUOD')
                out = fread(self.fid,1,'float64');
            elseif strcmp(t,'GNOL')
                out = fread(self.fid,1,'uint32');
            elseif strcmp(t,'LOOB')
                out = fread(self.fid,1,'uint32')>0;
            elseif strcmp(t,'GRTS')
                out = self.read_string();
            else
                out = t;
            end
        end
    end
    
end

