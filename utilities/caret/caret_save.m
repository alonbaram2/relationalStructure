function caret_save(filename,M)
% function caret_save(filename,M)
% Save a M-file structure as a binary or ascii metric, coord, topo, pait file
% M.encoding should be {'ASCII','BINARY'}
% --------------------------------------------------------------
% v.1.0 Joern Diedrichsen 04/11/28
if nargin==1
    M=filename;
    filename='';
end;
if isempty(filename)
    [F,P]=uiputfile('*.*','Save Structure as');
    filename = [P,F];
end
fid=fopen(filename,'w','ieee-be');
if (fid==-1)
    fprintf('Error opening file %s\n',filename);
    return;
end;
linefeed=sprintf('\n');

% Figure out the type of file we are loading:
s=strfind(filename,'.');
type=filename(s(end)+1:end);

% --------------------------------------------------------------------
% Now switch the format and make the header depending on the type of file
switch (type)
    % --------------------------------------------------------------------
    % Coordinate file
    case 'coord'
        if (~isfield(M,'header'))
            M.header{1}='BeginHeader';
            M.header{2}=['encoding ' M.encoding{1}];
            M.header{3}=['configuration_id ' M.configuration_id{1}];
            M.header{4}='EndHeader';
        end;
        M.num_nodes=size(M.data,1);
        for line=1:length(M.header)
            fprintf(fid,'%s\n',M.header{line});
        end;
        % Format string:
        if (strcmp(M.encoding,'ASCII'))
            fprintf(fid,'%d\n',M.num_nodes);
            M.data=[M.index M.data];
            format_str=['%i %3.3f %3.3f %3.3f\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.num_nodes,'int32');
            fwrite(fid,M.data','float32');
        end;
        % --------------------------------------------------------------------
        % Coordinate file
    case 'topo'
        if (~isfield(M,'header'))
            M.header{1}='BeginHeader';
            M.header{2}=['encoding ' M.encoding{1}];
            M.header{3}='filetype topo';
            M.header{4}=['perimeter_id ' M.perimeter_id{1}];
            if (isfield(M,'topo_file'))
                M.header{end+1}=['topo_file' M.topo_file];
            end;
            M.header{end+1}='EndHeader';
        end;
        M.num_tiles=size(M.data,1);
        for line=1:length(M.header)
            fprintf(fid,'%s\n',M.header{line});
        end;
        fprintf(fid,'tag-version 1\n');
        M.data=M.data-1;
        if (strcmp(M.encoding,'ASCII'))
            format_str=['%d %d %d\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.num_tiles,'int32');
            fwrite(fid,M.data','int32');
        end;
        % --------------------------------------------------------------
        % Metric file
    case {'surface_shape','metric'}
        M.header={};
        M.header{1}='BeginHeader';
        M.header{2}=['encoding ' M.encoding{1}];
        M.header{3}='EndHeader';
        M.header{4}='tag-version 2';
        M.header{5}=sprintf('tag-number-of-nodes %d',size(M.data,1));
        M.header{6}=sprintf('tag-number-of-columns %d',size(M.data,2));
        M.header{7}=sprintf('tag-title');
        where=7;
        for i=1:length(M.column_name)
            M.header{where+i}=sprintf('tag-column-name %d %s',i-1,M.column_name{i});
        end;
        where=where+length(M.column_name);
        for i=1:length(M.column_name)
            M.header{where+i}=sprintf('tag-column-color-mapping %d %8.8f %8.8f',i-1,M.column_color_mapping(i,1),M.column_color_mapping(i,2));
        end;
        where=where+length(M.column_name);
        M.header{where+1}='tag-BEGIN-DATA';
        
        for line=1:length(M.header)
            fprintf(fid,'%s\n',M.header{line});
        end;
        % Format string:
        if (strcmp(M.encoding,'ASCII'))
            M.data=[M.index M.data];
            format_str=['%i' repmat([' %6.6f'],1,size(M.data,2)-1) '\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.data','float32');
        end;
        % --------------------------------------------------------------
        % RGB_paint
    case 'RGB_paint'
        M.header={};
        M.header{1}='BeginHeader';
        if (~isfield(M,'encoding'))
            M.encoding={'BINARY'};
        end;
        
        M.header{2}=['encoding ' M.encoding{1}];
        M.header{3}='EndHeader';
        M.header{4}='tag-version 2';
        M.header{5}=sprintf('tag-number-of-nodes %d',size(M.data,1));
        M.header{6}=sprintf('tag-number-of-columns %d',size(M.data,2)/3);
        where=7;
        for i=1:length(M.scales)
            M.header{where}=sprintf('tag-scale-red %d %8.8f %8.8f',i-1,M.scales{i}(1,1),M.scales{i}(1,2));
            M.header{where+1}=sprintf('tag-scale-green %d %8.8f %8.8f',i-1,M.scales{i}(2,1),M.scales{i}(2,2));
            M.header{where+2}=sprintf('tag-scale-blue %d %8.8f %8.8f',i-1,M.scales{i}(3,1),M.scales{i}(3,2));
            where=where+3;
        end;
        M.header{where}='tag-BEGIN-DATA';
        
        for line=1:length(M.header)
            fprintf(fid,'%s\n',M.header{line});
        end;
        % Format string:
        if (strcmp(M.encoding,'ASCII'))
            M.data=[M.index M.data];
            format_str=['%i' repmat([' %6.6f'],1,size(M.data,2)-1) '\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.data','float32');
        end;
    case 'paint'
        M.header={};
        M.header{1}='BeginHeader';
        if (~isfield(M,'encoding'))
            M.encoding={'BINARY'};
        end;
        M.header{2}=['encoding ' M.encoding{1}];
        M.header{3}='EndHeader';
        M.header{4}='tag-version 1';
        M.header{5}=sprintf('tag-number-of-nodes %d',size(M.data,1));
        M.header{6}=sprintf('tag-number-of-columns %d',size(M.data,2));
        M.header{7}=sprintf('tag-title');
        M.num_paintnames=length(M.paintnames);
        M.header{8}=sprintf('tag-number-of-paint-names %d',M.num_paintnames);
        where=8;
        for i=1:length(M.column_name)
            M.header{where+i}=sprintf('tag-column-name %d %s',i-1,M.column_name{i});
        end;
        where=where+length(M.column_name);
        M.header{where+1}='tag-BEGIN-DATA';
        
        for line=1:length(M.header)
            fprintf(fid,'%s\n',M.header{line});
        end;
        
        
        % read number of paintnames
        for pn=1:M.num_paintnames
            fprintf(fid,'%d %s\n',pn-1,M.paintnames{pn});
        end;
        
        % Format string:
        if (strcmp(M.encoding,'ASCII'))
            M.data=[M.index M.data];
            format_str=['%i' repmat([' %6.6f'],1,size(M.data,2)-1) '\n'];
            fprintf(fid,format_str,M.data');
        else
            fwrite(fid,M.data','int32');
        end;
    case 'spec'
        for line=1:length(M.header)
            fprintf(fid,'%s\n',M.header{line});
        end;
        for line=1:length(M.data)
            fprintf(fid,'%s\n',M.data{line});
        end;
    case 'areacolor'
        % Make the area-color file (comma-seperated)
        
        AC{1}='CSVF-FILE,0,,,,,,,';
        AC{2}='csvf-section-start,header,2,,,,,,';
        AC{3}='tag,value,,,,,,,';
        AC{4}='Caret-Version,5.65,,,,,,,';
        AC{5}='Date,2012-08-19T11:19:34,,,,,,,';
        AC{6}='comment,,,,,,,,';
        AC{7}='encoding,COMMA_SEPARATED_VALUE_FILE,,,,,,,';
        AC{8}='pubmed_id,,,,,,,,';
        AC{9}='csvf-section-end,header,,,,,,,';
        AC{10}='csvf-section-start,Colors,9,,,,,,';
        AC{11}='Name,Red,Green,Blue,Alpha,Point-Size,Line-Size,Symbol,SuMSColorID';
        for i=1:size(M.data,1)
            AC{end+1}=sprintf('%s,%d,%d,%d,255,2.000,1.000,POINT,',M.paintnames{i},M.data(i,1),M.data(i,2),M.data(i,3));
        end;
        AC{end+1}='csvf-section-end,Colors,,,,,,,';
        
        for i=1:length(AC);
            fprintf(fid,'%s\n',AC{i});
        end;
end;
fclose(fid);
