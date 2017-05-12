function [ buf, len ] = Read_PCM( filename )
    %∂¡»°”Ô“Ù ˝æ›
    fp =fopen(filename,'rb');
    fseek(fp,0,'eof');
    len = ftell(fp)/2;
    fseek(fp,0,'bof');
    buf = fread(fp,len,'int16');
    fclose(fp);
end

