function volumeViewBin(fileName, imgSize)
    arguments
        fileName string
        imgSize {mustBeInteger}
    end

    fid = fopen(fileName);
    if (fid == 0)
        error("Unable to open %s.\n", fileName);
    end

    data = fread(fid, 'uint8');
    fclose(fid);

    data = reshape(data, imgSize);
    
    volumeViewer(data);
end