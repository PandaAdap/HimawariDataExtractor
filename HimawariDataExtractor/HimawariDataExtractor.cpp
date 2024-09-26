#include "HimawariDataExtractor.h"

HimawariDataExtractor::HimawariDataExtractor()
{
}

std::string HimawariDataExtractor::GetDateStrByTime(double time)
{
    // MJD to days since Unix epoch
    double daysSinceEpoch = time - 40587.0;

    // Days to seconds
    double secondsSinceEpoch = daysSinceEpoch * 86400.0;

    // Convert seconds to FILETIME
    ULARGE_INTEGER largeInt;
    largeInt.QuadPart = static_cast<unsigned long long>(secondsSinceEpoch * 10000000) + 116444736000000000LL;

    FILETIME fileTime;
    fileTime.dwLowDateTime = largeInt.LowPart;
    fileTime.dwHighDateTime = largeInt.HighPart;

    // Convert FILETIME to SYSTEMTIME
    SYSTEMTIME utcTime;
    FileTimeToSystemTime(&fileTime, &utcTime);

    //Format
    char dateUTC[32] = { 0 };
    sprintf_s(dateUTC, "%4d-%02d-%02d %02d:%02d:%02d", utcTime.wYear, utcTime.wMonth, utcTime.wDay, utcTime.wHour, utcTime.wMinute, utcTime.wSecond);

    return std::string(dateUTC);
}

BOOL HimawariDataExtractor::LoadHimawariDatas(const std::vector<std::string>& himawariDataFile)
{
    //Release all resource before load new data
    if (HimawariDataMode != -1)
    {
        ReleaseAll();
    }

    BOOL ret = -99;

    if (himawariDataFile.size() == 1)
    {
        //Check if netCDF format data
        try
        {
            ncInstance.open(himawariDataFile.front(), netCDF::NcFile::read);
            
            //netCDF success
            ret = HDE_SUCCESS;
            HimawariDataMode = 0;

            //Load all layers
            auto map_vars = ncInstance.getVars();
            for (auto it = map_vars.begin(), end = map_vars.end(); it != end; it = map_vars.upper_bound(it->first))
            {
                //Only phys value layer
                if ((int)it->first.find("albedo") >= 0 || (int)it->first.find("tbb") >= 0)
                    ncInfo.layers.emplace_back(it->first);
            }

            //Get satellite name
            auto att_satname = ncInstance.getAtt("title");
            char* _sat = new char[att_satname.getAttLength() + 1];
            memset(_sat, 0, att_satname.getAttLength() + 1);
            att_satname.getValues(_sat);
            ncInfo.satellite = std::string(_sat);
            ncInfo.satellite = ncInfo.satellite.substr(ncInfo.satellite.find("Himawari-", 0), ncInfo.satellite.find(" ", 0));
            delete[] _sat; _sat = nullptr;

            //Get date
            auto val_st = ncInstance.getVar("start_time");
            double st = 0.0;
            val_st.getVar(&st);
            ncInfo.date = GetDateStrByTime(st);
            
        }
        catch (netCDF::exceptions::NcException e)
        {
            ret = HDE_INVALID_NETCDF;
        }
    }

    //HSD data
    if(ret != HDE_SUCCESS)
        ret = hsdOpenDAT(himawariDataFile);

    //Free
    if (ret != HDE_SUCCESS)
        ReleaseAll();

    return ret;
}

void HimawariDataExtractor::ReleaseAll()
{
    //Release netCDF
    ncFreeAllResources();

    //Release HSD
    hsdFreeAllResources();

    //Reset to UnLoad (Reset already)
    //HimawariDataMode = -1;
}

HimawariDataExtractor::~HimawariDataExtractor()
{
    ReleaseAll();
}

void HimawariDataExtractor::hsdFreeAllResources()
{
    //Free hsdList
    for (auto& DAT : hsdList)
    {
        if (DAT.hsdFile)
        {
            fclose(DAT.hsdFile);
            DAT.hsdFile = nullptr;
        }
        if (DAT.hsdDataRangeRaw)
        {
            delete[] DAT.hsdDataRangeRaw;
            DAT.hsdDataRangeRaw = nullptr;
        }

        hisd_free(&DAT.hsdHeader);
    }hsdList.clear();

    //Free RAW data
    if (hsdRawData.lat)
    { 
        delete[] hsdRawData.lat;
        hsdRawData.lat = nullptr;
    }

    if (hsdRawData.lon)
    {
        delete[] hsdRawData.lon;
        hsdRawData.lon = nullptr;
    }

    if (hsdRawData.phys) 
    {
        delete[] hsdRawData.phys;
        hsdRawData.phys = nullptr;
    }
        
    //Clear satellite name
    memset(hsdRawData.satName, 0, 32);

    //Clear data mode status
    HimawariDataMode = -1;
}

BOOL HimawariDataExtractor::hsdOpenDAT(const std::vector<std::string>& inputDATs)
{
    //No inputs
    if (inputDATs.empty())
        return HDE_NOINPUT;

    //Release before load
    if (HimawariDataMode != -1)
        ReleaseAll();

    //Load inputs
    for (const auto& DAT : inputDATs)
    {
        HSDINSTANCE hsd;
        fopen_s(&hsd.hsdFile, DAT.c_str(), "rb");

        //Open failed
        if (!hsd.hsdFile)
            return HDE_OPENFAILED;

        //Invalid HSD data
        if (NORMAL_END != hisd_read_header(&hsd.hsdHeader, hsd.hsdFile))
        {
            //Close all file handle & clear memory
            hsdFreeAllResources();
            return HDE_INVALID_HSD;
        }

        //Seek to data range
        fseek(hsd.hsdFile, hsd.hsdHeader.basic->totalHeaderLen, 0);

        //alloc memory
        hsd.hsdDataRangeRaw = new unsigned short[hsd.hsdHeader.basic->dataLen];

        //Read raw data
        fread(hsd.hsdDataRangeRaw, sizeof(unsigned short) * hsd.hsdHeader.basic->dataLen, 1, hsd.hsdFile);

        //Load to memory
        hsdList.emplace_back(hsd);

    }

    //Single input-----------------------------------------------------------------------------------

    //Single data
    if (hsdList.size() == 1 && hsdList[0].hsdHeader.seg->totalSegNum == 1)
    {
        //Process data
        hsdProcessInputData();
        return HDE_SUCCESS;
    }
    //Not HSD sequence
    else if (hsdList.size() != 10)
    {
        //Close all file handle & clear memory
        hsdFreeAllResources();
        return HDE_INVALID_HSDSEQ;
    }

    //Multiple input---------------------------------------------------------------------------------

    //Sort by hsdHeader.seg->segSeqNo
    std::sort(hsdList.begin(), hsdList.end(), [](const HSDINSTANCE& lhs, const HSDINSTANCE& rhs) {
        return lhs.hsdHeader.seg->segSeqNo < rhs.hsdHeader.seg->segSeqNo;
        });

    //Check sequence
    for (size_t i = 0; i < hsdList.size() - 1; ++i)
    {
        const auto& curr = hsdList.at(i);
        const auto& next = hsdList.at(i + 1);

        //Missing part of sequence
        if (next.hsdHeader.seg->segSeqNo != curr.hsdHeader.seg->segSeqNo + 1)
        {
            //Close all file handle & clear memory
            hsdFreeAllResources();
            return HDE_HSD_SEQMISSINGPART;
        }

        //NOT Part of sequence
        if (strncmp(curr.hsdHeader.basic->fileName, next.hsdHeader.basic->fileName, 33) != 0)
        {
            //Close all file handle & clear memory
            hsdFreeAllResources();
            return HDE_HSD_SEQINVALIDPART;
        }
    }

    //Process Data------------------------------------------------------------------------------------
    hsdProcessInputData();

    return HDE_SUCCESS;
}

BOOL HimawariDataExtractor::hsdProcessInputData()
{
    //Get basic Raw Data
    hsdRawData.height = hsdList.front().hsdHeader.data->nLin * (int)hsdList.size();
    hsdRawData.width = hsdList.front().hsdHeader.data->nPix;
    hsdRawData.band = hsdList.front().hsdHeader.calib->bandNo;
    hsdRawData.bandWaveLength = hsdList.front().hsdHeader.calib->waveLen;
    strcpy_s(hsdRawData.satName, hsdList.front().hsdHeader.basic->satName);

    if (hsdRawData.band <= 6 && hsdRawData.band > 0)
    {
        hsdRawData.layerName = "albedo_" + std::to_string(hsdRawData.band);
    }
    else if (hsdRawData.band <= 16 && hsdRawData.band > 6)
    {
        hsdRawData.layerName = "tbb_" + std::to_string(hsdRawData.band);
    }

    //For calc scan time
    float minLine = 99999.0;
    float maxLine = -99999.0;

    //Letf top postition
    double ltlat = 0, ltlon = 0;

    //Get precision
    double dt = 0.02;//GetDeltaValue(hsdRawData.band);

    //Total size
    unsigned long n_size = hsdRawData.width * hsdRawData.height;

    //Loop index
    int jj = 0, ii = 0;

    //To judge if FLDK or TARGET DAT
    pixlin_to_lonlat(&hsdList.front().hsdHeader, 0, 0, &ltlon, &ltlat);
    if (ltlon == -9999 || ltlat == -9999)
    {
        ltlon = 80, ltlat = 60; //FLDK
    }
    else
    {
        double rblon, rblat;
        pixlin_to_lonlat(&hsdList.front().hsdHeader, hsdRawData.width, hsdRawData.height, &rblon, &rblat);
        dt = (rblon - ltlon) / hsdRawData.width;
    }

    //Pre alloc memory for data
    hsdRawData.phys = new float[n_size];
    hsdRawData.lon = new float[hsdRawData.width];
    hsdRawData.lat = new float[hsdRawData.height];

    //Initlizate Lon&Lat base on precision
    for (ii = 0; ii < hsdRawData.height; ++ii)
    {
        hsdRawData.lat[ii] = ltlat - dt * ii;
    }
    for (ii = 0; ii < hsdRawData.width; ++ii)
    {
        hsdRawData.lon[ii] = ltlon + dt * ii;
    }

#pragma omp parallel for num_threads(6)//Multi-thread
    //Each line
    for (int j = 0; j < hsdRawData.height; ++j)
    {
        //Each column
        for (int i = 0; i < hsdRawData.width; ++i)
        {
            //Index of total phys data
            unsigned long index_phys = j * hsdRawData.width + i;

            //Temporary Pix&Lin to get data from DAT
            float Pix, Lin;

            //For vaild getting data
            unsigned short count = hsdList.back().hsdHeader.calib->outCount;

             //Get Pix&Lin from long&lat
            lonlat_to_pixlin(&hsdList.back().hsdHeader, hsdRawData.lon[i], hsdRawData.lat[j], &Pix, &Lin);

            //Invalid data
            if (Pix == -9999 || Lin == -9999)
            {
                hsdRawData.phys[index_phys] = -1;
                continue;
            }

            //For calc scan time
            if (minLine > Lin) minLine = Lin;
            if (maxLine < Lin) maxLine = Lin;

            //Get count from DAT by Pix&Lin
            //Optmized for multi-threads
            int pt = Lin / hsdList.back().hsdHeader.data->nLin;
            int dat_index = Lin / (hsdList.back().hsdHeader.data->nLin + 0.5 / pt);

            //Safe check
            if (dat_index >= hsdList.size())
            {
                hsdRawData.phys[index_phys] = -1;
                continue;
            }

            //Not thread safed
            //hisd_getdata_by_pixlin(&hsdList[dat_index].hsdHeader, hsdList[dat_index].hsdFile, Pix, Lin, &count);

            //Load from memory
            int dataindex = 0;
            hisd_getdataindex_by_pixlin(&hsdList[dat_index].hsdHeader, Pix, Lin, &dataindex);
            count = hsdList[dat_index].hsdDataRangeRaw[dataindex];
            
            //Invalid count value
            if (count == hsdList.back().hsdHeader.calib->outCount ||
                count == hsdList.back().hsdHeader.calib->errorCount)
            {
                hsdRawData.phys[index_phys] = -1;
                continue;
            }

            //Convert count value to radiance
            float radiance = (float)count * hsdList.back().hsdHeader.calib->gain_cnt2rad + hsdList.back().hsdHeader.calib->cnst_cnt2rad;

            //Convert radiance to physical value
            if ((hsdList.back().hsdHeader.calib->bandNo >= 7 && strstr(hsdList.back().hsdHeader.basic->satName, "Himawari") != NULL) ||
                (hsdList.back().hsdHeader.calib->bandNo >= 2 && strstr(hsdList.back().hsdHeader.basic->satName, "MTSAT-2") != NULL))
            {
                //Infrared band
                double phys = 0;
                hisd_radiance_to_tbb(&hsdList.back().hsdHeader, radiance, &phys);
                hsdRawData.phys[index_phys] = (float)phys;
            }
            else
            {
                //Visible or near infrared band
                hsdRawData.phys[index_phys] = hsdList.back().hsdHeader.calib->rad2albedo * radiance;
            }
        }
	}
#pragma omp barrier

	//Convert maxLine & minLine to scanTime
	for (ii = 1; ii < hsdList.front().hsdHeader.obstime->obsNum; ii++) {
		if (minLine < hsdList.front().hsdHeader.obstime->lineNo[ii]) {
			hsdRawData.startTime = hsdList.front().hsdHeader.obstime->obsMJD[ii - 1];
			break;
		}
		else if (minLine == hsdList.front().hsdHeader.obstime->lineNo[ii]) {
			hsdRawData.startTime = hsdList.front().hsdHeader.obstime->obsMJD[ii];
			break;
		}
	}
	//EndTime
	for (ii = 1; ii < hsdList.back().hsdHeader.obstime->obsNum; ii++) {
		if (maxLine < hsdList.back().hsdHeader.obstime->lineNo[ii]) {
			hsdRawData.endTime = hsdList.back().hsdHeader.obstime->obsMJD[ii - 1];
		}
		else if (maxLine == hsdList.back().hsdHeader.obstime->lineNo[ii]) {
			hsdRawData.endTime = hsdList.back().hsdHeader.obstime->obsMJD[ii];
		}
	}
	if (hsdRawData.endTime < hsdRawData.startTime)
		hsdRawData.endTime = hsdRawData.startTime;

	//For reconize data mode
	HimawariDataMode = 1;

	//Release temporary raw
	hsdReleaseTemporaryRaw();

	return TRUE;
}

BOOL HimawariDataExtractor::hsdGetRaw(HSDRAWDATA& hsdRawOut)
{
    //NOT HSD or NOT load
    if (HimawariDataMode != 1)
        return HDE_INVALID_HSD;

    //Clean first
    hsdReleaseRaw(&hsdRawOut);

    //Get values
    hsdRawOut.band = hsdRawData.band;
    hsdRawOut.bandWaveLength = hsdRawData.bandWaveLength;
    hsdRawOut.height = hsdRawData.height;
    hsdRawOut.width = hsdRawData.width;
    hsdRawOut.startTime = hsdRawData.startTime;
    hsdRawOut.endTime = hsdRawData.endTime;

    strcpy_s(hsdRawOut.satName, hsdRawData.satName);
    hsdRawOut.layerName = hsdRawData.layerName;

    //Alloc memory
    unsigned long n_size = hsdRawOut.width * hsdRawOut.height;
    hsdRawOut.lon = new float[hsdRawOut.width];
    hsdRawOut.lat = new float[hsdRawOut.height];
    hsdRawOut.phys = new float[n_size];

    //Copy memory
    memcpy_s(hsdRawOut.lon, hsdRawOut.width, hsdRawData.lon, hsdRawData.width);
    memcpy_s(hsdRawOut.lat, hsdRawOut.height, hsdRawData.lat, hsdRawData.height);
    memcpy_s(hsdRawOut.phys, n_size, hsdRawData.phys, n_size);

    return HDE_SUCCESS;
}

void HimawariDataExtractor::hsdReleaseRaw(HSDRAWDATA* hsdRaw)
{
    if (hsdRaw)
    {
        if (hsdRaw->lat)
            delete[] hsdRaw->lat;
        if (hsdRaw->lon)
            delete[] hsdRaw->lon;
        if (hsdRaw->phys)
            delete[] hsdRaw->phys;

        *hsdRaw = HSDRAWDATA();
    }
}

BOOL HimawariDataExtractor::hsdConvert2netCDF(const std::string& outPath) const
{
    //NOT HSD or NOT load
    if (HimawariDataMode != 1)
        return HDE_INVALID_HSD;

    //Att: title                satName
    //Att: time_created         start_time convert to YYYY-MM-DDTHH:MM:SSZ (UTC)
    //Var: albedo/tbb           phys

    std::string phys_unit_name = "Invalid", phys_unit = "Invalid";

    if (hsdRawData.band <= 6 && hsdRawData.band > 0)
    {
        phys_unit_name = "reflectivity";
        phys_unit = "1";
    }
    else if (hsdRawData.band <= 16 && hsdRawData.band > 6)
    {
        phys_unit_name = "brightness_temperature";
        phys_unit = "K";
    }

    //Def netCDF---------------------------------------------------------------------------------------------------
 
    //New netCDF file
    netCDF::NcFile hsd2nc;
    hsd2nc.create(outPath, netCDF::NcFile::replace);

    //Def X,Y Dimmensions
    auto dim_latitude = hsd2nc.addDim(std::string("latitude"), hsdRawData.height);//y height
    auto dim_longitude = hsd2nc.addDim(std::string("longitude"), hsdRawData.width);//x width

    //Def start_time Dimmension
    auto dim_stime = hsd2nc.addDim(std::string("start_time"), 1);//obsvervation start time

    //Def X,Y Values
    auto var_latitude = hsd2nc.addVar(std::string("latitude"), netCDF::NcEnumType::nc_FLOAT, dim_latitude);
    auto var_longitude = hsd2nc.addVar(std::string("longitude"), netCDF::NcEnumType::nc_FLOAT, dim_longitude);
    var_latitude.putAtt(std::string("units"), std::string("degrees_north"));
    var_longitude.putAtt(std::string("units"), std::string("degrees_east"));

    //Def start_time Value
    auto var_stime = hsd2nc.addVar(std::string("start_time"), netCDF::NcEnumType::nc_DOUBLE, dim_stime);

    //Def Phys value
    std::vector<netCDF::NcDim> dims_latlon = { dim_latitude ,dim_longitude };
    auto var_phys = hsd2nc.addVar(hsdRawData.layerName, netCDF::NcEnumType::nc_FLOAT, dims_latlon);
    var_phys.putAtt(phys_unit_name, phys_unit);

    //Set "title" Attr
    hsd2nc.putAtt(std::string("title"), std::string(hsdRawData.satName));
    
    //Set "date_created" Attr
    //std::string date_created = GetDateStrByTime(hsdRawData.startTime);
    //hsd2nc.putAtt(std::string("date_created"), date_created);

    //Set "Institution" Attr
    hsd2nc.putAtt(std::string("Institution"), std::string("MSC/JMA"));
    
    //Set "Conventions" Attr
    hsd2nc.putAtt(std::string("Conventions"), std::string("Himawari_Data_Extractor"));

    //Add data to netCDF--------------------------------------------------------------------------------------------

    //Add start_time var
    var_stime.putVar(&hsdRawData.startTime);

    //Add latitude vars
    var_latitude.putVar(hsdRawData.lat);

    //Add longitude vars
    var_longitude.putVar(hsdRawData.lon);

    //Add phys vars
    std::vector<size_t> startPos = { 0,0 };
    std::vector<size_t> gridSize = { (size_t)hsdRawData.height,(size_t)hsdRawData.width };
    var_phys.putVar(startPos, gridSize, hsdRawData.phys);

    //End-----------------------------------------------------------------------------------------------------------
    hsd2nc.close();

    return HDE_SUCCESS;
}

BOOL HimawariDataExtractor::hsdGenerateBMPFormCurrent(const std::string& output, bool RePhase)
{
    //NOT HSD or NOT load
    if (HimawariDataMode != 1)
        return HDE_INVALID_HSD;

    //Bitmap 8bpp grayscale array
    unsigned char* bmp_raw = nullptr;

    //alloc
    int n_size = hsdRawData.height * hsdRawData.width;
    bmp_raw = new unsigned char[n_size];
    memset(bmp_raw, 0, n_size);

    //normalization
    auto minVal = *std::min_element(hsdRawData.phys, hsdRawData.phys + n_size - 1);
    auto maxVal = *std::max_element(hsdRawData.phys, hsdRawData.phys + n_size - 1);

    // Normalize the pixel values
    for (int i = 0; i < n_size; i++)
    {
        if (hsdRawData.phys[i] == -1)
            continue;
        auto sub = RePhase ? hsdRawData.phys[i] - minVal : maxVal - hsdRawData.phys[i];
        bmp_raw[i] = static_cast<unsigned char>(sub * 255.0 / (maxVal - minVal));
    }

    //Output bmp
    BOOL ret = WriteBMPFile(output, bmp_raw, hsdRawData.width, hsdRawData.height, 8);

    //Free memory
    delete[] bmp_raw; bmp_raw = nullptr;

    return ret;
}

BOOL HimawariDataExtractor::WriteBMPFile(const std::string& filename, unsigned char* pData, int biWidth, int biHeight, int bit)
{
    FILE* pf = nullptr;

    fopen_s(&pf, filename.c_str(), "wb");
    //pf = fopen(filename, "wb");
    if (!pf)
    {
        std::cerr << "Open file failed on ncWriteBMPFile.\n";
        return HDE_OPENFAILED;
    }

    //each line Must be multiplie of 4
    long bt = (biWidth * bit / 8) % 4;
    if (bt != 0)
    {
        bt = 4 - bt;
    }

    long lineByte = (biWidth * bit / 8) + bt;

    //Grayscale 8bit
    int colorTablesize = 0;
    if (bit == 8)
        colorTablesize = 1024;

    //BMP file header
    BITMAPFILEHEADER bitMapFileHeader;
    bitMapFileHeader.bfType = 0x4D42;//Windows BMP
    bitMapFileHeader.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + colorTablesize + lineByte * biHeight;
    bitMapFileHeader.bfReserved1 = 0;
    bitMapFileHeader.bfReserved2 = 0;
    bitMapFileHeader.bfOffBits = 54 + colorTablesize;

    //BMP info header
    BITMAPINFOHEADER bitMapInfoHeader;
    bitMapInfoHeader.biBitCount = bit;
    bitMapInfoHeader.biClrImportant = 0;
    bitMapInfoHeader.biClrUsed = 0;
    bitMapInfoHeader.biCompression = 0;
    bitMapInfoHeader.biHeight = biHeight;
    bitMapInfoHeader.biPlanes = 1;
    bitMapInfoHeader.biSize = 40;
    bitMapInfoHeader.biSizeImage = lineByte * biHeight;
    bitMapInfoHeader.biWidth = biWidth;
    bitMapInfoHeader.biXPelsPerMeter = 0;
    bitMapInfoHeader.biYPelsPerMeter = 0;

    //Write header
    fwrite(&bitMapFileHeader, sizeof(BITMAPFILEHEADER), 1, pf);
    //Write info
    fwrite(&bitMapInfoHeader, sizeof(BITMAPINFOHEADER), 1, pf);

    //palette
    if (bit == 8) {
        // Write the color palette (grayscale palette for 8-bit BMP)
        for (int i = 0; i < 256; ++i) {
            uint8_t palette[4] = { static_cast<uint8_t>(i), static_cast<uint8_t>(i), static_cast<uint8_t>(i), 0 };
            fwrite(reinterpret_cast<const char*>(palette), 4, 1, pf);
        }
    }

    //for assignment
    unsigned char null[1] = { 0 };

    //Bit
    int b = bit / 8;

    //write raw
    for (int y = biHeight - 1; y >= 0; --y) {
        fwrite(&pData[y * biWidth * b], sizeof(unsigned char), biWidth * b, pf);
        for (int k = 0; k < bt; ++k)
        {
            fwrite(null, sizeof(unsigned char), 1, pf);
        }
    }

    fclose(pf);
    return HDE_SUCCESS;
}

BOOL HimawariDataExtractor::ncGenerateBMPFromLayer(const std::string& var, const std::string& outPath, bool RePhase)
{
    //NOT netCDF or NOT load
    if (HimawariDataMode != 0)
        return HDE_INVALID_NETCDF;

    //Layers will be converted
    std::vector<std::string> layers;

    //BGR ?
    int bit = 1; //8 bpp
    if ((int)var.find("RGB", 0) == 0)
    {
        layers.emplace_back(std::string("albedo_01")); //B
        layers.emplace_back(std::string("albedo_02")); //G
        layers.emplace_back(std::string("albedo_03")); //R
        bit = 3; //24 bpp
    }
    else
        layers.emplace_back(var);

    //raw array & bitmap array & size
    int* layer_raw = nullptr;
    unsigned char* bmp_raw = nullptr;
    size_t x = 0, y = 0;

    int index = 0; //start from B
    for (const auto& layer : layers)
    {
        //Get nc raw data
        auto data = ncInstance.getVar(layer);

        //nodata or invalid
        if (data.isNull())
        {
            return HDE_NETCDF_NODATA;
        }

        //first initlizate
        if (index == 0)
        {
            //get nc data size
            auto dims = data.getDims();
            y = dims[0].getSize();
            x = dims[1].getSize();

            //alloc memory
            layer_raw = new int[x * y];
            bmp_raw = new unsigned char[x * y * bit];
        }

        //get data
        data.getVar(layer_raw);

        //normalization
        int minVal = *std::min_element(layer_raw, layer_raw + x * y - 1);
        int maxVal = *std::max_element(layer_raw, layer_raw + x * y - 1);

        // Normalize the pixel values
        for (int i = 0; i < x * y; ++i)
        {
            int sub = RePhase ? layer_raw[i] - minVal : maxVal - layer_raw[i];
            bmp_raw[i * bit + index] = static_cast<unsigned char>(sub * 255.0 / (maxVal - minVal));
        }

        index++;//Next color
    }

    //Output bmp
    BOOL ret =  WriteBMPFile(outPath, bmp_raw, x, y, bit * 8);

    //Free memory
    delete[] layer_raw; layer_raw = nullptr;
    delete[] bmp_raw; bmp_raw = nullptr;

    return ret;
}

BOOL HimawariDataExtractor::ncGetLayerRAW(const std::string& layer, NCLAYERDATA& ncLayerRawOut)
{
    //NOT netCDF or NOT load
    if (HimawariDataMode != 0)
        return HDE_INVALID_NETCDF;

    //Clean first
    ncReleaseLayerData(&ncLayerRawOut);

    //Get nc raw data
    auto data = ncInstance.getVar(layer);

    //nodata or invalid
    if (data.isNull())
    {
        return HDE_NETCDF_NODATA;
    }

    ncLayerRawOut.layer = layer;

    //Get band & wave length by layer name
    GetBandInfoByLayerName(layer, ncLayerRawOut.band, ncLayerRawOut.bandWaveLength);

    //get nc data size
    auto dims = data.getDims();
    ncLayerRawOut.height = dims[0].getSize();
    ncLayerRawOut.width = dims[1].getSize();

    //alloc memory
    ncLayerRawOut.phys = new float[ncLayerRawOut.height * ncLayerRawOut.width];

    //get data
    data.getVar(ncLayerRawOut.phys);

    return HDE_SUCCESS;
}

void HimawariDataExtractor::ncReleaseLayerData(NCLAYERDATA* ncLayer)
{
    if (ncLayer)
    {
        //if (ncLayer->lat)
        //   delete[] ncLayer->lat;
        //if (ncLayer->lon)
        //    delete[] ncLayer->lon;
        if (ncLayer->phys)
            delete[] ncLayer->phys;

        *ncLayer = NCLAYERDATA();
    }
}

void HimawariDataExtractor::ncFreeAllResources()
{
    //Close current netCDF instance
    ncInstance.close();

    //Clear base info
    ncInfo.date = "";
    ncInfo.satellite = "";
    ncInfo.layers.clear();

    //Reset to UnLoad
    HimawariDataMode = -1;
}
