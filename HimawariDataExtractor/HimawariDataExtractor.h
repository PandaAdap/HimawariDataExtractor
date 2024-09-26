#pragma once

#include <fstream>

#include <windows.h>
#include <map>
#include <algorithm>

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include <netcdf>

extern "C"
{
#include "hisd/hisd.h"
#include "hisd/hisd2netcdf.h"
#include "hisd/date_utl.h"
}

//Global Return Code
constexpr auto HDE_FAILED = 0;
constexpr auto HDE_SUCCESS = 1;
constexpr auto HDE_INVALID_NETCDF = -1;
constexpr auto HDE_INVALID_HSD = -2;
constexpr auto HDE_OPENFAILED = -3;
constexpr auto HDE_NOINPUT = -4;
constexpr auto HDE_INVALID_HSDSEQ = -5;
constexpr auto HDE_HSD_SEQMISSINGPART = -6;
constexpr auto HDE_HSD_SEQINVALIDPART = -7;
constexpr auto HDE_NETCDF_NODATA = -8;

class HimawariDataExtractor
{
public:
	HimawariDataExtractor();
	~HimawariDataExtractor();

	//Load Data
	//Single netCDF4 format data or HSD sequence data
	BOOL LoadHimawariDatas(const std::vector<std::string>& himawariDataFile);

	//Current data type
	//Data type: -1: Invalid/UnLoad, 0: netCDF, 1: HSD(Himawari Standard Data)
	int GetCurrentDataType() const noexcept
	{
		return HimawariDataMode;
	}

	//Get UTC date by MJD (Modified Julian Date)
	std::string GetDateStrByTime(double time);

	//Release all resources
	void ReleaseAll();

	//Global const value---------------------------------------

	const double delta_1km = 0.01;
	const double delta_2km = 0.02;
	const double delta_5km = 0.05;


	//Part of HSD----------------------------------------------

	//Instance for HSD Data
	struct HSDINSTANCE {
		HisdHeader	hsdHeader;
		FILE		*hsdFile = nullptr;
		unsigned short* hsdDataRangeRaw = nullptr;
	};

	//Raw data from input HSDs
	struct HSDRAWDATA {
		int					width = 0;
		int					height = 0;
		float*				lon = nullptr;			//width array(x)
		float*				lat = nullptr;			//height array(y)
		float*				phys = nullptr;			//z
		int					band = 0;				//Observation band number
		float				bandWaveLength = 0;		//Band wave length
		double				startTime = 0;			//Observation starttime (days since 1858-11-17 0:0:0)(MJD)
		double				endTime = 0;			//Observation endtime
		char				satName[32] = { 0 };	//Satellite name
		std::string			layerName = "";			//Layer name
	};

	//Open HSD Data
	//Load direct if total seg num is 1
	//Throw if missing or mismatching seg file
	BOOL hsdOpenDAT(const std::vector<std::string>& inputDATs);

	//Get RAW HSD data
	//Will alloc memory for HSDRAWDATA
	//Remember to call hsdReleaseRaw(...) to free after used.
	BOOL hsdGetRaw(HSDRAWDATA& hsdRawOut);

	//Release HSDRAWDATA
	void hsdReleaseRaw(HSDRAWDATA* hsdRaw);

	//TODO: Convert current HSD data into netCDF format
	BOOL hsdConvert2netCDF(const std::string& outPath) const;

	//Output to BMP
	BOOL hsdGenerateBMPFormCurrent(const std::string& output, bool RePhase = false);

	//Free all HSD resources
	void hsdFreeAllResources();


	//Part of netCDF-------------------------------------------

	//Raw layer data from input netCDF
	struct NCLAYERDATA {

		std::string layer = "";				//Layer key from ncData
		int			band = 0;				//Observation band number
		float		bandWaveLength = 0;		//Band wave length

		int			width = 0;
		int			height = 0;
		//float*		lon = nullptr;			//width array(x)
		//float*		lat = nullptr;			//height array(y)
		float*		phys = nullptr;			//z
	};

	//Base info from nc data
	struct NCBASEINFO {
		std::string satellite = "";			//Origional: title, hsdConverted: satName
		std::string date = "";				//Origional: var_start_time(MJD) hsdConverted: start_time(double)
		std::vector<std::string> layers;	//All layers
	};

	//All layers data
	std::vector<std::string> ncGetAllLayer() const noexcept
	{
		return ncInfo.layers;
	}

	//Current nc satellite name
	std::string ncGetSatelliteName() const noexcept
	{
		return ncInfo.satellite;
	}

	//Current nc create date
	std::string ncGetCreateDateUTC() const noexcept
	{
		return ncInfo.date;
	}

	//Generate bmp from layer
	//Set var to RGB to convert color BMP from albedo_01~03
	//Otherwish grayscale from single layer
	//(Recommand set RePhase to true while generating RGB bmp)
	BOOL ncGenerateBMPFromLayer(const std::string& var, const std::string& outPath, bool RePhase = false);

	//Get RAW layer data from ncFile
	//Will alloc memory for NCLAYERDATA 
	//Remember to call ncReleaseLayer(...) to free after used.
	BOOL ncGetLayerRAW(const std::string& layer, NCLAYERDATA& ncLayerRawOut);

	//Release NCLAYERDATA
	void ncReleaseLayerData(NCLAYERDATA* ncLayer);

	//Free all netCDF resources
	void ncFreeAllResources();

protected:

	//Data type: -1: Invalid/UnLoad, 0: netCDF, 1: HSD(Himawari Standard Data)
	int HimawariDataMode = -1;

	//Himawari Standard Data input sequence
	std::vector<HSDINSTANCE> hsdList;

	//Himawari RAW data struct
	HSDRAWDATA hsdRawData;

	//Mid-file instance
	netCDF::NcFile ncInstance;

	//nc File Data Base info
	NCBASEINFO ncInfo;

	//Process HSD input to RAW data
	BOOL hsdProcessInputData();

	//Free hsd temporary raw
	void hsdReleaseTemporaryRaw()
	{
		for (auto& DAT : hsdList)
		{
			if (DAT.hsdDataRangeRaw)
			{
				delete[] DAT.hsdDataRangeRaw;
				DAT.hsdDataRangeRaw = nullptr;
			}
		}
	}

	//Get delta value for each band
	double GetDeltaValue(int band) const noexcept
	{
		double dt_val = 0.020000;
		switch (band)
		{
		case 1:case 2:case 4:
			dt_val = delta_2km; break;
		case 3:
			dt_val = delta_1km; break;
		default:
			dt_val = delta_5km; break;
		}
		return dt_val;
	}

	//Get layer band info
	BOOL GetBandInfoByLayerName(const std::string layerName, int& bandNumOut, float& bandWidthOut) noexcept
	{
		//Get back 2
		auto backNum = layerName.substr(layerName.length() - 2, 2);

		int band = 0;
		
		//Avoid throw exception
		try
		{
			band = std::stoi(backNum);
		}
		catch (std::exception e)
		{
			return HDE_FAILED;
		}

		bandNumOut = band;

		//Base on official datasheet
		switch (band)
		{
		case 1:
			bandWidthOut = 0.47f; break;
		case 2:
			bandWidthOut = 0.51f; break;
		case 3:
			bandWidthOut = 0.64f; break;
		case 4:
			bandWidthOut = 0.86f; break;
		case 5:
			bandWidthOut = 1.6f; break;
		case 6:
			bandWidthOut = 2.3f; break;
		case 7:
			bandWidthOut = 3.9f; break;
		case 8:
			bandWidthOut = 6.2f; break;
		case 9:
			bandWidthOut = 6.9f; break;
		case 10:
			bandWidthOut = 7.3f; break;
		case 11:
			bandWidthOut = 8.6f; break;
		case 12:
			bandWidthOut = 9.6f; break;
		case 13:
			bandWidthOut = 10.4f; break;
		case 14:
			bandWidthOut = 11.2f; break;
		case 15:
			bandWidthOut = 12.4f; break;
		case 16:
			bandWidthOut = 13.3f; break;
		}

		return HDE_SUCCESS;
	}

	//Convert raw bgr buffer to BMP file
	BOOL WriteBMPFile(const std::string& filename, unsigned char* pData, int biWidth, int biHeight, int bit);

};

