#include <iostream>
#include <io.h>
#include "HimawariDataExtractor.h"


static void help_info()
{
    std::cout << "Usage: -i <InFile> [-i <InFile2> ...] [-o <OutFileDir>] [...]\n";
    std::cout << "          Input max of 10 HSD and must be in sequence.(DAT files)\n";
    std::cout << "          Input just 1 netCDF data.(nc files)\n";
    std::cout << "       -o [OutFileDir]\n";
    std::cout << "          Default \".\\\"\n";
    std::cout << "       -p Convert inputs to BMP.\n";
    std::cout << "       -c Convert inputs to netCDF format.(HSD only)\n";
    std::cout << "       -h Show this help.\n";

}

static bool isExist(const std::string& in)
{
    return (_access(in.c_str(), 0) == 0);
}

int main(int argc, char** argv)
{

    system("title Himawari Data Extractor");

    std::cout << "Himawari Data Extractor V1.2\n";
    std::cout << "Convert Himawari data (HSD & netCDF) to BMP.\n\n";

	std::string outDir = ".\\";
	std::vector<std::string> inputs;

    bool bmp = false;
    bool convert = false;

    //Get args
	for (int ii = 1; ii < argc; ii++)
	{
		char* ptr = argv[ii];
		if (*ptr == '-')
		{
			ptr++;
			switch (*ptr)
			{
            case 'h':
                help_info();
                return true;
            case 'p':
                bmp = true;
                break;
            case 'c':
                convert = true;
                break;
			case 'o':
                if (ii + 1 >= argc)
                {
                    std::cerr << "[Error]Invalid input.\n"; return false;
                }

				outDir = std::string(argv[ii + 1]);
                if (!isExist(outDir))
                {
                    std::cerr << "[Error]Invalid output dir.\n";
                    return false;
                }
				ii++;
				break;
			case 'i':
                if (ii + 1 >= argc)
                {
                    std::cerr << "[Error]Invalid input.\n"; return false;
                }
                if (!isExist(std::string(argv[ii + 1])))
                {
                    std::cerr << "[Error]Invalid inputs.\n";
                    return false;
                }
				if ((int)inputs.size() < 10) {
					inputs.emplace_back(std::string(argv[ii + 1]));
				}
				else {
					std::cerr << "[Error]Inputs out of 10.\n";
					return false;
				}
				ii++;
				break;
			}
		}
	}

    if (outDir.empty() || inputs.empty())
    {
        help_info();
        return false;
    }

    HimawariDataExtractor hde;

    BOOL ret = hde.LoadHimawariDatas(inputs);

    if (ret != HDE_SUCCESS)
    {
        std::cerr << "Load input(s) failed.\n";
        return false;
    }

    switch (hde.GetCurrentDataType())
    {
    case 0://netCDF
    { 
        auto layers = hde.ncGetAllLayer();

        std::cout << "\nData mode     : netCDF\n";
        std::cout << "Layers count  : " << (int)layers.size() << "\n";
        std::cout << "Satellite     : " << hde.ncGetSatelliteName() << "\n";
        std::cout << "Obs time(UTC) : " << hde.ncGetCreateDateUTC() << "\n\n";

        HimawariDataExtractor::NCLAYERDATA nc;

        std::string outFile, tmpstr;

        for (const auto& layer : layers)
        {
            hde.ncGetLayerRAW(layer, nc);

            if (bmp)
            {
                if (nc.band == 0)
                    tmpstr = "";
                else
                    tmpstr = "_BAND" + std::to_string(nc.band);

                outFile = outDir + "\\netCDF" + tmpstr + "_" + nc.layer + ".bmp";

                std::cout << "Output: " << outFile << "\n";

                hde.ncGenerateBMPFromLayer(layer, outFile);
            }
            else
            {
                std::cout << "Layer: " << layer << " " << nc.bandWaveLength << " nm\n";
            }

            hde.ncReleaseLayerData(&nc);
        }

        if (bmp)
        {
            //Recognize if is himawari standard netCDF data to generate RGB file
            auto it = std::find(layers.begin(), layers.end(), std::string("sd_albedo_03"));
            if (it != layers.end())
            {
                outFile = outDir + "\\netCDF_RGB.bmp";
                std::cout << "Output: " << outFile << "\n";
                hde.ncGenerateBMPFromLayer("RGB", outFile, true);
            }
        }
        break;
    }
    case 1://HSD
    {
        HimawariDataExtractor::HSDRAWDATA hsd;

        hde.hsdGetRaw(hsd);

        std::cout << "\nData mode     : Himawari Standard Data(HSD)\n";
        std::cout << "Layer         : " << hsd.layerName << " " << hsd.bandWaveLength <<" nm\n";
        std::cout << "Satellite     : " << hsd.satName << "\n";
        std::cout << "Obs time(UTC) : " << hde.GetDateStrByTime(hsd.startTime) << "\n\n";

        std::string outFile;

        if (convert)
        {
            outFile = outDir + "\\HSD_BAND" + std::to_string(hsd.band) + "_" + std::string(hsd.satName) + ".nc";
            std::cout << "Convert netCDF: " << outFile << "\n";

            hde.hsdConvert2netCDF(outFile);
        }

        if (bmp)
        {
            outFile = outDir + "\\HSD_BAND" + std::to_string(hsd.band) + "_" + std::string(hsd.satName) + ".bmp";
            std::cout << "Output BMP: " << outFile << "\n";

            hde.hsdGenerateBMPFormCurrent(outFile);
        }

        hde.hsdReleaseRaw(&hsd);
        break;
    }
    }

    hde.ReleaseAll();
    return true;
}

