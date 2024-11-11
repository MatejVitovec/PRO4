#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>
#include <fenv.h>
#include <iomanip>
#include <algorithm>

class doublet
{
    public:
        doublet(/* args */);

        double x;
        double p;

    private:

};




int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

	std::ifstream stream;
    std::string line;

    stream.open("fileName", std::ios_base::in);

	std::vector<doublet> data;

    if (stream.is_open())
    {
		while (std::getline(stream, line))
        {
            if(line != "\r")
            {			
                std::stringstream ss(line);
				std::string segment;
				std::vector<std::string> seglist;

				while(std::getline(ss, segment, ' '))
				{
					seglist.push_back(segment);
				}

				x.push_back(std::stod(seglist[0]));
                p.push_back(std::stod(seglist[1]));
            }            
        }
    }

    std::sort()

    return 0;
}