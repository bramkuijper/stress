#ifndef AUXILIARY_H_ 
#define AUXILIARY_H_

#include <string>
#include <iostream>
#include <sstream>
#include <chrono>

// convert from integer to string
std::string itos(int j)
{
	std::stringstream s;
	s << j;
	return s.str();
}

// convert double to string (doesn't work though...)
std::string dtos(double j)
{
    std::stringstream s;
    s << j;

    return s.str();
}

int get_nanoseconds()
{
    using namespace std::chrono;

// https://stackoverflow.com/questions/31255486/c-how-do-i-convert-a-stdchronotime-point-to-long-and-back 
    auto now = time_point_cast<nanoseconds>(system_clock::now());

    auto integral_duration = now.time_since_epoch();

    int duration = integral_duration.count();

    return(abs(duration));
}


void create_filename(std::string &filename)
{
	time_t timep = time(NULL);

	tm *timestruct = localtime(&timep);

	filename.append("_");
	filename.append(itos(timestruct->tm_mday));
	filename.append("_");
	filename.append(itos(timestruct->tm_mon + 1));
	filename.append("_");
	filename.append(itos(timestruct->tm_year + 1900));
	filename.append("_");

	size_t hour = timestruct->tm_hour;

	if (hour < 10)
		filename.append("0");

	filename.append(itos(hour));

	size_t minutes = timestruct->tm_min;

	if (minutes < 10)
		filename.append("0");

	filename.append(itos(minutes));

	size_t seconds = timestruct->tm_sec;

	if (seconds < 10)
		filename.append("0");

	filename.append(itos(seconds));

	filename.append("_");

	// now we are going to insert microsecs
	// to obtain unique filename
	filename.append(itos(get_nanoseconds()));
}

int linear_search(const int list[], const int maxsize, int value)
{
    int index = 0;
    int position = -1;
    bool found = false;

    while (index < maxsize && !found)
    {
        if (list[index] == value)
        {
            found = true;
            position = index;
        }

        ++index;
    }

    return(position);
}

#endif
