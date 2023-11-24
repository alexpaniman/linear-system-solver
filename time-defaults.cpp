#include "time-defaults.h"

#include <iomanip>
#include <sstream>
#include <chrono>

std::string describe_duration(duration range, bool with_seconds) {
    auto time_passed_seconds =
        std::chrono::duration_cast<std::chrono::seconds>(range).count();

    std::ostringstream oss;

    int   hours = time_passed_seconds / 3600;
    int minutes = time_passed_seconds % 3600 / 60;

    oss << hours << "h " << minutes << "m";

    if (with_seconds) {
	int seconds = time_passed_seconds % 60;
	oss << " " << seconds << "s";
    }

    return oss.str();
}

std::string describe_passed_time_since(time_point point) {
    return describe_duration(std::chrono::system_clock::now() - point);
}

std::string describe_time_point(time_point point, const char *format) {
    std::time_t time = std::chrono::system_clock::to_time_t(point);
    auto tm = std::localtime(&time);

    std::ostringstream oss;

    oss << std::put_time(tm, format);
    return oss.str();
}

time_point parse_time_stamp(std::stringstream &ss, const char *format) {
    std::tm tm = {};
    ss >> std::get_time(&tm, format);

    return std::chrono::system_clock::from_time_t(std::mktime(&tm));
}

