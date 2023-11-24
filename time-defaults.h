#pragma once

#include <chrono>


using duration   = std::chrono::duration<double>;
using time_point = std::chrono::time_point<std::chrono::system_clock>;

// Yeah, yeah, I know using namespace in header is bold, but to my
// defence, I'm really tired of using it anyway in most .cpp files
// for writing out, like, a single time literal...
using namespace std::chrono_literals;

std::string describe_duration(duration range, bool with_seconds = false);
std::string describe_passed_time_since(time_point point); 
std::string describe_time_point(time_point point, const char *format = "%d-%m-%Y %H:%M:%S"); 

time_point parse_time_stamp(std::stringstream &ss, const char *format); 

