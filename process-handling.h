#pragma once

#include "time-defaults.h"

#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>


using pipe_fileno = int;

struct process {
    pid_t id;
    pipe_fileno stdin, stdout;
};

// Some processes output colorful text filled with ansi escape-codes,
// it's very hard to process such text, so to avoid the hassle of
// parsing colors this function can be used to remove them  at first:
std::string filter_out_ansi_colors(std::string colored_text);

process create_process(const char* path, std::vector<const char*> args);
std::string communicate_command(process proc, std::string command = "", duration timeout = 0.15s);

