#include "process-handling.h"

#include <assert.h>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <thread>
#include <unistd.h>


process create_process(const char *path, std::vector<const char*> args) {
    args.insert(args.begin(), path);

    int in_fd[2], out_fd[2];

    pipe(in_fd);  // For child's stdin
    pipe(out_fd); // For child's stdout

    pid_t pid = fork();

    if (pid == 0) { // We're in the child process:
	std::cout << "Launching:";
	for (const auto &arg: args)
	    std::cout << " " << arg;

	close(out_fd[0]);
	dup2(out_fd[1], STDOUT_FILENO);
	close(out_fd[1]);

	close(in_fd[1]);

	dup2(in_fd[0], STDIN_FILENO);
	close(in_fd[0]);

	execvp(path, const_cast<char**>(&*args.begin()));
	_exit(0);
    }

    else if (pid == -1) {
	perror("fork");
	_exit(-1);
    }
    else {
	// You're in the parent
	close(out_fd[1]);
	close(in_fd[0]);

	// Make interactions with pipe non-blocking. 
        int process_out = out_fd[0];

	int flags = fcntl(process_out, F_GETFL, 0);
	fcntl(process_out, F_SETFL, flags | O_NONBLOCK);

        return { pid, in_fd[1], process_out };
    }

    assert(false);
}

std::string filter_out_ansi_colors(std::string colored_text) {
    std::string filtered_text;

    for (int i = 0; i < colored_text.size(); ++ i) {
        if (colored_text[i] == '\033') {
	    while (colored_text[i] != 'm')
		++ i;

	    continue;
        }

	filtered_text += colored_text[i];
    }

    return filtered_text;
}

std::string communicate_command(process proc, std::string command, duration timeout) {
    if (write(proc.stdin, &*command.begin(), command.size()) < 0) {
	perror("write");
	_exit(-1);
    }

    const std::size_t READ_SIZE = 4096;
    char buffer[READ_SIZE + 1] = {};

    std::string response;
    time_point start_time = std::chrono::system_clock::now();

    const auto RESPONSE_WAIT = 0.001s;

    int last_read_bytes = 0;
    while (true) {
	time_point now = std::chrono::system_clock::now();
        if (now - start_time >= timeout)
	    break;

	int read_bytes = read(proc.stdout, buffer, READ_SIZE);
        if (read_bytes < 0) {
	    if (errno != EAGAIN) {
		perror("read");
		_exit(-1);
	    }

	    if (last_read_bytes > 0)
		memset(buffer, '\0', last_read_bytes);
        }

        last_read_bytes = read_bytes;
	buffer[read_bytes] = '\0';

	std::this_thread::sleep_for(RESPONSE_WAIT);
	response += buffer;
    }

    return response;
}
