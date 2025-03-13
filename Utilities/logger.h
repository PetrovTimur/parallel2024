#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <ostream>

enum class LogLevel {
    INFO, DEBUG, WARNING, ERROR
};

class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    template<typename T>
    Logger& operator<<(const T& data) {
        logFile_ << data;
        return *this;
    }

    Logger& operator<<(std::ostream& (*manip)(std::ostream&)) {
        logFile_ << manip;
        return *this;
    }

private:
    std::ofstream logFile_;

    Logger() {
        logFile_.open("logfile.txt", std::ios::out | std::ios::trunc);
        // logFile_ << "OPEN";
    }
    ~Logger() {
        logFile_.close();
    }

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
};

#define LOG_INFO Logger::getInstance() << "[INFO] "
#define LOG_DEBUG Logger::getInstance() << "[DEBUG] "

#endif // LOGGER_H