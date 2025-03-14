#ifndef LOGGER_H
#define LOGGER_H

#include <ctime>
#include <iomanip>
#include <sstream>
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

    // Add this conversion operator
    operator std::ostream&() {
        return logFile_;
    }

private:
    std::ofstream logFile_;

    Logger() {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << LOG_DIR << "/log_" << std::put_time(&tm, "%Y-%m-%d_%H:%M:%S") << ".txt";
        logFile_.open(oss.str(), std::ios::out | std::ios::trunc);
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