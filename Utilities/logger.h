#ifndef LOGGER_H
#define LOGGER_H

#include <ctime>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    // call before first getInstance()
    static void setLogDirectory(const std::string& dir) {
        logDirOverride_ = dir;
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

    operator std::ostream&() {
        return logFile_;
    }

private:
    std::ofstream logFile_;
    std::string logFilePath_;
    static std::string logDirOverride_;  // declared here

    Logger() {
        std::string dir = logDirOverride_.empty() ? LOG_DIR : logDirOverride_;
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << dir << "/log_" << std::put_time(&tm, "%Y-%m-%d_%H:%M:%S") << ".txt";
        logFilePath_ = oss.str();
        logFile_.open(logFilePath_, std::ios::out | std::ios::trunc);
    }
    ~Logger() {
        std::cout << "Log file saved to: " << logFilePath_ << std::endl;
        logFile_.close();
    }

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
};

#define LOG Logger::getInstance()
#define LOG_INFO Logger::getInstance() << "[INFO]  \t"
#define LOG_DEBUG Logger::getInstance() << "[DEBUG] \t"
#define LOG_WARNING Logger::getInstance() << "[WARNING]  \t"
#define LOG_ERROR Logger::getInstance() << "[ERROR] \t"

#endif // LOGGER_H