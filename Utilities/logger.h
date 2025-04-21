#ifndef LOGGER_H
#define LOGGER_H

#include <ctime>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <ostream>
#include <string>

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
    std::string logFilePath_;

    Logger() {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << LOG_DIR << "/log_" << std::put_time(&tm, "%Y-%m-%d_%H:%M:%S") << ".txt";
        logFilePath_ = oss.str();
        logFile_.open(oss.str(), std::ios::out | std::ios::trunc);
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