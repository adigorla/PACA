#ifndef PACA_BASE_UTILS_H
#define PACA_BASE_UTILS_H

#include "pacabase.h"

class Logger {
public:
    static void SetVerbosity(int level);

    template<typename... Args>
    static void LogWARN(Args... args);

    template<typename... Args>
    static void LogINFO(Args... args);

    template<typename... Args>
    static void LogLOG(Args... args);

    template<typename... Args>
    static void LogDEBUG(Args... args);

private:
    template<typename T>
    static void concatenate_args(std::ostringstream& oss, const T& arg);

    template<typename T, typename... Args>
    static void concatenate_args(std::ostringstream& oss, const T& arg, const Args&... args);

    template<typename... Args>
    static void Log(int level, const std::string &prefix, const Args&... args);

    static std::mutex mtx;
    static int verbosity;
};

class Timer {
public:
    explicit Timer(const std::string& fname);
    ~Timer();

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_startTimept;
    std::string m_fname;

    void Stop();
};

// Template function definitions
#include "utils.inl"

Eigen::MatrixXd normalizeCPP(Eigen::MatrixXd& x, bool inplace = true);

Eigen::MatrixXd centerCPP(const Eigen::MatrixXd& x);

Eigen::MatrixXd covCPP(const Eigen::MatrixXd& x); 

Eigen::MatrixXd covSymCPP(const Eigen::MatrixXd& x); 

Eigen::VectorXd colwiseVars(const Eigen::MatrixXd& matrix);

#endif // PACA_UTILS_H