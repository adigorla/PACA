// utils.inl

// Logger template function definitions
template<typename... Args>
void Logger::LogWARN(Args... args) {
    Log(0, "WARN  : ", args...);
}

template<typename... Args>
void Logger::LogINFO(Args... args) {
    Log(1, "INFO  : ", args...);
}

template<typename... Args>
void Logger::LogLOG(Args... args) {
    Log(2, "LOG   : ", args...);
}

template<typename... Args>
void Logger::LogDEBUG(Args... args) {
    Log(3, "DEBUG : ", args...);
}

template<typename T>
void Logger::concatenate_args(std::ostringstream& oss, const T& arg) {
    oss << arg << " ";
}

template<typename T, typename... Args>
void Logger::concatenate_args(std::ostringstream& oss, const T& arg, const Args&... args) {
    oss << arg << " ";
    concatenate_args(oss, args...);
}

template<typename... Args>
void Logger::Log(int level, const std::string &prefix, const Args&... args) {
    std::lock_guard<std::mutex> lock(mtx);
    if (level > verbosity) {
        return;
    }
    const auto in_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::ostringstream oss;
    oss << "[" << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X") << "] " << prefix;

    concatenate_args(oss, args...);

    switch (level) {
    case 0:
        Rcpp::Rcout << "\033[1;33m";  // yellow
        break;
    case 1:
        // No color for INFO
        break;
    case 2:
        Rcpp::Rcout << "\033[1m\033[36m";  // Bold Cyan
        break;
    case 3:
        Rcpp::Rcout << "\033[1m\033[35m";  // Bold Magenta
        break;
    default:
        return;
    }

    Rcpp::Rcout << oss.str() << "\033[0m\n";  // reset color
}