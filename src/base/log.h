#ifndef OPT_LOG_H
#define OPT_LOG_H

#if defined(DISABLE_LOGGING) || !defined(WITH_FMT)

    #define LOG(...)         do {} while(0)
    #define LOG_T(...)       do {} while(0)
    #define LOG_PREFIX(...)  do {} while(0)

    #define LOG_WARNING(...) do {} while(0)
    #define LOG_ERROR(...)   do {} while(0)

#else

    #include <fmt/core.h>

    #define LOG(...) \
        do { fmt::println(__VA_ARGS__); } while(0)

    #define LOG_T(tabs, ...) \
        do { \
            fmt::print("{0:{1}}", "", (tabs) * 4); \
            fmt::println(__VA_ARGS__); \
        } while(0)

    #define LOG_PREFIX(c, ...) \
        do { \
            fmt::print("{} ", c); \
            fmt::println(__VA_ARGS__); \
        } while(0)

    #define LOG_WARNING(...) \
        do { \
            fmt::print("WARNING - "); \
            fmt::println(__VA_ARGS__); \
        } while(0)

    #define LOG_ERROR(...) \
        do { \
            fmt::print("ERROR - "); \
            fmt::println(__VA_ARGS__); \
        } while(0)

#endif

#endif // OPT_LOG_H
