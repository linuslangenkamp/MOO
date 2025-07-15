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
    #include <fmt/color.h>

    #define LOG(...) \
        do { \
            fmt::print(__VA_ARGS__); \
            fmt::print("\n"); \
        } while(0)

    #define LOG_T(tabs, ...) \
        do { \
            fmt::print("{0:{1}}", "", (tabs) * 4); \
            fmt::print(__VA_ARGS__); \
        } while(0)

    #define LOG_PREFIX(c, ...) \
        do { \
            fmt::print("{} ", c); \
            fmt::print(__VA_ARGS__); \
        } while(0)

    #define LOG_SUCCESS(...) \
        do { \
            fmt::print(fmt::fg(fmt::rgb(133, 153, 0)) | fmt::emphasis::bold, "\nSUCCESS"); \
            fmt::print(" - "); \
            fmt::print(__VA_ARGS__); \
            fmt::print("\n"); \
        } while(0)

    #define LOG_WARNING(...) \
        do { \
            fmt::print(fmt::fg(fmt::rgb(181, 137, 0)) | fmt::emphasis::bold, "\nWARNING"); \
            fmt::print(" - "); \
            fmt::print(__VA_ARGS__); \
            fmt::print("\n"); \
        } while(0)

    #define LOG_ERROR(...) \
        do { \
            fmt::print(fmt::fg(fmt::rgb(220, 50, 47)) | fmt::emphasis::bold, "\nERROR"); \
            fmt::print(" - "); \
            fmt::print(__VA_ARGS__); \
            fmt::print("\n"); \
        } while(0)

    #define LOG_START_MODULE(...) \
        do { \
            fmt::print("\n=== "); \
            fmt::print(fmt::emphasis::bold, __VA_ARGS__); \
            fmt::print(" ===\n\n"); \
        } while(0)

#endif

#endif // OPT_LOG_H
