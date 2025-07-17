#ifndef OPT_LOG_H
#define OPT_LOG_H

#if defined(DISABLE_LOGGING) || !defined(WITH_FMT)

    #define LOG(...)              do {} while(0)
    #define LOG_T(...)            do {} while(0)
    #define LOG_PREFIX(...)       do {} while(0)
    #define LOG_SUCCESS(...)      do {} while(0)
    #define LOG_WARNING(...)      do {} while(0)
    #define LOG_ERROR(...)        do {} while(0)
    #define LOG_START_MODULE(...) do {} while(0)
    #define LOG_ROW(...)          do {} while(0)
    #define LOG_HEADER(...)       do {} while(0)
    #define LOG_DASHES(...)       do {} while(0)

#else

    #include <fmt/core.h>
    #include <fmt/format.h>
    #include <fmt/color.h>
    #include <string>
    #include <array>

    // ---------------------------
    // General Log Macros
    // ---------------------------

    #define LOG(...) \
        do { fmt::print(__VA_ARGS__); fmt::print("\n"); } while (0)

    #define LOG_T(tabs, ...) \
        do { fmt::print("{0:{1}}", "", (tabs) * 4); fmt::print(__VA_ARGS__); } while (0)

    #define LOG_PREFIX(c, ...) \
        do { fmt::print("{} ", c); fmt::print(__VA_ARGS__); } while (0)

    #define LOG_SUCCESS(...) \
        do { \
            fmt::print(fmt::fg(fmt::rgb(110, 168, 1)) | fmt::emphasis::bold, "\nSUCCESS"); \
            fmt::print(" - "); fmt::print(__VA_ARGS__); fmt::print("\n"); \
        } while (0)

    #define LOG_WARNING(...) \
        do { \
            fmt::print(fmt::fg(fmt::rgb(199, 171, 15)) | fmt::emphasis::bold, "\nWARNING"); \
            fmt::print(" - "); fmt::print(__VA_ARGS__); fmt::print("\n"); \
        } while (0)

    #define LOG_ERROR(...) \
        do { \
            fmt::print(fmt::fg(fmt::rgb(211, 32, 86)) | fmt::emphasis::bold, "\nERROR"); \
            fmt::print(" - "); fmt::print(__VA_ARGS__); fmt::print("\n"); \
        } while (0)

    #define LOG_START_MODULE(...) \
        do { \
            fmt::print("\n=== "); \
            fmt::print(fmt::emphasis::bold, __VA_ARGS__); \
            fmt::print(" ===\n\n"); \
        } while (0)

    // ---------------------------
    // Table Logging Support
    // ---------------------------

    template <size_t N>
    struct FixedTableFormat {
        std::array<int, N> col_widths;
        std::array<std::string, N> fmt_strings;  // Precompiled "{:>w}" strings

        constexpr FixedTableFormat(const std::array<int, N>& widths)
            : col_widths(widths) {
            for (size_t i = 0; i < N; ++i) {
                fmt_strings[i] = fmt::format("{{:>{}}}", col_widths[i]);
            }
        }

        constexpr int total_width() const {
            int total = -1;
            for (auto w : col_widths)
                total += w + 3;  // +3 for ' | '
            return total;
        }
    };

    template <size_t N>
    constexpr inline void log_dashes_fast(const FixedTableFormat<N>& fmt) {
        fmt::print("{:-<{}}\n", "", fmt.total_width());
    }

    template <size_t N, typename... Args>
    constexpr inline void log_row_fast(const FixedTableFormat<N>& fmt, Args&&... args) {
        static_assert(sizeof...(Args) == N, "Number of columns must match format definition.");
        size_t i = 0;
        const char* sep = "";
        ((fmt::print("{}{}", sep, fmt::format(fmt.fmt_strings[i++], args)), sep = " | "), ...);
        fmt::print("\n");
    }

    #define LOG_ROW(fmt, ...)       do { log_row_fast(fmt, __VA_ARGS__); } while (0)
    #define LOG_HEADER(fmt, ...)    do { log_row_fast(fmt, __VA_ARGS__); } while (0)
    #define LOG_DASHES(fmt)         do { log_dashes_fast(fmt); }           while (0)

#endif

#endif // OPT_LOG_H
