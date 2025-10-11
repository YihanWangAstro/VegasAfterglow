//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <stdexcept>
#include <string>

namespace afterglow {

    class ValidationError final : public std::invalid_argument {
      public:
        explicit ValidationError(const std::string& message) : std::invalid_argument(message) {}
        explicit ValidationError(const char* message) : std::invalid_argument(message) {}
    };

    class LogicError final : public std::logic_error {
      public:
        explicit LogicError(const std::string& message) : std::logic_error(message) {}
        explicit LogicError(const char* message) : std::logic_error(message) {}
    };

#ifdef NDEBUG
    #define AFTERGLOW_REQUIRE(condition, message)                                                                      \
        do {                                                                                                           \
            if (!(condition)) [[unlikely]] {                                                                           \
                throw ::afterglow::ValidationError(message);                                                           \
            }                                                                                                          \
        } while (0)

    #define AFTERGLOW_ENSURE(condition, message)                                                                       \
        do {                                                                                                           \
            if (!(condition)) [[unlikely]] {                                                                           \
                throw ::afterglow::LogicError(message);                                                                \
            }                                                                                                          \
        } while (0)
#else
    #define AFTERGLOW_REQUIRE(condition, message)                                                                      \
        do {                                                                                                           \
            if (!(condition)) [[unlikely]] {                                                                           \
                throw ::afterglow::ValidationError(message);                                                           \
            }                                                                                                          \
        } while (0)

    #define AFTERGLOW_ENSURE(condition, message)                                                                       \
        do {                                                                                                           \
            if (!(condition)) [[unlikely]] {                                                                           \
                throw ::afterglow::LogicError(message);                                                                \
            }                                                                                                          \
        } while (0)
#endif

} // namespace afterglow
