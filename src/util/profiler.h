//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#ifdef AFTERGLOW_PROFILE

#include <chrono>
#include <string>
#include <unordered_map>

namespace afterglow {
namespace profiling {

/**
 * @brief Simple profiler for per-stage CPU timing
 *
 * When AFTERGLOW_PROFILE is defined, this provides scoped timing of code sections.
 * Used for performance analysis during validation and benchmarking.
 */
class Profiler {
public:
    static Profiler& instance() {
        static Profiler inst;
        return inst;
    }

    void start(const std::string& name) {
        start_times_[name] = std::chrono::high_resolution_clock::now();
    }

    void stop(const std::string& name) {
        auto end = std::chrono::high_resolution_clock::now();
        auto it = start_times_.find(name);
        if (it != start_times_.end()) {
            auto duration = std::chrono::duration<double>(end - it->second).count();
            times_[name] += duration;
            start_times_.erase(it);
        }
    }

    std::unordered_map<std::string, double> results() const {
        return times_;
    }

    void reset() {
        times_.clear();
        start_times_.clear();
    }

private:
    Profiler() = default;
    std::unordered_map<std::string, double> times_;
    std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> start_times_;
};

/**
 * @brief RAII scoped profiler
 */
class ScopedProfile {
public:
    explicit ScopedProfile(const std::string& name) : name_(name) {
        Profiler::instance().start(name_);
    }

    ~ScopedProfile() {
        Profiler::instance().stop(name_);
    }

private:
    std::string name_;
};

} // namespace profiling
} // namespace afterglow

// Macros for profiling (active when AFTERGLOW_PROFILE is defined)
#define AFTERGLOW_PROFILE_SCOPE(name) afterglow::profiling::ScopedProfile _profile_##name(#name)
#define AFTERGLOW_PROFILE_RESULTS() afterglow::profiling::Profiler::instance().results()
#define AFTERGLOW_PROFILE_RESET() afterglow::profiling::Profiler::instance().reset()

#endif // AFTERGLOW_PROFILE
