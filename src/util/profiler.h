#pragma once

#ifdef AFTERGLOW_PROFILE

    #include <chrono>
    #include <string>
    #include <unordered_map>

namespace afterglow {

    class ProfileData {
      public:
        void reset() {
            timings_.clear();
            counters_.clear();
        }

        void add(const std::string& name, double ms) { timings_[name] += ms; }

        void increment(const std::string& name, size_t count = 1) { counters_[name] += count; }

        void set_max(const std::string& name, size_t value) {
            auto it = counters_.find(name);
            if (it == counters_.end() || value > it->second) {
                counters_[name] = value;
            }
        }

        auto results() const -> std::unordered_map<std::string, double> { return timings_; }

        auto counter_results() const -> std::unordered_map<std::string, size_t> { return counters_; }

      private:
        std::unordered_map<std::string, double> timings_;
        std::unordered_map<std::string, size_t> counters_;
    };

    inline auto profiler() -> ProfileData& {
        static ProfileData instance;
        return instance;
    }

    class ScopedTimer {
      public:
        explicit ScopedTimer(const char* name) : name_(name), start_(std::chrono::high_resolution_clock::now()) {}

        ~ScopedTimer() {
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start_).count();
            profiler().add(name_, ms);
        }

        ScopedTimer(const ScopedTimer&) = delete;
        ScopedTimer& operator=(const ScopedTimer&) = delete;

      private:
        const char* name_;
        std::chrono::high_resolution_clock::time_point start_;
    };

} // namespace afterglow

    #define AFTERGLOW_PROFILE_SCOPE(name) afterglow::ScopedTimer _profiler_timer_##name(#name)
    #define AFTERGLOW_PROFILE_RESET() afterglow::profiler().reset()
    #define AFTERGLOW_PROFILE_RESULTS() afterglow::profiler().results()
    #define AFTERGLOW_PROFILE_COUNT(name) afterglow::profiler().increment(#name)
    #define AFTERGLOW_PROFILE_COUNT_N(name, n) afterglow::profiler().increment(#name, n)
    #define AFTERGLOW_PROFILE_MAX(name, n) afterglow::profiler().set_max(#name, n)
    #define AFTERGLOW_PROFILE_COUNTERS() afterglow::profiler().counter_results()

#else

    #define AFTERGLOW_PROFILE_SCOPE(name)
    #define AFTERGLOW_PROFILE_RESET()
    #define AFTERGLOW_PROFILE_RESULTS()                                                                                \
        std::unordered_map<std::string, double> {}
    #define AFTERGLOW_PROFILE_COUNT(name)
    #define AFTERGLOW_PROFILE_COUNT_N(name, n)
    #define AFTERGLOW_PROFILE_MAX(name, n)
    #define AFTERGLOW_PROFILE_COUNTERS()                                                                               \
        std::unordered_map<std::string, size_t> {}

#endif
