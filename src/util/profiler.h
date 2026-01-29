#pragma once

#ifdef AFTERGLOW_PROFILE

    #include <chrono>
    #include <string>
    #include <unordered_map>

namespace afterglow {

    class ProfileData {
      public:
        void reset() { timings_.clear(); }

        void add(const std::string& name, double ms) { timings_[name] += ms; }

        auto results() const -> std::unordered_map<std::string, double> { return timings_; }

      private:
        std::unordered_map<std::string, double> timings_;
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

#else

    #define AFTERGLOW_PROFILE_SCOPE(name)
    #define AFTERGLOW_PROFILE_RESET()
    #define AFTERGLOW_PROFILE_RESULTS()                                                                                \
        std::unordered_map<std::string, double> {}

#endif
