#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <iostream>
#include <string>
#include <chrono>

class Progressbar {
    public:
        // debug param
        std::size_t propagated_steps;

        // Constructor and Destructor
        Progressbar();
        Progressbar(const std::size_t max_step, const std::string preface="\b");
        Progressbar(const std::size_t max_step, const std::size_t state);
        ~Progressbar();

        // methods
        //! increase current_step by one
        void next();
        //! finish progressbar
        void close();

//    private:
        // variables
        std::size_t max_step;      // maximal number of steps
        std::string preface;
        int progress;              // in precision
        std::size_t current_step;  // current position
        double precision = 1000000.;  // 100. = 1%, 1000. = 0.1% and so on

        // times
        std::chrono::steady_clock::time_point t_start, t_curr;

        void _print_progress();
        //! return ETA time in seconds
        int _get_ETA_in_second();
        //! return elapsed time in seconds
        int _get_ET_in_second();
        //! return elapsed time in milliseconds
        float _get_ET_in_ms();
        std::string _print_time(int);

        void hide_cursor(){std::cout << "\e[?25l";}
        void show_cursor(){std::cout << "\e[?25h";}
};

class Progressspinner : public Progressbar {
    public:
        Progressspinner(const std::string preface, const int refreshing_rate=1);
        Progressspinner(const std::size_t state, const int refreshing_rate=1);
        ~Progressspinner();
        // methods
        //! increase current_step by one
        void next();
        //! finish progressbar
        void close();
    private:
        int refreshing_rate = 1;
};
#endif
