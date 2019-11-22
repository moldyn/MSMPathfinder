#include <iostream>
#include <sstream>
#include <iomanip>
#include "progressbar.hpp"
// Overloaded constructor
Progressbar::Progressbar()
    : max_step(0), preface(""), progress(0), current_step(0)
{}

Progressbar::Progressbar(const std::size_t max_step, const std::string preface)
    : max_step(max_step), preface(preface), progress(0), current_step(0)
{
    this->hide_cursor();
    t_start = std::chrono::steady_clock::now();
    this->_print_progress();
}
Progressbar::Progressbar(const std::size_t max_step, const std::size_t state)
    : max_step(max_step), progress(0), current_step(0)
{
    this->hide_cursor();
    std::ostringstream preface_ostring;
    preface_ostring << "state " << std::setw(3) << state;
    preface = preface_ostring.str();

    t_start = std::chrono::steady_clock::now();
    this->_print_progress();
}

// Destructor
Progressbar::~Progressbar()
{
    this->show_cursor();
}

void Progressbar::next() {
    current_step++;
    if (this->precision*current_step/max_step >= (unsigned) progress + 1) {
        progress = this->precision*current_step/max_step;
        t_curr = std::chrono::steady_clock::now();
        this->_print_progress();
    }
    return;
}


void Progressbar::_print_progress(){
    int barWidth = 20;
    // clear line
    std::cout << "\r                                                            "
              << "\r    " << this->preface << " [";
    int pos = barWidth * progress/this->precision;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else std::cout << "-";
    }
    {
        std::cout << "] "
                  << std::fixed << std::setprecision(1) << std::setw(5) << progress/10.
                  << " %";
        std::cout << std::setprecision(-1);
    }
    if ((progress > 0) && (progress < this->precision)){
        std::cout << " - ETA " << _print_time(this->_get_ETA_in_second());
    } else if (progress == this->precision) {
        std::cout << " - " << _print_time(this->_get_ET_in_second());
    }
    std::cout.flush();
    return;
}

int Progressbar::_get_ETA_in_second(){
    int t_curr = std::chrono::duration_cast<std::chrono::milliseconds>(this->t_curr -
                                                                       this->t_start).count();
    int remaining = this->precision - progress;
    return t_curr*remaining/(this->precision*progress);

}

int Progressbar::_get_ET_in_second(){
    int t_elapsed = std::chrono::duration_cast<std::chrono::seconds>(this->t_curr -
                                                                     this->t_start).count();
    return t_elapsed;
}

std::string Progressbar::_print_time(int time){
    std::ostringstream time_ostring;
    int min = time / 60;
    int sec = time % 60;
    if (min > 0) {
        time_ostring << min << "min ";
    }
    time_ostring << sec << "s";

    return time_ostring.str();
}

void Progressbar::close(){
    std::cout << std::endl;
    this->show_cursor();
    return;
}

// Progress Spinner
Progressspinner::Progressspinner(const std::string preface, const int refreshing_rate)
    : refreshing_rate(refreshing_rate)
{
    this->hide_cursor();
    this->preface = preface;
    t_start = std::chrono::steady_clock::now();
}

Progressspinner::Progressspinner(const std::size_t state, const int refreshing_rate)
    : refreshing_rate(refreshing_rate)
{
    this->hide_cursor();
    std::ostringstream preface_ostring;
    preface_ostring << "    state " << std::setw(3) << state;
    preface = preface_ostring.str();
    t_start = std::chrono::steady_clock::now();
}

Progressspinner::~Progressspinner()
{
    this->show_cursor();
}

void Progressspinner::next() {
    current_step++;
    if (current_step % refreshing_rate == 0) {
        std::cout << "\r" << preface << "... ";
        switch ((current_step/refreshing_rate)%8) {
            case 0:
                std::cout << "⡟";
                break;
            case 1:
                std::cout << "⣏";
                break;
            case 2:
                std::cout << "⣧";
                break;
            case 3:
                std::cout << "⣶";
                break;
            case 4:
                std::cout << "⣼";
                break;
            case 5:
                std::cout << "⣹";
                break;
            case 6:
                std::cout << "⢻";
                break;
            case 7:
                std::cout << "⠿";
                break;
        }
        std::cout.flush();
    }
    return;
}

void Progressspinner::close() {
    t_curr = std::chrono::steady_clock::now();
    std::cout << "\r" << preface << "... done in "
              << _print_time(this->_get_ET_in_second())
              << std::endl;
    this->show_cursor();
    return;
}
