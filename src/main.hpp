#ifndef TYPE_STATE
#define TYPE_STATE

#include <cstdint>
typedef std::uint16_t State;

#endif

#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include "pathfinder.hpp"
#include "transition_matrix.hpp"
#include "progressbar.hpp"

#include <string>
#include <vector>
#include <cmath>

//#include <exception>

#include <boost/program_options.hpp>
//#include <omp.h>

enum Mode {PATHS, MCMC, MCMCSINGLE};
typedef std::map<State, double> Weights;

//! check if states are legit values and unique
bool check_states(const std::vector<State> A, const std::string str_A,
                  const std::vector<State> B, const std::string str_B,
                  const std::vector<State> C, const std::string str_C,
                  const State no_of_states);

//! check if states in A are legit
bool check_states_A(const std::vector<State> A, const std::string str_A,
                    const std::vector<State> B, const std::string str_B,
                    const std::vector<State> C, const std::string str_C,
                    const State no_of_states);

#endif
