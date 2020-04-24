#include <iostream>
#include <fstream>
#include "main.hpp"
#include "transition_matrix.hpp"

#include <cmath>
#include <algorithm>
// Default constructor
Transition_matrix::Transition_matrix() {}
// Overloaded constructor
Transition_matrix::Transition_matrix(const std::string& filename,
                                     const bool keep_diag)
    : keep_diag(keep_diag)
{
    std::cout << "~~~ processing transition matrix\n"
              << "    open from: " << filename << std::endl;
    std::vector<double> T(0);
    std::ifstream ifs(filename);
    if (ifs.fail()) {
            std::cerr << "    ERROR: cannot open file" << std::endl;
            exit(EXIT_FAILURE);
    } else {
        while (!ifs.eof() && !ifs.bad()) {
            double buf;
            ifs >> buf;
            if ( ! ifs.fail()) {
                T.push_back(buf);
            } else {  // if conversion error, skip (comment) line
                ifs.clear();
                ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
    }
    if (T.empty()) {
        std::cerr << "    ERROR: opened empty file" << std::endl;
        exit(EXIT_FAILURE);
    }
    ifs.close();
    this->transition_matrix_input = T;

    // check if quadratic
    this->no_of_states = (State) sqrt(T.size());
    if (this->no_of_states != sqrt(T.size())) {
        std::cerr << "    ERROR: Marix is not quadratic" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // alloc memory for other vectors aswell and generate them
    this->transition_matrix = T;
//    this->accumulated_T = T;

    // check if normalized
    std::cout << "    normalize matrix" << std::endl;
    double sum;
    for (State i=0; i<this->no_of_states; ++i) {
        sum = 0.;
        for (State j=0; j<this->no_of_states; ++j) {
            sum += this->get_Tij(i, j, true);
//            this->accumulated_T[this->no_of_states*i + j] = sum;
        }
        if (std::abs(sum-1.) > 0) {
//            std::cerr << "   row " << i << " is not normalized with "
//                      << sum << std::endl;
            // normalize
            for (State j=0; j<this->no_of_states; ++j) {
                this->set_Tij(i, j,  this->get_Tij(i, j, true)/ sum, true);
//                this->accumulated_T[this->no_of_states*i + j] /= sum;
            }
        }
    }


    //TODO: check if diagonal element is equal to 1
    // generate matrix with T_ii=0
    if (!keep_diag) {
        std::cout << "    set T_ii=0" << std::endl;
    }
    this->remove_self_transition_rate();


    this->accumulated_T = T;
    this->accumulated_T_input = T;
    this->timescale_matrix = T;
    double sum_input;
    for (State i=0; i<this->no_of_states; ++i) {
        sum = 0.;
        sum_input = 0.;
        for (State j=0; j<this->no_of_states; ++j) {
            // fill accumulated matrices
            sum += this->get_Tij(i, j, false);
            sum_input += this->get_Tij(i, j, true);
            this->accumulated_T[this->no_of_states*i + j] = sum;
            this->accumulated_T_input[this->no_of_states*i + j] = sum_input;

            // fill timescale matrix
            this->timescale_matrix[this->no_of_states*i + j] = 1./this->get_Tij(i, j, true);
        }
    }

    /* TODO: print input and used matrix
    std::cout << "    Transition matrix (input):" << std::endl;
    for (std::size_t i=0; i<this->no_of_states; ++i) {
        std::cout << "    ";
        for (std::size_t j=0; j<this->no_of_states; ++j) {
            std::cout << this->transition_matrix_input[this->no_of_states*i + j]
                      << "\t";
        }
        std::cout << std::endl;
    }
    */
    std::cout << std::endl;
}

// copy constructor
// TODO: std::vector is always deep copy, default copy constructor is sufficient
Transition_matrix::Transition_matrix(const Transition_matrix& T)
    : no_of_states(T.no_of_states), keep_diag(T.keep_diag)
{
    this->transition_matrix = T.transition_matrix;
    this->transition_matrix_input = T.transition_matrix_input;
    this->accumulated_T = T.accumulated_T;
    this->accumulated_T_input = T.accumulated_T_input;
}

// Destructor
Transition_matrix::~Transition_matrix() {}

float Transition_matrix::get_Tij(const State i,
                                 const State j)
{
    return this->get_Tij(i, j, this->keep_diag);
}

float Transition_matrix::get_Tij(const State i,
                                 const State j,
                                 const bool keep_diag)
{
    std::vector<double>* T;
    if (keep_diag) {
        T = &this->transition_matrix_input;
    } else {
        T = &this->transition_matrix;
    }

    if (i < this->no_of_states && j < this->no_of_states) {
        return (*T)[this->no_of_states*i + j];
    } else {
        std::cerr << "    ERROR: Matrix is smaller than requested value"
                  << "    Dimension=" << this->no_of_states
                  << " vs. (" << i << ", " << j << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}


//TODO: modify Tij diagonal and accumulated aswell
void Transition_matrix::set_Tij(const State i,
                                const State j,
                                const float val)
{
    return this->set_Tij(i, j, val, this->keep_diag);
}
void Transition_matrix::set_Tij(const State i,
                                const State j,
                                const float val,
                                const bool keep_diag)
{
    std::vector<double>* T;
    if (keep_diag) {
        T = &this->transition_matrix_input;
    } else {
        T = &this->transition_matrix;
    }

    if (i < this->no_of_states && j < this->no_of_states) {
        (*T)[this->no_of_states*i + j] = val;
    } else {
        std::cerr << "    ERROR: Matrix is smaller than requested value"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

State Transition_matrix::get_mcmc(const State i,
                                  const double rand)
{
    return this->get_mcmc(i, rand, this->keep_diag);
}
State Transition_matrix::get_mcmc(const State i,
                                  const double rand,
                                  const bool keep_diag)
{
    std::vector<double>* T_acc;
    if (keep_diag) {
        T_acc = &this->accumulated_T_input;
    } else {
        T_acc = &this->accumulated_T;
    }
    // catch errors
    if (rand >= 1 || i >= this->no_of_states) {
        std::cerr << "    ERROR: in get_mcmc, rand=" << rand
                  << " i=" << i << std::endl;
        std::exit(EXIT_FAILURE);
    }
    for (State j=0; j<this->no_of_states; ++j) {
        if ((*T_acc)[this->no_of_states*i + j] > rand) {
            return j;
        }
    }
    // catch error if rand=1, this is due to float precision, see
    // http://open-std.org/JTC1/SC22/WG21/docs/lwg-active.html#2524
    // find highest transition rate which
    for (State j=this->no_of_states-1; j>1; --j) {
        if ((*T_acc)[this->no_of_states*i + j] >
            (*T_acc)[this->no_of_states*i + j - 1]) {
            return j;
        }
    }
    return 0;
}

double Transition_matrix::get_tauij(const State i, const State j)
{
    if (i < this->no_of_states && j < this->no_of_states) {
        return this->timescale_matrix[this->no_of_states*i + j];
    } else {
        std::cerr << "    ERROR: Matrix is smaller than requested value"
                  << "    Dimension=" << this->no_of_states
                  << " vs. (" << i << ", " << j << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void Transition_matrix::remove_self_transition_rate()
{
    float sum = 0.;
    for (State i=0; i<this->no_of_states; ++i) {
        sum = 1.-this->get_Tij(i, i, true);

        for (State j=0; j<this->no_of_states; ++j) {
            if (i!=j) {
                this->set_Tij(i, j, this->get_Tij(i, j, true)/sum, false);
            } else {
                this->set_Tij(i, j, 0., false);
            }
        }
    }
    return;
}

void Transition_matrix::clear()
{
    this->transition_matrix.clear();
    this->transition_matrix_input.clear();
    this->accumulated_T.clear();
    this->accumulated_T_input.clear();
    this->timescale_matrix.clear();
    this->no_of_states = 0;
    this->keep_diag = true;
    return;
}

State Transition_matrix::get_no_of_states()
{
    return this->no_of_states;
}
