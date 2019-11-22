#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

#include <string>
#include <vector>
#include <map>

/*
When declaring a class use the following order

    types: classes, enums, and aliases (using)
    constructors, assignments, destructor
    functions
    data

Use the public before protected before private order.

This is a recommendation for when you have no constraints or better ideas. This rule was added after many requests for guidance.
*/

class Transition_matrix {
    public:
        //! default constructor
        Transition_matrix();
        //! copy constructor
        Transition_matrix(const Transition_matrix& T);
        //! constructor from file,
        Transition_matrix(const std::string& filename,
                          const bool keep_diag);
        //! default destructor
        ~Transition_matrix();

        //! get entry (i,j) of transition matrix
        float get_Tij(const std::size_t i,
                      const std::size_t j);
        //! get entry (i,j) of diagonal or none-diagonal transition matrix
        float get_Tij(const std::size_t i,
                      const std::size_t j,
                      const bool keep_diag);
        //! get entry (i,j) of transition matrix
        void set_Tij(const std::size_t i,
                     const std::size_t j,
                     const float val);
        //! set entry (i,j) of diagonal or none-diagonal transition matrix
        void set_Tij(const std::size_t i,
                     const std::size_t j,
                     const float val,
                     const bool keep_diag);

        // class specific methods
        std::size_t get_mcmc(const std::size_t i,
                             const double rand);
        std::size_t get_mcmc(const std::size_t i,
                             const double rand,
                             const bool keep_diag);
        double get_tauij(const std::size_t i,
                         const std::size_t j);

        void set_keep_diag(const bool keep_diag) {this->keep_diag = keep_diag;}

        std::size_t get_no_of_states();

        //! setting the diagonal T_ii->0 and renormalize it
        void remove_self_transition_rate();
        //! reset all values to default
        void clear();
    private:
        // variables
        std::vector<double> transition_matrix;
        std::vector<double> transition_matrix_input;
        std::vector<double> accumulated_T;            // needed for MCMC
        std::vector<double> accumulated_T_input;
        std::vector<double> timescale_matrix;         // needed for paths method
        std::size_t no_of_states = 0;
        bool keep_diag = true;                        // if true, T_ii->0
};
#endif
