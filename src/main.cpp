#include <iostream>
#include <fstream>
#include "main.hpp"
#include "paths.hpp"
#include "pathfinder.hpp"
#include "transition_matrix.hpp"

#include <string>
#include <vector>
#include <cmath>

#include <chrono>

//#include <exception>

#include <boost/program_options.hpp>
//#include <omp.h>

//! check if states in A are legit
bool check_states(const std::vector<State> A, const std::string str_A,
                  const std::vector<State> B, const std::string str_B,
                  const std::vector<State> C, const std::string str_C,
                  const State no_of_states)
{
    return (check_states_A(A, str_A, B, str_B, C, str_C, no_of_states) &
            check_states_A(B, str_B, A, str_A, C, str_C, no_of_states) &
            check_states_A(C, str_C, A, str_A, B, str_B, no_of_states));
}

bool check_states_A(const std::vector<State> A, const std::string str_A,
                    const std::vector<State> B, const std::string str_B,
                    const std::vector<State> C, const std::string str_C,
                    const State no_of_states)
{
    // check if all values are legal and unique
    for (auto &state : A) {
        // because of size_t < 0 is not possible
        if (state >= no_of_states) {
            std::cout << "    ERROR: state " << int(state)+1
                      << " in " << str_A << " is larger/smaller than matrix"
                      << std::endl;
            return false;
        }
        // check if A has overlap with B or C
        if (std::count(B.begin(), B.end(), state)) {
            std::cout << "    ERROR: state " << state+1
                      << " is in " << str_A << " and " << str_B << std::endl;
            return false;
        }
        if (std::count(C.begin(), C.end(), state)) {
            std::cout << "    ERROR: state " << state+1
                      << " is in " << str_A << " and " << str_C << std::endl;
            return false;
        }
    }

    auto it_A = std::adjacent_find(A.begin(), A.end());
    if (!(it_A == A.end())) {
        std::cout << "    ERROR: state " << int(*it_A)+1
                  << " in " << str_A << " is not unique" << std::endl;
        return false;
    }

    return true;
}

//! MAIN ROUTINE
int main(int argc, char* argv[]){
    namespace b_po = boost::program_options;
    std::string version_number = "v0.8";
    // generate header string
    std::string leading_whitespace(25 - (22 + version_number.size())/2, ' ');
    std::ostringstream header_ostring;
    header_ostring << "\n" << leading_whitespace
                 << "~~~ MSMPathfinder " + version_number + " ~~~\n";
    if (argc > 2) {
        std::string leading_whitespace2nd(25 - (4 + strlen(argv[1]))/2, ' ');
        header_ostring << leading_whitespace2nd << "~ " << argv[1] << " ~\n";
    }
    std::string clustering_copyright =
        header_ostring.str() + ""
        "MSMPathfinder " + version_number +
        ": a framework for calculating pathways from states \n"
        "A={a_1, a2_,...} to states B={b_1, b2_,...} based on the provided Markov\n"
        "state model.\n"
        "Copyright (c) 2019, Daniel Nagel\n\n"
    ;
    std::string general_help =
        clustering_copyright +
        "modes:\n"
        "   paths: propagate transition matrix and weight results\n"
        "   mcmc: run many short markov chain monte carlo (mcmc) and weight them\n"
        "   mcmcsingle: run a single very long markov chain monte carlo (mcmc)\n"
        "\n"
        "usage:\n"
        "  msmpathfinder MODE --option1 --option2 ...\n"
        "\n"
        "for a list of available options per mode, run with '-h' option, e.g.\n"
        "  msmpathfinder paths -h\n"
    ;

    Mode mode;
    if (argc <= 2) {
        std::cerr << general_help;
        return EXIT_FAILURE;
    } else {
        std::string str_mode(argv[1]);
        if (str_mode.compare("paths") == 0) {
            mode = PATHS;
        } else if (str_mode.compare("mcmc") == 0) {
            mode = MCMC;
        } else if (str_mode.compare("mcmcsingle") == 0) {
            mode = MCMCSINGLE;
        }else {
            std::cerr << "\nerror: unrecognized mode '" << str_mode << "'\n\n";
            std::cerr << general_help;
            return EXIT_FAILURE;
        }
    }

    b_po::variables_map args;
    // header string for mode specific help
    std::string help_header = clustering_copyright + std::string(argv[1]) + ": \n";
    switch(mode){
        case PATHS:
            help_header +=
                "This program calculates the pathway distribution directly\n"
                "from the transition matrix."
                "\nThis method is suitable for up to 100 states.\n";
            break;
        case MCMC:
            help_header +=
                "This program estimates the pathway distribution by running\n"
                "multiple short markov chain monte carlo method runs.";
            break;
        case MCMCSINGLE:
            help_header +=
                "This program estimates the pathway distribution by running\n"
                "a single long markov chain monte carlo method run.";
            break;
        default:
            std::cerr << "    ERROR: unknown mode. this should never happen."
                      << std::endl;
            return EXIT_FAILURE;
    }
    help_header +=
        "\nTo keep the number of pathways feasible, loops are cut, so\n"
        "    1->2->5->3->2->7  => 1->2->7\n"
        "If multiple (nested) loops are found always the longest is\n"
        "cut first.\n\noptions";

    // parsing options
    b_po::options_description desc (help_header);
    desc.add_options()
        ("help,h", b_po::bool_switch()->default_value(false), "show this help.")
        ("input-file,i", b_po::value<std::string>()->required(),
         "input (required): transition matrix (multi-column ASCII, row normalized).")
        ("states-from,f", b_po::value<std::vector<State>>()->multitoken()->required(),
         "parameter: list of states (one-indexed) where pathways starts from.")
        ("states-to,t", b_po::value<std::vector<State>>()->multitoken()->required(),
         "parameter: list of states (one-indexed) where pathways ends.")
        ("states-forbidden", b_po::value<std::vector<State>>()->multitoken(),
         "parameter: list of states (one-indexed) which pathways should not enter.")
        ("steps", b_po::value<std::size_t>(),
         "parameter: number of steps which a transition can be maximal long.")
        ("threshold", b_po::value<std::size_t>()->default_value(100000),
         "parameter: maximal number of different pathways which are stored, -1 for all (default: 100000).")
    ;
    // add mode specific options
    switch(mode){
        case PATHS:
            desc.add_options()
                ("precision,p", b_po::value<double>()->required(),
                 "parameter: threshold up to which all paths are taken into account for the path search.")
            ;
            break;
        case MCMC:
            desc.add_options()
                ("iterations", b_po::value<std::size_t>()->required(),
                 "parameter: number of iterations per state")
            ;
            break;
        case MCMCSINGLE:
            desc.add_options()
                ("total-steps", b_po::value<std::size_t>()->required(),
                 "parameter: number of total steps the mcmc is generated")
            ;
            break;
        default:
            std::cerr << "    ERROR: unknown mode. this should never happen."
                      << std::endl;
            return EXIT_FAILURE;
    }
    // add optional and default options
    desc.add_options()
        // optional
        ("output,o", b_po::value<std::string>(), "output (optional): save pathways to file.")
        // defaults
        ("keep-diag", b_po::bool_switch()->default_value(false),
         "flag (optional): omits setting T_ii->0.")
        ("nthreads,n", b_po::value<int>()->default_value(0),
         "number of OpenMP threads. default: 0; i.e. use OMP_NUM_THREADS env-variable.")
        ("verbose,v", b_po::bool_switch()->default_value(false),
         "verbose mode: print runtime information to STDOUT.")
    ;

    // try to parse arguments
    try {
        b_po::store(b_po::parse_command_line(argc, argv, desc), args);
        if (args["help"].as<bool>()) {
            std::cout << desc << std::endl;
            std::exit(EXIT_SUCCESS);
        }
        b_po::notify(args);
    } catch (std::exception& e) {
        std::cerr << "\nerror parsing arguments:\n    " << e.what() << std::endl;
        std::cerr << "\nuse --help/-h for further usage information.\n" << std::endl;
//        std::cerr << desc << std::endl;
        std::exit(EXIT_FAILURE);
    } catch (...) { // This should never be reached!
        std::cerr << "\nunhandled error occured.\n" << std::endl;
    }
    // SET NUMBER OF THREADS
//    omp_set_num_threads(6);

    // get parsed input
    std::string output;
    if (args.count("output")) {
        output = args["output"].as<std::string>();
    } else {
        output = "not specified! only first 50 pathways are printed on console";
    }
    std::string file_name = args["input-file"].as<std::string>();
    std::vector<State> A = args["states-from"].as<std::vector<State>>();
    std::vector<State> B = args["states-to"].as<std::vector<State>>();
    std::vector<State> C;
    if (args.count("states-forbidden")) {
        C = args["states-forbidden"].as<std::vector<State>>();
    }

    // shift one-based states to zero based
    for (auto &val : A) {
        val -= 1;
    }
    for (auto &val : B) {
        val -= 1;
    }
    if (args.count("states-forbidden")) {
        for(auto &val : C) {
            val -= 1;
        }
    }
    // sort states
    std::sort(A.begin(), A.end());
    std::sort(B.begin(), B.end());
    std::sort(C.begin(), C.end());

    bool keep_diag = args["keep-diag"].as<bool>();

    // print header
    std::cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "          ~~~ Markov State Model PATHFINDER ~~~\n"
              << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "\n~~~ Input parameters:\n"
              << "    output: " << output << "\n"
              << std::endl;

    // read transition matrix from file
    // TODO: remove steps and iterations values to input default values
    std::size_t steps = 0;
    if (args.count("steps")) {
        steps = args["steps"].as<std::size_t>();
    }
    double cut_off = 0.;
    if (args.count("precision")) {
        cut_off = args["precision"].as<double>();
    }
    std::size_t iterations = 0;
    if (args.count("iterations")) {
        iterations = args["iterations"].as<std::size_t>();
    }
    std::size_t total_steps = 0;
    if (args.count("total-steps")) {
        total_steps = args["total-steps"].as<std::size_t>();
    }

    if (!(args.count("output"))) {
        output = "";
    }

    std::size_t threshold = args["threshold"].as<std::size_t>();

    // generate transition matrix
    Transition_matrix T(file_name, keep_diag);

    // check states
    bool legit_states = check_states(A, "states-from",
                                     B, "states-to",
                                     C, "states-forbidden",
                                     T.get_no_of_states());
    if (legit_states) {
        std::cout << "    from: ";
        for (auto &val : A) {
            std::cout << val+1 << " ";
        }
        std::cout << "\n    to: ";
        for (auto &val : B) {
            std::cout << val+1 << " ";
        }
        if (args.count("states-forbidden")) {
            std::cout << "\n    forbidden: ";
            for (auto &val : C) {
                std::cout << val+1 << " ";
            }
        }
        std::cout << std::endl;
    } else {
        std::cout << "    ERROR: Bad selection of states." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // call main routine
    run_paths(mode, A, B, C, cut_off, steps, iterations, total_steps, T, output, threshold, argc, argv);
}
