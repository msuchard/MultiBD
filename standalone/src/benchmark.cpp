
#include <chrono>
#include <random>
#include <iostream>

#include <boost/program_options.hpp>

//#include "Rcpp.h"
#include <RInside.h>
#include "AbstractEngine.h"

// template <typename T, typename PRNG, typename D>
// void generateLocation(T& locations, D& d, PRNG& prng) {
// 	for (auto& location : locations) {
// 		location = d(prng);
// 	}
// }

typedef std::shared_ptr<np_cluster::AbstractEngine> EnginePtr;

__attribute__((constructor))
static void initialize(void) {
	RInside R(0, nullptr); // This is a horrible hack to allow Rcpp memory objects as global variables
	// TODO Real solution is to move RInside() into main() and *not* use global variables
}

int main(int argc, char* argv[]) {

// Set-up CLI
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("sse", "use SSE vectorization")
		("unroll", "unroll loop")				
		("tbb", po::value<int>()->default_value(0), "use TBB with specified number of threads")		
		("iterations", po::value<int>()->default_value(1), "number of iterations")
	;
	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	} catch (std::exception& e) {
		std::cout << desc << std::endl;
		return 1;
	}

	if (vm.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}



// 	auto useSSE = false;
// 	auto useTBB = false;
	
// 	RInside R(0, nullptr); // This is a horrible hack to allow Rcpp memory objects as global variables
	
    using namespace Rcpp;
    
	unsigned int seed = 666;

	// Set the R RNG seed
	Rcpp::Environment baseEnv("package:base");
	Rcpp::Function setSeed = baseEnv["set.seed"];
	setSeed(seed);    
		
	auto N = 200000;
// 	auto P = 2000;
	auto K = 100;
	
	auto tau = 0.5;
	auto tau0 = 2.0;
// 	auto tauInt = 0.5;
	
	auto extraSort = true;
// 	auto colSum = false;
	
	auto specialMode = 0;
	int threads = vm["tbb"].as<int>();
	if (threads > 0) {
		specialMode |= static_cast<long>(np_cluster::SpecialComputeMode::TBB);
	}
	if (vm.count("sse")) {
		specialMode |= static_cast<long>(np_cluster::SpecialComputeMode::SSE);
	}
	if (vm.count("unroll")) {
		specialMode |= static_cast<long>(np_cluster::SpecialComputeMode::UNROLL);
	}
		
	auto prng = std::mt19937(seed);
	auto normal = std::normal_distribution<double>(0.0, 1.0);
	auto uniform = std::uniform_int_distribution<int>(0, 5);
// 	auto binomial = std::bernoulli_distribution(0.75);
// 	auto normalData = std::normal_distribution<double>(0.0, 1.0);
	
	EnginePtr engine = EnginePtr(
	    new np_cluster::AbstractEngine(extraSort, specialMode));
	
	std::cout << "Loading data" << std::endl;
	
// 	NumericMatrix X(N, P);
// 	NumericMatrix A(N, K);
// 	IntegerMatrix S(N, K);
	
	IntegerVector nVec(K);
	NumericVector Y(N);
	NumericVector Xsd(N);
	IntegerVector rowSubsetI(N);
	IntegerVector cmVec(1000);
	NumericVector phiV(K);
	
// 	for (auto& x : X) {
// 	    x = normal(prng);
// 	}
// 	
// 	for (auto it = std::begin(X); it != std::end(X); ++it) {
// 	    *it = normal(prng);
// 	}
// 	
// 	for (auto it = std::begin(A); it != std::end(A); ++it) {
// 	    *it = normal(prng);
// 	}
// 	
// 	for (auto it = std::begin(S); it != std::end(S); ++it) {
// 	    *it = uniform(prng) + 1;
// 	}
	
	for (auto& n : nVec) {
		n = uniform(prng) + 1;
	}
	
	for (auto& y : Y) {
		y = normal(prng);
	}
	
	for (auto& x : Xsd) {
		x = normal(prng);
	}
	
	for (auto& c : cmVec) {
		c = uniform(prng) + 1;
	}
	
	for (auto& phi : phiV) {
		phi = std::exp(normal(prng));
	}
	
	for (int i = 0; i < N; ++i) {
		rowSubsetI[i] = i + 1;
	}

	auto n0 = 0;
	auto n2 = 200;
	auto epsilon = 0.001;
	auto epsilon2 = 0.00001;
	auto maxNeighborhoodSize = 28;
	auto cutOff = 0.1;
	
	int iterations = vm["iterations"].as<int>();
	
    std::cout << "Running benchmark ..." << std::endl;
    auto startTime = std::chrono::steady_clock::now();
    
    // Run
    for (int i = 0; i < iterations; ++i) {
    
	    auto list = engine->computePmfAndNeighborhoods(
	    	n0, nVec, epsilon, epsilon2, K, N, 
    		Y, Xsd, rowSubsetI, cmVec, n2, phiV, 
    		tau, tau0, maxNeighborhoodSize, cutOff);
    
    }
    
//   Rcpp::List computePmfAndNeighborhoods(int n0, const IntegerVector& nVec,
//                                               const double epsilon, const double epsilon2,
//                                               int K, int N,  // NB: K -> G for column-space
//                                               const NumericVector& Y, const NumericMatrix& X,
//                                               const NumericMatrix& A, const IntegerMatrix& S,
//                                               const IntegerVector& rowSubsetI,
//                                               const IntegerVector& CmVec, const int P,
//                                               const NumericVector& phiV, const double tau, const double tau0, const double tauInt,
//                                               const int maxNeighborhoodSize,
//                                               const double cutOff,
//                                               const bool collectMax)     

	auto endTime = std::chrono::steady_clock::now();
	auto duration = endTime - startTime;
	std::cout << "End benchmark" << std::endl;
	std::cout << std::chrono::duration<double, std::milli> (duration).count() << " ms "
			  << std::endl;
			  
	engine->printTiming();

}
