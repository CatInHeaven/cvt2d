
// http://stackoverflow.com/questions/21012580/is-it-possible-to-write-data-to-file-using-only-javascript
// http://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html#details
// http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html

#include <Eigen/Sparse>

#ifdef EMSCRIPTEN

#include "emscripten.h"

extern "C" int process() {
	// Load matrix using JS global variables (easier than allocating data in Emscripten heap)
	std::vector<Eigen::Triplet<float> > triplets;
	int size = EM_ASM_INT({return inputSize;}, 0);
	int l = EM_ASM_INT({return inputI.length;}, 0);
	for (int k = 0; k < l; ++k) {
		int i = EM_ASM_INT({return inputI[$0];}, k);
		int j = EM_ASM_INT({return inputJ[$0];}, k);
		float v = EM_ASM_DOUBLE({return inputV[$0];}, k);
		triplets.push_back(Eigen::Triplet<float>(i, j, v));
	}
	Eigen::SparseMatrix<float> in(size, size);
	in.setFromTriplets(triplets.begin(), triplets.end());
	// Compute
	Eigen::SimplicialLLT<Eigen::SparseMatrix<float> > cholesky(in);
	Eigen::SparseMatrix<float> out = cholesky.matrixL();
	Eigen::SparseMatrix<float> outSqr = out * out.transpose();
	float norm = (outSqr - in).norm();
	// Write results in JS global variables
	int m = 0;
	for (int k = 0; k < out.outerSize(); ++k)
		for (Eigen::SparseMatrix<float>::InnerIterator it(out, k); it; ++it) {
			EM_ASM_ARGS({
				outputI[$0] = $1;
				outputJ[$0] = $2;
				outputV[$0] = $3;
			}, m, it.row(), it.col(), it.value());
			++m;
		}
	EM_ASM_ARGS({norm = $0;}, norm);
	return m;
}

#else

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <ratio>
typedef std::chrono::duration<double, std::ratio<1, 1000> > millis;

int main(int argc, char const * const * argv) {
	if (argc < 3) {
		std::cerr << "2 arguments required" << std::endl;
		return -1;
	}
	std::chrono::high_resolution_clock::time_point t1, t2, t3, t4;
	t1 = std::chrono::high_resolution_clock::now();
	// Load matrix
	std::ifstream inFile(argv[1]);
	std::vector<Eigen::Triplet<float> > triplets;
	int size;
	inFile >> size;
	if (!inFile) {
		std::cerr << "Failed to read size" << std::endl;
		return -1;
	}
	for (;;) {
		int i, j;
		float v;
		inFile >> i >> j >> v;
		if (!inFile)
			break;
		triplets.push_back(Eigen::Triplet<float>(i, j, v));
		if (i >= size || j >= size) {
			std::cerr << "File contains invalid triplet" << std::endl;
			return -1;
		}
	}
	Eigen::SparseMatrix<float> in(size, size);
	in.setFromTriplets(triplets.begin(), triplets.end());
	//std::cout << "Input " << size << "x" << size<< std::endl << in << std::endl;
	// Compute Cholesky factorization
	t2 = std::chrono::high_resolution_clock::now();
	Eigen::SimplicialLLT<Eigen::SparseMatrix<float> > cholesky(in);
	Eigen::SparseMatrix<float> out = cholesky.matrixL();
	//std::cout << "Output" << std::endl << out << std::endl;
	Eigen::SparseMatrix<float> outSqr = out * out.transpose();
	float norm = (outSqr - in).norm();
	t3 = std::chrono::high_resolution_clock::now();
	//std::cout << "Input from output" << std::endl << outSqr << std::endl;
	std::cout << "Frobenius norm: " << norm << std::endl;
	// Store result
	std::ofstream outFile(argv[2]);
	outFile << size << std::endl;
	for (int k = 0; k < out.outerSize(); ++k)
		for (Eigen::SparseMatrix<float>::InnerIterator it(out, k); it; ++it)
			outFile << it.row() << ' ' << it.col() << ' ' << it.value() << std::endl;
	t4 = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<millis>(t4 - t1).count() << "ms (" << std::chrono::duration_cast<millis>(t3 - t2).count() << "ms in C++ code)" << std::endl;
	return 0;
}

#endif
