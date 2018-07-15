
#include <Eigen/Sparse>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>

int main(int argc, char const * const * argv) {
	// Get params
	srand(time(NULL));
	int size = 0, count = 0;
	if (argc > 1)
		size = atoi(argv[1]);
	if (size < 3)
		size = 3;
	if (argc > 2)
		count = atoi(argv[2]);
	if (count < 0)
		count = 0;
	// Generate L matrix
	std::vector<Eigen::Triplet<float> > triplets;
	for (int k = 0; k < count; ++k) {
		int i = rand() % size;
		int j = rand() % size;
		if (j > i) {
			int t = i;
			i = j;
			j = t;
		}
		float v = (rand() % 100) * 0.1f;
		triplets.push_back(Eigen::Triplet<float>(i, j, v));
	}
	Eigen::SparseMatrix<float> L(size, size);
	L.setFromTriplets(triplets.begin(), triplets.end());
	// Compute final LL* matrix and add I to garantee positive eigen values
	Eigen::SparseMatrix<float> out(size, size);
	out.setIdentity();
	out += L * L.transpose();
	// Check if Cholesky yield something
	Eigen::SimplicialLLT<Eigen::SparseMatrix<float> > cholesky(out);
	Eigen::SparseMatrix<float> LX = cholesky.matrixL();
	// Output triplets
	std::cout << size << std::endl;
	for (int k = 0; k < out.outerSize(); ++k)
		for (Eigen::SparseMatrix<float>::InnerIterator it(out, k); it; ++it)
			std::cout << it.row() << ' ' << it.col() << ' ' << it.value() << std::endl;
	return 0;
}
