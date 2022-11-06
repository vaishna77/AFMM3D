#include "kernel.hpp"
#include "ACA.hpp"
#include "FMM3DTreeRAMeff2.hpp"
// #include <fstream>
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <filesystem>
// #include <cereal/archives/portable_binary.hpp>
// #include <cereal/types/vector.hpp>
//
// namespace storedata
// {
// Eigen::VectorXd inline load_vec(std::string fname)
// {
//         Eigen::VectorXd X;
//         std::vector<double> K;
//         std::ifstream infile (fname, std::ios::binary);
//         cereal::PortableBinaryInputArchive iarchive(infile);      // Create an input archive
//         iarchive(K);  // Read the data from the archive
//         X = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(K.data(), K.size());
//         return X;
// }
//
// void inline save_vec(std::string fname,Eigen::VectorXd X)
// {
//         std::vector<double> K(X.data(), X.data() + X.size());
//         std::ofstream outfile (fname,std::ios::binary);
//         cereal::PortableBinaryOutputArchive oarchive(outfile);    // Create an output archive
//         oarchive(CEREAL_NVP(K));                                  // Write the data to the archive
// }
// }

int main(int argc, char* argv[]) {
	int cubeRootN		=	atoi(argv[1]);
	int nParticlesInLeafAlong1D	=	atoi(argv[2]); // assuming the particles are located at tensor product chebyNodes
	double L			=	atof(argv[3]);
	int TOL_POW = atoi(argv[4]);
	int Qchoice = atoi(argv[5]);
	double start, end;
	int nLevels		=	ceil(3*log(double(cubeRootN)/nParticlesInLeafAlong1D)/log(8));
	// std::cout << "nLevels: " << nLevels << std::endl;
	// std::cout << "nLevels: " << nLevels << std::endl;
  // omp_set_dynamic(0);
  // omp_set_num_threads(20);
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	std::vector<pts3D> particles;
	userkernel* mykernel		=	new userkernel(particles, Qchoice);
	FMM3DTree<userkernel>* A	=	new FMM3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

	// A->set_Uniform_Nodes();
	A->set_Standard_Cheb_Nodes();
	A->createTree();
	A->assign_Tree_Interactions();
	A->assign_Center_Location();
	// A->assignLeafChargeLocations();
	A->assignChargeLocations();
	A->assignNonLeafChargeLocations();
	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	// std::cout << std::endl << "Number of particles is: " << A->N << std::endl;
	// std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	// A->getUVtTree();
	A->getNodes();
	A->assemble_M2L();
	// A->assemble_NearField();
	end		=	omp_get_wtime();

	double timeAssemble =	(end-start);
	// std::cout << std::endl << "Time taken to assemble is: " << timeAssemble << std::endl;
	/////////////////////////////////////////////////////////////////////////
	int N = A->N;
	// Eigen::VectorXd b=Eigen::VectorXd::Ones(N);

	Eigen::VectorXd b=Eigen::VectorXd::Zero(N);
  int n = N/500;
  srand(time(NULL));
  std::set<int> s;
  while(s.size() < n) {
    int index	=	rand()%N;
    s.insert(index);
  }
  std::set<int>::iterator it;
  for (it = s.begin(); it != s.end(); it++) {
    b(*it) = 1.0;
  }

	A->assignCharges(b);
	// A->check();
	// A->assignLeafCharges(b);
	// A->assignNonLeafCharges();

	start	=	omp_get_wtime();
	// A->evaluateFarField();
	// A->evaluate_NearField();

	A->evaluate_M2M();
	A->evaluate_M2L();
	A->evaluate_L2L();
	A->evaluate_NearField();
	Eigen::VectorXd AFMM_Ab;
	A->collectPotential(AFMM_Ab);
	A->reorder(AFMM_Ab);

	end		=	omp_get_wtime();
	// timeHODLR3DRepres += A->assTime;

	double timeMatVecProduct = end-start;
	// std::cout << std::endl << "Time taken to do Mat-Vec product is: " << timeMatVecProduct << std::endl;
	//////////////////
	double err;
	// Eigen::VectorXd true_Ab = Afull*b; //comment this
  std::string fname; //uncomment this
  Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
  for (it = s.begin(); it != s.end(); it++) {
    true_Ab = true_Ab + mykernel->getCol(*it);
  }

	// Eigen::VectorXd true_Ab = Eigen::VectorXd::Zero(N);
	// #pragma omp parallel for
	// for (size_t i = 0; i < N; i++) {
	// 	#pragma omp parallel for
	// 	for (size_t j = 0; j < N; j++) {
	// 		true_Ab(i) += A->K->getMatrixEntry(i,j)*b(j);
	// 	}
	// }


  // #pragma omp parallel for
  // for (size_t i = 0; i < N; i++) {
  //   true_Ab(i) = A->K->getMatrixEntry(i,0);
  // }
  // ------------------- write b to file--------------------------------
	// #pragma omp parallel for
	// for (size_t i = 0; i < N; i++) {
	// 	#pragma omp parallel for
	// 	for (size_t j = 0; j < N; j++) {
	// 		true_Ab(i) += A->K->getMatrixEntry(i,j)*b(j);
	// 	}
	// }
  // fname = "trueb/trueb_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)N) + ".bin";
  // storedata::save_vec(fname, true_Ab);
  // ---------------------------------------------------

  // ------------------- read b from file--------------------------------
  // fname = "../trueb/trueb_" + std::to_string((int)Qchoice) + "_" + std::to_string((int)N) + ".bin";
  // true_Ab = storedata::load_vec(fname);
  // ---------------------------------------------------

	err = (true_Ab - AFMM_Ab).norm()/true_Ab.norm();
	// std::cout << "err: " << err << std::endl;
	// std::cout << "true_Ab - AFMM_Ab: " << std::endl << true_Ab - AFMM_Ab << std::endl;
	// std::cout << "AFMM_Ab.norm(): " << AFMM_Ab.norm() << std::endl;
	// std::cout << "true_Ab.norm(): " << true_Ab.norm() << std::endl;
	// std::cout << "AFMM_Ab: " << AFMM_Ab << std::endl;
	// std::cout << "true_Ab: " << true_Ab << std::endl;

	double sum;
	A->findMemory2(sum);
	// std::cout << "Memory in GB: " << sum/8*pow(10,-9) << std::endl;
	// std::cout << "CR: " << double(sum)/N/N << std::endl;
	// std::cout << "max rank: " << A->getAvgRank() << std::endl;
  std::cout << N << " " << TOL_POW << " " << Qchoice << " " << sum/8*pow(10,-9) << " " << timeAssemble << " " << timeMatVecProduct << " " << A->getAvgRank() << " " << err << std::endl;
	delete A;
	delete mykernel;
}
