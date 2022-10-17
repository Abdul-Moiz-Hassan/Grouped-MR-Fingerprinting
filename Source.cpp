// Grouped MRF Pattern Matching C++ Implementation
// @author: Abdul Moiz Hassan - MIPRG, COMSATS University Islamabad

#include <fstream>
#include <iostream> 
#include <complex>  
#include <cstdlib>
#include <chrono>

using namespace std;

/**
 * Compute inner product of complex matrix A with complex matrix B
 * (where matrices are in row major format)
 *
 * Inner products are used to compute correlations of signals
 * arranged in matrices.
 */
complex<float>* innerProduct(std::complex<float>* A, std::complex<float>* B, int ARows, int ACols, int BRows, int BCols)
{
	int CRows = ARows; int CCols = BCols;
	std::complex<float>* C = (std::complex<float> *)malloc(CRows * CCols * sizeof(std::complex<float>));

	// Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
	for (int i = 0; i < ARows; i++)
	{
		for (int j = 0; j < BCols; j++)
		{
			for (int k = 0; k < ACols; k++)
			{
				*(C + i * CCols + j) += (*(A + i * ACols + k)) * (*(B + k * BCols + j));
			}
		}
	}
	return C;
}

/**
 * Find the maximum number in an array and return it's index/position
 *
 */
int Max(float* arr, int size)
{
	int index;
	float max_n;
	max_n = 0;

	for (int i = 0; i < size; i++)
	{
		if (arr[i] > max_n)
		{
			max_n = arr[i];
			index = i;
		}
	}
	return index;
}

int main()
{
	//Grouped dictionary parameters are defined when grouping data in MATLAB
	const int t = 1000;
	const int n_groups = 23;
	const int dict_sz_3D = 5865;
	const int imadapt_sz = 57600;

	// Input File Streams to read data (all data files are in row major format)
	// Acquired Signals to be pattern matched with grouped dictionary
	std::ifstream acuired_real("Acquired Signals/imadapt_r_VUMC.txt");
	std::ifstream acuired_imag("Acquired Signals/imadapt_i_VUMC.txt");

	// Mean Representatives of all Signals
	std::ifstream Means_rr("Grouped Dictionary/Mean_Shift_Clustering/G_Means_r_22_mnsh.txt");
	std::ifstream Means_ii("Grouped Dictionary/Mean_Shift_Clustering/G_Means_i_22_mnsh.txt");

	// Grouped Signals
	std::ifstream Grouped_dict_i("Grouped Dictionary/Mean_Shift_Clustering/Grouped_i_22_mnsh.txt");
	std::ifstream Grouped_dict_r("Grouped Dictionary/Mean_Shift_Clustering/Grouped_r_22_mnsh.txt");

	// Grouped Signals Sizes (individual size of each group)	
	std::ifstream Grouped_labels("Grouped Dictionary/Mean_Shift_Clustering/G_lables_22_mnsh.txt");
	std::ifstream Grouped_Sizes("Grouped Dictionary/Mean_Shift_Clustering/G_Sizes_22_mnsh.txt");

	// Output file stream to store resultant/reconstructed images
	std::ofstream T1_cpu("T1_CPU.txt");
	std::ofstream T2_cpu("T2_CPU.txt");

	// Allocate memory to store/load data
	float* G_lables = (float*)malloc(2 * dict_sz_3D * n_groups * sizeof(float));

	std::complex<float>* Grouped_dict = (std::complex<float> *)malloc(t * dict_sz_3D * n_groups * sizeof(std::complex<float>));
	std::complex<float>* Grouped_Means = (std::complex<float> *)malloc(t * n_groups * sizeof(std::complex<float>));

	std::complex<float>* acuired_signals = (std::complex<float> *)malloc(imadapt_sz * t * sizeof(std::complex<float>));

	float G_sizes[n_groups];
	float Reald = 0;
	float Imagd = 0;

	// Data loading can be imrpoved signnificantly but is not the scope if this project
	// this program is used to accelerate patten matching time (not data loading)
	cout << "Loading Grouped Dictionary : .....  ";
	for (int k = 0; k < n_groups; k++)
	{
		for (int j = 0; j < dict_sz_3D; j++)    // Reading real part of dict
		{
			for (int i = 0; i < t; i++)
			{
				Grouped_dict_r >> Reald;
				Grouped_dict_i >> Imagd;

				*(Grouped_dict + ((i * dict_sz_3D + j) * n_groups) + k) = { Reald ,Imagd };
			}
		}
	}
	cout << " Completed !! " << endl;

	cout << endl << "Loading Grouped Dictionary Labels : .....  ";
	for (int k = 0; k < n_groups; k++)
	{
		for (int j = 0; j < dict_sz_3D; j++)    // Reading real part of dict
		{
			for (int i = 0; i < 2; i++)
			{
				Grouped_labels >> Reald;

				*(G_lables + ((i * dict_sz_3D + j) * n_groups) + k) = Reald;
			}
		}
	}
	cout << " Completed !! " << endl;

	cout << endl << "Loading Grouped Dictionary Sizes : .....  ";
	for (int k = 0; k < n_groups; k++)
	{
		Grouped_Sizes >> Reald;

		G_sizes[k] = Reald;

	}
	cout << " Completed !! " << endl;

	cout << endl << "Loading Grouped Means : ..... ";
	for (int i = 0; i < t; i++)
	{
		for (int j = 0; j < n_groups; j++)
		{
			Means_rr >> Reald;
			Means_ii >> Imagd;

			*(Grouped_Means + i * n_groups + j) = { Reald ,Imagd };
		}
	}
	cout << " Completed !! " << endl;

	cout << endl << "Loading Acquired Scan/Signals : ..... ";
	for (int i = 0; i < imadapt_sz; i++)    // Reading real part of dict
	{

		for (int j = 0; j < t; j++)
		{
			acuired_real >> Reald;
			acuired_imag >> Imagd;

			*(acuired_signals + i * t + j) = { Reald,Imagd };
		}
	}
	cout << " Completed !! " << endl;

	//////////////////////////////////////////////////////////////
	//              Pattern Matching Starts            //
	//////////////////////////////////////////////////////////////

	// Allocate memory for inner products of mean signal
	std::complex<float>* IPs_m = (std::complex<float> *)malloc(imadapt_sz * n_groups * sizeof(std::complex<float>));

	float RIPs_m[n_groups];

	float T1[imadapt_sz] = {};
	float T2[imadapt_sz] = {};
	int index_m[imadapt_sz] = {};
	int groups_IP_sz[n_groups] = {};

	auto t1 = std::chrono::high_resolution_clock::now();

	// Correlate Acquired Signals with Mean representative signals
	IPs_m = innerProduct(acuired_signals, Grouped_Means, imadapt_sz, t, t, n_groups);

	// Find maximum correlaction index against each acquired signal
	for (int ii = 0; ii < imadapt_sz; ii++)
	{
		for (int j = 0; j < n_groups; j++)
		{
			RIPs_m[j] = abs(*(IPs_m + ii * n_groups + j));
		}
		index_m[ii] = Max(RIPs_m, n_groups);
		*RIPs_m = {};
		groups_IP_sz[(index_m[ii])]++;

	}

	int count = 0;
	int idx = 0;

	// Correlating acquired signal with its respective group 
	// (selected using max correlated mean representative)
	for (int g = 0; g < n_groups; g++)
	{
		int Gsz = G_sizes[g];
		float* RIPs = (float*)malloc(Gsz * sizeof(float));
		std::complex<float>* Signals = (std::complex<float> *)malloc(groups_IP_sz[g] * t * sizeof(std::complex<float>));
		std::complex<float>* tmp_G = (std::complex<float> *)malloc(t * Gsz * sizeof(std::complex<float>));
		int* Pixels = (int*)malloc(groups_IP_sz[g] * sizeof(int));
		*Pixels = 0;
		std::complex<float>* IPs = (std::complex<float> *)malloc(groups_IP_sz[g] * dict_sz_3D * sizeof(std::complex<float>));
		count = 0;

		// Group acuired signals based on coorelation index of mean  representatives
		for (int m = 0; m < imadapt_sz; m++)   //iterate over whole acquired Scan
		{
			if (g == index_m[m])
			{
				for (int k = 0; k < t; k++)
				{
					*(Signals + count * t + k) = *(acuired_signals + m * t + k);
				}

				*(Pixels + count) = m;
				count++;
			}
		}

		// Extract individual group from gouped dictionary
		for (int i = 0; i < t; i++)   //iterate over whole acquired Scan
		{
			for (int j = 0; j < Gsz; j++)
			{
				*(tmp_G + i * Gsz + j) = *(Grouped_dict + ((i * dict_sz_3D + j) * n_groups) + g);
			}
		}

		// Correlate aquired_signals_group with respective selected dictionary_group
		IPs = innerProduct(Signals, tmp_G, groups_IP_sz[g], t, t, Gsz);

		// Find maximum correlation index in the dictionary_group
		// and assign to respective pixel location
		for (int k = 0; k < groups_IP_sz[g]; k++)
		{
			for (int l = 0; l < Gsz; l++)
			{
				*(RIPs + l) = abs(*(IPs + k * Gsz + l));
			}
			idx = Max(RIPs, Gsz);

			T1[*(Pixels + k)] = *(G_lables + ((0 * dict_sz_3D + idx) * n_groups) + g);
			T2[*(Pixels + k)] = *(G_lables + ((1 * dict_sz_3D + idx) * n_groups) + g);

		}
	}

	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	std::cout << endl << "Time Elapsed :: " << duration << " micro seconds" << endl;

	// Store data to the output stream file for T1 and T2 maps
	for (int i = 0; i < imadapt_sz; i++)
	{

		if (T1_cpu.is_open() && T2_cpu.is_open())
		{
			T1_cpu << T1[i] << endl;
			T2_cpu << T2[i] << endl;
		}
		else cout << "Unable to open file";
	}

	T1_cpu.close();
	T2_cpu.close();

	return 0;
}
