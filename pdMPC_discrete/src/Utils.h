#pragma once
#include <string>
#include "DefineSettings.h"


class Utils{
	public:

	static void LoadSparseMatrixRowPlain(const char* str,int_t*& row_ptr, int_t*& col_idx, real_t*& val,
								int_t& nr, int_t& nc, int_t& nnz);
	
	static void LoadVec(const char* str, real_t** vec, int_t& nv);
	
	// for loading integer vectors
	static void LoadVec(const char* str, int_t** vec, int_t& nv);

	static void PrintVec(const real_t* v, int_t nv);

	static void PrintVec(const int_t* v, int_t nv);

	static void PrintMat(const real_t* v, int_t nc, int_t nr);
	
	// Important; The function assumes that col_idx is sorted (ascending)!
	// Moreover, col_ptr must be an array of lenght nCols+1 (already initialized).
	static void ConvertCOSToCCS(const int_t* col_idx, const int_t& nnZ, int_t* col_ptr);

	static void DotProduct(const real_t* a, const real_t* b, const int& n, real_t& res);

	// implements a+b=c, where a,b,c are n-dimensional
	static void VectorAdd(const real_t* a, const real_t* b, real_t* c, int_t n);

	// implements a-b=c, where a,b,c, are n-dimensional
	static void VectorSubstract(const real_t* a, const real_t* b, real_t* c, int_t n);

	// implements a*alpha+b*beta=c, where a,b,c are n-dimensional, and alpha and beta are scalars.
	static void VectorAddMultiply(const real_t* a, real_t alpha, const real_t* b, real_t beta, real_t* c, int_t n);

	// implements b=a
	static void VectorCopy(const real_t* a, real_t* b, int_t n);

	// returns ||a||_2
	static real_t VectorNorm(const real_t* a, int_t n);

	// returns ||a-b||_2
	static real_t VectorNormDiff(const real_t* a, const real_t* b, int_t n);
	
	// component wise multiplication of the vectors a and b
	static void VectorMult(const real_t* a, const real_t* b, real_t* c, int_t n);
	
	static int_t factorial(int_t n);

	static int_t isZero(const real_t& a,const real_t& tresh){
		return abs(a)<=tresh? 1:0;
	}

	static real_t getSqrt(const real_t& a){
		return a>=0? sqrt(a): 0;
	}

	// read matrices
	static int_t readFromFile(real_t* data, int_t nrow, int_t ncol, const char* datafilename);
	static int_t readFromFile(int_t* data, int_t n, const char* datafilename);
	static int_t readFromFile(real_t* data, int_t n, const char* datafilename);

	// used for timing the code
	static real_t getCPUtime();
	
	// multiplication of the matrices A and B
	static void MatrixMult(const real_t* matA, const real_t* matB, real_t* matC, \
								int_t rowsA, int_t colsA, int_t colsB);

	// Multiplies matA with vec1
	static void MatVecMult(const real_t* matA, const real_t* vec1, real_t* vec2, \
						const int_t rowsA, const int_t colsA);
	
	// Multiplies transpose of matA with vec1
	static void MatTVecMult(const real_t* const matA, const real_t* const vec1, real_t* const vec2,
						 const int_t rowsA, const int_t colsA);
	
	// B = AT (matrices stored row wise)
	static void MatrixTranspose(const real_t* matA, real_t* matB,const int_t rowsA,const int_t colsA);

	// solves x = inv(A)*b, where A is lower triangular
	static void ForwardSubstitution(const real_t* matA,  real_t* v1,const int_t rowsA, const int_t colsA);

	// solves x = inv(A)*b, where A is upper triangular
	static void BackwardSubstitution(const real_t* matA, real_t* v1, const int_t colsA, const int_t activeCols);

	//  min value in a vector
	static real_t min_value(const real_t* v1, const int_t n);
	
	//  max value in a vector
	static real_t max_value(const real_t* v1, const int_t n);

	//  min value and index in a vector
	static void min_value_idx(const real_t* v1, real_t &val, int_t &idx, const int_t n);

	//  max value and index in a vector
	static void max_value_idx(const real_t* v1, real_t &val, int_t &idx, const int_t n);

	// copy MatrixA into MatrixB
	static void MatrixCopy(const real_t* MatA, real_t* MatB, int_t rowsA, int_t colsA);

	// find infinity norm of vector
	static real_t VectorInfNorm(const real_t* v1,  const int_t n );

	// perform sequential search in a vector for val
	static bool SeqSearch(const int_t *v1, const int_t val, const int_t n);

	// multiply vector by a scalar
	static void ScalarVectorMult(real_t *vec, const real_t a, const int_t n);

	// check if there is a positive element
	static bool anyPositive(const real_t* const vec1, const int_t n);

	static int_t absolute(const int_t a) {
		if (a > 0)	{return a;}
		else		{return -a;}
	};
};