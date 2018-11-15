#pragma once
#include <string>
#include "DefineSettings.h"

/*!
 * \brief This class contains all the utilites used for the toolbox.
 */
class Utils{
	public:
	
	/*!
	 * \brief function to load data from .txt files
	 *
	 * \param str is the path to the file
	 * \param vec is the address of the pointer where the data must be stored
	 * \param nv is the size of the vector which is being loaded
	 */
	static void LoadVec(const char* str, real_t** vec, int_t& nv);
		
	/*!
	 * \brief function to load data from .txt files
	 *
	 * \param str is the path to the file
	 * \param vec is the address of the pointer where the data must be stored
	 * \param nv is the size of the vector which is being loaded
	 */
	static void LoadVec(const char* str, int_t** vec, int_t& nv);


	/*!
	 * \brief performs dot product a.b = res, where a and b are n-dimensional
	 */
	static void DotProduct(const real_t* a, const real_t* b, const int& n, real_t& res);

	
	/*!
	 * \brief implements a+b=c, where a,b,c are n-dimensional
	 */
	static void VectorAdd(const real_t* a, const real_t* b, real_t* c, int_t n);

	/*!
	 * \brief implements a-b=c, where a,b,c, are n-dimensional
	 */
	static void VectorSubstract(const real_t* a, const real_t* b, real_t* c, int_t n);

	/*!
	 * \brief implements a*alpha+b*beta=c, where a,b,c are n-dimensional, and alpha and beta are scalars.
	 */
	static void VectorAddMultiply(const real_t* a, real_t alpha, const real_t* b, real_t beta, real_t* c, int_t n);

	/*!
	 * \brief implements b=a
	 */
	static void VectorCopy(const real_t* a, real_t* b, int_t n);

	/// returns the 2-norm of the vector a which is n dimensional
	static real_t VectorNorm(const real_t* a, int_t n);

	/// returns ||a-b||_2
	static real_t VectorNormDiff(const real_t* a, const real_t* b, int_t n);
	
	/// component wise multiplication of the vectors a and b
	static void VectorMult(const real_t* a, const real_t* b, real_t* c, int_t n);	

	/// returns the square root of a
	static real_t getSqrt(const real_t& a){
		return a>=0? sqrt(a): 0;
	}

	/// read matrix from file
	static int_t readFromFile(int_t* data, int_t n, const char* datafilename);

	/// read matrix from file
	static int_t readFromFile(real_t* data, int_t n, const char* datafilename);
	
	/// multiplication of the matrices A and B
	static void MatrixMult(const real_t* matA, const real_t* matB, real_t* matC, \
								int_t rowsA, int_t colsA, int_t colsB);

	/// performs vec2 = matA*vec1
	static void MatVecMult(const real_t* matA, const real_t* vec1, real_t* vec2, \
						const int_t rowsA, const int_t colsA);
	
	/// performs vec2 = matA^T*vec1
	static void MatTVecMult(const real_t* const matA, const real_t* const vec1, real_t* const vec2,
						 const int_t rowsA, const int_t colsA);
	
	/// B = A^T (matrices stored row wise)
	static void MatrixTranspose(const real_t* matA, real_t* matB,const int_t rowsA,const int_t colsA);

	/// solves x = inv(A)*b, where A is lower triangular
	static void ForwardSubstitution(const real_t* matA,  real_t* v1,const int_t rowsA, const int_t colsA);

	/// solves x = inv(A)*b, where A is upper triangular
	static void BackwardSubstitution(const real_t* matA, real_t* v1, const int_t colsA, const int_t activeCols);

	///  returns the minimum value in the vector v1
	static real_t min_value(const real_t* v1, const int_t n);
	
	///  returns the maximum value in the vector v1
	static real_t max_value(const real_t* v1, const int_t n);

	///  returns the minimum value and index in the vector v1
	static void min_value_idx(const real_t* v1, real_t &val, int_t &idx, const int_t n);

	///  returns the maximum value and index in the vector v1
	static void max_value_idx(const real_t* v1, real_t &val, int_t &idx, const int_t n);

	/// copy MatrixA into MatrixB
	static void MatrixCopy(const real_t* MatA, real_t* MatB, int_t rowsA, int_t colsA);

	/// return infinity norm of vector v1
	static real_t VectorInfNorm(const real_t* v1,  const int_t n );

	/// perform sequential search in a vector for val
	static bool SeqSearch(const int_t *v1, const int_t val, const int_t n);

	/// multiply vector by a scalar
	static void ScalarVectorMult(real_t *vec, const real_t a, const int_t n);

	/// check existence of a positive element in the vector vec1
	static bool anyPositive(const real_t* const vec1, const int_t n);

	/// return the absolute value of a
	static int_t absolute(const int_t a) {
		if (a > 0)	{return a;}
		else		{return -a;}
	};
};