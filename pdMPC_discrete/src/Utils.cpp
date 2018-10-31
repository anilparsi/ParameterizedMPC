#include "Utils.h"
#include <stdio.h>
#include <cassert>


void Utils::LoadSparseMatrixRowPlain(const char* str,int_t*& row_ptr, int_t*& col_idx, real_t*& val,
								int_t& nr, int_t& nc, int_t& nnz){
		std::string filename_base=str;
		filename_base.append("_n.txt");

		int_t size[4];
		readFromFile(size,4,filename_base.c_str());
		
		nr=size[0];
		nc=size[1];
		int_t row_ptr_length=size[2];
		nnz=size[3];

		// read column pointer
		filename_base=str;
		filename_base.append("_ia.txt");
		row_ptr=new int_t[row_ptr_length];
		readFromFile(row_ptr,row_ptr_length,filename_base.c_str());


		// read row indices
		filename_base=str;
		filename_base.append("_ja.txt");
		col_idx=new int_t[nnz];
		readFromFile(col_idx,nnz,filename_base.c_str());


		filename_base=str;
		filename_base.append("_ra.txt");
		val=new real_t[nnz];
		readFromFile(val,nnz,filename_base.c_str());

	}

void Utils::LoadVec(const char* str, real_t** vec, int_t& nv){
		std::string filename_base=str;
		filename_base.append(".txt");

		int_t size[1];
		readFromFile(size,1,filename_base.c_str());

		nv=size[0];
		real_t* tmp=new real_t[nv+1];
		
		readFromFile(tmp,nv+1,filename_base.c_str());
		*vec=new real_t[nv];

		for(int_t k=1;k<nv+1;++k){
			(*vec)[k-1]=tmp[k];
		}

		delete [] tmp;
	}

void Utils::LoadVec(const char* str, int_t** vec, int_t& nv){
		std::string filename_base=str;
		filename_base.append(".txt");

		int_t size[1];
		readFromFile(size,1,filename_base.c_str());

		nv=size[0];
		int_t* tmp=new int_t[nv+1];
		
		readFromFile(tmp,nv+1,filename_base.c_str());
		*vec=new int_t[nv];

		for(int_t k=1;k<nv+1;++k){
			(*vec)[k-1]=tmp[k];
		}

		delete [] tmp;
	}

void Utils::PrintVec(const real_t* v, int_t nv){
		for(int_t k=0;k<nv;++k){
			printf("%f\n",v[k]);
		}
	}

void Utils::PrintVec(const int_t* v, int_t nv){
		for(int_t k=0;k<nv;++k){
			printf("%i\n",v[k]);
		}
	}

// modified definition and parameter order: Anil Parsi
void Utils::PrintMat(const real_t* val, int_t nr, int_t nc){
		for(int_t j=0;j<nr;++j){
			for(int_t k=0;k<nc;++k){
				printf("%f\t",val[k*nr+j]);
			}
			printf("\n");
		}
	}


	// Important; The function assumes that col_idx is sorted (ascending)!
	// Moreover, col_ptr must be an array of lenght nCols+1 (already initialized).
void Utils::ConvertCOSToCCS(const int_t* col_idx, const int_t& nnZ,int_t* col_ptr){
		int_t oldCol=0;
		col_ptr[0]=0;
		for(int_t k=0;k<nnZ;++k){
			while(oldCol!=col_idx[k]){
				col_ptr[++oldCol]=k;
			}
		}
		col_ptr[oldCol+1]=nnZ;
	}

void Utils::DotProduct(const real_t* a, const real_t* b, const int& n, real_t& res){
		res=0;
		for(int_t i=0;i<n;++i){
			res+=a[i]*b[i];
		}
	}


void Utils::VectorAdd(const real_t* a, const real_t* b, real_t* c, int_t n){
		for(int_t k=0;k<n;++k){
			c[k]=a[k]+b[k];
		}
	}

void Utils::VectorSubstract(const real_t* a, const real_t* b, real_t* c, int_t n){
		for(int_t k=0;k<n;++k){
			c[k]=a[k]-b[k];
		}
	}

void Utils::VectorAddMultiply(const real_t* a, real_t alpha, const real_t* b, real_t beta, real_t* c, int_t n){
		for(int_t k=0;k<n;++k){
			c[k]=alpha*a[k]+beta*b[k];
		}
	}

void Utils::VectorCopy(const real_t* a, real_t* b, int_t n){
		for(int_t k=0;k<n;++k){
			b[k]=a[k];
		}
	}


real_t Utils::VectorNorm(const real_t* a, int_t n){
		real_t res=0;
		for(int_t k=0;k<n;++k){
			res+=a[k]*a[k];
		}
		return getSqrt(res);
	}

	// returns ||a-b||_2
real_t Utils::VectorNormDiff(const real_t* a, const real_t* b, int_t n){
		real_t res=0;
		for(int_t k=0;k<n;++k){
			res+=(a[k]-b[k])*(a[k]-b[k]);
		}
		return Utils::getSqrt(res);
	}

	// component wise multiplication of the vectors a and b
void Utils::VectorMult(const real_t* a, const real_t* b, real_t* c, int_t n){
		for(int_t k=0;k<n;++k){
			c[k]=a[k]*b[k];
		}
	}


int_t Utils::factorial(int_t n){
		return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
	}



	/*
	 *	r e a d F r o m F i l e
	 */
int_t Utils::readFromFile(real_t* data, int_t n, const char* datafilename){
		return readFromFile( data, n, 1, datafilename );
	}


	
int_t Utils::readFromFile(int_t* data, int_t n, const char* datafilename){
		
		int_t i;
		FILE* datafile;

		if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
		{
			printf("\n\runable to read file %s\n",datafilename);
			return -1;
		}

		for(i=0; i<n; ++i)
		{
			if(fscanf( datafile, "%d\n", &(data[i]) ) == 0)
			{
				fclose( datafile );
				printf("\n\runable to read file %s\n",datafilename);
				return -1;
			}
		}

		fclose( datafile );

		return 1;
	}

	
	/*
	*	This definition is never directly used: called with one column from above. 
	*/
int_t Utils::readFromFile(real_t* data, int_t nrow, int_t ncol, const char* datafilename){
		int_t i, j;
		real_t float_data;
		FILE* datafile;

		if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
		{
			printf("\n\runable to read file %s\n",datafilename);
			return -1;
		}

		for( i=0; i<nrow; ++i )
		{
			for( j=0; j<ncol; ++j )
			{
				
				if ( fscanf( datafile, "%lf ", &float_data ) == 0 )
				{
					fclose( datafile );
					printf("\n\runable to read file %s\n",datafilename);
					return -1;
				}
				data[i*ncol + j] = ( (real_t) float_data );
			}
		}

		fclose( datafile );

		return -1;
}

	// Multiplies two matrices A and B (matrices stored row wise)
void Utils::MatrixMult(const real_t* matA, const real_t* matB, real_t* matC, \
						const int_t rowsA, const int_t colsA, const int_t colsB){
	
	for (int_t i = 0; i<rowsA; ++i){
		for (int_t j = 0; j<colsB; ++j){
			matC[i*colsB + j] = 0;
			for (int_t k = 0; k<colsA; ++k){
				matC[i*colsB +j] += matA[i*colsA+k]*matB[k*colsB+j];				
			}
		}
	}
	
}
	// Multiplies matA with vec1, with rowsA rows
void Utils::MatVecMult(const real_t* matA, const real_t* vec1, real_t* vec2, \
						const int_t rowsA, const int_t colsA){
	for (int_t i = 0; i<rowsA; ++i){
		vec2[i] = 0;
		for (int_t k = 0; k<colsA; ++k){
			vec2[i] += matA[i*colsA+k]*vec1[k];
		}
	}
	
}

	// B = AT (matrices stored row wise)
void Utils::MatrixTranspose(const real_t* matA, real_t* matB,const int_t rowsA,const int_t colsA){
	
	for (int_t j = 0; j<colsA; ++j){
		for (int_t i = 0; i<rowsA; ++i){
			matB[j*rowsA+i]=matA[i*colsA+j];
		}
	}
}

void Utils::ForwardSubstitution(const real_t* matA, real_t* v1, const int_t rowsA, const int_t colsA){
	// This function assumes matA is lower triangular, and # of rows <= # of cols
	// It updates the same vector v1 with new values
	for (int i=0; i<rowsA; ++i){
		for (int j=0;j<i; ++j){
			v1[i] -= matA[i*colsA+j]*v1[j];
		}
		v1[i] = v1[i]/matA[i*colsA+i];
	}
}
void Utils::BackwardSubstitution(const real_t* matA, real_t* v1, const int_t colsA, const int_t activeCols){
	// This function assumes matA is upper triangular, and # of cols <= # of rows
	// It updates the same vector v1 with new values
	for (int i=activeCols-1; i>=0; --i){
		for (int j=activeCols-1;j>i; --j){
			v1[i] -= matA[i*colsA+j]*v1[j];
		}
		v1[i] = v1[i]/matA[i*colsA+i];
	}
}

real_t Utils::min_value(const real_t* v1, const int_t n){
	real_t val = INFVAL;
	for (int i=0; i<n; ++i){
		(v1[i]<val)?val = v1[i]:0;
	}
	return val;
}

real_t Utils::max_value(const real_t* v1, const int_t n){
	real_t val = -INFVAL;
	for (int i=0; i<n; ++i){
		(v1[i]>val)?val = v1[i]:0;
	}
	return val;
}

void Utils::min_value_idx(const real_t* v1, real_t &val, int_t &idx, const int_t n){
	val = INFVAL;
	idx = 0;
	for (int i=0; i<n; ++i){
		if(v1[i]<val){
			val = v1[i];
			idx = i;
		}
	}
}

void Utils::max_value_idx(const real_t* v1, real_t &val, int_t &idx, const int_t n){
	val = -INFVAL;
	idx = 0;
	for (int i=0; i<n; ++i){
		if(v1[i]>val){
			val = v1[i];
			idx = i;
		}
	}
}
void Utils::MatrixCopy(const real_t* MatA, real_t* MatB, int_t rowsA, int_t colsA){
	for(int_t i=0;i<rowsA;++i){
		for(int_t j=0; j<colsA; ++j){
			MatB[i*colsA+j] = MatA[i*colsA+j];
		}
	}
}

real_t Utils::VectorInfNorm(const real_t* v1,  const int_t n ){
	real_t abs;
	real_t norm  = 0;
	for(int i = 0; i<n; ++i){
		abs = (v1[i]>0)?v1[i]:-v1[i];
		norm = (norm>abs)?norm:abs;
	}
	return norm;

}

bool Utils::SeqSearch(const int_t *v1, const int_t val, const int_t n){
	for(int_t i=0; i<n; ++i){
		if(val==v1[i]){
			return true;
			break;
		}
	}
	return false;
}


void Utils::ScalarVectorMult(real_t *vec, const real_t a, const int_t n){
	for(int_t i = 0; i<n; ++i){
		vec[i] = a*vec[i];
	}
}

bool Utils::anyPositive(const real_t* const vec1, const int_t n){
	for(int_t i = 0; i<n; ++i){
		if(vec1[i]>0){
			return true;		
		}
	}
	return false;
}

void  Utils::MatTVecMult(const real_t* const matA, const real_t* const vec1, real_t* const vec2,
						 const int_t rowsA, const int_t colsA) {
	
	// vec2 = matA^T * vec1;
	// vec2 must be a zero vector
	for (int_t k = 0; k<colsA; ++k){
		vec2[k] = 0.0;
	}
	
	for (int_t i = 0; i<rowsA; ++i){		
		for (int_t k = 0; k<colsA; ++k){	
			vec2[k] += matA[i*colsA+k]*vec1[i]; 
		}

	}

}