#pragma once
#include "shared.h"
#include "gridmath.h"

class MatCell {
public:
	Float diag, plusi, plusj, plusk;
	MatCell() :diag(0), plusi(0), plusj(0), plusk(0) {}
};

class PressureSolver {
public:
	aryi idx;
	vector<int> cells;
	vector<double> div, pressure;
	vector<MatCell> coefs;
	void init(int n, int m, int w) {
		idx.init(n, m, w);
	}
	~PressureSolver() {
		idx.clear();
	}
	void Send_Back_To(aryf &pa);//send pressure back to a 3d-field
	void Add_Scaled_Vector(vector<double>& v1, const vector<double>& v2, double scale);
	void Add_Scaled_Vectors(vector<double>& v1, double s1, vector<double>& v2, double s2, vector<double>& result);
	int Get_Cell_Idx(int i, int j, int k);
	void Build_Water_Index(const aryi &mask);
	void Calc_Neg_Divergence(const aryf & vx, const aryf & vy, const aryf & vz, const aryi & mask);
	//void Calc_Neg_Divergence(const aryf &vx, const aryf &vy, const aryf &vz);
	void Calc_Matrix_Coeffs(vector<MatCell>& A, const aryi & mask);
	void Calc_Preconditioner(vector<MatCell>& A, vector<double>& precon);
	void Apply_Preconditioner(vector<MatCell>& A, vector<double>& precon, vector<double>& residual, vector<double>& vect);
	void Apply_Matrix(vector<MatCell>& A, vector<double>& x, vector<double>& result);
	void Solver_LA_System(vector<MatCell>& A, vector<double>& b, vector<double>& precon, vector<double>& pressure);
	void Solve_Pressure(const aryf & vx, const aryf & vy, const aryf & vz, const aryi & mask);
};