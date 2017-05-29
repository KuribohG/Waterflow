#include "pressuresolver.h"

double Maximum_Abs(vector<double> &d) {
	double ret = -INF;
	for (double a : d) ret = max(ret, fabs(a));
	return ret;
}

double Dot(const vector<double> &a, const vector<double> &b) {
	assert(a.size() == b.size());
	double ret = 0;
	for (int i = 0; i < a.size(); i++) {
		ret += a[i] * b[i];
	}
	/*if (ret == 0) {
		for (int i = 0; i < 10; i++) cout << a[i] << " "; cout << endl;
		for (int i = 0; i < 10; i++) cout << b[i] << " "; cout << endl;
	}*/
	assert(ret != 0);
	return ret;
}

// v1 += v2*scale
void PressureSolver::Add_Scaled_Vector(vector<double> &v1, vector<double> &v2, double scale) {
	assert(v1.size() == v2.size());
	assert(scale == scale);
	for (unsigned int idx = 0; idx < v1.size(); idx++) {
		v1[idx] += v2[idx] * scale;
	}
}

// result = v1*s1 + v2*s2
void PressureSolver::Add_Scaled_Vectors(vector<double> &v1, double s1, vector<double> &v2, double s2, vector<double> &result) {
	assert(v1.size() == v2.size() && v2.size() == result.size());
	assert(s1 == s1);
	assert(s2 == s2);
	for (unsigned int idx = 0; idx < v1.size(); idx++) {
		result[idx] = v1[idx] * s1 + v2[idx] * s2;
	}
}

int PressureSolver::Get_Cell_Idx(int i, int j, int k){
	if (!idx.inside(i, j, k)) return -1;
	return idx.get(i, j, k);
}

void PressureSolver::Send_Back_To(aryf & pa) {
	// todo: consider solid pressure
	pa.set(0);
	for (int t = 0; t < pressure.size(); t++) {
		int i = cells[t * 3 + 0], j = cells[t * 3 + 1], k = cells[t * 3 + 2];
		pa(i, j, k) = pressure[t];
	}
}

void PressureSolver::Build_Water_Index(const aryi & mask){
	assert(mask.n == idx.n);
	assert(mask.m == idx.m);
	assert(mask.w == idx.w);
	idx.set(-1);
	cells.clear();
	int p = 0;
	for (int i = 0; i < mask.n; i++) {
		for (int j = 0; j < mask.m; j++) {
			for (int k = 0; k < mask.w; k++) {
				if (mask.get(i, j, k) == WATER) {
					idx(i, j, k) = p++;
					cells.push_back(i);
					cells.push_back(j);
					cells.push_back(k);
				}
			}
		}
	}
}

void PressureSolver::Calc_Neg_Divergence(const aryf & vx, const aryf & vy, const aryf & vz, const aryi &mask){
	for (int t = 0; t < div.size(); t++) {
		int i = cells[t * 3 + 0], j = cells[t * 3 + 1], k = cells[t * 3 + 2];
		div[t]=-(
			+ vx.get(i + 1, j, k) - vx.get(i, j, k)
			+ vy.get(i, j + 1, k) - vy.get(i, j, k)
			+ vz.get(i, j, k + 1) - vz.get(i, j, k)
		);
		if (mask.is(i - 1, j, k, SOLID)) div[t] -= vx.get(i, j, k);
		if (mask.is(i + 1, j, k, SOLID)) div[t] += vx.get(i + 1, j, k);
		if (mask.is(i, j - 1, k, SOLID)) div[t] -= vy.get(i, j, k);
		if (mask.is(i, j + 1, k, SOLID)) div[t] += vy.get(i, j + 1, k);
		if (mask.is(i, j, k - 1, SOLID)) div[t] -= vz.get(i, j, k);
		if (mask.is(i, j, k + 1, SOLID)) div[t] += vz.get(i, j, k + 1);
	}
}

int Nonsolid_Neighbors(const aryi &mask, int i, int j, int k) {
	int n = 0;
	if (!mask.is(i - 1, j, k, SOLID)) n++;
	if (!mask.is(i + 1, j, k, SOLID)) n++;
	if (!mask.is(i, j - 1, k, SOLID)) n++;
	if (!mask.is(i, j + 1, k, SOLID)) n++;
	if (!mask.is(i, j, k - 1, SOLID)) n++;
	if (!mask.is(i, j, k + 1, SOLID)) n++;
	return n;
}

void PressureSolver::Calc_Matrix_Coeffs(vector<MatCell> &A, const aryi &mask) {
	for (int t = 0; t < pressure.size(); t++) {
		int i = cells[t * 3 + 0], j = cells[t * 3 + 1], k = cells[t * 3 + 2];
		A[t].diag = Nonsolid_Neighbors(mask, i, j, k);
		if (mask.is(i + 1, j, k, WATER)) A[t].plusi = 1;
		if (mask.is(i, j + 1, k, WATER)) A[t].plusj = 1;
		if (mask.is(i, j, k + 1, WATER)) A[t].plusk = 1;
	}
}

void PressureSolver::Calc_Preconditioner(vector<MatCell> &A, vector<double> &precon) {
	assert(A.size() == precon.size());
	assert(pressure.size() == precon.size());

	double scale = TIME_DELTA / (DENSITY);
	double negscale = -scale;

	double tau = 0.97;      // Tuning constant
	double sigma = 0.25;    // safety constant
	for (int t = 0; t < pressure.size(); t++) {
		int i = cells[t * 3 + 0], j = cells[t * 3 + 1], k = cells[t * 3 + 2];
		int vidx = t;
		int vidx_im1 = Get_Cell_Idx(i - 1, j, k);
		int vidx_jm1 = Get_Cell_Idx(i, j - 1, k);
		int vidx_km1 = Get_Cell_Idx(i, j, k - 1);

		double diag = (double)A[vidx].diag*scale;

		double plusi_im1 = vidx_im1 != -1 ? (double)A[vidx_im1].plusi * negscale : 0.0;
		double plusi_jm1 = vidx_jm1 != -1 ? (double)A[vidx_jm1].plusi * negscale : 0.0;
		double plusi_km1 = vidx_km1 != -1 ? (double)A[vidx_km1].plusi * negscale : 0.0;

		double plusj_im1 = vidx_im1 != -1 ? (double)A[vidx_im1].plusj * negscale : 0.0;
		double plusj_jm1 = vidx_jm1 != -1 ? (double)A[vidx_jm1].plusj * negscale : 0.0;
		double plusj_km1 = vidx_km1 != -1 ? (double)A[vidx_km1].plusj * negscale : 0.0;

		double plusk_im1 = vidx_im1 != -1 ? (double)A[vidx_im1].plusk * negscale : 0.0;
		double plusk_jm1 = vidx_jm1 != -1 ? (double)A[vidx_jm1].plusk * negscale : 0.0;
		double plusk_km1 = vidx_km1 != -1 ? (double)A[vidx_km1].plusk * negscale : 0.0;

		double precon_im1 = vidx_im1 != -1 ? precon[vidx_im1] : 0.0;
		double precon_jm1 = vidx_jm1 != -1 ? precon[vidx_jm1] : 0.0;
		double precon_km1 = vidx_km1 != -1 ? precon[vidx_km1] : 0.0;

		double v1 = plusi_im1 * precon_im1;
		double v2 = plusj_jm1 * precon_jm1;
		double v3 = plusk_km1 * precon_km1;
		double v4 = precon_im1 * precon_im1;
		double v5 = precon_jm1 * precon_jm1;
		double v6 = precon_km1 * precon_km1;

		double e = diag - v1*v1 - v2*v2 - v3*v3 -
			tau*(plusi_im1*(plusj_im1 + plusk_im1)*v4 +
				plusj_jm1*(plusi_jm1 + plusk_jm1)*v5 +
				plusk_km1*(plusi_km1 + plusj_km1)*v6);

		if (e < sigma*diag) {
			e = diag;
		}

		if (fabs(e) > 10e-9) {
			precon[vidx] = 1.0 / sqrt(e);
		}
	}
}

void PressureSolver::Apply_Preconditioner(vector<MatCell> &A, vector<double> &precon, vector<double> &residual, vector<double> &vect) {

	double scale = TIME_DELTA / (DENSITY);
	double negscale = -scale;

	// Solve A*q = residual
	vector<double> q(A.size());
	for (int idx = 0; idx < pressure.size(); idx++) {
		int i = cells[idx * 3 + 0], j = cells[idx * 3 + 1], k = cells[idx * 3 + 2];
		int vidx = idx;

		int vidx_im1 = Get_Cell_Idx(i - 1, j, k);
		int vidx_jm1 = Get_Cell_Idx(i, j - 1, k);
		int vidx_km1 = Get_Cell_Idx(i, j, k - 1);

		double plusi_im1 = 0.0;
		double precon_im1 = 0.0;
		double q_im1 = 0.0;
		if (vidx_im1 != -1) {
			plusi_im1 = (double)A[vidx_im1].plusi * negscale;
			precon_im1 = precon[vidx_im1];
			q_im1 = q[vidx_im1];
		}

		double plusj_jm1 = 0.0;
		double precon_jm1 = 0.0;
		double q_jm1 = 0.0;
		if (vidx_jm1 != -1) {
			plusj_jm1 = (double)A[vidx_jm1].plusj * negscale;
			precon_jm1 = precon[vidx_jm1];
			q_jm1 = q[vidx_jm1];
		}

		double plusk_km1 = 0.0;
		double precon_km1 = 0.0;
		double q_km1 = 0.0;
		if (vidx_km1 != -1) {
			plusk_km1 = (double)A[vidx_km1].plusk * negscale;
			precon_km1 = precon[vidx_km1];
			q_km1 = q[vidx_km1];
		}

		double t = residual[vidx] - plusi_im1 * precon_im1 * q_im1 -
			plusj_jm1 * precon_jm1 * q_jm1 -
			plusk_km1 * precon_km1 * q_km1;

		t = t*precon[vidx];
		q[vidx] = t;
	}

	// Solve transpose(A)*z = q
	for (int idx = (int)pressure.size() - 1; idx >= 0; idx--) {
		int i = cells[idx * 3 + 0], j = cells[idx * 3 + 1], k = cells[idx * 3 + 2];
		int vidx = idx;


		int vidx_ip1 = Get_Cell_Idx(i + 1, j, k);
		int vidx_jp1 = Get_Cell_Idx(i, j + 1, k);
		int vidx_kp1 = Get_Cell_Idx(i, j, k + 1);

		double vect_ip1 = vidx_ip1 != -1 ? vect[vidx_ip1] : 0.0;
		double vect_jp1 = vidx_jp1 != -1 ? vect[vidx_jp1] : 0.0;
		double vect_kp1 = vidx_kp1 != -1 ? vect[vidx_kp1] : 0.0;

		double plusi = (double)A[vidx].plusi * negscale;
		double plusj = (double)A[vidx].plusj * negscale;
		double plusk = (double)A[vidx].plusk * negscale;

		double preconval = precon[vidx];
		double t = q[vidx] - plusi * preconval * vect_ip1 -
			plusj * preconval * vect_jp1 -
			plusk * preconval * vect_kp1;

		t = t*preconval;
		vect[vidx] = t;
	}
}

void PressureSolver::Apply_Matrix(vector<MatCell> &A, vector<double> &x, vector<double> &result) {
	assert(A.size() == x.size() && x.size() == result.size());

	double scale = TIME_DELTA / (DENSITY);
	double negscale = -scale;

	for (unsigned int idx = 0; idx < A.size(); idx++) {
		int i = cells[idx * 3 + 0], j = cells[idx * 3 + 1], k = cells[idx * 3 + 2];

		// val = dot product of column vector x and idxth row of matrix A
		double val = 0.0;
		int vidx = Get_Cell_Idx(i - 1, j, k);
		if (vidx != -1) { val += x[vidx]; }

		vidx = Get_Cell_Idx(i + 1, j, k);
		if (vidx != -1) { val += x[vidx]; }

		vidx = Get_Cell_Idx(i, j - 1, k);
		if (vidx != -1) { val += x[vidx]; }

		vidx = Get_Cell_Idx(i, j + 1, k);
		if (vidx != -1) { val += x[vidx]; }

		vidx = Get_Cell_Idx(i, j, k - 1);
		if (vidx != -1) { val += x[vidx]; }

		vidx = Get_Cell_Idx(i, j, k + 1);
		if (vidx != -1) { val += x[vidx]; }

		val *= negscale;

		vidx = Get_Cell_Idx(i, j, k);
		val += (double)A[vidx].diag * scale * x[vidx];

		result[vidx] = val;
	}
}

// Solve (A*pressure = b) with Modified Incomplete Cholesky 
// Conjugate Gradient method (MICCG(0))
void PressureSolver::Solver_LA_System(vector<MatCell> &A, vector<double> &b, vector<double> &precon, vector<double> &pressure) {
	double tol = TOLERANCE;
	if (Maximum_Abs(b) < tol) {
		return;
	}
	vector<double> residual(b);
	vector<double> auxillary(pressure.size());

	Apply_Preconditioner(A, precon, residual, auxillary);

	vector<double> search(auxillary);

	double alpha = 0.0;
	double beta = 0.0;
	double sigma = Dot(auxillary,residual);
	double sigmaNew = 0.0;
	int iterationNumber = 0;

	while (iterationNumber < LINSOLVER_ITER) {
		Apply_Matrix(A, search, auxillary);
		alpha = sigma / Dot(auxillary, search);
		Add_Scaled_Vector(pressure, search, alpha);
		Add_Scaled_Vector(residual, auxillary, -alpha);

		if (Maximum_Abs(residual) < tol) {
			printf("CG Iterations: %d\n", iterationNumber);
			return;
		}

		Apply_Preconditioner(A, precon, residual, auxillary);
		sigmaNew = Dot(auxillary, residual);
		beta = sigmaNew / sigma;
		Add_Scaled_Vectors(auxillary, 1.0, search, beta, search);
		sigma = sigmaNew;

		iterationNumber++;

		/*if (iterationNumber % 10 == 0) {
			cerr << "\tIteration #: " << iterationNumber <<
				"\tEstimated Error: " << Maximum_Abs(residual) << std::endl;
		}*/
	}
}

void PressureSolver::Solve_Pressure(const aryf &vx, const aryf &vy, const aryf &vz, const aryi &mask) {
	int t0 = clock();
	Build_Water_Index(mask);
	pressure = vector<double>(cells.size() / 3);
	div = vector<double>(pressure.size());
	Calc_Neg_Divergence(vx, vy, vz, mask);
	if (Maximum_Abs(div) < TOLERANCE) return; //with pressure=0
	vector<MatCell> A(pressure.size());
	vector<double> precon(pressure.size());
	Calc_Matrix_Coeffs(A, mask);
	Calc_Preconditioner(A, precon);
	Solver_LA_System(A, div, precon, pressure);
	printf("pressure solve done, time cost: %.2fs\n", (clock() - t0 + 0.0) / CLOCKS_PER_SEC);
}