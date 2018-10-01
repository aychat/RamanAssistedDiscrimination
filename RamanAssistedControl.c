#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <nlopt.h>

typedef double complex cmplx;

typedef struct parameters{
    double* time;
    cmplx* rho_0;

    double A_R;
    double width_R;
    double t0_R;

    double A_EE;
    double width_EE;
    double t0_EE;

    double w_R;
    double w_v;
    double w_EE;

    int nDIM;
    int timeDIM;

    cmplx* field_out;

    double* lower_bounds;
    double* upper_bounds;
    double* guess;

    int MAX_EVAL;

} parameters;

typedef struct molecule{
    int nDIM;
    double* energies;
    double* gamma_decay;
    double* gamma_pure_dephasing;
    cmplx* mu;

    cmplx* rho;
    cmplx* dyn_rho;
    cmplx* rho_0;

    double* abs_spectra;
    double* vib_spectra;
} molecule;

typedef struct mol_system{
    molecule* moleculeA;
    molecule* moleculeB;
    parameters* params;
} mol_system;


//====================================================================================================================//
//                                                                                                                    //
//                                        AUXILIARY FUNCTIONS FOR MATRIX OPERATIONS                                   //
//                                                                                                                    //
//====================================================================================================================//


void print_complex_mat(cmplx *A, int nDIM)
//----------------------------------------------------//
// 	          PRINTS A COMPLEX MATRIX                 //
//----------------------------------------------------//
{
	int i,j;
	for(i=0; i<nDIM; i++)
	{
		for(j=0; j<nDIM; j++)
		{
			printf("%3.3e + %3.3eJ  ", creal(A[i * nDIM + j]), cimag(A[i * nDIM + j]));
		}
	    printf("\n");
	}
	printf("\n\n");
}

void print_complex_vec(cmplx *A, int vecDIM)
//----------------------------------------------------//
// 	          PRINTS A COMPLEX MATRIX                 //
//----------------------------------------------------//
{
	int i;
	for(i=0; i<vecDIM; i++)
	{
		printf("%3.3e + %3.3eJ  ", creal(A[i]), cimag(A[i]));
	}
	printf("\n");
}

void print_double_mat(double *A, int nDIM)
//----------------------------------------------------//
// 	            PRINTS A REAL MATRIX                  //
//----------------------------------------------------//
{
	int i,j;
	for(i=0; i<nDIM; i++)
	{
		for(j=0; j<nDIM; j++)
		{
			printf("%3.3e  ", A[i * nDIM + j]);
		}
	    printf("\n");
	}
	printf("\n\n");
}

void print_double_vec(double *A, int vecDIM)
//----------------------------------------------------//
// 	          PRINTS A COMPLEX MATRIX                 //
//----------------------------------------------------//
{
	int i;
	for(i=0; i<vecDIM; i++)
	{
		printf("%3.3e  ", A[i]);
	}
	printf("\n");
}


void copy_mat(const cmplx *A, cmplx *B, int nDIM)
//----------------------------------------------------//
// 	        COPIES MATRIX A ----> MATRIX B            //
//----------------------------------------------------//
{
    int i, j = 0;
    for(i=0; i<nDIM; i++)
    {
        for(j=0; j<nDIM; j++)
        {
            B[i * nDIM + j] = A[i * nDIM + j];
        }
    }
}

void add_mat(const cmplx *A, cmplx *B, int nDIM)
//----------------------------------------------------//
// 	        ADDS A to B ----> MATRIX B = A + B        //
//----------------------------------------------------//
{
    int i, j = 0;
    for(i=0; i<nDIM; i++)
    {
        for(j=0; j<nDIM; j++)
        {
            B[i * nDIM + j] += A[i * nDIM + j];
        }
    }
}


void scale_mat(cmplx *A, double factor, int nDIM)
//----------------------------------------------------//
// 	     SCALES A BY factor ----> MATRIX B = A + B    //
//----------------------------------------------------//
{
    for(int i=0; i<nDIM; i++)
    {
        for(int j=0; j<nDIM; j++)
        {
            A[i * nDIM + j] *= factor;
        }
    }
}


double complex_abs(cmplx z)
//----------------------------------------------------//
// 	            RETURNS ABSOLUTE VALUE OF Z           //
//----------------------------------------------------//
{

    return sqrt((creal(z)*creal(z) + cimag(z)*cimag(z)));
}


cmplx complex_trace(cmplx *A, int nDIM)
//----------------------------------------------------//
// 	                 RETURNS TRACE[A]                 //
//----------------------------------------------------//
{
    cmplx trace = 0.0 + I * 0.0;
    for(int i=0; i<nDIM; i++)
    {
        trace += A[i*nDIM + i];
    }
    printf("Trace = %3.3e + %3.3eJ  \n", creal(trace), cimag(trace));

    return trace;
}


void multiply_complex_mat(cmplx *product, const cmplx *A, const cmplx *B, const int nDIM)
//----------------------------------------------------//
// 	             RETURNS A*B MATRIX PRODUCT           //
//----------------------------------------------------//
{
    for (int i=0; i<nDIM; i++)
    {
        for (int j=0; j<nDIM; j++)
        {
            for (int k=0; k<nDIM; k++)
            {
                product[i*nDIM + j] += A[i*nDIM + k]*B[k*nDIM + j];
            }
        }
    }
}


void commute_complex_mat(cmplx *commutator, const cmplx *A, const cmplx *B, const int nDIM)
//-----------------------------------------------------------------//
// 	          RETURNS commutator = [A, B] MATRIX COMMUTATOR        //
//-----------------------------------------------------------------//
{
    for (int i=0; i<nDIM; i++)
    {
        for (int j=0; j<nDIM; j++)
        {
            commutator[i*nDIM + j] = 0. + I * 0.;
            for (int k=0; k<nDIM; k++)
            {
                commutator[i*nDIM + j] += A[i*nDIM + k]*B[k*nDIM + j] - B[i*nDIM + k]*A[k*nDIM + j];
            }
        }
    }
}


double complex_max_element(cmplx *A, int nDIM)
//----------------------------------------------------//
// 	   RETURNS ELEMENT WITH MAX ABSOLUTE VALUE        //
//----------------------------------------------------//
{
    double max_el = 0.0;
    int i, j = 0;
    for(i=0; i<nDIM; i++)
    {
        for(j=0; j<nDIM; j++)
        {
            if(complex_abs(A[i * nDIM + j]) > max_el)
            {
                max_el = complex_abs(A[i * nDIM + j]);
            }
        }
    }
    return max_el;
}


double integrate_Simpsons(double *f, double Tmin, double Tmax, int Tsteps)
{
    double dt = (Tmax - Tmin) / Tsteps;
    double integral = 0.0;
    int k;

    for(int tau=0; tau<Tsteps; tau++)
    {
        if(tau == 0) k=1;
	    else if(tau == (Tsteps-1)) k=1;
	    else if( tau % 2 == 0) k=2;
        else if( tau % 2 == 1) k=4;

        integral += f[tau] * k;
    }

    integral *= dt/3.;
    return integral;
}


//====================================================================================================================//
//                                                                                                                    //
//                                         FUNCTIONS FOR PROPAGATION STEP                                             //
//                                                                                                                    //
//====================================================================================================================//

void CalculateField(cmplx* field, parameters* params)
//----------------------------------------------------//
//   RETURNS THE ENTIRE FIELD AS A FUNCTION OF TIME   //
//----------------------------------------------------//
{
    int i;
    int nDIM = params->nDIM;
    int timeDIM = params->timeDIM;

    double* t = params->time;

    double A_R = params->A_R;
    double width_R = params->width_R;
    double t0_R = params->t0_R;

    double A_EE = params->A_EE;
    double width_EE = params->width_EE;
    double t0_EE = params->t0_EE;

    double w_R = params->w_R;
    double w_v = params->w_v;
    double w_EE = params->w_EE;


    for(i=0; i<timeDIM; i++)
    {
        field[i] = A_R * exp(-pow(t[i] - t0_R, 2) / (2. * pow(width_R, 2))) * (cos((w_R + w_v) * t[i]) + cos(w_R * t[i]))
        + A_EE * exp(-pow(t[i] - t0_EE, 2) / (2. * pow(width_EE, 2))) * cos(w_EE * t[i]);
        field[i] = A_EE * exp(-pow(t[i] - t0_EE, 2) / (2. * pow(width_EE, 2))) * cos(w_EE * t[i]);
    }
}

void L_operate(cmplx* Qmat, const cmplx field_ti, molecule* mol, parameters* params)
//----------------------------------------------------//
// 	    RETURNS Q <-- L[Q] AT A PARTICULAR TIME (t)   //
//----------------------------------------------------//
{
    int m, n, k;
    int nDIM = params->nDIM;
    double* gamma_pure_dephasing = mol->gamma_pure_dephasing;
    double* gamma_decay = mol->gamma_decay;
    cmplx* mu = mol->mu;
    double* energies = mol->energies;

    cmplx* Lmat = (cmplx*)calloc(nDIM * nDIM,  sizeof(cmplx));

    for(m = 0; m < nDIM; m++)
        {
        for(n = 0; n < nDIM; n++)
            {
                Lmat[m * nDIM + n] += - I * (energies[m] - energies[n]) * Qmat[m * nDIM + n];
                for(k = 0; k < nDIM; k++)
                {
                    Lmat[m * nDIM + n] -= 0.5 * (gamma_decay[n * nDIM + k] + gamma_decay[m * nDIM + k]) * Qmat[m * nDIM + n];
                    Lmat[m * nDIM + n] += I * field_ti * (mu[m * nDIM + k] * Qmat[k * nDIM + n] - Qmat[m * nDIM + k] * mu[k * nDIM + n]);

                    if (m == n)
                    {
//                        Lmat[m * nDIM + m] += gamma_decay[k * nDIM + m] * Qmat[k * nDIM + k];
                    }
                    else
                    {
                        Lmat[m * nDIM + n] -= gamma_pure_dephasing[m * nDIM + n] * Qmat[m * nDIM + n];
                    }


                }

            }

        }

    for(m = 0; m < nDIM; m++)
        {
        for(n = 0; n < nDIM; n++)
            {
                Qmat[m * nDIM + n] = Lmat[m * nDIM + n];
            }
        }
    free(Lmat);

}


void Propagate(molecule* mol, parameters* params)
//----------------------------------------------------------------------//
//    GETTING rho(T)_{k=[3,4]} FROM rho(0) USING PROPAGATE FUNCTION     //
//----------------------------------------------------------------------//
{
    int i, j, k;
    int tau_index, t_index;
    int nDIM = params->nDIM;
    int timeDIM = params->timeDIM;

    cmplx *rho_0 = mol->rho_0;
    double *time = params->time;

    cmplx* field = params->field_out;

    double dt = time[1] - time[0];

    cmplx* L_rho_func = (cmplx*)calloc(nDIM * nDIM, sizeof(cmplx));
    copy_mat(rho_0, L_rho_func, nDIM);
    copy_mat(rho_0, mol->rho, nDIM);

    for(t_index=0; t_index<timeDIM; t_index++)
    {
        k=1;
        do
        {
            L_operate(L_rho_func, field[t_index], mol, params);
            scale_mat(L_rho_func, dt/k, nDIM);
            add_mat(L_rho_func, mol->rho, nDIM);
            k+=1;
        }while(complex_max_element(L_rho_func, nDIM) > 1.0E-8);

        for(i=0; i<nDIM; i++)
        {
            mol->dyn_rho[i * timeDIM + t_index] = mol->rho[i * nDIM + i];
        }

        mol->dyn_rho[nDIM * timeDIM + t_index] = 0.0*I;

        for(i=0; i<nDIM; i++)
        {
            for(j=0; j<nDIM; j++)
            {
                mol->dyn_rho[nDIM * timeDIM + t_index] += mol->rho[i * nDIM + j]*mol->rho[j * nDIM + i];
            }
        }


        copy_mat(mol->rho, L_rho_func, nDIM);

    }

    free(L_rho_func);
}


double calculateJ(molecule* molA, molecule* molB, int nDIM)
{
    double molA_excited_pop = creal(molA->rho[2*nDIM + 2]);
//    double molA_excited_pop = creal(molA->rho[1*nDIM + 1]) + creal(molA->rho[2*nDIM + 2]);
    double molB_excited_pop = creal(molB->rho[2*nDIM + 2]);
//    double molB_excited_pop = creal(molB->rho[1*nDIM + 1]) + creal(molB->rho[2*nDIM + 2]);

    return molA_excited_pop - molB_excited_pop;
}


double nloptJ(unsigned N, const double *opt_params, double *grad_J, void *nloptJ_params)
{

    mol_system** Ensemble = (mol_system**)nloptJ_params;

    parameters* params = (*Ensemble)->params;
    molecule* moleculeA = (*Ensemble)->moleculeA;
    molecule* moleculeB = (*Ensemble)->moleculeB;
    double J;

//    params->A_R = opt_params[0];
    params->A_EE = opt_params[0];
//    params->width_R = -params->time[0] / opt_params[2];
    params->width_EE = -params->time[0] / opt_params[1];
//    params->w_v = opt_params[4];
    params->w_EE = opt_params[2];
//    params->w_R = opt_params[6];

    CalculateField(params->field_out, params);
    Propagate(moleculeA, params);
    Propagate(moleculeB, params);

    int nDIM = params->nDIM;

    J = calculateJ(moleculeA, moleculeB, nDIM);
    printf("%g %g %g %g \n", params->A_EE, -params->time[0] / params->width_EE, params->w_EE * 27.211385, J);
    return J;
}


cmplx* RamanControlFunction(molecule* molA, molecule* molB, parameters* func_params)
//------------------------------------------------------------//
//    GETTING rho(T) FROM rho(0) USING PROPAGATE FUNCTION     //
//------------------------------------------------------------//
{
    mol_system* Ensemble;
    Ensemble->moleculeA = molA;
    Ensemble->moleculeB = molB;
    Ensemble->params = func_params;

    nlopt_opt opt;

    double *lower_bounds = func_params->lower_bounds;
    double *upper_bounds = func_params->upper_bounds;

    opt = nlopt_create(NLOPT_LN_COBYLA, Ensemble->params->nDIM);
//    opt = nlopt_create(NLOPT_GN_DIRECT_L, Ensemble->params->nDIM);
    nlopt_set_lower_bounds(opt, func_params->lower_bounds);
    nlopt_set_upper_bounds(opt, func_params->upper_bounds);
    nlopt_set_max_objective(opt, nloptJ, (void*)&Ensemble);
    nlopt_set_xtol_rel(opt, 1.E-6);
    nlopt_set_maxeval(opt, func_params->MAX_EVAL);

    double *x =  func_params->guess;
    double maxf;

    if (nlopt_optimize(opt, x, &maxf) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found maximum at f(%g, %g, %g) = %0.10g\n", x[0], x[1], x[2], maxf);
    }

    nlopt_destroy(opt);
}
