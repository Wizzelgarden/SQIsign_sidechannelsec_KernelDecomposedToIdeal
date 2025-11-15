#include <quaternion.h>
#include <ec.h>
#include <mp.h>
#include <endomorphism_action.h>
#include <id2iso.h>
#include <inttypes.h>
#include <locale.h>
#include <biextension.h>
#include <torsion_constants.h>

#include "../../../quaternion/ref/generic/internal_quaternion_headers/internal.h"
#include "intbig_secure.h"


// Scalar multiplication [x]P + [y]Q where x and y are stored
// inside an ibz_vec_2_t [x, y] and P, Q \in E[2^f]
void
ec_biscalar_mul_ibz_vec(ec_point_t *res,
                        const ibz_vec_2_t *scalar_vec,
                        const int f,
                        const ec_basis_t *PQ,
                        const ec_curve_t *curve)
{
    digit_t scalars[2][NWORDS_ORDER];
    ibz_to_digit_array(scalars[0], &(*scalar_vec)[0]);
    ibz_to_digit_array(scalars[1], &(*scalar_vec)[1]);
    ec_biscalar_mul(res, scalars[0], scalars[1], f, PQ, curve);
}

// Given an ideal, computes the scalars s0, s1 which determine the kernel generator
// of the equivalent isogeny
void
id2iso_ideal_to_kernel_dlogs_even(ibz_vec_2_t *vec, const quat_left_ideal_t *lideal)
{
    ibz_t tmp;
    ibz_init(&tmp);

    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    // construct the matrix of the dual of alpha on the 2^f-torsion
    {
        quat_alg_elem_t alpha;
        quat_alg_elem_init(&alpha);

        int lideal_generator_ok UNUSED = quat_lideal_generator(&alpha, lideal, &QUATALG_PINFTY);
        assert(lideal_generator_ok);
        quat_alg_conj(&alpha, &alpha);

        ibz_vec_4_t coeffs;
        ibz_vec_4_init(&coeffs);
        quat_change_to_O0_basis(&coeffs, &alpha);

        for (unsigned i = 0; i < 2; ++i) {
            ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
            for (unsigned j = 0; j < 2; ++j) {
                ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
            }
        }

        ibz_vec_4_finalize(&coeffs);
        quat_alg_elem_finalize(&alpha);
    }

    // find the kernel of alpha modulo the norm of the ideal
    {
        const ibz_t *const norm = &lideal->norm;

        ibz_mod(&(*vec)[0], &mat[0][0], norm);
        ibz_mod(&(*vec)[1], &mat[1][0], norm);
        ibz_gcd(&tmp, &(*vec)[0], &(*vec)[1]);
        if (ibz_is_even(&tmp)) {
            ibz_mod(&(*vec)[0], &mat[0][1], norm);
            ibz_mod(&(*vec)[1], &mat[1][1], norm);
        }
#ifndef NDEBUG
        ibz_gcd(&tmp, &(*vec)[0], norm);
        ibz_gcd(&tmp, &(*vec)[1], &tmp);
        assert(!ibz_cmp(&tmp, &ibz_const_one));
#endif
    }

    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&tmp);
}

// helper function to apply a matrix to a basis of E[2^f]
// works in place
int
matrix_application_even_basis(ec_basis_t *bas, const ec_curve_t *E, ibz_mat_2x2_t *mat, int f)
{
    digit_t scalars[2][NWORDS_ORDER] = { 0 };
    int ret;

    ibz_t tmp, pow_two;
    ibz_init(&tmp);
    ibz_init(&pow_two);
    ibz_pow(&pow_two, &ibz_const_two, f);

    ec_basis_t tmp_bas;
    copy_basis(&tmp_bas, bas);

    // reduction mod 2f
    ibz_mod(&(*mat)[0][0], &(*mat)[0][0], &pow_two);
    ibz_mod(&(*mat)[0][1], &(*mat)[0][1], &pow_two);
    ibz_mod(&(*mat)[1][0], &(*mat)[1][0], &pow_two);
    ibz_mod(&(*mat)[1][1], &(*mat)[1][1], &pow_two);

    // For a matrix [[a, c], [b, d]] we compute:
    //
    // first basis element R = [a]P + [b]Q
    ibz_to_digit_array(scalars[0], &(*mat)[0][0]);
    ibz_to_digit_array(scalars[1], &(*mat)[1][0]);
    ec_biscalar_mul(&bas->P, scalars[0], scalars[1], f, &tmp_bas, E);

    // second basis element S = [c]P + [d]Q
    ibz_to_digit_array(scalars[0], &(*mat)[0][1]);
    ibz_to_digit_array(scalars[1], &(*mat)[1][1]);
    ec_biscalar_mul(&bas->Q, scalars[0], scalars[1], f, &tmp_bas, E);

    // Their difference R - S = [a - c]P + [b - d]Q
    ibz_sub(&tmp, &(*mat)[0][0], &(*mat)[0][1]);
    ibz_mod(&tmp, &tmp, &pow_two);
    ibz_to_digit_array(scalars[0], &tmp);
    ibz_sub(&tmp, &(*mat)[1][0], &(*mat)[1][1]);
    ibz_mod(&tmp, &tmp, &pow_two);
    ibz_to_digit_array(scalars[1], &tmp);
    ret = ec_biscalar_mul(&bas->PmQ, scalars[0], scalars[1], f, &tmp_bas, E);

    ibz_finalize(&tmp);
    ibz_finalize(&pow_two);

    return ret;
}

// helper function to apply some endomorphism of E0 on the precomputed basis of E[2^f]
// works in place
void
endomorphism_application_even_basis(ec_basis_t *bas,
                                    const int index_alternate_curve,
                                    const ec_curve_t *E,
                                    const quat_alg_elem_t *theta,
                                    int f)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_vec_4_t coeffs;
    ibz_vec_4_init(&coeffs);
    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    ibz_t content;
    ibz_init(&content);

    // decomposing theta on the basis
    quat_alg_make_primitive(&coeffs, &content, theta, &EXTREMAL_ORDERS[index_alternate_curve].order);
    assert(ibz_is_odd(&content));

    ibz_set(&mat[0][0], 0);
    ibz_set(&mat[0][1], 0);
    ibz_set(&mat[1][0], 0);
    ibz_set(&mat[1][1], 0);

    // computing the matrix

    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &CURVES_WITH_ENDOMORPHISMS[index_alternate_curve].action_gen2[i][j], &coeffs[1]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &CURVES_WITH_ENDOMORPHISMS[index_alternate_curve].action_gen3[i][j], &coeffs[2]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &CURVES_WITH_ENDOMORPHISMS[index_alternate_curve].action_gen4[i][j], &coeffs[3]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&mat[i][j], &mat[i][j], &content);
        }
    }

    // and now we apply it
    matrix_application_even_basis(bas, E, &mat, f);

    ibz_vec_4_finalize(&coeffs);
    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&content);

    ibz_finalize(&tmp);
}






// compute the ideal whose kernel is generated by vec2[0]*BO[0] + vec2[1]*B0[1] where B0 is the
// canonical basis of E0
void
id2iso_kernel_dlogs_to_ideal_even(quat_left_ideal_t *lideal, const ibz_vec_2_t *vec2, int f)
{

    // algorithm: apply endomorphisms 1 and j+(1+k)/2 to the kernel point,
    // the result should form a basis of the respective torsion subgroup.
    // then apply i to the kernel point and decompose over said basis.
    // hence we have an equation a*P + b*[j+(1+k)/2]P == [i]P, which will
    // easily reveal an endomorphism that kills P.

    struct timespec start, end;

    ibz_t two_pow;
    ibz_init_secure(&two_pow);
    // even though this is public, we require a fixed number of allocated limbs for ct algorithms

    ibz_vec_2_t vec;
    ibz_vec_2_init_secure(&vec);

if (f == TORSION_EVEN_POWER) {
        ibz_copy_secure(&two_pow, &TORSION_PLUS_2POWER); 
    } else {
        ibz_pow(&two_pow, &ibz_const_two, f); 
    }
      
    clock_gettime(CLOCK_MONOTONIC, &start);


    // In this block we compute vector vec = [a,b] as in the specification document

    {
        ibz_mat_2x2_t mat;
        ibz_mat_2x2_init_secure(&mat);

        ibz_mat_2x2_t inv;
        ibz_mat_2x2_init_secure(&inv);

        ibz_copy_secure(&mat[0][0], &(*vec2)[0]);
        ibz_copy_secure(&mat[1][0], &(*vec2)[1]);

        ibz_mat_2x2_eval_secure(&vec, &ACTION_J, vec2, &two_pow);
        ibz_copy_secure(&mat[0][1], &vec[0]);
        ibz_copy_secure(&mat[1][1], &vec[1]);

        ibz_mat_2x2_eval_secure(&vec, &ACTION_GEN4, vec2, &two_pow);
        ibz_add_secure_modpow2(&mat[0][1], &mat[0][1], &vec[0], &two_pow);
        ibz_add_secure_modpow2(&mat[1][1], &mat[1][1], &vec[1], &two_pow);

        {
            int inv_ok UNUSED = ibz_mat_2x2_inv_mod_secure(&inv, &mat, &two_pow);
            assert(inv_ok);
        }

        ibz_mat_2x2_eval_secure(&vec, &ACTION_I, vec2, &two_pow);
        ibz_mat_2x2_eval_secure(&vec, &inv, &vec, &two_pow);


        ibz_mat_2x2_finalize_secure(&mat);
        ibz_mat_2x2_finalize_secure(&inv);
    }


        //gmp_printf("a= %Zd\n", &vec[0]);
        //gmp_printf("b= %Zd\n", &vec[1]);

    // Next, we have to compute the corresponding ideal. In this special instance, we do not need to use the general 
    // functions because we can simplify most calculations


    /* 
    As this function is only called in context of the challenge translation, 
    one can rewrite the ideal creating for this special case
    Now the correspodning ideal is computed. However, we can greatly simplify it in our case.
    instead of building  alpha = a - i + b*(j+(1+k)/2) and computing the ideal using standard functions, we can hardcode it.
    Further, we can avoid the redudction of the denom as we know all entries are even and the denom is 4. Therefore just
    divide every entry by to when generating and set the denom also to two.
    Additionally, we can compute the norm of the ideal by hand using special form of alpha
    */ 

    ibz_t temp;
    ibz_init_secure2(&temp);

    ibz_t temp2;
    ibz_init_secure2(&temp2);

    ibz_vec_4_t generators[8];
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init_secure2(&(generators[i]));

    // secure init for fixed number of limbs allocated
    ibz_t secure_const_one;
    ibz_init_secure2(&secure_const_one);
    ibz_set_secure(&secure_const_one, 1);

    ibz_t secure_const_two;
    ibz_init_secure2(&secure_const_two);
    ibz_set_secure(&secure_const_two, 2);

/*
    In generators we copy the rows of the lattices corresponding to O*N and O*alpha, where N=2^f. Therefore, the matrix 
    corresponding to maxorder is scaled by this factor (note that the denominator is one by default and later calculations)
    are optimized such that the denominator is not explictly required. 
*/

    ibz_mul_secure(&temp, &two_pow, &ibz_const_two);
    ibz_copy_secure(&(generators[0][0]), &temp);
    ibz_copy_secure(&(generators[1][1]), &temp);
    ibz_copy_secure(&(generators[2][1]), &two_pow);
    ibz_copy_secure(&(generators[2][2]), &two_pow);
    ibz_copy_secure(&(generators[3][0]), &two_pow);
    ibz_copy_secure(&(generators[3][3]), &two_pow);


/* 
    The lattice for maxorder is given by 1/2 *
    2 0 0 1
    0 2 1 0
    0 0 1 0 
    0 0 0 1
    and in the standard basis alpha can be given by
    1/2 * (2a+b, -2, 2b, b)
    The product of both can explicitely given in terms of a,b,p, and this is how we now compute lideal->lattice.basis
    However, every entry will be a multiple of two, and to avoid later denominator reduction, we can directly ignore it
*/
 
    // All entries are bounded by 2^f*(p+3). 65 is the largest cofactor of p in all parametersets. 
    // Thus 2^(2f+10) will be a sufficient upper bound in all cases

    // 0 0 
    ibz_add_secure(&temp, &vec[0], &vec[0]);
    ibz_add_secure(&temp, &temp, &vec[1]);
    ibz_copy_secure(&generators[4][0], &temp);
    ibz_set_secure(&temp, 0);

    // 1 0
    ibz_set_secure(&generators[4][1], -2);

    // 2 0
    ibz_mul_secure(&temp, &vec[1], &secure_const_two);
    ibz_copy_secure(&generators[4][2], &temp);
    ibz_set_secure(&temp, 0);

    // 3 0   
    ibz_copy_secure(&temp, &vec[1]);
    ibz_copy_secure(&generators[4][3], &temp);
    ibz_set_secure(&temp, 0);

    // 0 1
    ibz_set_secure(&generators[5][0], 2);

    // 1 1
    ibz_add_secure(&temp, &vec[0], &vec[0]);
    ibz_add_secure(&temp, &temp, &vec[1]);
    ibz_copy_secure(&generators[5][1], &temp);
    ibz_set_secure(&temp, 0);  
    
    // 2 1
    ibz_set_secure(&temp, -1);
    ibz_mul_secure(&temp, &temp, &vec[1]);
    ibz_copy_secure(&generators[5][2], &temp);
    ibz_set_secure(&temp, 0);

    // 3 1
    ibz_add_secure(&temp, &vec[1], &vec[1]);
    ibz_copy_secure(&generators[5][3], &temp);
    ibz_set_secure(&temp, 0);

    // 0 2
    ibz_set_secure(&temp, -1);
    ibz_mul_secure(&temp, &temp, &vec[1]);
    ibz_mul_secure(&temp, &temp, &(QUATALG_PINFTY.p));
    ibz_add_secure(&temp, &temp, &secure_const_one);
    ibz_copy_secure(&generators[6][0], &temp);
    ibz_set_secure(&temp, 0);

    // 1 2
    ibz_copy_secure(&temp, &vec[1]);
    ibz_mul_secure(&temp, &temp, &(QUATALG_PINFTY.p));
    ibz_add_secure(&temp2, &vec[0], &vec[0]);
    ibz_add_secure(&temp2, &temp2, &vec[1]);
    ibz_add_secure(&temp, &temp, &temp2);
    ibz_div2_secure(&temp, &temp);
    ibz_copy_secure(&generators[6][1], &temp);
    ibz_set_secure(&temp, 0);
    ibz_set_secure(&temp2, 0);

    // 2 2
    ibz_copy_secure(&temp, &vec[0]);
    ibz_copy_secure(&generators[6][2], &temp);
    ibz_set_secure(&temp, 0);

    // 3 2
    ibz_copy_secure(&temp, &vec[1]);
    ibz_add_secure(&temp, &temp, &secure_const_one);
    ibz_copy_secure(&generators[6][3], &temp);
    ibz_set_secure(&temp, 0);

    // 0 3
    ibz_copy_secure(&temp, &vec[1]);
    ibz_mul_secure(&temp, &temp, &(QUATALG_PINFTY.p));
    ibz_add_secure(&temp2, &vec[0], &vec[0]);
    ibz_add_secure(&temp2, &temp2, &vec[1]);
    ibz_sub_secure(&temp, &temp2, &temp);
    ibz_div2_secure(&temp, &temp);
    ibz_copy_secure(&generators[7][0], &temp);
    ibz_set_secure(&temp, 0);
    ibz_set_secure(&temp2, 0);

    // 1 3
    ibz_set_secure(&temp, -1);
    ibz_mul_secure(&temp, &temp, &vec[1]);
    ibz_mul_secure(&temp, &temp, &(QUATALG_PINFTY.p));
    ibz_sub_secure(&temp, &temp, &secure_const_one);
    ibz_copy_secure(&generators[7][1], &temp);
    ibz_set_secure(&temp, 0);

    // 2 3
    ibz_copy_secure(&temp, &vec[1]);
    ibz_sub_secure(&temp, &temp, &secure_const_one);
    ibz_copy_secure(&generators[7][2], &temp);
    ibz_set_secure(&temp, 0);

    // 3 3
    ibz_add_secure(&temp, &vec[0], &vec[1]);
    ibz_copy_secure(&generators[7][3], &temp);
    ibz_set_secure(&temp, 0);
  
    // would be nice to get rid of this hnf function call

    ibz_set(&temp2, 8);
    ibz_mul(&temp, &two_pow, &two_pow);
    ibz_mul(&temp, &temp, &temp2); // determinant of hnf
    ibz_mat_4xn_hnf_mod_core(&(lideal->lattice.basis), 8, generators, &temp);
    ibz_copy_secure(&(lideal->lattice.denom), &ibz_const_two);

    /*
    gmp_printf("hnf output\n");
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            gmp_printf(" %Zd ",&(lideal->lattice.basis[i][j]));
        }
        gmp_printf("\n");
    }*/

    lideal->parent_order = &MAXORD_O0;
    ibz_copy(&lideal->norm, &two_pow);
    // both is necessary for further computations


    clock_gettime(CLOCK_MONOTONIC, &end);

    //double elapsed_sec = (end.tv_sec - start.tv_sec)
    //                   + (end.tv_nsec - start.tv_nsec) / 1e9;

    //printf("Elapsed time: %f seconds\n", elapsed_sec);

    ibz_finalize_secure(&two_pow);
    ibz_vec_2_finalize_secure(&vec);
    for (int i = 0; i < 8; i++){
        ibz_vec_4_finalize_secure(&(generators[i]));
    }

    ibz_finalize_secure(&temp);
    ibz_finalize_secure(&temp2);




}

// finds mat such that:
// (mat*v).B2 = v.B1
// where "." is the dot product, defined as (v1,v2).(P,Q) = v1*P + v2*Q
// mat encodes the coordinates of the points of B1 in the basis B2
// specifically requires B1 or B2 to be "full" w.r.t to the 2^n torsion, so that we use tate
// full = 0 assumes B2 is "full" so the easier case.
// if we want to switch the role of B2 and B1, we invert the matrix, e.g. set full = 1
static void
_change_of_basis_matrix_tate(ibz_mat_2x2_t *mat,
                             const ec_basis_t *B1,
                             const ec_basis_t *B2,
                             ec_curve_t *E,
                             int f,
                             bool invert)
{
    digit_t x1[NWORDS_ORDER] = { 0 }, x2[NWORDS_ORDER] = { 0 }, x3[NWORDS_ORDER] = { 0 }, x4[NWORDS_ORDER] = { 0 };

#ifndef NDEBUG
    int e_full = TORSION_EVEN_POWER;
    int e_diff = e_full - f;
#endif

    // Ensure the input basis has points of order 2^f
    if (invert) {
        assert(test_basis_order_twof(B1, E, e_full));
        ec_dlog_2_tate(x1, x2, x3, x4, B1, B2, E, f);
        mp_invert_matrix(x1, x2, x3, x4, f, NWORDS_ORDER);
    } else {
        assert(test_basis_order_twof(B2, E, e_full));
        ec_dlog_2_tate(x1, x2, x3, x4, B2, B1, E, f);
    }

#ifndef NDEBUG
    {
        if (invert) {
            ec_point_t test, test2;
            ec_biscalar_mul(&test, x1, x2, f, B2, E);
            ec_dbl_iter(&test2, e_diff, &B1->P, E);
            assert(ec_is_equal(&test, &test2));

            ec_biscalar_mul(&test, x3, x4, f, B2, E);
            ec_dbl_iter(&test2, e_diff, &B1->Q, E);
            assert(ec_is_equal(&test, &test2));
        } else {
            ec_point_t test;
            ec_biscalar_mul(&test, x1, x2, f, B2, E);
            ec_dbl_iter(&test, e_diff, &test, E);
            assert(ec_is_equal(&test, &(B1->P)));

            ec_biscalar_mul(&test, x3, x4, f, B2, E);
            ec_dbl_iter(&test, e_diff, &test, E);
            assert(ec_is_equal(&test, &(B1->Q)));
        }
    }
#endif

    // Copy the results into the matrix
    ibz_copy_digit_array(&((*mat)[0][0]), x1);
    ibz_copy_digit_array(&((*mat)[1][0]), x2);
    ibz_copy_digit_array(&((*mat)[0][1]), x3);
    ibz_copy_digit_array(&((*mat)[1][1]), x4);
}

void
change_of_basis_matrix_tate(ibz_mat_2x2_t *mat, const ec_basis_t *B1, const ec_basis_t *B2, ec_curve_t *E, int f)
{
    _change_of_basis_matrix_tate(mat, B1, B2, E, f, false);
}

void
change_of_basis_matrix_tate_invert(ibz_mat_2x2_t *mat, const ec_basis_t *B1, const ec_basis_t *B2, ec_curve_t *E, int f)
{
    _change_of_basis_matrix_tate(mat, B1, B2, E, f, true);
}
