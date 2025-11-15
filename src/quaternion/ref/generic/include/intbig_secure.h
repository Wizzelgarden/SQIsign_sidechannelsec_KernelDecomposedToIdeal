#include "intbig.h"  // damit ibz_t bekannt ist
#include "quaternion.h" // für ibz_mat_2x2_t

#define DESIRED_BITS 512
#define DESIRED_BITS2 1048
#define LIMBS ((DESIRED_BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)
#define LIMBS2 ((DESIRED_BITS2 + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

void ibz_init_secure(ibz_t *x);
void ibz_init_secure2(ibz_t *x);
void ibz_finalize_secure(ibz_t *x);
static inline void ibz_reduce_mod2k_secure(mp_limb_t *rp, const ibz_t *mod);
void ibz_set_secure(ibz_t *x, long v);
void ibz_copy_secure(ibz_t *target, const ibz_t *value);
void ibz_add_secure(ibz_t *sum, const ibz_t *a, const ibz_t *b);
void ibz_add_secure_modpow2(ibz_t *sum, const ibz_t *a, const ibz_t *b, const ibz_t *mod);
void ibz_mul_secure(ibz_t *prod, const ibz_t *a, const ibz_t *b);
void ibz_mul_secure_modpow2(ibz_t *prod, const ibz_t *a, const ibz_t *b, const ibz_t *mod);
void ibz_sub_secure(ibz_t *res, const ibz_t *a, const ibz_t *b);
void ibz_sub_secure_modpow2(ibz_t *res, const ibz_t *a, const ibz_t *b, const ibz_t *mod);
void ibz_div2_secure(ibz_t *res, const ibz_t *x);
void ibz_neg_secure_modpow2(ibz_t *res, const ibz_t *x, const ibz_t *mod);
int ibz_mat_2x2_inv_mod_secure(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *mod_pow2);


int ibz_invmod2k_secure(ibz_t *inv, const ibz_t *a, const ibz_t *mod);



void ibz_vec_2_init_secure(ibz_vec_2_t *x);
void ibz_vec_2_finalize_secure(ibz_vec_2_t *x);


void ibz_vec_4_init_secure(ibz_vec_4_t *x);
void ibz_vec_4_init_secure2(ibz_vec_4_t *x);
void ibz_vec_4_finalize_secure(ibz_vec_4_t *x);



// Matrizen

void ibz_mat_2x2_init_secure(ibz_mat_2x2_t *mat);
void ibz_mat_2x2_finalize_secure(ibz_mat_2x2_t *mat);
void ibz_mat_2x2_copy_secure(ibz_mat_2x2_t *dest, const ibz_mat_2x2_t *src);
void ibz_mat_2x2_eval_secure(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec, const ibz_t *mod);

void ibz_mat_hnf_mod_modpow2();

