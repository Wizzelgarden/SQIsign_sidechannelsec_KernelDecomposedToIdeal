#include "intbig_secure.h"


// initialize with always same amount of space allocated
void ibz_init_secure(ibz_t *x) {
    mpz_init2(*x, LIMBS * GMP_NUMB_BITS);
}

void ibz_init_secure2(ibz_t *x) {
    mpz_init2(*x, LIMBS2 * GMP_NUMB_BITS);
}




void ibz_finalize_secure(ibz_t *x) {
    if (!x || !(*x)) return;

    mp_limb_t *xp = (*x)->_mp_d;

    size_t alloc_limbs = (size_t)(*x)->_mp_alloc;
    for (size_t i = 0; i < alloc_limbs; i++) {
        ((volatile mp_limb_t*)xp)[i] = 0;
    }

    mpz_clear(*x);
}



static inline void ibz_reduce_mod2k_secure(mp_limb_t *rp, const ibz_t *mod)
{
    size_t mod_bits = mpz_sizeinbase(*mod, 2)-1;

    size_t limb_index = mod_bits / GMP_NUMB_BITS;  
    size_t top_bits   = mod_bits % GMP_NUMB_BITS;  

    // 1. alles über dem höchsten benötigten Limb nullen
    for (size_t i = limb_index + (top_bits != 0); i < LIMBS; ++i)
        rp[i] = 0;

    // 2. Falls genau auf Limbgrenze -> nichts maskieren
    if (top_bits == 0) {
        // volle Limbs bleiben unverändert
        // rp[limb_index] bleibt ebenfalls voll erhalten,
        // weil mod_bits genau ein Vielfaches von limbsize ist
        return;
    }

    // 3. Im obersten Limb sauber maskieren:
    // exakt die unteren top_bits behalten
    mp_limb_t mask = ((mp_limb_t)1 << top_bits) - 1;

    rp[limb_index] &= mask;
}




void ibz_set_secure(ibz_t *x, long v) {
    mp_limb_t *dst = (*x)->_mp_d;

    mp_limb_t low = (mp_limb_t)v;
    mp_limb_t sign_mask = (mp_limb_t)0 - (mp_limb_t)(v < 0); //fills other limbs depending from the sign

    dst[0] = low;

    size_t alloc_limbs = (size_t)(*x)->_mp_alloc;
    for (size_t i = 1; i < alloc_limbs; ++i) {
        dst[i] = sign_mask;
    }

    (*x)->_mp_size = LIMBS;
}



void ibz_copy_secure(ibz_t *target, const ibz_t *value) {
    
    mp_limb_t *dst = (*target)->_mp_d;
    const mp_limb_t *src = (*value)->_mp_d;
    
    size_t src_limbs = (size_t)((*value)->_mp_size);
    size_t dst_limbs = (size_t)(*target)->_mp_alloc;

    for (size_t i = 0; i < dst_limbs; i++) {
        mp_limb_t limb = (i < src_limbs) ? src[i] : 0;
        dst[i] = limb;
    }

    (*target)->_mp_size = (*value)->_mp_size;
}

void ibz_add_secure(ibz_t *sum, const ibz_t *a, const ibz_t *b) {
    mp_limb_t rp[LIMBS2];

    (void) mpn_add_n(rp, (*a)->_mp_d, (*b)->_mp_d, LIMBS2); //mpn_add_n is ct

    mpz_import(*sum, LIMBS2, -1, sizeof(mp_limb_t), 0, 0, rp);

    volatile mp_limb_t *v = rp;
    for (size_t i = 0; i < LIMBS2; ++i) v[i] = 0;
}

void ibz_add_secure_modpow2(ibz_t *sum, const ibz_t *a, const ibz_t *b, const ibz_t *mod) {
    mp_limb_t rp[LIMBS];

    (void) mpn_add_n(rp, (*a)->_mp_d, (*b)->_mp_d, LIMBS); //mpn_add_n is ct

    ibz_reduce_mod2k_secure(rp, mod);

    mpz_import(*sum, LIMBS, -1, sizeof(mp_limb_t), 0, 0, rp);

    volatile mp_limb_t *v = rp;
    for (size_t i = 0; i < LIMBS; ++i) v[i] = 0;
}


void ibz_mul_secure(ibz_t *prod, const ibz_t *a, const ibz_t *b) {
    mp_limb_t tmp[2 * LIMBS2];
    mp_limb_t rp[LIMBS2];

    (void) mpn_mul(tmp, (*a)->_mp_d, LIMBS2, (*b)->_mp_d, LIMBS2); //mpn_mul is ct

    for (size_t i = 0; i < LIMBS; ++i) {
        rp[i] = tmp[i];
    }

    mpz_import(*prod, LIMBS2, -1, sizeof(mp_limb_t), 0, 0, rp);

    volatile mp_limb_t *wipe_tmp = tmp;
    for (size_t i = 0; i < 2 * LIMBS2; ++i) wipe_tmp[i] = 0;

    volatile mp_limb_t *wipe_rp = rp;
    for (size_t i = 0; i < LIMBS2; ++i) wipe_rp[i] = 0;
}


void ibz_mul_secure_modpow2(ibz_t *prod, const ibz_t *a, const ibz_t *b, const ibz_t *mod) {
    mp_limb_t tmp[2 * LIMBS];
    mp_limb_t rp[LIMBS];

    (void) mpn_mul(tmp, (*a)->_mp_d, LIMBS, (*b)->_mp_d, LIMBS); //mpn_mul is ct

    for (size_t i = 0; i < LIMBS; ++i) {
        rp[i] = tmp[i];
    }

    ibz_reduce_mod2k_secure(rp, mod);

    mpz_import(*prod, LIMBS, -1, sizeof(mp_limb_t), 0, 0, rp);

    volatile mp_limb_t *wipe_tmp = tmp;
    for (size_t i = 0; i < 2 * LIMBS; ++i) wipe_tmp[i] = 0;

    volatile mp_limb_t *wipe_rp = rp;
    for (size_t i = 0; i < LIMBS; ++i) wipe_rp[i] = 0;
}

#include <gmp.h>
#include <stdint.h>

/* res := arithmetische Division von x durch 2 (Right shift, sign-extended)
 * - arbeitet über die tatsächliche Ziel-Alloc-Größe (r->_mp_alloc)
 * - liest fehlende Quell-Limbs als sign-extended (sign_mask)
 * - keine datenabhängigen Branches auf geheimen Daten
 * - volatile Writes um Zeroing/Stores gegen Optimizer zu schützen
 */
void ibz_div2_secure(ibz_t *res, const ibz_t *x) {
    mpz_ptr r = *res;
    mpz_srcptr a = *x;

    size_t r_alloc = (size_t) r->_mp_alloc;    
    size_t a_alloc = (size_t) a->_mp_alloc;                    

    volatile mp_limb_t *rd = (volatile mp_limb_t*) r->_mp_d;
    const mp_limb_t *ad = a->_mp_d;

    mp_limb_t top = (a_alloc > 0) ? ad[a_alloc - 1] : 0;
    mp_limb_t signbit = (top >> (GMP_NUMB_BITS - 1)) & 1u;
    mp_limb_t sign_mask = (mp_limb_t)0 - signbit; 

    mp_limb_t carry = signbit; 
    for (size_t ii = 0; ii < r_alloc; ++ii) {
        size_t i = r_alloc - 1 - ii;


        mp_limb_t limb = (i < a_alloc) ? ad[i] : sign_mask;

        mp_limb_t new_carry = limb & 1u; 
        rd[i] = (limb >> 1) | (carry << (GMP_NUMB_BITS - 1));
        carry = new_carry;
    }

    r->_mp_size = signbit ? -(mp_size_t)r_alloc : (mp_size_t)r_alloc;
}



void ibz_sub_secure(ibz_t *res, const ibz_t *a, const ibz_t *b) {
    mp_limb_t rp[LIMBS2];

    (void) mpn_sub_n(rp, (*a)->_mp_d, (*b)->_mp_d, LIMBS2);

    mpz_import(*res, LIMBS, -1, sizeof(mp_limb_t), 0, 0, rp);
 
    volatile mp_limb_t *wipe = rp;
    for (size_t i = 0; i < LIMBS2; ++i) wipe[i] = 0;
}


void ibz_sub_secure_modpow2(ibz_t *res, const ibz_t *a, const ibz_t *b, const ibz_t *mod) {
    mp_limb_t rp[LIMBS];

    (void) mpn_sub_n(rp, (*a)->_mp_d, (*b)->_mp_d, LIMBS);

    ibz_reduce_mod2k_secure(rp, mod);

    mpz_import(*res, LIMBS, -1, sizeof(mp_limb_t), 0, 0, rp);
 
    volatile mp_limb_t *wipe = rp;
    for (size_t i = 0; i < LIMBS; ++i) wipe[i] = 0;
}



// negate x modulo 2^k in ct: -x = 2^x - x mod 2^x -> inverts all bits and add one leading 1
void ibz_neg_secure_modpow2(ibz_t *res, const ibz_t *x, const ibz_t *mod)
{
    mp_limb_t *rp = (*res)->_mp_d;
    const mp_limb_t *xp = (*x)->_mp_d;

    size_t mod_bits = mpz_sizeinbase(*mod, 2);
    size_t full_limbs = mod_bits / GMP_NUMB_BITS;
    size_t top_bits = mod_bits % GMP_NUMB_BITS;

    mp_limb_t borrow = 1; 

    for (size_t i = 0; i < full_limbs; i++) {
        mp_limb_t xi = xp[i];
        mp_limb_t tmp = (~xi) + borrow;
        rp[i] = tmp;
        borrow = (tmp < borrow);
    }

    if (top_bits > 0) {
        mp_limb_t xi = xp[full_limbs];
        mp_limb_t tmp = (~xi) + borrow;

        mp_limb_t mask = ((mp_limb_t)1 << top_bits) - 1;
        rp[full_limbs] = tmp & mask;
    }

    for (size_t i = full_limbs + (top_bits>0); i < LIMBS; i++) {
        rp[i] = 0;
    }

    (*res)->_mp_size = full_limbs + (top_bits>0);
}









// 2d
void ibz_vec_2_init_secure(ibz_vec_2_t *x){
    for (int i=0; i<2; i++){
        ibz_init_secure(&(*x)[i]);
    }
}

void ibz_vec_2_finalize_secure(ibz_vec_2_t *x){
    for (int i=0; i<2; i++){
        ibz_finalize_secure(&(*x)[i]);
    }
}



// 4d
void ibz_vec_4_init_secure(ibz_vec_4_t *x){
    for (int i=0; i<4; i++){
        ibz_init_secure(&(*x)[i]);
    }
}

void ibz_vec_4_init_secure2(ibz_vec_4_t *x){
    for (int i=0; i<4; i++){
        ibz_init_secure2(&(*x)[i]);
    }
}

void ibz_vec_4_finalize_secure(ibz_vec_4_t *x){
    for (int i=0; i<4; i++){
        ibz_finalize_secure(&(*x)[i]);
    }
}








// 2x2 matrices
void ibz_mat_2x2_init_secure(ibz_mat_2x2_t *mat) {
    for (int i=0;i<2;i++)
        for (int j=0;j<2;j++)
            ibz_init_secure(&(*mat)[i][j]);
}


void ibz_mat_2x2_finalize_secure(ibz_mat_2x2_t *mat) {
    for (int i=0;i<2;i++)
        for (int j=0;j<2;j++)
            ibz_finalize_secure(&(*mat)[i][j]);
}

void
ibz_mat_2x2_eval_secure(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec, const ibz_t *mod)
{
    ibz_t prod;
    ibz_vec_2_t matvec;
    ibz_init_secure(&prod);
    ibz_vec_2_init_secure(&matvec);
    ibz_mul_secure_modpow2(&prod, &((*mat)[0][0]), &((*vec)[0]), mod);
    ibz_copy_secure(&(matvec[0]), &prod);
    ibz_mul_secure_modpow2(&prod, &((*mat)[0][1]), &((*vec)[1]), mod);
    ibz_add_secure_modpow2(&(matvec[0]), &(matvec[0]), &prod, mod);
    ibz_mul_secure_modpow2(&prod, &((*mat)[1][0]), &((*vec)[0]), mod);
    ibz_copy_secure(&(matvec[1]), &prod);
    ibz_mul_secure_modpow2(&prod, &((*mat)[1][1]), &((*vec)[1]), mod);
    ibz_add_secure_modpow2(&(matvec[1]), &(matvec[1]), &prod, mod);
    ibz_copy_secure(&((*res)[0]), &(matvec[0]));
    ibz_copy_secure(&((*res)[1]), &(matvec[1]));
    ibz_finalize(&prod);
    ibz_vec_2_finalize(&matvec);
}


// ------------------------
// Newton-Raphson für mod 2^k
// see arxiv.org/pdf/1209.6626 for reference of this implementation
// ------------------------

int ibz_invmod2k_secure(ibz_t *U, const ibz_t *a, const ibz_t *mod)
{
    // temporaries
    ibz_t b;         // holds (a - 1) mod 2^m, then its successive squares
    ibz_t tmp;       // temporary for squaring
    ibz_t tmp2;      // temporary for (b^2 + 1)
    ibz_t secure_one;// const one with required 4 limbs init
    ibz_t secure_two;// same for two
    ibz_init_secure(&b);
    ibz_init_secure(&tmp);
    ibz_init_secure(&tmp2);

    ibz_init_secure(&secure_one);
    ibz_init_secure(&secure_two);

    // check invertibility: a must be odd
    mp_limb_t a0 = ((*a)->_mp_d)[0];
    unsigned int invertible_flag = (unsigned int)(a0 & 1UL);
    if (!invertible_flag) {
        // define U = 0 if not invertible
        ibz_set_secure(U, 0);
        ibz_finalize_secure(&b);
        ibz_finalize_secure(&tmp);
        ibz_finalize_secure(&tmp2);
        return 0;
    }
        ibz_set_secure(&secure_two, 2);
        ibz_set_secure(&secure_one, 1);
        
    //  U = 2 - a  (mod 2^m)
    ibz_sub_secure_modpow2(U, &secure_two, a, mod); 

    // b = a - 1  (mod 2^m)
    ibz_sub_secure_modpow2(&b, a, &secure_one, mod); 

    //    we need to cover log(mod) many bits
    size_t m_bits = mpz_sizeinbase(*mod, 2);
    //printf("m_bits = %zu\n", m_bits);
    size_t i = 1;

    while (i < m_bits) {
        // tmp = b^2 mod 2^m
        ibz_mul_secure_modpow2(&tmp, &b, &b, mod);

        // tmp2 = tmp + 1  (mod 2^m)
        ibz_add_secure_modpow2(&tmp2, &tmp, &secure_one, mod);

        // U = U * tmp2  (mod 2^m)
        ibz_mul_secure_modpow2(U, U, &tmp2, mod);

        // set b := tmp (the squared value) for next iteration
        ibz_copy_secure(&b, &tmp);

        // double i
        i <<= 1;
    }

    ibz_finalize_secure(&b);
    ibz_finalize_secure(&tmp);
    ibz_finalize_secure(&tmp2);

    return 1;
}


int ibz_mat_2x2_inv_mod_secure(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *mod)
{
    ibz_t det, temp;

    ibz_init_secure(&det);
    ibz_init_secure(&temp);

    // det = a*d - b*c  mod 2^k
    ibz_mul_secure_modpow2(&det, &((*mat)[0][0]), &((*mat)[1][1]), mod);
    ibz_mul_secure_modpow2(&temp, &((*mat)[0][1]), &((*mat)[1][0]), mod);
    ibz_sub_secure_modpow2(&det, &det, &temp, mod);

    // det^-1  mod 2^k (ct)
    int invertible = ibz_invmod2k_secure(&det, &det, mod);

    // adjugate matrix
    ibz_copy_secure(&((*inv)[0][0]), &((*mat)[1][1]));
    ibz_copy_secure(&((*inv)[1][1]), &((*mat)[0][0]));
    ibz_neg_secure_modpow2(&((*inv)[0][1]), &((*mat)[0][1]), mod);
    ibz_neg_secure_modpow2(&((*inv)[1][0]), &((*mat)[1][0]), mod);

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            ibz_mul_secure_modpow2(&((*inv)[i][j]), &((*inv)[i][j]), &det, mod);

    ibz_finalize_secure(&det);
    ibz_finalize_secure(&temp);

    return invertible;
}


void ibz_mat_hnf_mod_modpow2(){
    
}