#include <iostream>
#include <sys/time.h>
#include <gmpxx.h>

using namespace std;

typedef mpz_class ZZZ;

// L function
ZZZ L_function(const ZZZ x, const ZZZ n) 
{ 
    return (x -1) / n;
}

void keyGeneration(ZZZ &p, ZZZ &q, ZZZ &p_phi, ZZZ &q_phi, ZZZ &n, ZZZ &n2, ZZZ &g, const unsigned long securityPara) // ZZZ &lambda, ZZZ &mu,
{
    //generate p,q
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    mpz_urandomb(p.get_mpz_t(), state, securityPara/2);
    mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());

    mpz_urandomb(q.get_mpz_t(), state, securityPara/2);
    mpz_nextprime(q.get_mpz_t(), q.get_mpz_t()); 

    //compute n,n^2, g=n+1
    n = p * q;
    n2 = n * n;
    p_phi = p - 1;
    q_phi = q - 1;
    g = n + 1;

    //compute lambda,mu
    // mpz_lcm(lambda.get_mpz_t(), p_phi.get_mpz_t(), q_phi.get_mpz_t());
    // mpz_invert(mu.get_mpz_t(), lambda.get_mpz_t(), n.get_mpz_t());
    
}

ZZZ encrypt(const ZZZ m, const ZZZ n, const ZZZ n2)
{
    //generate random r in [0, n]
    ZZZ r;
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    mpz_urandomm(r.get_mpz_t(), state, n.get_mpz_t());


    ZZZ t1, t2, t3, t4, c;
    //compute m*n+1 mod n^2
    t1 = m * n + 1;
    mpz_mod(t2.get_mpz_t(), t1.get_mpz_t(), n2.get_mpz_t());

    //compute r^n mod n^2
    mpz_powm(t3.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t(), n2.get_mpz_t());

    //compute (m*n+1)*r^n mod n^2
    t4 = t2 * t3;
    mpz_mod(c.get_mpz_t(), t4.get_mpz_t(), n2.get_mpz_t());

    return c;
}

ZZZ decrypt(const ZZZ p, const ZZZ q, const ZZZ p_phi, const ZZZ q_phi, const ZZZ n, const ZZZ g, const ZZZ c)
{
    ZZZ x_p, m_p, h_p, p2, p_inv;
    ZZZ x_q, m_q, h_q, q2;
    p2 = p * p;
    q2 = q * q;
    mpz_invert(p_inv.get_mpz_t(), p.get_mpz_t(), q.get_mpz_t());
    
    mpz_powm(m_p.get_mpz_t(), c.get_mpz_t(), p_phi.get_mpz_t(), p2.get_mpz_t());
    m_p = L_function(m_p, p);

    mpz_powm(h_p.get_mpz_t(), g.get_mpz_t(), p_phi.get_mpz_t(), p2.get_mpz_t());
    h_p = L_function(h_p, p);
    mpz_invert(h_p.get_mpz_t(), h_p.get_mpz_t(), p.get_mpz_t());
    mpz_mod(h_p.get_mpz_t(), h_p.get_mpz_t(), p.get_mpz_t());
    m_p = m_p * h_p;
    mpz_mod(x_p.get_mpz_t(), m_p.get_mpz_t(), p.get_mpz_t());

    mpz_powm(m_q.get_mpz_t(), c.get_mpz_t(), q_phi.get_mpz_t(), q2.get_mpz_t());
    m_q = L_function(m_q, q);

    mpz_powm(h_q.get_mpz_t(), g.get_mpz_t(), q_phi.get_mpz_t(), q2.get_mpz_t());
    h_q = L_function(h_q, q);
    mpz_invert(h_q.get_mpz_t(), h_q.get_mpz_t(), q.get_mpz_t());
    mpz_mod(h_q.get_mpz_t(), h_q.get_mpz_t(), q.get_mpz_t());
    m_q = m_q * h_q;
    mpz_mod(x_q.get_mpz_t(), m_q.get_mpz_t(), q.get_mpz_t());

    ZZZ t, m;
    t = x_q - x_p;
    t = t * p_inv;
    mpz_mod(t.get_mpz_t(), t.get_mpz_t(), q.get_mpz_t());
    t = t * p;
    m = t + x_p;

    return m;
}

ZZZ homoAdd(const ZZZ c1, const ZZZ c2, const ZZZ n2)
{
    ZZZ sum;
    sum = c1 * c2;
    mpz_mod(sum.get_mpz_t(), sum.get_mpz_t(), n2.get_mpz_t());

    return sum;
}

//test
int main()
{
    unsigned long securityPara = 1024;
    ZZZ p, q, p_phi, q_phi, n, n2, g; //, lambda, mu
    ZZZ m1, m2, c1, c2, d1, d2;
    keyGeneration(p, q, p_phi, q_phi, n, n2, g, securityPara);
    cout << "Please input plaintext m1 : ";
    cin >> m1;
    cout << "Please input plaintext m2 : ";
    cin >> m2;

    c1 = encrypt(m1, n, n2);
    cout << "ciphertext c1 : " << c1 << endl;
    c2 = encrypt(m2, n, n2);
    cout << "ciphertext c2 : " << c2 << endl;

    d1 = decrypt(p, q, p_phi, q_phi, n, g, c1);
    cout << "decrypt c1 : " << d1 << endl;
    d2 = decrypt(p, q, p_phi, q_phi, n, g, c2);
    cout << "decrypt c2 : " << d2 << endl;

    ZZZ c_sum = homoAdd(c1, c2, n2);
    cout << "ciphertext c1 + c2 : " << c_sum << endl;
    ZZZ m_sum = decrypt(p, q, p_phi, q_phi, n, g, c_sum);
    cout << "decrypt c1 + c2 : " << m_sum << endl;

    return 0;
}
