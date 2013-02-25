#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

typedef unsigned char uchar;
typedef unsigned long long ull;
typedef long long ll;
typedef unsigned int uint;

void error(char *s) {
	puts(s);
	exit(1);
}

/* base 2 logarithm */
double log2(double a) {
	const double z=1.44269504088896340736; /* 1/log(2) */
	return log(a)*z;
}

#define BIGDEG 10
#define MAXDEG 10

/* input parameters */

mpz_t opt_n;    /* number to factorize */
ull opt_Ba;     /* bound for algebraic factor base */
ull opt_Br;     /* bound for rational factor base */
int opt_Bq;     /* number of quadratic characters */
int opt_deg;    /* degree of polynomial (must be odd, >=3) */
mpz_t opt_m;    /* m value for base-m algorithm */
ull opt_sievew; /* width of line sieve */
int opt_thr;    /* threshold for accepting a number in the sieve */
int opt_skip;   /* skip this amount of smallest primes in the sieve */
int opt_extra;  /* number of extra relations wanted for linear algebra */
int opt_signb;  /* -1: a-bm, 1: a+bm */

void getnextline(char *s) {
	int l;
loop:
	if(!fgets(s,1048570,stdin)) { s[0]=0; return; }
	if(s[0]=='\n' || s[0]=='\r') goto loop;
	if(s[0]==';' || s[0]=='%' || s[0]=='#') goto loop;
	l=strlen(s);
	while(l && (s[l]=='\n' || s[l]=='\r')) s[l--]=0;
}

gmp_randstate_t gmpseed;

/* get a d-digit random number */
void getmpzrandom(mpz_t r,int d) {
	static char t[1024];
	mpz_t a,b;
	int i;
	if(d>1022) error("too many digits");
	mpz_init(a); mpz_init(b);
	t[i]='1';
	for(i=1;i<=d;i++) t[i]='0';
	t[i]=0;
	mpz_set_str(a,t,10);
	t[i-1]=0;
	mpz_set_str(b,t,10);
	do mpz_urandomm(r,gmpseed,a); while(mpz_cmp(r,b)<0);
	mpz_clear(b); mpz_clear(a);
}

/* return a mod p where a is mpz and p is int */
int mpz_mod_int(mpz_t a,int p) {
	mpz_t b;
	int r;
	mpz_init(b);
	r=mpz_mod_ui(b,a,p);
	mpz_clear(b);
	return r;
}

void readoptions() {
	static char s[1048576],t[4096];
	int z,i;
	/* read n */
	mpz_init(opt_n);
	mpz_init(opt_m);
	getnextline(s);
	sscanf(s,"%4090s",t);
	if(t[0]=='c') {
		/* generate a z-digit composite number without small prime factors */
		z=strtol(t+1,NULL,10);
		do {
			getmpzrandom(opt_n,z);
			if(!mpz_mod_int(opt_n,2)) continue;
			if(mpz_probab_prime_p(opt_n,25)) continue;
			if(z>10) for(i=3;i<20000;i+=2) if(!mpz_mod_int(opt_n,i)) continue;
		} while(0);
	} else if(t[0]=='r') {
		/* generate RSA number: a z-digit number that is the product
		   of two similarly sized primes */
		z=strtol(t+1,NULL,10);
		/* TODO, pick two primes of z/2 digits and multiply */
		error("not implemented yet");
	} else {
		/* take literal number */
		mpz_set_str(opt_n,t,10);
	}
	getnextline(s); sscanf(s,"%I64d",&opt_Ba);
	getnextline(s); sscanf(s,"%I64d",&opt_Br);
	getnextline(s); sscanf(s,"%d",&opt_Bq);
	getnextline(s); sscanf(s,"%d",&opt_deg);
	getnextline(s); mpz_set_str(opt_m,s,10);
	getnextline(s); sscanf(s,"%I64d",&opt_sievew);
	getnextline(s); sscanf(s,"%d",&opt_thr);
	/* TODO support percentage for skip (that is, skip x percent of the
	   primes) */
	getnextline(s); sscanf(s,"%d",&opt_skip);
	getnextline(s); sscanf(s,"%d",&opt_extra);
	if(opt_deg>MAXDEG) error("too high degree");
	if(!(opt_deg&1)) error("degree must be odd");
	getnextline(s); sscanf(s,"%d",&opt_signb);
	if(opt_signb!=1 && opt_signb!=-1) error("wrong sign");
}

/* all polynomials have the following format:
   coefficients in f[i], f[0]=a_0, f[1]=a_1, ...,  f[i]=a_i
   size of f[] is MAXDEG+1 */

/* auxilliary routines */

/* need these since gmp doesn't support long long */
ull mpz_get_ull(mpz_t a) {
	static char s[1048576];
	ull ret;
	mpz_get_str(s,10,a);
	sscanf(s,"%I64d",&ret);
	return ret;
}

void mpz_set_ull(mpz_t b,ull a) {
  mpz_import(b, 1, 1, sizeof(a), 0, 0, &a);
}

/* return a mod p where a is mpz */
ull mpz_mod_ull(mpz_t a,ull p) {
	ull r;
	mpz_t b;
	mpz_init(b);
	mpz_set_ull(b,p);
	mpz_mod(b,a,b);
	r=mpz_get_ull(b);
	mpz_clear(b);
	return r;
}

ull gcd(ull a,ull b) {
	return b?gcd(b,a%b):a;
}

/* nfs init stage: create polynomials, determine bounds, create factor base */
/* includes many subroutines for polynomials */

/* calculate asymptotically optimal d (degree of polynomial)
   warning, doesn't work for n with more than 307 digits or so
   (n must fit in double) */
int findhighestdegree(mpz_t n) {
	double N=mpz_get_d(n);
	return pow(3*log(N)/log(log(N)),1./3);
}

/* calculate upper bound for factor base
   warning, doesn't work for n with more than 307 digits or so
   (n must fit in double) */
ull findB(mpz_t n) {
	double N=mpz_get_d(n),z=1./3;
	return exp(pow(8./9,z)*pow(log(N),z)*pow(log(log(N)),z+z));
}

/* calculate number of quadratic characters to obtain
   warning, doesn't work for n with more than 307 digits or so
   (n must fit in double) */
ull findK(mpz_t n) {
	double N=mpz_get_d(n);
	return 3*log(N)/log(2.728182818);
}

/* given n, m and d, return polynomial of degree d which is the
   base-m expansion of n. return 0 if something went wrong (degree doesn't 
   match expansion, polynomial isn't monic etc) */
/* assume that *f is allocated with (d+1) uninitialized elements.
   f[0]=a_0, f[1]=a_1, ..., f[d]=1,
   polynomial is f(x)=a_d x^d + ... + a_1 x + a_0 */
int getpolynomial(mpz_t n,mpz_t m,int d,mpz_t *f) {
	mpz_t N;
	int i,r=1;
	mpz_init_set(N,n);
	for(i=0;i<=d;i++) {
		mpz_init(f[i]);
		mpz_fdiv_qr(N,f[i],N,m);
	}
	/* error if base-m expansion of n requires degree!=d or f isn't monic */
	if(mpz_cmp_si(N,0) || mpz_cmp_si(f[d],1)) r=0;
	mpz_clear(N);
	return r;
}

void printmpzpoly(mpz_t *f,int d) {
	for(;d>1;d--) gmp_printf("%Zd x^%d + ",f[d],d);
	gmp_printf("%Zd x + %Zd\n",f[1],f[0]);
}

void printullpoly(ull *f,int d) {
	printf("(%d) ",d);
	for(;d>-1;d--) printf("%I64u ",f[d]);
}

/* calculate the norm of a-b*alpha without using division.
   needs the minimal polynomial f (must be monic!) and its
   degree d (f[] has d+1 elements, where f[0]=a_0, f[i]=a_i and f[d]=1).
   put the answer in r */
void calcnorm(mpz_t r,mpz_t a,mpz_t b,mpz_t *f,int d) {
	static mpz_t x[MAXDEG+1];
	mpz_t y,temp;
	int i;
	mpz_init(temp);
	for(i=0;i<=d;i++) mpz_init(x[i]);
	mpz_set_si(x[0],1);
	mpz_set(x[1],a);
	mpz_init_set_si(y,1);
	mpz_set_si(r,0);
	for(i=2;i<=d;i++) mpz_mul(x[i],x[i-1],a);
	for(i=d;i>=0;i--) {
		mpz_mul(temp,y,x[i]);
		mpz_addmul(r,temp,f[i]);
		mpz_mul(y,y,b);
	}
	mpz_clear(temp);
	mpz_clear(y);
	for(i=0;i<=d;i++) mpz_clear(x[i]);
}

/* factor base routines */
uchar *sieve;
#define SETBIT(p) sieve[(p)>>3]|=1<<((p)&7);
#define CLEARBIT(p) sieve[(p)>>3]&=~(1<<((p)&7));
#define CHECKBIT(p) (sieve[(p)>>3]&(1<<((p)&7)))

/* allocate and generate bit-packed sieve up to (not including) N */
void createsieve(ull N) {
	ull i,j;
	sieve=malloc((N+7)>>3);
	memset(sieve,0xaa,(N+7)>>3);
	sieve[0]=172;
	for(i=2;i*i<N;i++) if(CHECKBIT(i)) for(j=i*i;j<N;j+=i) CLEARBIT(j);
}

/* algebraic factor base */
ull *p1,*r1;
ull bn1;
/* rational factor base */
ull *p2;
ull bn2;
/* quadratic characters */
ull *p3,*r3;
ull bn3;

/* evaluate f(x), assume f monic */
void evalpoly(mpz_t *f,int deg,mpz_t x,mpz_t ret) {
	int i;
	mpz_set(ret,f[deg]);
	for(i=deg-1;i>=0;i--) {
		mpz_mul(ret,ret,x);
		mpz_add(ret,ret,f[i]);
	}
}

/* warning, requires 64-bit compiler, i think */
typedef __uint128_t ulll;
ull ullmulmod2(ull a,ull b,ull mod) { return (ulll)a*b%mod; }

/* evaluate f(x)%p, assume f monic. requires p<2^63 */
ull evalpolymod(ull *f,int df,ull x,ull p) {
	ull r=f[df];
	int i;
	for(i=df-1;i>=0;i--) r=(ullmulmod2(r,x,p)+f[i])%p;
	return r;
}

/* start of routine that finds all roots (aka linear factors) of a polynomial
   modulo a prime */
/* begins with various routines for doing polynomial arithmetic over Z_p */
/* in general, all routines that do stuff modulo m should be fed numbers in
   0, 1, ..., m-1 */

/* calculate inverse of a mod m (m can be composite, but 0 will be
   returned if an inverse doesn't exist). warning, don't use if m>=2^63 */
ll inverse(ll a,ll m) {
	ll b=m,x=0,y=1,t,q,lastx=1,lasty=0;
	while(b) {
		q=a/b;
		t=a,a=b,b=t%b;
		t=x,x=lastx-q*x,lastx=t;
		t=y,y=lasty-q*y,lasty=t;
	}
	return a==1?(lastx%m+m)%m:0;
}

/* modular square root! */

/* calculate the jacobi symbol, returns 0, 1 or -1 */
/* 1: a is quadratic residue mod m, -1: a is not, 0: a mod m=0 */
/* based on algorithm 2.3.5 in "prime numbers" (crandall, pomerance) */
/* WARNING, m must be an odd positive number */
int jacobi(ll a,ll m) {
	int t=1;
	ll z;
	a%=m;
	while(a) {
		while(!(a&1)) {
			a>>=1;
			if((m&7)==3 || (m&7)==5) t=-t;
		}
		z=a,a=m,m=z;
		if((a&3)==3 && (m&3)==3) t=-t;
		a%=m;
	}
	if(m==1) return t;
	return 0;
}

ull ullpowmod(ull n,ull k,ull mod) {
	int i,j;
	ull v=n,ans=1;
	if(!k) return 1;
	/* find topmost set bit */
	for(i=63;!(k&(1ULL<<i));i--);
	for(j=0;j<=i;j++) {
		if(k&(1ULL<<j)) ans=ullmulmod2(ans,v,mod);
		v=ullmulmod2(v,v,mod);
	}
	return ans;
}

/* calculate legendre symbol, returns 0, 1 or -1 */
/* 1: a is quadratic residue mod p, -1: a is not, 0: a mod p=0 */
/* WARNING, p must be an odd prime */
int legendre(ll a,ll p) {
	a%=p;
	if(a<0) a+=p;
	int z=ullpowmod(a,(p-1)>>1,p);
	return z==p-1?-1:z;
}

ull rand64() {
	return (rand()&32767) +
	       ((rand()&32767)<<15) +
	       ((rand()&32767ULL)<<30) +
	       ((rand()&32767ULL)<<45) +
	       ((rand()&15ULL)<<60);
}

/* find square root of a modulo p (p prime) using tonelli-shanks */
/* runtime O(ln^4 p) */
/* mod 3,5,7: algorithm 2.3.8 from "prime numbers" (crandall, pomerance) */
/* mod 1: from http://www.mast.queensu.ca/~math418/m418oh/m418oh11.pdf */
ull sqrtmod(ull a,ull p) {
	int p8,alpha,i;
	ull x,c,s,n,b,J,r2a,r;
	if(p==2) return a&1;
	a%=p;
	if(legendre(a,p)!=1) return 0; /* no square root */
	p8=p&7;
	if(p8==3 || p8==5 || p8==7) {
		if((p8&3)==3) return ullpowmod(a,(p+1)/4,p);
		x=ullpowmod(a,(p+3)/8,p);
		c=ullmulmod2(x,x,p);
		return c==a?x:ullmulmod2(x,ullpowmod(2,(p-1)/4,p),p);
	}
	alpha=0;
	s=p-1;
	while(!(s&1)) s>>=1,alpha++;
	r=ullpowmod(a,(s+1)/2,p);
	r2a=ullmulmod2(r,ullpowmod(a,(s+1)/2-1,p),p);
	do n=rand64()%(p-2)+2; while(legendre(n,p)!=-1);
	b=ullpowmod(n,s,p);
	J=0;
	for(i=0;i<alpha-1;i++) {
		c=ullpowmod(b,2*J,p);
		c=ullmulmod2(r2a,c,p);
		c=ullpowmod(c,1ULL<<(alpha-i-2),p);
		if(c==p-1) J+=1ULL<<i;
	}
	return ullmulmod2(r,ullpowmod(b,J,p),p);
}

/* set b(x)=a(x) */
void polyset(ull *a,int da,ull *b,int *db) {
	int i;
	for(*db=da,i=0;i<=*db;i++) b[i]=a[i];
}

/* set c(x)=a(x)+b(x) */
void polyaddmod(ull *a,int da,ull *b,int db,ull *c,int *dc,ull p) {
	static ull r[MAXDEG+1];
	int i,dr;
	dr=da>db?da:db;
	for(i=da+1;i<=dr;i++) r[i]=0;
	for(i=0;i<=da;i++) r[i]=a[i];
	for(i=0;i<=db;i++) {
		r[i]+=b[i];
		if(r[i]>=p) r[i]-=p;
	}
	while(dr>-1 && !r[dr]) dr--;
	*dc=dr;
	for(i=0;i<=dr;i++) c[i]=r[i];
}

/* negates a (modifies a) */
void polynegmod(ull *a,int da,ull p) {
	for(;da>-1;da--) a[da]=((ll)p-(ll)a[da])%p;
}

/* given polynomials a(x) and b(x), calculate quotient and
   remainder of a(x)/b(x) (mod p)
   dega, degb are the degrees of a and b, respectively. assume that *c, *d
   has enough pre-allocated memory to hold the results. don't assume that
   any of a,b,c,d are non-overlapping memory areas.
   if c is non-NULL, return quotient.
   if d is non-NULL, return remainder. remainder==0 has degree -1.
*/
void polydivmod(ull *a,int dega,ull *b,int degb,ull *c,int *degc,ull *d,int *degd,ull p) {
	static ull u[MAXDEG+1],q[MAXDEG+1];
	ull inv=inverse(b[degb],p);
	int k,j;
	for(k=0;k<=dega;k++) u[k]=a[k];
	for(k=dega-degb;k>-1;k--) {
		q[k]=ullmulmod2(u[degb+k],inv,p);
		for(j=degb+k-1;j>=k;j--) {
			u[j]=u[j]-ullmulmod2(q[k],b[j-k],p);
			if(u[j]>=p) u[j]+=p;
		}
	}
	if(c) for(*degc=dega-degb,k=*degc;k>-1;k--) c[k]=q[k];
	if(d) for(*degd=-1,k=0;k<degb && k<=dega;k++) if((d[k]=u[k])) *degd=k;
}

/* make polynomial monic, destroy input polynomial */
void polymonic(ull *a,int da,ull p) {
	ull z;
	int i;
	if(da<0 || a[da]==1) return;
	z=inverse(a[da],p);
	for(i=0;i<da;i++) a[i]=ullmulmod2(a[i],z,p);
	a[da]=1;
}

/* return a(x)*b(x) over Z_p */
void polymulmod(ull *a,int dega,ull *b,int degb,ull *c,int *degc,ull p) {
	static ull r[2*MAXDEG+1];
	int i,j;
	*degc=dega+degb;
	for(i=0;i<=*degc;i++) r[i]=0;
	for(i=0;i<=dega;i++) {
		for(j=0;j<=degb;j++) {
			r[i+j]=r[i+j]+ullmulmod2(a[i],b[j],p);
			if(r[i+j]>=p) r[i+j]-=p;
		}
	}
	for(i=0;i<=*degc;i++) c[i]=r[i];
	while(*degc>-1 && !c[*degc]) (*degc)--;
}

/* reduce a(x) mod v(x) over Z_p */
/* runtime: O(degree^2) */
void polyreduce(ull *a,int da,ull *v,int dv,ull *c,int *dc,ull p) {
	static ull w[2*MAXDEG+1];
	ull t;
	int i,j,z;
	for(i=0;i<=da;i++) w[i]=a[i];
	for(i=da+1;i<=dv;i++) w[i]=0;
	/* for each i=da, da-1, ..., dv, subtract a(i)*v(x)*x^(i-dv) */
	for(i=da;i>=dv;i--) for(j=0;j<=dv;j++) {
		z=i-dv; t=w[i];
		w[z+j]=(w[z+j]+p-ullmulmod2(t,v[j],p))%p;
	}
	/* tighten dc */
	for(*dc=-1,i=0;i<dv;i++) if((c[i]=w[i])) *dc=i;
}

/* given f, return g=f' (mod p) */
void polyderivemod(ull *f,int df,ull *g,int *dg,ull p) {
	int i;
	*dg=df-1;
	for(i=1;i<=df;i++) g[i-1]=ullmulmod2(f[i],i,p);
	while(*dg>-1 && !g[*dg]) (*dg)--;
}

/* return a(x)*b(x) mod v(x) over Z_p */
/* this can probably also be used to multiply two elements
   in the quotient ring Z_p/<v(x)> */
/* WARNING, not efficient. integrate mulmod and reduce more tightly */
void polymulmodmod(ull *a,int da,ull *b,int db,ull *v,int dv,ull *c,int *dc,ull p) {
	static ull d[2*MAXDEG+1];
	int dd;
	polymulmod(a,da,b,db,d,&dd,p);
	polyreduce(d,dd,v,dv,c,dc,p);
}

/* return a(x)^n mod v(x) over Z_p, put result in c */
/* warning, not very efficient, really, but care about that later */
void polypowmodmod(ull *a,int da,ull n,ull *v,int dv,ull *c,int *dc,ull p) {
	ull z[MAXDEG+1],y[MAXDEG+1]={1};
	int dz,dy=0,i;
	polyset(a,da,z,&dz);
	while(n) {
		if(n&1) {
			n>>=1;
			polymulmodmod(y,dy,z,dz,v,dv,y,&dy,p);
			if(!n) break;
		} else n>>=1;
		polymulmodmod(z,dz,z,dz,v,dv,z,&dz,p);
	}
	for(*dc=dy,i=0;i<=*dc;i++) c[i]=y[i];
}

/* return a(x)^n mod v(x) over Z_p, put result in c, exponent is mpz */
/* warning, not very efficient, really, but care about that later */
void polypowmodmodmpz(ull *a,int da,mpz_t N,ull *v,int dv,ull *c,int *dc,ull p) {
	ull z[MAXDEG+1],y[MAXDEG+1]={1};
	int dz=da,dy=0,i;
	mpz_t t,n;
	mpz_init(t);
	mpz_init_set(n,N);
	for(i=0;i<=dz;i++) z[i]=a[i];
	while(mpz_cmp_si(n,0)>0) {
		if(mpz_mod_ui(t,n,2)) {
			mpz_fdiv_q_2exp(n,n,1);
			polymulmodmod(y,dy,z,dz,v,dv,y,&dy,p);
			if(!mpz_cmp_si(n,0)) break;
		} else mpz_fdiv_q_2exp(n,n,1);
		polymulmodmod(z,dz,z,dz,v,dv,z,&dz,p);
	}
	for(*dc=dy,i=0;i<=*dc;i++) c[i]=y[i];
	mpz_clear(n);
	mpz_clear(t);
}

/* return a(x)^n over Z_p */
void polypowmod(ull *a,int da,ull n,ull *c,int *dc,ull p) {
	ull z[MAXDEG+1],y[MAXDEG+1]={1};
	int dz=da,dy=0,i;
	for(i=0;i<=dz;i++) z[i]=a[i];
	while(n) {
		if(n&1) {
			n>>=1;
			polymulmod(y,dy,z,dz,y,&dy,p);
			if(!n) break;
		} else n>>=1;
		polymulmod(z,dz,z,dz,z,&dz,p);
	}
	for(*dc=dy,i=0;i<=*dc;i++) c[i]=y[i];
}

/* given polynomials a(x), b(x), calculate g(x)=gcd(a(x),b(x)) mod p. */
void polygcdmod(ull *a,int da,ull *b,int db,ull *g,int *dg,ull p) {
	static ull c[MAXDEG+1],d[MAXDEG+1],e[MAXDEG+1];
	int dc,dd,de,i;
	polyset(a,da,c,&dc);
	polyset(b,db,d,&dd);
	/* sanity check: a==0 */
	if(da<0) {
		for(*dg=dd,i=0;i<=*dg;i++) g[i]=d[i];
		goto end;
	}
	while(dd>-1) {
		polydivmod(c,dc,d,dd,NULL,NULL,e,&de,p);
		polyset(d,dd,c,&dc);
		polyset(e,de,d,&dd);
	}
	polyset(c,dc,g,dg);
end:
	/* make output monic */
	polymonic(g,*dg,p);
}

/* calculate inverse of a mod m (m can be composite, but 0 will be
   returned if an inverse doesn't exist). warning, don't use if m>=2^63 */
ll inversemal(ll a,ll m) {
	ll b=m,x=0,y=1,t,q,lastx=1,lasty=0;
	while(b) {
		q=a/b;
		t=a,a=b,b=t%b;
		t=x,x=lastx-q*x,lastx=t;
		t=y,y=lasty-q*y,lasty=t;
	}
	return a==1?(lastx%m+m)%m:0;
}

/* find the inverse g(x)=a^1(x) of a(x) mod f(x) mod p using
   the extended euclid algorithm */
/* f(x) is assumed to be monic. if an inverse doesn't exist, return g=0 */
void polyinversemodmod(ull *in,int din,ull *f,int df,ull *g,int *dg,ull p) {
	ull b[MAXDEG+1],x[MAXDEG+1],y[MAXDEG+1],lastx[MAXDEG+1],lasty[MAXDEG+1];
	ull t[MAXDEG+1],q[MAXDEG+1],a[MAXDEG+1],z[MAXDEG+1],v;
	int db,dx,dy,lastdx,lastdy,dt,dq,da,dz,i;
	if(din<0) { *dg=-1; return; }
	polyset(f,df,b,&db);
	polyset(in,din,a,&da);
	dx=-1; y[0]=1; dy=0;
	lastx[0]=1; lastdx=0; lastdy=-1;
	while(db>-1) {
		/* set a=b, b=a%b, q=a/b */
		polyset(a,da,t,&dt);
		polyset(b,db,a,&da);
		polydivmod(t,dt,a,da,q,&dq,b,&db,p);
		/* set x=lastx-q*x, lastx=x */
		polyset(x,dx,t,&dt);
		polymulmod(q,dq,x,dx,z,&dz,p);
		polyset(lastx,lastdx,x,&dx);
		polynegmod(z,dz,p);
		polyaddmod(x,dx,z,dz,x,&dx,p);
		polyset(t,dt,lastx,&lastdx);
		/* set y=lasty-q*y, lasty=y */
		polyset(y,dy,t,&dt);
		polymulmod(q,dq,y,dy,z,&dz,p);
		polyset(lasty,lastdy,y,&dy);
		polynegmod(z,dz,p);
		polyaddmod(y,dy,z,dz,y,&dy,p);
		polyset(t,dt,lasty,&lastdy);
	}
	/* now a is gcd(a,f). if !=1 return failure */
	if(da>0) { *dg=-1; return; }
	/* lastx is inverse, multiply with inverse of a[0] */
	if(a[0]!=1) {
		v=inverse(a[0],p);
		for(i=0;i<=lastdx;i++) lastx[i]=ullmulmod2(lastx[i],v,p);
	}
	for(*dg=lastdx,i=0;i<=lastdx;i++) g[i]=lastx[i];
}

/* return 1 if u(x) is squarefree. u is squarefree iff gcd(u,u')==1.
   u(x) must be monic. unpredictable results if deg u <= p */
int ispolymodsquarefree(ull *u,int du,ull p) {
	static ull ud[MAXDEG+1],g[MAXDEG+1];
	int dud,dg,i;
	for(dud=du-1,i=0;i<du;i++) ud[i]=ullmulmod2(u[i+1],i+1,p);
	while(dud>-1 && !ud[dud]) dud--;
	polygcdmod(u,du,ud,dud,g,&dg,p);
	return dg==0;
}

/* find all roots by naive method (evaluate in(x) for all x), inefficient */
void polylinmodnaive(ull *in,int dv,ull p,ull *f,int *fn) {
	ll x;
	for(*fn=x=0;x<p;x++) if(!evalpolymod(in,dv,x,p)) f[(*fn)++]=x;
}

/* find roots of u(x) mod p, p must be an odd prime larger than
   the degree of u(x) */
/* based on algorithm 1.6.1 in cohen */
void polyfindrootmod(ull *z,int dz,ull p,ull *f,int *fn) {
	/* cast out gcd(f',f) */
	static ull g[MAXDEG+1],ud[MAXDEG+1],u[MAXDEG+1],m1[MAXDEG+1];
	static ull q[MAXDEG+1][MAXDEG+1];
	ull d,e;
	int du,dg,dud,qn=1,done,i,dm1;
	static int dq[MAXDEG+1];
	*fn=0;
	polyset(z,dz,u,&du);
	/* force u monic */
	polymonic(u,du,p);
	polyderivemod(u,du,ud,&dud,p);
	polygcdmod(u,du,ud,dud,g,&dg,p);
	/* force gcd monic */
	polymonic(g,dg,p);
	/* divide out squares */
	polydivmod(u,du,g,dg,u,&du,NULL,NULL,p);
	/* cast out 0-factor */
	if(!u[0]) {
		g[0]=0; g[1]=1; dg=1;
		polydivmod(u,du,g,dg,u,&du,NULL,NULL,p);
		f[(*fn)++]=0;
	}
	/* m1(x)=-1 (p-1) */
	m1[0]=p-1; dm1=0;
	/* take gcd(x^(p-1)-1, u(x)) and isolate roots */
	/* first take d=x^(p-1) mod u, then take gcd(d-1,u) */
	g[0]=0; g[1]=1; dg=1;
	polypowmodmod(g,dg,p-1,u,du,g,&dg,p);
	polyaddmod(g,dg,m1,dm1,g,&dg,p);
	polygcdmod(g,dg,u,du,q[0],&dq[0],p);
	do {
		done=1;
		/* if deg>2, try to split polynomial. benchmarking shows it's faster
		   to split down to deg 2 rather than deg 1. */
		for(i=0;i<qn;i++) if(dq[i]>2) {
			do {
				g[0]=rand64()%p; g[1]=1; dg=1;
				polypowmodmod(g,dg,p>>1,q[i],dq[i],g,&dg,p);
				polyaddmod(g,dg,m1,dm1,g,&dg,p);
				polygcdmod(g,dg,q[i],dq[i],g,&dg,p);
			} while(!dg || dg==dq[i]);
			polydivmod(q[i],dq[i],g,dg,q[i],&dq[i],NULL,NULL,p);
			polyset(g,dg,q[qn],&dq[qn]);
			qn++;
			done=0;
		}
	} while(!done);
	/* go through each item in the list, and output roots */
	for(i=0;i<qn;i++) {
		if(dq[i]==1) {
			if(q[i][1]==1) f[(*fn)++]=(p-q[i][0])%p;
			else f[(*fn)++]=ullmulmod2((p-q[i][0])%p,inverse(q[i][1],p),p);
		} else if(dq[i]==2) {
			d=ullmulmod2(q[i][1],q[i][1],p);
			e=ullmulmod2(q[i][0],q[i][2],p);
			e=sqrtmod((d+p-ullmulmod2(e,4,p))%p,p);
			d=ullmulmod2(inverse(2,p),q[i][2],p);
			f[(*fn)++]=ullmulmod2((p+e-q[i][1])%p,d,p);
			f[(*fn)++]=ullmulmod2((p+p-e-q[i][1])%p,d,p);
		}
	}
}

/* entry point for new routine */
void findideals2(ull *u,int du,ull p,ull *f,int *fn) {
	/* naive algorithm for small enough p: evaluate f(r) for all 0<=r<p */
	if(p<200 || p<=du) return polylinmodnaive(u,du,p,f,fn);
	polyfindrootmod(u,du,p,f,fn);
}

/* determinant using stupid and slow O(n!) algorith, but n will never
   be huge (say, never larger than 6 and in practice it will always be 3).
   generate permutations using fancy loop-free algorithm by knuth
   where successively generated permutations have alternating parity
   [an easy O(n^3) algorithm: gauss-jordan and return product of diagonal
   times the numbers we divided the rows with] */
ull calcdet(ull A[MAXDEG+1][MAXDEG+1],int n,ull p) {
	ull res=0,r;
  int o[100],c[100],j,s,q,a[100],sign=1;
  char t;
  for(j=0;j<n;j++) c[j]=0,o[j]=1,a[j]=j;
p2:
	/* visit permutation */
	r=sign?1:p-1;
	for(j=0;j<n;j++) r=ullmulmod2(r,A[j][a[j]],p);
	res+=r;
	if(res>=p) res-=p;
	sign^=1;
	/* end visit */
  j=n; s=0;
p4:
  q=c[j-1]+o[j-1];
  if(q<0) goto p7;
  if(q==j) goto p6;
  t=a[j-c[j-1]+s-1]; a[j-c[j-1]+s-1]=a[j-q+s-1]; a[j-q+s-1]=t;
  c[j-1]=q;
  goto p2;
p6:
  if(j==1) return res;
  s++;
p7:
  o[j-1]=-o[j-1]; j--;
  goto p4;
}

/* calculate norm mod p of general element a(x) in field with minimal
   polynomial f(x). uses determinant method */
/* tested against calcnorm() with tens of millions of numbers of the form
   a+b*alpha with degrees 3-6, with a,b huge modulo a huge prime */
ull calcnormmod(ull *a,int da,ull *f,int df,ull p) {
	static ull A[MAXDEG+1][MAXDEG+1];
	ull b[MAXDEG+1]={0,1},c[MAXDEG+1];
	int i,j,db=1,dc;
	polyset(a,da,c,&dc);
	for(i=0;i<=dc;i++) A[i][0]=a[i];
	for(;i<df;i++) A[i][0]=0;
	for(j=1;j<df;j++) {
		polymulmodmod(c,dc,b,db,f,df,c,&dc,p);
		for(i=0;i<=dc;i++) A[i][j]=c[i];
		for(;i<df;i++) A[i][j]=0;
	}
	return calcdet(A,df,p);
}

/* B1 and B2 are upper bound for primes (algebraic and rational)
   f is polynomial, deg is degree
   p1,r1 is algebraic factor base, bn1 is number of primes
   p2 is rational factor base, bn2 is number of primes */
void createfactorbases(ull B1,ull B2,ull Bk,mpz_t *f,int deg,ull **_p1,ull **_r1,ull *bn1,ull **_p2,ull *bn2,ull **_p3,ull **_r3,ull *bn3) {
	static ull b[MAXDEG+1];
	static ull root[MAXDEG+1];
	ull B=B1>B2?B1:B2,i,j,q;
	ull *p1,*r1,*p2,*p3,*r3;
	int fn;
	int db,k;
	char *sieve=malloc(B+1);
	memset(sieve,1,B+1);
	for(i=2;i*i<=B;i++) if(sieve[i]) for(j=i*i;j<=B;j+=i) sieve[j]=0;
	/* generate rational factor base */
	for(*bn2=0,i=2;i<=B2;i++) if(sieve[i]) (*bn2)++;
	if(!(p2=malloc(*bn2*sizeof(ull)))) error("couldn't allocate rational factor base");
	for(*bn2=0,i=2;i<=B2;i++) if(sieve[i]) p2[(*bn2)++]=i;

	/* generate algebraic factor base */
	for(*bn1=0,i=2;i<=B1;i++) if(sieve[i]) {
		/* find all eligible r: r such that f(r)=0 (mod p) using factorization */
		db=deg;
		for(k=0;k<=db;k++) b[k]=mpz_mod_ull(f[k],i);
		findideals2(b,db,i,root,&fn);
		*bn1+=fn;
	}
	if(!(p1=malloc(*bn1*sizeof(ull)))) error("couldn't allocate algebraic factor base");
	if(!(r1=malloc(*bn1*sizeof(ull)))) error("couldn't allocate algebraic factor base");
	for(*bn1=0,i=2;i<=B1;i++) if(sieve[i]) {
		/* find all roots again. we happily waste some computing resources since
		   the sieve stage will dominate the runtime anyway */
		/* slow method again TODO replace with factorization */
		db=deg;
		for(k=0;k<=db;k++) b[k]=mpz_mod_ull(f[k],i);
		findideals2(b,db,i,root,&fn);
		for(j=0;j<fn;j++) p1[*bn1]=i,r1[(*bn1)++]=root[j];
	}

	/* generate quadratic characters */
	*bn3=Bk;
	if(!(p3=malloc(*bn3*sizeof(ull)))) error("couldn't allocate quadratic characters");
	if(!(r3=malloc(*bn3*sizeof(ull)))) error("couldn't allocate quadratic characters");
	for(i=0,q=B1+1;i<*bn3;q++) {
		/* check if q is prime */
		for(j=0;j<*bn2 && p2[j]*p2[j]<=q;j++) if(q%p2[j]==0) goto noprime;
		db=deg;
		for(k=0;k<=db;k++) b[k]=mpz_mod_ull(f[k],q);
		findideals2(b,db,q,root,&fn);
		if(!fn) continue;
		/* find value from root such that f'(value)!=0 mod q */
		polyderivemod(b,db,b,&db,q);
		for(k=0;k<fn;k++) if(evalpolymod(b,db,root[k],q)) {
			p3[i]=q;
			r3[i]=root[k];
			i++;
			break;
		}
	noprime:;
	}

	// TEMP output factor bases
//	for(i=0;i<*bn1;i++) printf("(%I64d %I64d)\n",p1[i],r1[i]);printf("\n");
//	for(i=0;i<*bn2;i++) printf("%I64d ",p2[i]);printf("\n");
//	for(i=0;i<*bn3;i++) printf("(%I64d %I64d)\n",p3[i],r3[i]);printf("\n");
	free(sieve);
	*_p1=p1; *_r1=r1; *_p2=p2; *_p3=p3; *_r3=r3;
}

/* return index of v in p, or -1 if it doesn't exist */
ull bs(ull *p,ull bn,ull v) {
	ull lo=0,hi=bn,mid;
	while(lo<hi) {
		mid=(lo+hi)>>1;
		if(v>p[mid]) lo=mid+1;
		else hi=mid;
	}
	return lo<bn && p[lo]==v?lo:-1;
}

/* matrix (global) */
uint **M;
int notsmooth,missed,smooth;

/* gaussian elimination mod 2 on bitmasks, A is n*m, b is n*o */
/* a is a malloced array of pointers, each a[i] is of size
   sizeof(uint)*(m+o+31)/32 */
/* return 0: no solutions, 1: one solution, 2: free variables */
#define ISSET(a,row,col) (a[(row)][(col)>>5]&(1U<<((col)&31)))
#define MSETBIT(a,row,col) a[(row)][(col)>>5]|=(1U<<((col)&31))
#define MTOGGLEBIT(a,row,col) a[(row)][(col)>>5]^=(1U<<((col)&31))
int bitgauss32(uint **a,int n,int m,int o) {
	int i,j,k,z=m+o,c=0,fri=0,bz=(z+31)>>5;
	uint t;
	/* process each column */
	for(i=0;i<m;i++) {
		/* TODO check words instead of bits */
		for(j=c;j<n;j++) if(ISSET(a,j,i)) break;
		if(j==n) { fri=1; continue; }
		/* swap? */
		if(j>c)  for(k=0;k<bz;k++) {
			t=a[j][k],a[j][k]=a[c][k],a[c][k]=t;
		}
		/* subtract multiples of this row */
		for(j=0;j<n;j++) if(j!=c && ISSET(a,j,i)) {
			for(k=0;k<bz;k++) a[j][k]^=a[c][k];
		}
		c++;
	}
	/* detect no solution: rows with 0=b */
	for(i=0;i<n;i++) {
		/* TODO make bit-efficient solution later */
		for(j=0;j<m;j++) if(ISSET(a,i,j)) goto ok;
		for(;j<z;j++) if(ISSET(a,i,j)) return 0;
	ok:;
	}
	return 1+fri;
}

/* find all free variables: variable i is free if there is no row having its first
   1-element in column i */
int findfreevars(uint **a,int rows,int cols,uchar *freevar) {
	int i,j,r=cols;
	memset(freevar,1,cols);
	for(i=0;i<rows;i++) {
		for(j=0;j<cols;j++) if(ISSET(a,i,j)) {
			freevar[j]=0;
			r--;
			break;
		}
	}
	return r;
}

/* find exponents of square. id is a bitmask that contains the values of the
   free variables from the linear algebra
	 rows: factor base
	 cols: relations */
void getsquare(uint **a,int rows,int cols,uchar *freevar,int id,uchar *v) {
	int i,j,k;
	memset(v,0,cols);
	/* set id-th free variable */
	for(j=i=0;i<cols;i++) if(freevar[i]) {
		if(id==j) { v[i]=1; break; }
		j++;
	}
	/* get solution vector by back substitution! set the first 1-element to the
	   xor of the others. */
	for(i=rows-1;i>=0;i--) {
		for(j=0;j<cols;j++) if(ISSET(a,i,j)) goto ok;
		continue;
	ok:
		for(k=j++;j<cols;j++) if(ISSET(a,i,j) && v[j]) v[k]^=1;
	}
}

/* store rational factors for pairs (a,b) */
ull **faclist;
int *facn;
/* store algebraic factors for pairs (a,b) */
ull **alglist;
int *algn;

/* get rational square root! */
void getratroot(mpz_t n,uchar *v,int cols,mpz_t *f,int df,mpz_t m,mpz_t root,int *aval,int *bval) {
	mpz_t t;
	static mpz_t fd[MAXDEG+1];
	static int *ev;
	int dfd;
	mpz_init(t);
	mpz_set_si(root,1);
	int i,j;
	ev=calloc(bn2,sizeof(int));
	if(!ev) error("out of memory");
	for(i=0;i<cols;i++) if(v[i]) for(j=0;j<facn[i];j++) ev[faclist[i][j]]++;
	/* sanity */
	for(i=0;i<bn2;i++) if(ev[i]&1) error("odd exponent in rat");
	for(i=0;i<bn2;i++) if(ev[i]) {
		mpz_set_ull(t,p2[i]);
		for(j=0;j+j<ev[i];j++) mpz_mul(root,root,t);
		mpz_mod(root,root,n);
	}
	/* multiply value with f'(m)^2 */
	dfd=df-1;
	for(i=0;i<=dfd;i++) {
		mpz_init_set(fd[i],f[i+1]);
		mpz_mul_ui(fd[i],fd[i],i+1);
	}
	evalpoly(fd,dfd,m,t);
	mpz_mod(t,t,n); /* t = f'(m) mod n */
	mpz_mul(root,root,t); /* multiply in f'(m) */
	mpz_mod(root,root,n); /* and reduce mod n */
	mpz_mul(t,root,root);
	mpz_mod(t,t,n);
	gmp_printf("rational root: %Zd, square %Zd\n",root,t);
	for(i=0;i<=dfd;i++) mpz_clear(fd[i]);
	free(ev);
	mpz_clear(t);
}

/* start of routines for algebraic square root */

/* return 1 if f(x) is irreducible mod p */
int polyirredmod(mpz_t *in,int df,ull p) {
	/* check if gcd(x^(p^d)-x,f) is a non-constant
	   polynomial for 1<=d<=df/2 */
	static ull g[MAXDEG+1],h[MAXDEG+1],f[MAXDEG+1];
	int dg,i,j,dh;
	for(i=0;i<=df;i++) f[i]=mpz_mod_ull(in[i],p);
	for(i=1;i+i<=df;i++) {
		/* form x^p^i - x */
		/* use that x^p^i = ((x^p)^p) ... ^p (i times) */
		g[0]=0; g[1]=1; dg=1;
		for(j=0;j<i;j++) polypowmodmod(g,dg,p,f,df,g,&dg,p);
		h[0]=0; h[1]=p-1; dh=1;
		polyaddmod(g,dg,h,dh,g,&dg,p);
		polygcdmod(g,dg,f,df,g,&dg,p);
		if(dg>0) return 0;
	}
	return 1;
}

/* calculate the legendre symbol of the element a (in polynomial format)
   in the field F_p^df:
   1 if element is a quadratic residue, -1 if not.
	 p must be an odd prime! */
int polylegendre(ull *a,int da,ull *f,int df,ull p) {
	ull b[MAXDEG+1];
	mpz_t n,P;
	int db,i;
	for(i=0;i<=da;i++) if(a[i]) goto notzero;
	return 0;
notzero:
	mpz_init(n);
	mpz_init(P);
	mpz_set_ull(P,p);
	mpz_pow_ui(n,P,df);
	mpz_sub_ui(n,n,1);
	mpz_divexact_ui(n,n,2);
	polypowmodmodmpz(a,da,n,f,df,b,&db,p);
	mpz_clear(n);
	mpz_clear(P);
	if(b[0]==p-1) return -1;
	if(b[0]==1) return 1;
	error("error in polylegendre, res not 1 or -1");
	return 0;
}

int findexpdiv2(mpz_t P,int df) {
	mpz_t s;
	int r=0;
	mpz_init(s);
	mpz_pow_ui(s,P,df);
	mpz_sub_ui(s,s,1);
	while(!mpz_tstbit(s,0)) {
		r++;
		mpz_fdiv_q_2exp(s,s,1);
	}
	mpz_clear(s);
	return r;
}

/* given a, find b such that b^2=a in the field F_{p^df} given by the
   minimal polynomial f with degree df */
/* based on description in briggs */
/* algorithm is pretty much tonelli-shanks, adapted to F_{p^df} */
/* warning, i took a dubious short cut when implementing. p^df-1 should
   not have a divisor 2^s for a large s. this was circumvented by avoiding
   finite fields with this property */
void polysqrtmod(ull *a,int da,ull *f,int df,ull *b,int *db,ull p) {
	mpz_t s,z;
	ull j,c[MAXDEG+1],d[MAXDEG+1],e[MAXDEG+1];
	int r=0,dc,i,dd,t,de;
	/* does the square root exist? */
	if(1!=polylegendre(a,da,f,df,p)) { *db=-1; printf("not a square\n"); return; }
	mpz_init(s);
	mpz_init(z);
	/* write p^df-1 as 2^r * s for s odd */
	mpz_set_ull(s,p);
	mpz_pow_ui(s,s,df);
	mpz_sub_ui(s,s,1);
	while(!mpz_tstbit(s,0)) {
		r++;
		mpz_fdiv_q_2exp(s,s,1);
	}
	if(r>10) error("error, unsuitable r");
	/* find an element in F_{p^df} which is a non-residue */
	for(j=1;;j++) {
		for(dc=df-1,i=0;i<=dc;i++) c[i]=j;
//		for(dc=1,i=0;i<=dc;i++) c[i]=j;
		if(-1==polylegendre(c,dc,f,df,p)) break;
	}
	/* d=a^s */
	polypowmodmodmpz(a,da,s,f,df,d,&dd,p);
	/* find t such that c^2st = d. guaranteed to be <2^r */
	for(t=0;t<(1<<r);t++) {
		mpz_mul_ui(z,s,2*t);
		polypowmodmodmpz(c,dc,z,f,df,e,&de,p);
		/* c^2st == d? */
		if(de==dd) {
			for(i=0;i<=de;i++) if(e[i]!=d[i]) goto noteq;
			goto eq;
		}
	noteq:;
	}
	error("didn't find t in sqrt");
eq:;
	mpz_mul_ui(z,s,t);
	polypowmodmodmpz(c,dc,z,f,df,e,&de,p);
	/* calculate the inverse of e */
	polyinversemodmod(e,de,f,df,e,&de,p);
	/* the root is a^(s+1)/2 * e^-1 */
	mpz_add_ui(s,s,1);
	mpz_fdiv_q_2exp(s,s,1);
	polypowmodmodmpz(a,da,s,f,df,c,&dc,p);
	polymulmodmod(c,dc,e,de,f,df,b,db,p);
	mpz_clear(z);
	mpz_clear(s);
}

/* here follows some subroutines for polynomial arithmetic over Z */

/* multiply two polynomials, c(x)=a(x)*b(x) */
void polymulmpz(mpz_t *a,int da,mpz_t *b,int db,mpz_t *c,int *dc) {
	static mpz_t r[2*BIGDEG+2];
	int i,j;
	for(i=0;i<=da+db;i++) mpz_init_set_ui(r[i],0);
	for(i=0;i<=da;i++) for(j=0;j<=db;j++) mpz_addmul(r[i+j],a[i],b[j]);
	for(*dc=da+db,i=0;i<=*dc;i++) mpz_set(c[i],r[i]);
	for(i=0;i<=da+db;i++) mpz_clear(r[i]);
}

/* reduce a(x) mod f(x), return result in b(x) */
void polyreducempz(mpz_t *a,int da,mpz_t *f,int df,mpz_t *b,int *db) {
	mpz_t w[2*BIGDEG+2],t;
	int i,j,z;
	mpz_init(t);
	for(i=0;i<=da;i++) mpz_init_set(w[i],a[i]);
	for(;i<=df;i++) mpz_init_set_ui(w[i],0);
	/* for each i=da, da-1, ..., dv, subtract a(i)*v(x)*x^(i-dv) */
	for(i=da;i>=df;i--) for(j=0;j<=df;j++) {
		z=i-df;
		mpz_set(t,w[i]);
		mpz_submul(w[z+j],t,f[j]);
	}
	for(i=0;i<df;i++) mpz_set(b[i],w[i]);
	*db=df-1;
	/* tighten db */
	while(*db>-1 && !mpz_cmp_si(b[*db],0)) (*db)--;
	for(i=0;i<=da;i++) mpz_clear(w[i]);
	for(;i<=df;i++) mpz_clear(w[i]);
	mpz_clear(t);
}

/* given f, return g=f' */
void polyderivempz(mpz_t *f,int df,mpz_t *g,int *dg) {
	int i;
	*dg=df-1;
	for(i=1;i<=df;i++) mpz_mul_si(g[i-1],f[i],i);
	while(*dg>-1 && !mpz_cmp_si(g[*dg],0)) (*dg)--;
}

/* calculate the algebraic number and display it */
void printalgnum(mpz_t n,uchar *v,int cols,mpz_t *f,int df,mpz_t m,int *aval,int *bval) {
	mpz_t a[2*BIGDEG+2],b[BIGDEG+1];
	int da,db,i;
	for(i=0;i<2*BIGDEG+2;i++) mpz_init_set_ui(a[i],i==0);
	for(i=0;i<BIGDEG+1;i++) mpz_init(b[i]);
	da=0;
	/* multiply with f'(alpha)^2 */
	polyderivempz(f,df,b,&db);
	polymulmpz(a,da,b,db,a,&da);
	polyreducempz(a,da,f,df,a,&da);
	polymulmpz(a,da,b,db,a,&da);
	polyreducempz(a,da,f,df,a,&da);
	for(i=0;i<cols;i++) if(v[i]) {
		mpz_set_si(b[0],aval[i]);
		mpz_set_si(b[1],-bval[i]);
		db=1;
		polymulmpz(a,da,b,db,a,&da);
		polyreducempz(a,da,f,df,a,&da);
	}
	printf("algebraic square:\n");
	printmpzpoly(a,da);
	for(i=0;i<BIGDEG+1;i++) mpz_clear(b[i]);
	for(i=0;i<2*BIGDEG+2;i++) mpz_clear(a[i]);
}

/* get algebraic square root! v is the subset of (a,b) pairs */
/* use couveignes' algorithm */
int getalgroot(mpz_t n,uchar *v,int cols,mpz_t *in,int df,mpz_t m,mpz_t root,int *aval,int *bval) {
	double logest=0,b;
	mpz_t P,M,temp,ans;
	ull *q,pp,*ai,f[MAXDEG+1],fd[MAXDEG+1],g[MAXDEG+1],h[MAXDEG+1];
	ull n1,n2,xi;
	const ull MAX=(1ULL<<61)-1; /* start here to check for primes */
	int i,s,maxu,qn,dfd,j,dg,dh,k,ret=0;
	double zp;
	static int *ev;
	/* populate exponent vector */
	ev=calloc(bn1,sizeof(int));
	if(!ev) error("out of memory");
	for(i=0;i<cols;i++) if(v[i]) for(j=0;j<algn[i];j++) {
		if(alglist[i][j]<0 || alglist[i][j]>=bn1) error("error");
		ev[alglist[i][j]]++;
	}
	mpz_init(P);
	mpz_init_set_si(M,1);
	mpz_init(temp);
	mpz_init(ans);
	mpz_set_ui(ans,0);
	/* rough estimate:
	   d^(d+5)/2 * n * (2*u*sqrt(d)*m)^(s/2)
	   calculate log2 of this since it's huge */
	/* if this turns out to be bad, check the paper of couveignes for a
	   tighter bound using complex roots and direct evaluation of stuff */
	logest=log2(df)*(df+5)*.5;
	logest+=mpz_sizeinbase(n,2);
	/* get u and s */
	maxu=0;
	for(i=0;i<cols;i++) {
		if(maxu<-aval[i]) maxu=-aval[i];
		if(maxu<aval[i]) maxu=aval[i];
		if(maxu<-bval[i]) maxu=-bval[i];
		if(maxu<bval[i]) maxu=bval[i];
	}
	for(s=i=0;i<cols;i++) s+=v[i];
	b=2*maxu*sqrt(df)*mpz_get_d(m);
	logest+=s*.5*log2(b);
	printf("estimate: %f bits\n",logest);
	/* find multiple q such that their product has >= logest digits */
	qn=(int)(1+logest/log2(MAX));
	q=malloc(qn*sizeof(ull));
	if(!q) error("out of memory in algroot");
	ai=malloc(qn*sizeof(ull));
	if(!ai) error("out of memory in algroot");
	/* don't be super duper tight and take primes just below 2^63.
	   it seems there are overflow issues in some of the subroutines,
	   the suspects are polyderivemod and polymulmodmod (and their callees) */
	for(pp=MAX,i=0;i<qn;pp+=2) {
		mpz_set_ull(P,pp);
		/* P must be prime and f(x) mod P must be irreducible */
		if(!mpz_probab_prime_p(P,30)) continue;
		if(!polyirredmod(in,df,pp)) continue;
		/* we also want to avoid P such that 2^r for large r divides P^df-1 */
		if(findexpdiv2(P,df)>5) continue;
		q[i++]=pp;
		mpz_mul(M,M,P);
	}
	/* for each i, compute a_i */
	for(zp=i=0;i<qn;i++) {
		mpz_set_ull(P,q[i]);
		mpz_fdiv_q(temp,M,P);
		pp=mpz_mod_ull(temp,q[i]);
		ai[i]=inverse(pp,q[i]);
	}
	/* for each q_i, calculate f'^2 * prod(a-bx) mod f, mod q_i
	   and calculate its square root in Z_p/<f> */
	dfd=df-1;
	for(i=0;i<qn;i++) {
//		printf("prime %I64d\n",q[i]);
		for(j=0;j<=df;j++) f[j]=mpz_mod_ull(in[j],q[i]);
		polyderivemod(f,df,fd,&dfd,q[i]);
		/* form f'^2 * prod_{(a,b)} (a-b*alpha) mod q[i] */
		polymulmodmod(fd,dfd,fd,dfd,f,df,g,&dg,q[i]);
		for(j=0;j<cols;j++) if(v[j]) {
			h[0]=(aval[j]%(ll)q[i]+(ll)q[i])%(ll)q[i];
			h[1]=((-(ll)bval[j])%(ll)q[i]+(ll)q[i])%(ll)q[i];
			dh=1;
			polymulmodmod(g,dg,h,dh,f,df,g,&dg,q[i]);
		}
//		printf("%d: take sqrt of ",i);printullpoly(g,dg);printf("mod %I64d\n",q[i]);
//		printf("norm before root: %I64d\n",calcnormmod(g,dg,f,df,q[i]));
		/* take square root of g */
		polysqrtmod(g,dg,f,df,g,&dg,q[i]);
		/* sanity, not a square */
		if(dg<0) {
			printf("failed in %d of %d\n",i+1,qn);
			puts("error!");
			printf("p %I64d, f(x) mod p = ",q[i]);
			printullpoly(f,df);printf("\n");
			printf("g(x) mod p is not square: ");
			printullpoly(g,dg);printf("\n");
			goto quit;
		}
		/* norm of root (in g,dg) */
		n1=calcnormmod(g,dg,f,df,q[i]);
		/* norm of f'(alpha) */
		n2=calcnormmod(fd,dfd,f,df,q[i]);
		/* norm of square root of all prime ideals */
		for(j=0;j<bn1;j++) if(ev[j]) {
			/* norm of prime factor represented by the pair (p,r) is p */
			for(k=0;k+k<ev[j];k++) n2=ullmulmod2(n2,p1[j],q[i]);
		}
		/* if the norms are different, negate the root */
//		printf("before neg %I64d %I64d, n^2=%I64d\n",n1,n2,ullmulmod2(n2,n2,q[i]));
		if(n1!=n2) for(j=0;j<=dg;j++) g[j]=(q[i]-g[j])%q[i];
//		printf("sqrt is ");printullpoly(g,dg);printf("\n");
		n1=calcnormmod(g,dg,f,df,q[i]);
//		printf("after  neg %I64d %I64d\n",n1,n2);
		if(n1!=n2) { printf("error %d/%d, norms are not equal!\n",i+1,qn); goto quit; }
		/* calculate a_i*x_i*P_i mod n and add it to result */
		mpz_set_ull(P,q[i]);
		mpz_fdiv_q(temp,M,P);
		mpz_set_ull(P,ai[i]);
		mpz_mul(temp,temp,P);
		xi=evalpolymod(g,dg,mpz_mod_ull(m,q[i]),q[i]);
		mpz_set_ull(P,xi);
		mpz_mul(temp,temp,P);
		mpz_add(ans,ans,temp);
		mpz_fdiv_r(ans,ans,M);
	}
	ret=1;
	mpz_set(root,ans);
	mpz_fdiv_r(root,root,n);
	mpz_mul(temp,root,root);
	mpz_fdiv_r(temp,temp,n);
	gmp_printf("root %Zd root^2 %Zd\n",root,temp);
quit:	
	free(q);
	mpz_clear(ans);
	mpz_clear(temp);
	mpz_clear(M);
	mpz_clear(P);
	free(ev);
	return ret;
}

/* use trial division to check that a-bm (rational) and a-b*alpha (algebraic)
   are smooth with regard to our factor base.	return 1 if smooth and also
	 return the indexes of the factors in *f1,*f2,*f3. also set f0 to 1 if
	 a-bm is negative. f3 will contain list of indexes where legendre
	 symbol=-1. */
int trialsmooth(mpz_t a,mpz_t b,mpz_t *f,int deg,mpz_t m,int *f0,ull *f1,int *fn1,ull *f2,int *fn2,ull *f3,int *fn3) {
	mpz_t rat,alg,t,u,div;
	ull i,j,r,A,B;
	int ret=0;
	mpz_init(t);
	mpz_init(u);
	mpz_set(t,a);
	mpz_set(u,b);
	mpz_abs(t,t);
	mpz_abs(u,u);
	mpz_gcd(t,t,u);
	if(mpz_cmp_si(t,1)) goto cleanupgcd;
	mpz_init(div);
	/* rat = a-bm */
	mpz_init(rat);
	mpz_mul(rat,b,m);
	mpz_sub(rat,a,rat);
	// gmp_printf("\nrat %Zd: ",rat);
	/* check for negative a-bm */
	if(mpz_cmp_si(rat,0)<0) *f0=1,mpz_abs(rat,rat);
	else *f0=0;
	*fn2=0;
	/* trial division on a-bm */
	for(i=0;i<bn2;i++) {
		/* break if p2[i]^2 > rat */
		mpz_set_ull(div,p2[i]);
		mpz_mul(t,div,div);
		if(mpz_cmp(t,rat)>0) break;
		/* factor out div from rat and keep count */
		mpz_fdiv_qr(t,u,rat,div);
		if(mpz_cmp_si(u,0)) continue;
		mpz_set(rat,t);
		f2[(*fn2)++]=i;
//		printf("%I64d ",p2[i]);
		while(1) {
			mpz_fdiv_qr(t,u,rat,div);
			if(mpz_cmp_si(u,0)) break;
			mpz_set(rat,t);
			f2[(*fn2)++]=i;
//			printf("%I64d ",p2[i]);
		}
	}
	/* if remainder of rat > largest prime in factor base, number isn't smooth */
	mpz_set_ull(div,p2[bn2-1]);
	if(mpz_cmp(div,rat)<0) goto cleanuprat;
	if(mpz_cmp_si(rat,1)>0) {
		/* add remainder to primes */
		f2[(*fn2)++]=bs(p2,bn2,mpz_get_ull(rat));
		if(mpz_get_ull(rat)!=p2[bs(p2,bn2,mpz_get_ull(rat))]) error("sanity test failed, rational remainder is not equal to prime found");
//		printf("%I64d ",mpz_get_ull(rat));
	}
	/* alg = norm(a-b*alpha) */
	mpz_init(alg);
	calcnorm(alg,a,b,f,deg);
	mpz_abs(alg,alg);
	// gmp_printf("\nalg %Zd: \n",alg);
	*fn1=0;
	/* trial division on norm(a-b*alpha) */
	for(i=0;i<bn1;i++)  {
		/* break if p1[i]^2 > alg */
		mpz_set_ull(div,p1[i]);
		mpz_mul(t,div,div);
		if(mpz_cmp(t,alg)>0) break;
		/* check if p1[i] divides alg */
		mpz_fdiv_r(t,alg,div);
		if(mpz_cmp_si(t,0)) continue;
		/* if a-br=0 mod p this is the prime we want */
		mpz_set_ull(t,r1[i]);
		mpz_mul(t,b,t);
		mpz_sub(t,a,t);
		mpz_fdiv_r(t,t,div);
		if(mpz_cmp_si(t,0)) continue;
		mpz_fdiv_q(alg,alg,div);
//		gmp_printf("[%I64d %I64d] ",p1[i],r1[i]);
		/* factor out div from alg and keep count */
		f1[(*fn1)++]=i;
		while(1) {
			mpz_fdiv_qr(t,u,alg,div);
			if(mpz_cmp_si(u,0)) break;
			mpz_set(alg,t);
//			gmp_printf("[%I64d %I64d] ",p1[i],r1[i]);
			f1[(*fn1)++]=i;
		}
	}
	/* check if alg>largest prime in factor base */
	mpz_set_ull(div,p1[bn1-1]);
//	gmp_printf("rest %Zd vs %Zd\n",alg,div);
	if(mpz_cmp(div,alg)<0) goto cleanupalg;
	if(mpz_cmp_si(alg,1)>0) {
		/* add reminder to primes */
		/* find index of first eligible pair (p,r) */
		i=bs(p1,bn1,mpz_get_ull(alg));
		if(i==-1) {
			gmp_printf("a = %Zd, b = %Zd\n",a,b);
			printf("tried to find %I64d, not in factor base\n",mpz_get_ull(alg));
			r=mpz_get_ull(alg);
			for(i=0;i<bn1;i++) if(p1[i]>r-1000 && p1[i]<r+1000) printf("[%I64d %I64d] ",p1[i],r1[i]);
			error("\n");
		}
		/* find r such that a-br=0 (mod p) which is a*inverse(b) mod p */
		mpz_fdiv_r(u,a,alg);
		A=mpz_get_ull(u);
		mpz_fdiv_r(u,b,alg);
		B=mpz_get_ull(u);
		r=ullmulmod2(inverse(B,p1[i]),A,p1[i]);
		for(j=i;j<bn1;j++) {
			if(p1[j]!=p1[i]) break;
			if(r1[j]==r) goto ok;
		}
		error("(p,r) not found, shouldn't happen!");
	ok:
		f1[(*fn1)++]=j;
//		printf("[%I64d %I64d]\n",p1[i],r1[i]);
	}
	/* we won, (a,b) is smooth. now get the quadratic characters */
	*fn3=0;
	for(i=0;i<bn3;i++) {
		/* if legendre(a-br/p)==-1, then add this (p,r) */
		mpz_set_ull(t,p3[i]);
		mpz_set_ull(u,r3[i]);
		mpz_mul(u,b,u);
		mpz_sub(u,a,u);
		if(mpz_legendre(u,t)<0) f3[(*fn3)++]=i;
	}
	ret=1;
cleanupalg:
	mpz_clear(alg);
cleanuprat:
	mpz_clear(rat);
	mpz_clear(div);
cleanupgcd:
	mpz_clear(u);
	mpz_clear(t);
	return ret;
}

/* sieve from a1,b to a2,b, inclusive. restriction: a1 and a2 are int */
/* return 1 whenever enough relations are found */
int linesieve(int a1,int a2,int b,mpz_t n,mpz_t *f,int fn,mpz_t m,int extra,int *aval,int *bval) {
	int *sieve;
	mpz_t rat,norm,A,B,t,u;
	double invl2=1./log(2);
	ull j,z;
	int size=a2-a1+1,i,a,v,lgp,ret=0;
	int flog[MAXDEG+1];
	int blog[MAXDEG+1];
	double temp;
	if(!(sieve=malloc(size*sizeof(int)))) error("out of memory in line sieve");
	mpz_init(rat);
	mpz_init(norm);
	mpz_init(t);
	mpz_init(u);
	mpz_init(A);
	mpz_init(B);
	/* initialize rat=a-bm */
	mpz_set_si(B,b);
	mpz_mul(t,B,m);
	mpz_set_si(A,a1);
	mpz_sub(rat,A,t);
	mpz_set(t,rat);
	/* precalculate values for fast log_2(norm) */
	for(i=0;i<=fn;i++) flog[i]=mpz_sizeinbase(f[i],2);
	for(temp=0,i=0;i<=fn;i++,temp+=log(temp)*invl2) blog[i]=(int)(0.5+temp);
	for(a=a1,i=0;i<size;i++) {
		calcnorm(norm,A,B,f,fn);
		/* store lg norm + lg rat in sieve */
		v=mpz_sizeinbase(t,2)+mpz_sizeinbase(norm,2);
		/* fast version! approximate log2(norm) faster than calculating
		   the full norm every time */
		/* TODO */
		/* v+=mpz_sizeinbase(t,2); */
		sieve[i]=v;
		mpz_add_ui(t,t,1);
		mpz_add_ui(A,A,1);
	}
	/* process each rational prime */
	for(j=opt_skip;j<bn2;j++) {
		lgp=.5+log(p2[j])*invl2;
		/* find starting point: first smallest i>=0 such that a+i-bm=0 mod p */
		z=mpz_mod_ull(rat,p2[j]);
		/* subtract lg(prime) for each eligible element in sieve */
//		gmp_printf("start for %d (%Zd, mod=%d): %d\n",(int)p2[j],rat,(int)z,(int)(z?p2[j]-z:0));
		for(i=z?p2[j]-z:0;i<size;i+=p2[j]) sieve[i]-=lgp;
	}
	/* process each algebraic prime */
	mpz_set_si(A,a1);
	for(j=opt_skip;j<bn1;j++) {
		lgp=.5+log(p1[j])*invl2;
//		printf("process %I64d, subtract %d\n",p1[j],lgp);
		/* find starting point: find smallest i>=0 such that a+i-br=0 mod p */
		mpz_set_ull(t,r1[j]);
		mpz_mul(t,t,B);
		mpz_sub(t,A,t);
		z=mpz_mod_ull(t,p1[j]);
		for(i=z?p1[j]-z:0;i<size;i+=p1[j]) sieve[i]-=lgp;
	}
	/* find candidates for smooth numbers by taking the ones with small
	   remaining log values. only taking 0-values is too strict, since
	   sieve doesn't subtract powers of primes, and all logs are rounded
	   to int */
	for(i=0;i<size;i++) {
		/* WARNING, magic constants */
		static ull f1[100000],f2[100000],f3[100000];
		int fn1=0,fn2=0,fn3=0,f0;
		a=a1+i;
		if(a==0) continue;
		if(gcd(a>0?a:-a,b>0?b:-b)>1) continue;
		if(sieve[i]<=opt_thr) {
			mpz_add_ui(t,A,i);
//			printf("check %d %d\n",a1+i,b);
			if(trialsmooth(t,B,f,fn,m,&f0,f1,&fn1,f2,&fn2,f3,&fn3)) {
				/* insert in transposed matrix:
				   column i is the ith relation we find
				   row corresponds to -1, prime or quadratic character */
				if(f0) MSETBIT(M,0,smooth);
				for(j=0;j<fn1;j++) MTOGGLEBIT(M,1+f1[j],smooth);
				for(j=0;j<fn2;j++) MTOGGLEBIT(M,1+bn1+f2[j],smooth);
				for(j=0;j<fn3;j++) MSETBIT(M,1+bn1+bn2+f3[j],smooth);
				/* store the rational divisors */
				faclist[smooth]=malloc(fn2*sizeof(ull));
				if(!faclist[smooth]) error("out of memory trialsmooth");
				alglist[smooth]=malloc(fn1*sizeof(ull));
				if(!alglist[smooth]) error("out of memory trialsmooth");
				memcpy(faclist[smooth],f2,sizeof(ull)*fn2);
				facn[smooth]=fn2;
				memcpy(alglist[smooth],f1,sizeof(ull)*fn1);
				algn[smooth]=fn1;
				/* store the actual a,b pair */
				aval[smooth]=a1+i;
				bval[smooth]=b;
				smooth++;
				if(smooth%100==0) {
//					printf("(%d, %d)\n",a1+i,b);
					printf("%d/%I64d found: (%d, %d) is smooth, log %d\n",smooth,extra+1+bn1+bn2+bn3,a1+i,b,sieve[i]);
				}
				if(smooth==extra+1+bn1+bn2+bn3) {
					puts("==> enough relations gathered!");
					ret=1;
					goto end;
				}
			} else notsmooth++;
		} else {
			/* remove the continue if you want to benchmark 
			   smooth numbers not found by the sieving */
			continue;
			mpz_add_ui(t,A,i);
			if(trialsmooth(t,B,f,fn,m,&f0,f1,&fn1,f2,&fn2,f3,&fn3)) {
				printf("%d - %d*alpha is smooth, log %d MISSED\n",a1+i,b,sieve[i]);
				missed++;
			}
		}
	}
end:
	mpz_clear(B);
	mpz_clear(A);
	mpz_clear(u);
	mpz_clear(t);
	mpz_clear(norm);
	mpz_clear(rat);
	free(sieve);
	return ret;
}

void testsieve(mpz_t n,mpz_t *f,int fn,mpz_t m,int extra,int *aval,int *bval) {
	int B;
	/* factor lists */
	puts("start sieve");
	notsmooth=missed=smooth=0;
	faclist=malloc((1+bn1+bn2+bn3+extra)*sizeof(ull*));
	if(!faclist) error("out of memory");
	facn=malloc((1+bn1+bn2+bn3+extra)*sizeof(int));
	if(!facn) error("out of memory");
	alglist=malloc((1+bn1+bn2+bn3+extra)*sizeof(ull*));
	if(!alglist) error("out of memory");
	algn=malloc((1+bn1+bn2+bn3+extra)*sizeof(int));
	if(!algn) error("out of memory");
	for(B=1;;B++) if(linesieve(-opt_sievew,opt_sievew,-1*opt_signb*B,n,f,fn,m,extra,aval,bval)) break;
	printf("smooth numbers found: %d\n",smooth);
	printf("nonsmooth numbers trial-divided: %d\n",notsmooth);
	printf("smooth numbers missed: %d\n",missed);
	puts("end sievetest");
}

void takegcd(mpz_t ans,mpz_t alg,mpz_t rat,mpz_t n) {
	mpz_t sub;
	mpz_init(sub);
	mpz_sub(sub,alg,rat);
	mpz_gcd(ans,n,sub);
	mpz_clear(sub);
}

/* takes a number n and returns a factor p, if found
   return values:
   1: factor found
   0: factor not found
   -1: n is even
   -2: n is a perfect power
   -3: n is probably prime
   -4: mysterious error */
int donfs(mpz_t n) {
	mpz_t m,f[MAXDEG+1],r,temp;
	mpz_t ratrot,algrot;
	ull Br=opt_Br,Ba=opt_Ba,rows,k;
	int *aval,*bval;
	int deg=opt_deg;
	int err,retval=0,i,Bk=opt_Bq,j;
	int extra=opt_extra,zero;
	uchar *v;
	uchar *freevar;
	mpz_init(r); mpz_init(m); mpz_init(temp);
	mpz_init(ratrot); mpz_init(algrot);
	for(i=0;i<=MAXDEG;i++) mpz_init(f[i]);
	/* check prerequisites: n cannot be even, prime or perfect power */
	/* (if n is a perfect power, try running nfs again on the root */
	mpz_fdiv_r_ui(r,n,2);
	if(!mpz_cmp_si(r,0)) { retval=-1; goto end; }
	if(mpz_perfect_power_p(n)) { retval=-2; goto end; }
	/* we want to be really, REALLY sure that n is composite */
	if(mpz_probab_prime_p(n,100)) { retval=-3; goto end; }
	mpz_set(m,opt_m);
	if(!mpz_cmp_si(m,0)) mpz_root(m,n,deg); /* deg-th root of n, get our base m */
	gmp_printf("m = %Zd\n",m);
	err=getpolynomial(n,m,deg,f);
	if(!err) error("polynomial isn't monic or is otherwise wrong");
	printmpzpoly(f,deg);
	/* for now, only try to find linear factors when
	   the a_0 coefficient is small enough */
	/* TODO replace with better way to find all linear factors.
	   fully factorize f[0] (possibly by pollard rho or even qs) and
	   generate all divisors by generating all exponent tuples. in this way,
	   the program will have full degree 3 support */
	/* TODO move this to a function */
	if(mpz_cmp_si(f[0],2000000000)<0) {
		if(!mpz_cmp_si(f[0],0)) {
			printmpzpoly(f,deg);
			gmp_printf("f(x) factored, found factor %Zd\n",m);
			retval=1;
			goto end;
		}
		j=mpz_get_si(f[0]);
		for(i=1;i*i<=j;i++) if(j%i==0) {
			if(i>1) {
				mpz_set_si(r,-i);
				evalpoly(f,deg,r,temp);
				if(!mpz_cmp_si(temp,0)) {
					printmpzpoly(f,deg);
					gmp_printf("f(x) factored, found factor %d\n",i);
					retval=1;
					goto end;
				}
			}
			mpz_set_si(r,-j/i);
			evalpoly(f,deg,r,temp);
			if(!mpz_cmp_si(temp,0)) {
				printmpzpoly(f,deg);
				gmp_printf("f(x) factored, found factor %d\n",i);
				retval=1;
				goto end;
			}
		}
	}
	/* TODO try to factorize polynomial properly and terminate early */

	/* factor base */
	if(!Ba) Ba=findB(n)*1;
	if(!Br) Br=findB(n)*1;
	if(!Bk) Bk=findK(n)*0.25; /* number of quadratic characters */
	puts("factor base info:");
	printf("  bound %I64d\n",Ba);
	createfactorbases(Ba,Br,Bk,f,deg,&p1,&r1,&bn1,&p2,&bn2,&p3,&r3,&bn3);
	printf("  %I64d rational primes\n",bn2);
	printf("  %I64d algebraic primes\n",bn1);
	printf("  %d quadratic characters\n",Bk);
	printf("  total size %I64d\n",bn1+bn2+Bk);
	if(!extra) extra=3+(bn1+bn2+bn3+1)/1000;

	puts("continue with factorization, altsaa");
	/* allocate memory for matrix, uncompressed */
	rows=1+bn1+bn2+bn3;
	M=malloc(sizeof(uint *)*rows);
	for(i=0;i<rows;i++) {
		M[i]=calloc(((rows+31+extra)/32),sizeof(uint));
		if(!M[i]) error("out of memory while allocating matrix");
	}
	aval=malloc(sizeof(int)*(rows+extra));
	if(!aval) error("out of memory");
	bval=malloc(sizeof(int)*(rows+extra));
	if(!bval) error("out of memory");

	testsieve(n,f,deg,m,extra,aval,bval);
	puts("start gauss");
//	printf("before solve\n");for(i=0;i<rows;i++) { for(j=0;j<rows+extra;j++) printf("%d",ISSET(M,i,j)?1:0); printf("\n");}
	bitgauss32(M,rows,rows+extra,0);
//	printf("after solve\n");for(i=0;i<rows;i++) { for(j=0;j<rows+extra;j++) printf("%d",ISSET(M,i,j)?1:0); printf("\n");}
	v=malloc(rows+extra);
	if(!v) error("out of memory");
	freevar=malloc(rows+extra);
	if(!freevar) error("out of memory");
	zero=findfreevars(M,rows,rows+extra,freevar);
	printf("gauss done, %d free variables found\n",zero);
	for(k=0;k<zero;k++) {
		puts("-------------------------------------");
		getsquare(M,rows,rows+extra,freevar,k,v);
//		printalgnum(n,v,rows+extra,f,deg,m,aval,bval);
		if(!getalgroot(n,v,rows+extra,f,deg,m,algrot,aval,bval)) continue;
//		printf("takevector:\n");for(i=0;i<rows+extra;i++) printf("%d",v[i]);printf("\n");
		getratroot(n,v,rows+extra,f,deg,m,ratrot,aval,bval);
		gmp_printf("algroot %Zd ratroot %Zd\n",algrot,ratrot);
		takegcd(temp,algrot,ratrot,n);
		/* trivial result, try next linear combination */
		if(!mpz_cmp_si(temp,1) || !mpz_cmp(temp,n)) continue;
		gmp_printf("found factor %Zd after %d tries\n",temp,k+1);
		retval=1;
		break;
	}
	if(!retval) puts("no factor found");
	free(v);
end:
	for(i=0;i<=MAXDEG;i++) mpz_clear(f[i]);
	mpz_clear(ratrot); mpz_clear(algrot);
	mpz_clear(m); mpz_clear(r); mpz_clear(temp);
	return retval;
}

int main() {
	gmp_randinit_mt(gmpseed);
	gmp_randseed_ui(gmpseed,time(0));
	readoptions();
	gmp_printf("try to factor %Zd\n",opt_n);
	printf("return %d\n",donfs(opt_n));
	return 0;
}
