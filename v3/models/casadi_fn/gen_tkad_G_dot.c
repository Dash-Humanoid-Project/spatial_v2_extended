/* This file was automatically generated by CasADi 3.6.2.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) gen_tkad_G_dot_ ## ID
#endif

#include <math.h>
#include <string.h>
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_fill CASADI_PREFIX(fill)
#define casadi_from_mex CASADI_PREFIX(from_mex)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_to_mex CASADI_PREFIX(to_mex)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

void casadi_fill(casadi_real* x, casadi_int n, casadi_real alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

#ifdef MATLAB_MEX_FILE
casadi_real* casadi_from_mex(const mxArray* p, casadi_real* y, const casadi_int* sp, casadi_real* w) {
  casadi_int nrow, ncol, is_sparse, c, k, p_nrow, p_ncol;
  const casadi_int *colind, *row;
  mwIndex *Jc, *Ir;
  const double* p_data;
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError",
      "\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
  p_nrow = mxGetM(p);
  p_ncol = mxGetN(p);
  is_sparse = mxIsSparse(p);
  Jc = 0;
  Ir = 0;
  if (is_sparse) {
    Jc = mxGetJc(p);
    Ir = mxGetIr(p);
  }
  p_data = (const double*)mxGetData(p);
  if (p_nrow==1 && p_ncol==1) {
    casadi_int nnz;
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    nnz = sp[ncol];
    casadi_fill(y, nnz, v);
  } else {
    casadi_int tr = 0;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError",
                                 "\"from_mex\" failed: Dimension mismatch. "
                                 "Expected %d-by-%d, got %d-by-%d instead.",
                                 nrow, ncol, p_nrow, p_ncol);
    }
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=colind[c]; k<colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}

#endif

#define casadi_to_double(x) ((double) x)

#ifdef MATLAB_MEX_FILE
mxArray* casadi_to_mex(const casadi_int* sp, const casadi_real* x) {
  casadi_int nrow, ncol, c, k;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int nnz;
#endif
  const casadi_int *colind, *row;
  mxArray *p;
  double *d;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int i;
  mwIndex *j;
#endif /* CASADI_MEX_NO_SPARSE */
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
#ifndef CASADI_MEX_NO_SPARSE
  nnz = sp[ncol];
  if (nnz!=nrow*ncol) {
    p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *colind++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *row++;
    if (x) {
      d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = casadi_to_double(*x++);
    }
    return p;
  }
#endif /* CASADI_MEX_NO_SPARSE */
  p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  if (x) {
    d = (double*)mxGetData(p);
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        d[row[k]+c*nrow] = casadi_to_double(*x++);
      }
    }
  }
  return p;
}

#endif

static const casadi_int casadi_s0[6] = {2, 1, 0, 2, 0, 1};
static const casadi_int casadi_s1[13] = {4, 2, 0, 4, 8, 0, 1, 2, 3, 0, 1, 2, 3};

/* tkad_G_dot:(i0[2],i1[2],i2[2],i3[2])->(o0[4x2]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a3, a4, a5, a6, a7, a8, a9;
  a0=0.;
  if (res[0]!=0) res[0][0]=a0;
  if (res[0]!=0) res[0][1]=a0;
  if (res[0]!=0) res[0][2]=a0;
  a1=6.;
  a2=1.6000000000000001e-03;
  a3=13.;
  a4=8.3333333333333329e-02;
  a5=arg[0]? arg[0][0] : 0;
  a5=(a5/a1);
  a6=2.;
  a7=(a5/a6);
  a8=arg[0]? arg[0][1] : 0;
  a8=(a8/a1);
  a9=(a8/a6);
  a7=(a7-a9);
  a9=1.0325367854798453e+00;
  a7=(a7+a9);
  a9=cos(a7);
  a10=5.0000000000000000e-01;
  a11=1.6666666666666666e-01;
  a12=arg[2]? arg[2][0] : 0;
  a12=(a11*a12);
  a13=(a10*a12);
  a14=arg[2]? arg[2][1] : 0;
  a11=(a11*a14);
  a14=(a10*a11);
  a13=(a13-a14);
  a9=(a9*a13);
  a13=(a4*a9);
  a13=(a3*a13);
  a13=(a2*a13);
  a14=1.6000000000000001e-04;
  a15=21.;
  a16=(a5/a6);
  a17=(a8/a6);
  a16=(a16-a17);
  a17=1.3816026358787112e+00;
  a16=(a16+a17);
  a17=cos(a16);
  a18=(a10*a12);
  a19=(a10*a11);
  a18=(a18-a19);
  a17=(a17*a18);
  a18=(a4*a17);
  a18=(a15*a18);
  a18=(a14*a18);
  a13=(a13-a18);
  a18=4.0000000000000002e-04;
  a19=7.;
  a5=(a5/a6);
  a8=(a8/a6);
  a5=(a5-a8);
  a8=arg[1]? arg[1][1] : 0;
  a5=(a5+a8);
  a6=1.4514158059584845e+00;
  a5=(a5+a6);
  a6=sin(a5);
  a12=(a10*a12);
  a10=(a10*a11);
  a12=(a12-a10);
  a6=(a6*a12);
  a10=(a4*a6);
  a10=(a19*a10);
  a10=(a18*a10);
  a13=(a13+a10);
  a10=2.0000000000000001e-04;
  a11=91.;
  a20=4.1887902047863906e-01;
  a20=(a8+a20);
  a21=cos(a20);
  a21=(a11*a21);
  a21=(a10*a21);
  a22=cos(a5);
  a22=(a19*a22);
  a22=(a18*a22);
  a21=(a21-a22);
  a22=2.0000000000000002e-05;
  a23=147.;
  a24=6.9813170079773182e-02;
  a8=(a8+a24);
  a24=cos(a8);
  a24=(a23*a24);
  a24=(a22*a24);
  a21=(a21-a24);
  a13=(a13/a21);
  a7=sin(a7);
  a24=(a4*a7);
  a24=(a3*a24);
  a24=(a2*a24);
  a16=sin(a16);
  a25=(a4*a16);
  a25=(a15*a25);
  a25=(a14*a25);
  a24=(a24-a25);
  a25=cos(a5);
  a26=(a4*a25);
  a26=(a19*a26);
  a26=(a18*a26);
  a24=(a24-a26);
  a24=(a24/a21);
  a26=(a24/a21);
  a27=sin(a5);
  a27=(a27*a12);
  a27=(a19*a27);
  a27=(a18*a27);
  a26=(a26*a27);
  a13=(a13-a26);
  a13=(a1*a13);
  a26=sin(a5);
  a12=arg[3]? arg[3][1] : 0;
  a26=(a26*a12);
  a4=(a4*a26);
  a4=(a19*a4);
  a4=(a18*a4);
  a4=(a4/a21);
  a24=(a24/a21);
  a5=sin(a5);
  a5=(a5*a12);
  a5=(a19*a5);
  a5=(a18*a5);
  a20=sin(a20);
  a20=(a20*a12);
  a11=(a11*a20);
  a10=(a10*a11);
  a5=(a5-a10);
  a8=sin(a8);
  a8=(a8*a12);
  a23=(a23*a8);
  a22=(a22*a23);
  a5=(a5+a22);
  a24=(a24*a5);
  a4=(a4-a24);
  a4=(a1*a4);
  a13=(a13+a4);
  a13=(-a13);
  if (res[0]!=0) res[0][3]=a13;
  if (res[0]!=0) res[0][4]=a0;
  if (res[0]!=0) res[0][5]=a0;
  if (res[0]!=0) res[0][6]=a0;
  a0=-8.3333333333333329e-02;
  a9=(a0*a9);
  a9=(a3*a9);
  a9=(a2*a9);
  a17=(a0*a17);
  a17=(a15*a17);
  a17=(a14*a17);
  a9=(a9-a17);
  a6=(a0*a6);
  a6=(a19*a6);
  a6=(a18*a6);
  a9=(a9+a6);
  a9=(a9/a21);
  a7=(a0*a7);
  a3=(a3*a7);
  a2=(a2*a3);
  a16=(a0*a16);
  a15=(a15*a16);
  a14=(a14*a15);
  a2=(a2-a14);
  a25=(a0*a25);
  a25=(a19*a25);
  a25=(a18*a25);
  a2=(a2-a25);
  a2=(a2/a21);
  a25=(a2/a21);
  a25=(a25*a27);
  a9=(a9-a25);
  a9=(a1*a9);
  a0=(a0*a26);
  a19=(a19*a0);
  a18=(a18*a19);
  a18=(a18/a21);
  a2=(a2/a21);
  a2=(a2*a5);
  a18=(a18-a2);
  a1=(a1*a18);
  a9=(a9+a1);
  a9=(-a9);
  if (res[0]!=0) res[0][7]=a9;
  return 0;
}

CASADI_SYMBOL_EXPORT int tkad_G_dot(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int tkad_G_dot_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int tkad_G_dot_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tkad_G_dot_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int tkad_G_dot_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void tkad_G_dot_release(int mem) {
}

CASADI_SYMBOL_EXPORT void tkad_G_dot_incref(void) {
}

CASADI_SYMBOL_EXPORT void tkad_G_dot_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int tkad_G_dot_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_int tkad_G_dot_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real tkad_G_dot_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tkad_G_dot_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* tkad_G_dot_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tkad_G_dot_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s0;
    case 3: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* tkad_G_dot_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int tkad_G_dot_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_tkad_G_dot(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  int mem;
  casadi_real w[44];
  casadi_int *iw = 0;
  const casadi_real* arg[4] = {0};
  casadi_real* res[1] = {0};
  if (argc>4) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"tkad_G_dot\" failed. Too many input arguments (%d, max 4)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"tkad_G_dot\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s0, w+16);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+2, casadi_s0, w+16);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+4, casadi_s0, w+16);
  if (--argc>=0) arg[3] = casadi_from_mex(argv[3], w+6, casadi_s0, w+16);
  --resc;
  res[0] = w+8;
  tkad_G_dot_incref();
  mem = tkad_G_dot_checkout();
  i = tkad_G_dot(arg, res, iw, w+16, mem);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"tkad_G_dot\" failed.");
  tkad_G_dot_release(mem);
  tkad_G_dot_decref();
  if (res[0]) resv[0] = casadi_to_mex(casadi_s1, res[0]);
}
#endif


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[11];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_tkad_G_dot(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "tkad_G_dot")==0) {
    mex_tkad_G_dot(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'tkad_G_dot'");
}
#endif
#ifdef __cplusplus
} /* extern "C" */
#endif