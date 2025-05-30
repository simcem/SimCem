#include <memory>

#include <sundials/sundials_config.h>
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>

namespace sundials {
  namespace detail {
    class SundialsContext {
    public:
      static const SUNContext& getCtx() {
	static SundialsContext ctx;
	return ctx.sunctx;
      }

    private:
      SundialsContext() {	
	SUNContext_Create(SUN_COMM_NULL, &sunctx);
      }

      SUNContext sunctx;
    };
    
    int check_flag(void *flagvalue, const char *funcname, int opt)
    {
      int *errflag;
	
      if (opt == 0 && flagvalue == NULL) {
	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	fprintf(stderr, 
		"\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
		funcname);
	return(1);
      } else if (opt == 1) {
	/* Check if flag < 0 */
	errflag = (int *) flagvalue;
	if (*errflag < 0) {
	  fprintf(stderr, 
		  "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
		  funcname, *errflag);
	  return(1); 
	}
      } else if (opt == 2 && flagvalue == NULL) {
	/* Check if function returned NULL pointer - no memory allocated */
	fprintf(stderr, 
		"\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
		funcname);
	return(1); 
      }

      return(0);
    }    
  }

  struct Serial_Vector {
    Serial_Vector(const long int N) {
      _vec = N_VNew_Serial(N, detail::SundialsContext::getCtx());
      if (detail::check_flag((void *)_vec, "N_VNew_Serial", 0))
	stator_throw() << "Exit!";
      _owner = true;
    }

    Serial_Vector(const N_Vector N) {
      _owner = false;
      _vec = N;
    }

    ~Serial_Vector() {
      if (_owner)
	N_VDestroy_Serial(_vec);
    }
    
    sunrealtype& operator[](const size_t i) {
      return N_VGetArrayPointer_Serial(_vec)[i];
    }

    const sunrealtype& operator[](const size_t i) const {
      return N_VGetArrayPointer_Serial(_vec)[i];
    }
    
    operator N_Vector() { return _vec; }

    long int size() const { return NV_LENGTH_S(_vec); }
    
    N_Vector _vec;
    bool _owner;
  };

  struct IDA {
    SUNMatrix A;
    SUNLinearSolver LS;

    IDA(sunrealtype t0, IDAResFn resrob, N_Vector y0, N_Vector yprime0)
    {
      std::sort(_stop_points.begin(), _stop_points.end());
      
      _mem = IDACreate(detail::SundialsContext::getCtx());
      if(sundials::detail::check_flag((void *)_mem, "IDACreate", 0))
	stator_throw() << "Exit!";

      if (NV_LENGTH_S(y0) != NV_LENGTH_S(yprime0))
	stator_throw() << "Mismatched dimensions of y and y'";

      _N = NV_LENGTH_S(y0);
      
      int retval = IDAInit(_mem, resrob, t0, y0, yprime0);
      if(detail::check_flag(&retval, "IDAInit", 1))
	stator_throw() << "Exit!";
    }

    void denseSolver(N_Vector y0) {
      A = SUNDenseMatrix(_N, _N, detail::SundialsContext::getCtx());
      if (sundials::detail::check_flag((void *)A, "SUNDenseMatrix", 0))
	stator_throw() << "Exit!";
 
      LS = SUNLinSol_Dense(y0, A, detail::SundialsContext::getCtx());
      if(sundials::detail::check_flag((void *)LS, "SUNDenseLinearSolver", 0)) 
	stator_throw() << "Exit!";

      int retval = IDASetLinearSolver(_mem, LS, A);
      if(sundials::detail::check_flag(&retval, "IDADlsSetLinearSolver", 1))
	stator_throw() << "Exit!";
    }

    ~IDA() {
      IDAFree(&_mem);
      SUNLinSolFree(LS);
      SUNMatDestroy(A);
    }

    void reinit(sunrealtype t0, N_Vector y0, N_Vector yprime0) {
      _stop_points.clear();
      int retval = IDAReInit(_mem, t0, y0, yprime0);
      if(detail::check_flag(&retval, "IDAReInit", 1))
	stator_throw() << "Exit!";
    }

    void setTol(sunrealtype rel, N_Vector abs) {
      int retval = IDASVtolerances(_mem, rel, abs);
      if (detail::check_flag(&retval, "IDASVtolerances", 1))
	stator_throw() << "Exit!";
    }

    void setTol(sunrealtype rel, sunrealtype abs) {
      int retval = IDASStolerances(_mem, rel, abs);
      if (detail::check_flag(&retval, "IDASStolerances", 1))
	stator_throw() << "Exit!";
    }
    
    template<class T>
    void setUserData(T& data) {
      int retval = IDASetUserData(_mem, (void*)(&data));
      if (detail::check_flag(&retval, "IDASetUserData", 1))
	stator_throw() << "Exit!";
    }
    
    void setId(N_Vector Id) {
      int retval = IDASetId(_mem, Id);
      if (detail::check_flag(&retval, "IDASetId", 1))
	stator_throw() << "Exit!";
      IDASetId(_mem, Id);
    }

    void calcIC(int icopt, sunrealtype tout1) {
      int retval = IDACalcIC(_mem, icopt, tout1);
      if (detail::check_flag(&retval, "IDACalcIC", 1))
	stator_throw() << "Exit!";
    }

    void add_stop_point(sunrealtype tstop) {
      _stop_points.push_back(tstop);
    }
    
    sunrealtype solve(sunrealtype tfinal, N_Vector y, N_Vector yp, int mode) {
      sunrealtype t;
      int retval = IDASolve(_mem, tfinal, &t, y, yp, mode);
      if(detail::check_flag(&retval, "IDASolve", 1))
	stator_throw() << "Exit!";      
      return t;
    }
    
    void interpolate(sunrealtype tinterp, N_Vector y, long int order = 0) {
      int retval = IDAGetDky(_mem, tinterp, order, y);
      if(detail::check_flag(&retval, "IDAGetDky", 1))
	stator_throw() << "Exit!";
    }
  

    
    void* getMem() const {
      return _mem;
    }
    
    void* _mem;
    long int _N;

    std::vector<sunrealtype> _stop_points;
  };
  
}
