#include <Geex/symbolics/function.h>
namespace Geex {

   class Ltheta20: public Function {
      public:
      Ltheta20();
      virtual void eval(bool do_f, bool do_g, bool do_H) ;
   };
}
