#include <ma.h>

#ifndef AW_EXPLF_H_
#define AW_EXPLF_H_

class ExpLengthFunc {
  protected:
  enum class Axis { X, Y, Z } axis;
  double midpoint, h0; // h0 = longest_edge/20
  ma::Mesh *m;

  double ldist (const ma::Vector &v) const;
  double eval (const ma::Vector &v) const;

  public:
  ExpLengthFunc(ma::Mesh *mesh);
};

#endif // AW_EXPLF_H_

