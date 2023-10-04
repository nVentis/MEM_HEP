#include "array"
#include "cmath"

using vec3 = std::array<double,3>;

std::array<double, 2> get_angles(double x, double y, double z){
  if (x == 0. && y == 0. && z == 0.) {
    // to handle use_transer_funcs=false case
    return std::array<double, 2> {
      0,
      0
    };
  } else {
    return std::array<double, 2> {
      std::acos(z/std::sqrt(std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.))),
      std::atan2(y, x)
    };
  }
}

std::array<std::array<double, 2>, 4> get_reco_angles(double reco_kin[]) {
  std::array<std::array<double, 2>, 4> res;
  // reco_kin[0-3]: mu-
  // reco_kin[4-7]: mu+

  int i;
  for (i = 0; i < 4; i++)
    res[i] = get_angles(reco_kin[9+(i*4)], reco_kin[10+(i*4)], reco_kin[11+(i*4)]);

  return res;
}

const vec3 vec3_sph(double theta, double phi, double len) {
  return vec3 {
    len*std::sin(theta)*std::cos(phi),
    len*std::sin(theta)*std::sin(phi),
    len*std::cos(theta)
  };
};

vec3 vec3_add(vec3 a, vec3 b) {
  return vec3{
    a[0] + b[0],
    a[1] + b[1],
    a[2] + b[2]
  };
};

const double vec3_dot(vec3 a, vec3 b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}