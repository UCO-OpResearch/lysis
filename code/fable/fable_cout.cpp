// fingerprint 202a11e5c11e50b6a1202da76277416c
#include <fem.hpp> // Fortran EMulation library of fable module

namespace example {

using namespace fem::major_types;

using fem::common;

//C cctbx_project/compcomm/newsletter09/sf_times.py, svn rev. 9983
//C cos_wrapper, exp_wrapper replaced with COS, SIN, EXP
//C
void
sf(
  arr_cref<float> abcss,
  int const& n_scatt,
  arr_cref<float, 2> xyz,
  arr_cref<float> b_iso,
  int const& n_refl,
  arr_cref<int, 2> hkl,
  arr_ref<float, 2> f_calc)
{
  abcss(dimension(3));
  xyz(dimension(3, star));
  b_iso(dimension(star));
  hkl(dimension(3, star));
  f_calc(dimension(2, star));
  int i_refl = fem::int0;
  float a = fem::float0;
  float b = fem::float0;
  int i_scatt = fem::int0;
  float phi = fem::float0;
  int j = fem::int0;
  float cphi = fem::float0;
  float sphi = fem::float0;
  float dss = fem::float0;
  int h = fem::int0;
  float ldw = fem::float0;
  float dw = fem::float0;
  FEM_DO_SAFE(i_refl, 1, n_refl) {
    a = 0;
    b = 0;
    FEM_DO_SAFE(i_scatt, 1, n_scatt) {
      phi = 0;
      FEM_DO_SAFE(j, 1, 3) {
        phi += hkl(j, i_refl) * xyz(j, i_scatt);
      }
      phi = phi * 2 * 3.1415926535897931f;
      cphi = fem::cos(phi);
      sphi = fem::sin(phi);
      dss = 0;
      FEM_DO_SAFE(j, 1, 3) {
        h = hkl(j, i_refl);
        dss += h * h * abcss(j);
      }
      ldw = -0.25f * dss * b_iso(i_scatt);
      dw = fem::exp(ldw);
      a += dw * cphi;
      b += dw * sphi;
    }
    f_calc(1, i_refl) = a;
    f_calc(2, i_refl) = b;
  }
}

void
program_run(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  common_write write(cmn);
  arr_1d<3, float> abcss(fem::fill0);
  abcss(1) = 1 / (11.0f * 11.0f);
  abcss(2) = 1 / (12.0f * 12.0f);
  abcss(3) = 1 / (13.0f * 13.0f);
  int jr = 0;
  int i = fem::int0;
  const int n_scatt = 10;
  int j = fem::int0;
  arr_2d<3, n_scatt, float> xyz(fem::fill0);
  FEM_DO_SAFE(i, 1, n_scatt) {
    FEM_DO_SAFE(j, 1, 3) {
      jr = fem::mod(jr * 1366 + 150889, 714025);
      xyz(j, i) = (fem::mod(jr, 20000) - 10000) / 10000.0f;
    }
  }
  arr_1d<n_scatt, float> b_iso(fem::fill0);
  FEM_DO_SAFE(i, 1, n_scatt) {
    jr = fem::mod(jr * 1366 + 150889, 714025);
    b_iso(i) = fem::mod(jr, 10000) / 100.0f;
  }
  if (n_scatt <= 10) {
    FEM_DO_SAFE(i, 1, n_scatt) {
      write(6, "(4(1x,f9.6))"), xyz(1, i), xyz(2, i), xyz(3, i), b_iso(i);
    }
  }
  const int n_refl = 10;
  arr_2d<3, n_refl, int> hkl(fem::fill0);
  FEM_DO_SAFE(i, 1, n_refl) {
    FEM_DO_SAFE(j, 1, 3) {
      jr = fem::mod(jr * 1366 + 150889, 714025);
      hkl(j, i) = fem::mod(jr, 10) - 5;
    }
  }
  arr_2d<2, n_refl, float> f_calc(fem::fill0);
  sf(abcss, n_scatt, xyz, b_iso, n_refl, hkl, f_calc);
  float max_a = fem::float0;
  float max_b = fem::float0;
  float a = fem::float0;
  float b = fem::float0;
  if (n_refl <= 100) {
    FEM_DO_SAFE(i, 1, n_refl) {
      write(6, "(3(1x,i3),1x,f12.6,1x,f12.6)"), hkl(1, i), hkl(2, i),
        hkl(3, i), f_calc(1, i), f_calc(2, i);
    }
  }
  else {
    max_a = 0;
    max_b = 0;
    FEM_DO_SAFE(i, 1, n_refl) {
      a = f_calc(1, i);
      b = f_calc(2, i);
      if (max_a < a) {
        max_a = a;
      }
      if (max_b < b) {
        max_b = b;
      }
    }
    write(6, "(2(1x,f12.6))"), max_a, max_b;
  }
}

} // namespace example

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    example::program_run);
}
