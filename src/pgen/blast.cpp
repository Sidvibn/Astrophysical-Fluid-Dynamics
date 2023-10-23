//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
// #include"cnpy.h"

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

using namespace std;
// int file(){
//   string density;
//   vector<float>den_val;
//   string fname;
//   fname = "/home/sid/athena/vis/python/Newton_Raphson_numerical.txt";
//   int i = 0;
//   ifstream coeff(fname);
//   if (coeff.is_open()){
//     string line;
//     getline(coeff,line);
//     while (!coeff.eof()){
//       getline(coeff,density,'\n');
//       den_val.push_back(stof(density));
//       i += 1;
//     }
//     coeff.close();
//     cout << "Number of entries: " << i-1 << endl;
//   }
//   else cout << "Unable to open file";
//   return 0;
// }

// int file()
// {
//   string fname;
//   fname = "/home/sid/athena/vis/python/Newton_Raphson_numerical.csv";
//   vector<vector<string>> content;
//   vector<string> row;
//   string line, word;
//   fstream file (fname, ios::in);
//   if(file.is_open())
//   {
//     while(getline(file, line))
//     {
//       row.clear();
//       stringstream str(line);
//       while(getline(str, word, ','))
//         row.push_back(word);
//       content.push_back(row);
//     }
//   }
//   else
//     cout<<"Could not open the file\n";
//   for(int i=0;i<content.size();i++)
//   {
//     for(int j=0;j<content[i].size();j++)
//     {
//       cout<<content[i][j]<<" ";
//     }
//     cout<<"\n";
//   }
//   return 0;
// }

Real threshold;
// cnpy::NpyArray arr = cnpy::npy_load("~/athena/vis/python/arr_analytical.npy");

void MyBoundary_inner(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void MyBoundary_outer(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

int RefinementCondition(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
    threshold = pin->GetReal("problem","thr");
  }
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, MyBoundary_inner);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, MyBoundary_outer);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rout = pin->GetReal("problem", "radius");
  Real rin  = rout - pin->GetOrAddReal("problem", "ramp", 0.0);
  Real pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  Real da   = pin->GetOrAddReal("problem", "damb", 1.0);
  Real prat = pin->GetReal("problem", "prat");
  Real drat = pin->GetOrAddReal("problem", "drat", 1.0);
  Real b0, angle;
  if (MAGNETIC_FIELDS_ENABLED) {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI/180.0)*pin->GetReal("problem", "angle");
  }
  // Real gamma = peos->GetGamma();
  // Real gm1 = gamma - 1.0;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  string fname;
  fname = "/home/sid/athena/vis/python/iso.csv";
  vector<vector<string>> content;
  vector<string> row;
  string line, word;
  fstream file (fname, ios::in);
  if(file.is_open())
  {
    while(getline(file, line))
    {
      row.clear();
      stringstream str(line);
      while(getline(str, word, ','))
        row.push_back(word);
      content.push_back(row);
    }
  }
  else
    cout<<"Could not open the file\n";

  // for(int i=0;i<content.size();i++)
  // {
  //   for(int j=0;j<content[i].size();j++)
  //   {
  //     cout<<content[i][j]<<" ";
  //   }
  //   cout<<"\n";
  // }

  // for(int i=0;i<content.size();i++){
  //   cout<<stod(content[i][0])<<endl;
  // } 

  // Real density[128] = {1.0, 0.9871693909736833, 0.9733976795140618, 0.959893075261259, 
  // 0.9468504722780056, 0.9343025558692525, 0.922243306924561, 0.9106549657302999, 0.8995160912557445, 0.8888044914041097, 
  // 0.8784983362558971, 0.8685766840380443, 0.8590196667153225, 0.8498085693078063, 0.8409257202559285, 0.8323545990789284, 
  // 0.8240796573697433, 0.8160863069162922, 0.8083608825342887, 0.8008905384522612, 0.7936631801569554, 0.7866674998625202, 
  // 0.7798928288938325, 0.7733291068169059, 0.7669668877722642, 0.7607972636608807, 0.7548118119430746, 0.7490026150277986, 
  // 0.743362162907678, 0.737883383359354, 0.7325595671644427, 0.7273843750510526, 0.7223518064690125, 0.717456152810441, 
  // 0.7126920134372352, 0.708054254899285, 0.7035380576569171, 0.6991387568968548, 0.694851960691487, 0.6906735070204089, 
  // 0.6865993846962263, 0.6826258256362943, 0.6787492114648364, 0.6749661168713943, 0.6712732200461977, 0.6676674256036895, 
  // 0.6641457048762273, 0.6607052184699671, 0.6573432368089472, 0.6540571351107541, 0.6508444272534968, 0.6477027092809696, 
  // 0.6446296837009517, 0.6416231632258648, 0.6386810513107097, 0.6358013302494365, 0.6329820468902326, 0.6302213628154943, 
  // 0.6275174953861032, 0.6248687243264305, 0.6222734088563884, 0.6197299751939657, 0.6172368879206027, 0.6147927035048099, 
  // 0.6123960035791889, 0.6100454443940616, 0.607739737719825, 0.6054775849414863, 0.6032578210447919, 0.6010792663901576, 
  // 0.5989407707207999, 0.596841292360886, 0.5947797857398626, 0.5927552162661999, 0.5907666240914967, 0.5888130572083545, 
  // 0.5868936209535329, 0.5850074505668694, 0.5831536864379953, 0.5813314864967164, 0.5795401055523645, 0.577778743078465, 
  // 0.5760466587650793, 0.5743431578574653, 0.572667524622651, 0.5710191191477827, 0.5693972383374124, 0.5678013041530606, 
  // 0.5662306751202254, 0.5646847796311995, 0.5631630339906399, 0.5616648630684186, 0.5601897659166418, 0.5587371816418331, 
  // 0.5573066398920229, 0.5558976115077531, 0.5545096544957786, 0.5531422691905034, 0.5517950227932582, 0.5504674859991062, 
  // 0.5491592248905658, 0.5478698262962383, 0.5465988806066412, 0.5453460058265519, 0.5441108303441391, 0.5428929693638507, 
  // 0.5416920796101777, 0.5405077966193075, 0.5393397884837879, 0.5381877320107138, 0.5370512835436118, 0.5359301519011452, 
  // 0.5348240318995932, 0.5337326124779957, 0.5326556047423818, 0.5315927475430485, 0.5305437392189766, 0.5295083326073062, 
  // 0.5284862668047939, 0.5274772682475128, 0.5264811087402399, 0.5254975340115511, 0.524526327295383, 0.5235672339972397, 
  // 0.522620048544062, 0.5216845461110134, 0.5207605074637244, 0.5198477362051701};

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }

        // Real den = da;
        // if (rad < rout) {
        //   if (rad < rin) {
        //     den = drat*da;
        //   } else {   // add smooth ramp in density
        //     Real f = (rad-rin) / (rout-rin);
        //     Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
        //     den = std::exp(log_den);
        //   }
        // }

        // std::cout<<std::stold(content[i][0])<<std::endl;

        // std::cout<<std::setprecision(8)<<pcoord->x1v(is-2)<<std::endl;

        // std::cout<<content.size()<<std::endl;

        Real den = std::stold(content[i][0]); //density[i];
        // std::cout << den << std::endl;

        phydro->u(IDN,k,j,i) = den;
        // phydro->u(IPR,k,j,i) = 0.011*std::pow(den,1.666666666667); 
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        // phydro->u(IEN,k,j,i) = 0.011*std::pow(den,1.666666666667)/(1.666666666667-1);

        // std::cout << phydro->u(IDN,k,j,i) << std::endl;
        // if (NON_BAROTROPIC_EOS) {
        //   Real pres = pa;
        //   if (rad < rout) {
        //     if (rad < rin) {
        //       pres = prat*pa;
        //     } else {  // add smooth ramp in pressure
        //       Real f = (rad-rin) / (rout-rin);
        //       Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
        //       pres = std::exp(log_pres);
        //     }
        //   }
          // phydro->u(IEN,k,j,i) = pres/gm1;
        //   if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
        //     phydro->u(IEN,k,j,i) += den;
        // }
      }
    }
  }

  // initialize interface B and total energy
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            pfield->b.x1f(k,j,i) = b0 * std::cos(angle);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real phi = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) =
                b0 * (std::cos(angle) * std::cos(phi) + std::sin(angle) * std::sin(phi));
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0 * std::abs(std::sin(theta))
                                   * (std::cos(angle) * std::cos(phi)
                                      + std::sin(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
            pfield->b.x2f(k,j,i) = b0 * std::sin(angle);
          } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            Real phi = pcoord->x2v(j);
            pfield->b.x2f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real theta = pcoord->x2v(j);
            Real phi = pcoord->x3v(k);
            pfield->b.x2f(k,j,i) = b0 * std::cos(theta)
                                   * (std::cos(angle) * std::cos(phi)
                                      + std::sin(angle) * std::sin(phi));
            if (std::sin(theta) < 0.0)
              pfield->b.x2f(k,j,i) *= -1.0;
          }
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0
              || std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            pfield->b.x3f(k,j,i) = 0.0;
          } else { //if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            Real phi = pcoord->x3v(k);
            pfield->b.x3f(k,j,i) =
                b0 * (std::sin(angle) * std::cos(phi) - std::cos(angle) * std::sin(phi));
          }
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IEN,k,j,i) += 0.5*b0*b0;
        }
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;
  MeshBlock *pmb = my_blocks(0);

  // analysis - check shape of the spherical blast wave
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pmb->phydro->w, 4, IPR, 1);

  // get coordinate location of the center, convert to Cartesian
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // find indices of the center
  int ic, jc, kc;
  for (ic=is; ic<=ie; ic++)
    if (pmb->pcoord->x1f(ic) > x1_0) break;
  ic--;
  for (jc=pmb->js; jc<=pmb->je; jc++)
    if (pmb->pcoord->x2f(jc) > x2_0) break;
  jc--;
  for (kc=pmb->ks; kc<=pmb->ke; kc++)
    if (pmb->pcoord->x3f(kc) > x3_0) break;
  kc--;

  // search pressure maximum in each direction
  Real rmax = 0.0, rmin = 100.0, rave = 0.0;
  int nr = 0;
  for (int o=0; o<=6; o++) {
    int ios = 0, jos = 0, kos = 0;
    if (o == 1) ios=-10;
    else if (o == 2) ios =  10;
    else if (o == 3) jos = -10;
    else if (o == 4) jos =  10;
    else if (o == 5) kos = -10;
    else if (o == 6) kos =  10;
    for (int d=0; d<6; d++) {
      Real pmax = 0.0;
      int imax(0), jmax(0), kmax(0);
      if (d == 0) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i>=is; i--) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 1) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i<=ie; i++) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 2) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j>=js; j--) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 3) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j<=je; j++) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 4) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k>=ks; k--) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      } else { // if (d == 5) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k<=ke; k++) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      }

      Real xm, ym, zm;
      Real x1m = pmb->pcoord->x1v(imax);
      Real x2m = pmb->pcoord->x2v(jmax);
      Real x3m = pmb->pcoord->x3v(kmax);
      if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
        xm = x1m*std::cos(x2m);
        ym = x1m*std::sin(x2m);
        zm = x3m;
      } else {  // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        xm = x1m*std::sin(x2m)*std::cos(x3m);
        ym = x1m*std::sin(x2m)*std::sin(x3m);
        zm = x1m*std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm-x0)+SQR(ym-y0)+SQR(zm-z0));
      if (rad > rmax) rmax = rad;
      if (rad < rmin) rmin = rad;
      rave += rad;
      nr++;
    }
  }
  rave /= static_cast<Real>(nr);

  // use physical grid spacing at center of blast
  Real dr_max;
  Real  x1c = pmb->pcoord->x1v(ic);
  Real dx1c = pmb->pcoord->dx1f(ic);
  Real  x2c = pmb->pcoord->x2v(jc);
  Real dx2c = pmb->pcoord->dx2f(jc);
  Real dx3c = pmb->pcoord->dx3f(kc);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), dx3c);
  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), x1c*std::sin(x2c)*dx3c);
  }
  Real deform=(rmax-rmin)/dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
    }
    std::fprintf(pfile,"# Offset blast wave test in %s coordinates:\n",COORDINATE_SYSTEM);
    std::fprintf(pfile,"# Rmax       Rmin       Rave        Deformation\n");
    std::fprintf(pfile,"%e  %e  %e  %e \n",rmax,rmin,rave,deform);
    std::fclose(pfile);
  }
  return;
}


// refinement condition: check the maximum pressure gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;
  if (pmb->pmy_mesh->f3) {
    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                               +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                               +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
          maxeps = std::max(maxeps, eps);
        }
      }
    }
  } else if (pmb->pmy_mesh->f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                             + SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
        maxeps = std::max(maxeps, eps);
      }
    }
  } else {
    return 0;
  }

  if (maxeps > threshold) return 1;
  if (maxeps < 0.25*threshold) return -1;
  return 0;
}

void MyBoundary_inner(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
  AthenaArray<Real> &w = pmb->phydro->w;
  string fname;
  fname = "/home/sid/athena/vis/python/iso.csv";
  vector<vector<string>> content;
  vector<string> row;
  string line, word;
  fstream file (fname, ios::in);
  if(file.is_open())
  {
    while(getline(file, line))
    {
      row.clear();
      stringstream str(line);
      while(getline(str, word, ','))
        row.push_back(word);
      content.push_back(row);
    }
  }
  else
    cout<<"Could not open the file\n";
  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      for (int i=il-ngh; i<=il-1; ++i){
        prim(IM1,k,j,i) = 0;
        prim(IDN,k,j,i) = std::stold(content[i][0]);
        // prim(IPR,k,j,i) = 0.011*std::pow(std::stold(content[i][0]),1.666666666667);
        // std::cout << std::stold(content[i][0]) << std::endl;
      }
    }
  }
}

void MyBoundary_outer(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh){
  AthenaArray<Real> &w = pmb->phydro->w;
  string fname;
  fname = "/home/sid/athena/vis/python/iso.csv";
  vector<vector<string>> content;
  vector<string> row;
  string line, word;
  fstream file (fname, ios::in);
  if(file.is_open())
  {
    while(getline(file, line))
    {
      row.clear();
      stringstream str(line);
      while(getline(str, word, ','))
        row.push_back(word);
      content.push_back(row);
    }
  }
  else
    cout<<"Could not open the file\n";
  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      for (int i=iu+1; i<=iu+ngh; ++i){
        prim(IM1,k,j,i) = 0;
        prim(IDN,k,j,i) = std::stold(content[i][0]);
        // prim(IPR,k,j,i) = 0.011*std::pow(std::stold(content[i][0]),1.666666666667);
        // std::cout << w(IDN,k,j,iu) << std::endl;
      }
    }
  }
}
