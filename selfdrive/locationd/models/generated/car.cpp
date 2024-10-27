#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8122324459395638877) {
   out_8122324459395638877[0] = delta_x[0] + nom_x[0];
   out_8122324459395638877[1] = delta_x[1] + nom_x[1];
   out_8122324459395638877[2] = delta_x[2] + nom_x[2];
   out_8122324459395638877[3] = delta_x[3] + nom_x[3];
   out_8122324459395638877[4] = delta_x[4] + nom_x[4];
   out_8122324459395638877[5] = delta_x[5] + nom_x[5];
   out_8122324459395638877[6] = delta_x[6] + nom_x[6];
   out_8122324459395638877[7] = delta_x[7] + nom_x[7];
   out_8122324459395638877[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2909140782474590728) {
   out_2909140782474590728[0] = -nom_x[0] + true_x[0];
   out_2909140782474590728[1] = -nom_x[1] + true_x[1];
   out_2909140782474590728[2] = -nom_x[2] + true_x[2];
   out_2909140782474590728[3] = -nom_x[3] + true_x[3];
   out_2909140782474590728[4] = -nom_x[4] + true_x[4];
   out_2909140782474590728[5] = -nom_x[5] + true_x[5];
   out_2909140782474590728[6] = -nom_x[6] + true_x[6];
   out_2909140782474590728[7] = -nom_x[7] + true_x[7];
   out_2909140782474590728[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7515354840926240222) {
   out_7515354840926240222[0] = 1.0;
   out_7515354840926240222[1] = 0;
   out_7515354840926240222[2] = 0;
   out_7515354840926240222[3] = 0;
   out_7515354840926240222[4] = 0;
   out_7515354840926240222[5] = 0;
   out_7515354840926240222[6] = 0;
   out_7515354840926240222[7] = 0;
   out_7515354840926240222[8] = 0;
   out_7515354840926240222[9] = 0;
   out_7515354840926240222[10] = 1.0;
   out_7515354840926240222[11] = 0;
   out_7515354840926240222[12] = 0;
   out_7515354840926240222[13] = 0;
   out_7515354840926240222[14] = 0;
   out_7515354840926240222[15] = 0;
   out_7515354840926240222[16] = 0;
   out_7515354840926240222[17] = 0;
   out_7515354840926240222[18] = 0;
   out_7515354840926240222[19] = 0;
   out_7515354840926240222[20] = 1.0;
   out_7515354840926240222[21] = 0;
   out_7515354840926240222[22] = 0;
   out_7515354840926240222[23] = 0;
   out_7515354840926240222[24] = 0;
   out_7515354840926240222[25] = 0;
   out_7515354840926240222[26] = 0;
   out_7515354840926240222[27] = 0;
   out_7515354840926240222[28] = 0;
   out_7515354840926240222[29] = 0;
   out_7515354840926240222[30] = 1.0;
   out_7515354840926240222[31] = 0;
   out_7515354840926240222[32] = 0;
   out_7515354840926240222[33] = 0;
   out_7515354840926240222[34] = 0;
   out_7515354840926240222[35] = 0;
   out_7515354840926240222[36] = 0;
   out_7515354840926240222[37] = 0;
   out_7515354840926240222[38] = 0;
   out_7515354840926240222[39] = 0;
   out_7515354840926240222[40] = 1.0;
   out_7515354840926240222[41] = 0;
   out_7515354840926240222[42] = 0;
   out_7515354840926240222[43] = 0;
   out_7515354840926240222[44] = 0;
   out_7515354840926240222[45] = 0;
   out_7515354840926240222[46] = 0;
   out_7515354840926240222[47] = 0;
   out_7515354840926240222[48] = 0;
   out_7515354840926240222[49] = 0;
   out_7515354840926240222[50] = 1.0;
   out_7515354840926240222[51] = 0;
   out_7515354840926240222[52] = 0;
   out_7515354840926240222[53] = 0;
   out_7515354840926240222[54] = 0;
   out_7515354840926240222[55] = 0;
   out_7515354840926240222[56] = 0;
   out_7515354840926240222[57] = 0;
   out_7515354840926240222[58] = 0;
   out_7515354840926240222[59] = 0;
   out_7515354840926240222[60] = 1.0;
   out_7515354840926240222[61] = 0;
   out_7515354840926240222[62] = 0;
   out_7515354840926240222[63] = 0;
   out_7515354840926240222[64] = 0;
   out_7515354840926240222[65] = 0;
   out_7515354840926240222[66] = 0;
   out_7515354840926240222[67] = 0;
   out_7515354840926240222[68] = 0;
   out_7515354840926240222[69] = 0;
   out_7515354840926240222[70] = 1.0;
   out_7515354840926240222[71] = 0;
   out_7515354840926240222[72] = 0;
   out_7515354840926240222[73] = 0;
   out_7515354840926240222[74] = 0;
   out_7515354840926240222[75] = 0;
   out_7515354840926240222[76] = 0;
   out_7515354840926240222[77] = 0;
   out_7515354840926240222[78] = 0;
   out_7515354840926240222[79] = 0;
   out_7515354840926240222[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_7622769248808891409) {
   out_7622769248808891409[0] = state[0];
   out_7622769248808891409[1] = state[1];
   out_7622769248808891409[2] = state[2];
   out_7622769248808891409[3] = state[3];
   out_7622769248808891409[4] = state[4];
   out_7622769248808891409[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7622769248808891409[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7622769248808891409[7] = state[7];
   out_7622769248808891409[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3328511493680243564) {
   out_3328511493680243564[0] = 1;
   out_3328511493680243564[1] = 0;
   out_3328511493680243564[2] = 0;
   out_3328511493680243564[3] = 0;
   out_3328511493680243564[4] = 0;
   out_3328511493680243564[5] = 0;
   out_3328511493680243564[6] = 0;
   out_3328511493680243564[7] = 0;
   out_3328511493680243564[8] = 0;
   out_3328511493680243564[9] = 0;
   out_3328511493680243564[10] = 1;
   out_3328511493680243564[11] = 0;
   out_3328511493680243564[12] = 0;
   out_3328511493680243564[13] = 0;
   out_3328511493680243564[14] = 0;
   out_3328511493680243564[15] = 0;
   out_3328511493680243564[16] = 0;
   out_3328511493680243564[17] = 0;
   out_3328511493680243564[18] = 0;
   out_3328511493680243564[19] = 0;
   out_3328511493680243564[20] = 1;
   out_3328511493680243564[21] = 0;
   out_3328511493680243564[22] = 0;
   out_3328511493680243564[23] = 0;
   out_3328511493680243564[24] = 0;
   out_3328511493680243564[25] = 0;
   out_3328511493680243564[26] = 0;
   out_3328511493680243564[27] = 0;
   out_3328511493680243564[28] = 0;
   out_3328511493680243564[29] = 0;
   out_3328511493680243564[30] = 1;
   out_3328511493680243564[31] = 0;
   out_3328511493680243564[32] = 0;
   out_3328511493680243564[33] = 0;
   out_3328511493680243564[34] = 0;
   out_3328511493680243564[35] = 0;
   out_3328511493680243564[36] = 0;
   out_3328511493680243564[37] = 0;
   out_3328511493680243564[38] = 0;
   out_3328511493680243564[39] = 0;
   out_3328511493680243564[40] = 1;
   out_3328511493680243564[41] = 0;
   out_3328511493680243564[42] = 0;
   out_3328511493680243564[43] = 0;
   out_3328511493680243564[44] = 0;
   out_3328511493680243564[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3328511493680243564[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3328511493680243564[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3328511493680243564[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3328511493680243564[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3328511493680243564[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3328511493680243564[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3328511493680243564[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3328511493680243564[53] = -9.8000000000000007*dt;
   out_3328511493680243564[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3328511493680243564[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3328511493680243564[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3328511493680243564[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3328511493680243564[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3328511493680243564[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3328511493680243564[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3328511493680243564[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3328511493680243564[62] = 0;
   out_3328511493680243564[63] = 0;
   out_3328511493680243564[64] = 0;
   out_3328511493680243564[65] = 0;
   out_3328511493680243564[66] = 0;
   out_3328511493680243564[67] = 0;
   out_3328511493680243564[68] = 0;
   out_3328511493680243564[69] = 0;
   out_3328511493680243564[70] = 1;
   out_3328511493680243564[71] = 0;
   out_3328511493680243564[72] = 0;
   out_3328511493680243564[73] = 0;
   out_3328511493680243564[74] = 0;
   out_3328511493680243564[75] = 0;
   out_3328511493680243564[76] = 0;
   out_3328511493680243564[77] = 0;
   out_3328511493680243564[78] = 0;
   out_3328511493680243564[79] = 0;
   out_3328511493680243564[80] = 1;
}
void h_25(double *state, double *unused, double *out_2975096576012810940) {
   out_2975096576012810940[0] = state[6];
}
void H_25(double *state, double *unused, double *out_43261314303816022) {
   out_43261314303816022[0] = 0;
   out_43261314303816022[1] = 0;
   out_43261314303816022[2] = 0;
   out_43261314303816022[3] = 0;
   out_43261314303816022[4] = 0;
   out_43261314303816022[5] = 0;
   out_43261314303816022[6] = 1;
   out_43261314303816022[7] = 0;
   out_43261314303816022[8] = 0;
}
void h_24(double *state, double *unused, double *out_2322939789125350196) {
   out_2322939789125350196[0] = state[4];
   out_2322939789125350196[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4912076179331522874) {
   out_4912076179331522874[0] = 0;
   out_4912076179331522874[1] = 0;
   out_4912076179331522874[2] = 0;
   out_4912076179331522874[3] = 0;
   out_4912076179331522874[4] = 1;
   out_4912076179331522874[5] = 0;
   out_4912076179331522874[6] = 0;
   out_4912076179331522874[7] = 0;
   out_4912076179331522874[8] = 0;
   out_4912076179331522874[9] = 0;
   out_4912076179331522874[10] = 0;
   out_4912076179331522874[11] = 0;
   out_4912076179331522874[12] = 0;
   out_4912076179331522874[13] = 0;
   out_4912076179331522874[14] = 1;
   out_4912076179331522874[15] = 0;
   out_4912076179331522874[16] = 0;
   out_4912076179331522874[17] = 0;
}
void h_30(double *state, double *unused, double *out_4788168680807861501) {
   out_4788168680807861501[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4570957644431424220) {
   out_4570957644431424220[0] = 0;
   out_4570957644431424220[1] = 0;
   out_4570957644431424220[2] = 0;
   out_4570957644431424220[3] = 0;
   out_4570957644431424220[4] = 1;
   out_4570957644431424220[5] = 0;
   out_4570957644431424220[6] = 0;
   out_4570957644431424220[7] = 0;
   out_4570957644431424220[8] = 0;
}
void h_26(double *state, double *unused, double *out_1975616577538008118) {
   out_1975616577538008118[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3784764633177872246) {
   out_3784764633177872246[0] = 0;
   out_3784764633177872246[1] = 0;
   out_3784764633177872246[2] = 0;
   out_3784764633177872246[3] = 0;
   out_3784764633177872246[4] = 0;
   out_3784764633177872246[5] = 0;
   out_3784764633177872246[6] = 0;
   out_3784764633177872246[7] = 1;
   out_3784764633177872246[8] = 0;
}
void h_27(double *state, double *unused, double *out_370691255204113289) {
   out_370691255204113289[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6745720956231849131) {
   out_6745720956231849131[0] = 0;
   out_6745720956231849131[1] = 0;
   out_6745720956231849131[2] = 0;
   out_6745720956231849131[3] = 1;
   out_6745720956231849131[4] = 0;
   out_6745720956231849131[5] = 0;
   out_6745720956231849131[6] = 0;
   out_6745720956231849131[7] = 0;
   out_6745720956231849131[8] = 0;
}
void h_29(double *state, double *unused, double *out_213504586852984144) {
   out_213504586852984144[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4060726300117032036) {
   out_4060726300117032036[0] = 0;
   out_4060726300117032036[1] = 1;
   out_4060726300117032036[2] = 0;
   out_4060726300117032036[3] = 0;
   out_4060726300117032036[4] = 0;
   out_4060726300117032036[5] = 0;
   out_4060726300117032036[6] = 0;
   out_4060726300117032036[7] = 0;
   out_4060726300117032036[8] = 0;
}
void h_28(double *state, double *unused, double *out_8611277916802727476) {
   out_8611277916802727476[0] = state[0];
}
void H_28(double *state, double *unused, double *out_9143125317186562610) {
   out_9143125317186562610[0] = 1;
   out_9143125317186562610[1] = 0;
   out_9143125317186562610[2] = 0;
   out_9143125317186562610[3] = 0;
   out_9143125317186562610[4] = 0;
   out_9143125317186562610[5] = 0;
   out_9143125317186562610[6] = 0;
   out_9143125317186562610[7] = 0;
   out_9143125317186562610[8] = 0;
}
void h_31(double *state, double *unused, double *out_8492665488885234008) {
   out_8492665488885234008[0] = state[8];
}
void H_31(double *state, double *unused, double *out_4410972735411223722) {
   out_4410972735411223722[0] = 0;
   out_4410972735411223722[1] = 0;
   out_4410972735411223722[2] = 0;
   out_4410972735411223722[3] = 0;
   out_4410972735411223722[4] = 0;
   out_4410972735411223722[5] = 0;
   out_4410972735411223722[6] = 0;
   out_4410972735411223722[7] = 0;
   out_4410972735411223722[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_8122324459395638877) {
  err_fun(nom_x, delta_x, out_8122324459395638877);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2909140782474590728) {
  inv_err_fun(nom_x, true_x, out_2909140782474590728);
}
void car_H_mod_fun(double *state, double *out_7515354840926240222) {
  H_mod_fun(state, out_7515354840926240222);
}
void car_f_fun(double *state, double dt, double *out_7622769248808891409) {
  f_fun(state,  dt, out_7622769248808891409);
}
void car_F_fun(double *state, double dt, double *out_3328511493680243564) {
  F_fun(state,  dt, out_3328511493680243564);
}
void car_h_25(double *state, double *unused, double *out_2975096576012810940) {
  h_25(state, unused, out_2975096576012810940);
}
void car_H_25(double *state, double *unused, double *out_43261314303816022) {
  H_25(state, unused, out_43261314303816022);
}
void car_h_24(double *state, double *unused, double *out_2322939789125350196) {
  h_24(state, unused, out_2322939789125350196);
}
void car_H_24(double *state, double *unused, double *out_4912076179331522874) {
  H_24(state, unused, out_4912076179331522874);
}
void car_h_30(double *state, double *unused, double *out_4788168680807861501) {
  h_30(state, unused, out_4788168680807861501);
}
void car_H_30(double *state, double *unused, double *out_4570957644431424220) {
  H_30(state, unused, out_4570957644431424220);
}
void car_h_26(double *state, double *unused, double *out_1975616577538008118) {
  h_26(state, unused, out_1975616577538008118);
}
void car_H_26(double *state, double *unused, double *out_3784764633177872246) {
  H_26(state, unused, out_3784764633177872246);
}
void car_h_27(double *state, double *unused, double *out_370691255204113289) {
  h_27(state, unused, out_370691255204113289);
}
void car_H_27(double *state, double *unused, double *out_6745720956231849131) {
  H_27(state, unused, out_6745720956231849131);
}
void car_h_29(double *state, double *unused, double *out_213504586852984144) {
  h_29(state, unused, out_213504586852984144);
}
void car_H_29(double *state, double *unused, double *out_4060726300117032036) {
  H_29(state, unused, out_4060726300117032036);
}
void car_h_28(double *state, double *unused, double *out_8611277916802727476) {
  h_28(state, unused, out_8611277916802727476);
}
void car_H_28(double *state, double *unused, double *out_9143125317186562610) {
  H_28(state, unused, out_9143125317186562610);
}
void car_h_31(double *state, double *unused, double *out_8492665488885234008) {
  h_31(state, unused, out_8492665488885234008);
}
void car_H_31(double *state, double *unused, double *out_4410972735411223722) {
  H_31(state, unused, out_4410972735411223722);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
