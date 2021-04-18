#include "IK.h"
#include "FK.h"
#include "minivectorTemplate.h"
#include <Eigen/Dense>
#include <adolc/adolc.h>
#include <cassert>
#if defined(_WIN32) || defined(WIN32)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

namespace
{

// Converts degrees to radians.
template<typename real>
inline real deg2rad(real deg) { return deg * M_PI / 180.0; }

template<typename real>
Mat3<real> Euler2Rotation(const real angle[3], RotateOrder order)
{
  Mat3<real> RX = Mat3<real>::getElementRotationMatrix(0, deg2rad(angle[0]));
  Mat3<real> RY = Mat3<real>::getElementRotationMatrix(1, deg2rad(angle[1]));
  Mat3<real> RZ = Mat3<real>::getElementRotationMatrix(2, deg2rad(angle[2]));

  switch(order)
  {
    case RotateOrder::XYZ:
      return RZ * RY * RX;
    case RotateOrder::YZX:
      return RX * RZ * RY;
    case RotateOrder::ZXY:
      return RY * RX * RZ;
    case RotateOrder::XZY:
      return RY * RZ * RX;
    case RotateOrder::YXZ:
      return RZ * RX * RY;
    case RotateOrder::ZYX:
      return RX * RY * RZ;
  }
  assert(0);
}

// Performs forward kinematics, using the provided "fk" class.
// This is the function whose Jacobian matrix will be computed using adolc.
// numIKJoints and IKJointIDs specify which joints serve as handles for IK:
//   IKJointIDs is an array of integers of length "numIKJoints"
// Input: numIKJoints, IKJointIDs, fk, eulerAngles (of all joints)
// Output: handlePositions (world-coordinate positions of all the IK joints; length is 3 * numIKJoints)
template<typename real>
void forwardKinematicsFunction(
    int numIKJoints, const int * IKJointIDs, const FK & fk,
    const std::vector<real> & eulerAngles, std::vector<real> & handlePositions)
{
  int numJoints = fk.getNumJoints();

  // compute local transforms of all joints
  vector<Mat3<real>> localTransformsR;
  vector<Vec3<real>> localTransformsT;

  for (int i = 0; i < numJoints; i++)
  {
      // rotation
      Vec3<real> eulerAngle = {eulerAngles[i * 3 + 0], eulerAngles[i * 3 + 1], eulerAngles[i * 3 + 2]};
      Mat3<real> RJoint = Euler2Rotation(eulerAngle.data(), fk.getJointRotateOrder(i));

      Vec3<real> orient = {fk.getJointOrient(i)[0], fk.getJointOrient(i)[1], fk.getJointOrient(i)[2]};
      Mat3<real> ROrient = Euler2Rotation(orient.data(), fk.getJointRotateOrder(i));

      localTransformsR.push_back(ROrient * RJoint);

      // translation
      Vec3<real> T = {fk.getJointRestTranslation(i)[0], fk.getJointRestTranslation(i)[1], fk.getJointRestTranslation(i)[2]};
      localTransformsT.push_back(T);
  }

  // compute global transforms of all joints (root -> leaves)
  vector<Mat3<real>> globalTransformsR;
  vector<Vec3<real>> globalTransformsT;

  for (int i = 0; i < numJoints; i++)
  {
      int jointIdx = fk.getJointUpdateOrder(i);
      int parentIdx = fk.getJointParent(jointIdx);

      if (parentIdx == -1) // root
      {
          globalTransformsR.push_back(localTransformsR[jointIdx]);
          globalTransformsT.push_back(localTransformsT[jointIdx]);
      } else    // leaves
      {
          Mat3<real> Rout;
          Vec3<real> Tout;
          multiplyAffineTransform4ds(globalTransformsR[parentIdx], globalTransformsT[parentIdx],
                                     localTransformsR[jointIdx], localTransformsT[jointIdx],
                                     Rout, Tout);
          globalTransformsR.push_back(Rout);
          globalTransformsT.push_back(Tout);
      }
  }

  // compute world position of all IK joints
  for (int i = 0; i < numIKJoints; i++)
      for (int j = 0; j < 3; j++)
        handlePositions[i * 3 + j] = globalTransformsT[IKJointIDs[i]][j];
}

} // end anonymous namespaces

IK::IK(int numIKJoints, const int * IKJointIDs, FK * inputFK, int adolc_tagID)
{
  this->numIKJoints = numIKJoints;
  this->IKJointIDs = IKJointIDs;
  this->fk = inputFK;
  this->adolc_tagID = adolc_tagID;

  FKInputDim = fk->getNumJoints() * 3;
  FKOutputDim = numIKJoints * 3;

  train_adolc();
}

void IK::train_adolc()
{
  // first, call trace_on to ask ADOL-C to begin recording how function f is implemented
  trace_on(adolc_tagID); // start tracking computation with ADOL-C

  vector<adouble> eulerAngles(FKInputDim); // define the input of the function f
  for(int i = 0; i < FKInputDim; i++)
      eulerAngles[i] <<= 0.0; // The <<= syntax tells ADOL-C that these are the input variables.

  vector<adouble> handlePositions(FKOutputDim); // define the output of the function f

  // The computation of f goes here:
  //   Use the "forwardKinematicsFunction" as the function that will be computed by adol_c.
  //   This will later make it possible for you to compute the gradient of this function in IK::doIK
  //   (in other words, compute the "Jacobian matrix" J).
  // See ADOLCExample.cpp .
  forwardKinematicsFunction(numIKJoints, IKJointIDs, *fk, eulerAngles, handlePositions);

  vector<double> output(FKOutputDim);
  for(int i = 0; i < FKOutputDim; i++)
      handlePositions[i] >>= output[i]; // Use >>= to tell ADOL-C that y[i] are the output variables

  // Finally, call trace_off to stop recording the function f.
  trace_off(); // ADOL-C tracking finished

}

void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles, bool useDLS)
{
  int numJoints = fk->getNumJoints(); // Note that is NOT the same as numIKJoints!

  // evaluate forwardKinematicsFunction to get handlePositions(n)
  Vec3d handlePositions[numIKJoints];
  for (int i = 0; i < numIKJoints; i++)
      handlePositions[i] = Vec3d(0.0);
  ::function(adolc_tagID, FKOutputDim, FKInputDim, jointEulerAngles->data(), handlePositions->data());

  // evaluate jocobian matrix (mxn)
  double jacobianMatrix[FKOutputDim*FKInputDim]; // row-major order.
  double * jacobianMatrixEachRow[FKOutputDim]; // pointer array where each pointer points to one row of the jacobian matrix
  for (int i = 0; i < FKOutputDim; i++)
      jacobianMatrixEachRow[i] = &jacobianMatrix[i * FKInputDim];
  ::jacobian(adolc_tagID, FKOutputDim, FKInputDim, jointEulerAngles->data(), jacobianMatrixEachRow); // each row is the gradient of one output component of the function

  // convert it to Eigen types
  Eigen::MatrixXd J(FKOutputDim, FKInputDim); // column-major
  for (int r = 0; r < FKOutputDim; r++)
      for (int c = 0; c < FKInputDim; c++)
          J(r, c) = jacobianMatrix[r * FKInputDim + c];

  // compute db (nx1)
  Eigen::VectorXd db(FKOutputDim);
  for (int i = 0; i < numIKJoints; i++)
      for (int j = 0; j < 3; j++)
          db(i * 3 + j) = targetHandlePositions[i][j] - handlePositions[i][j];

  // compute dtheta (mx1)
  Eigen::VectorXd dtheta(FKInputDim);
  solveIK(J, db, dtheta, useDLS);

  // update targetEulerAngles
  for (int i = 0; i < numJoints; i++)
      for (int j = 0; j < 3; j++)
          jointEulerAngles[i][j] += dtheta[i * 3 + j];
}

void IK::solveIK(const Eigen::MatrixXd & J, Eigen::VectorXd & db, Eigen::VectorXd & dtheta, bool useDLS) {

    // subdivide if db is too large
    bool useSubdivision = false;
    double maxDistance = 0.8;
    for (int i = 0; i < db.size(); i++)
        if (db(i) > maxDistance) useSubdivision = true;

    if (useSubdivision)
    {
        db *= 0.5;
        solveIK(J, db, dtheta, useDLS);
        dtheta *= 2.0;
    }
    else if (useDLS)
    {
        // Damped least squares, same as solving the following:
        // (J^T * J + alpha * I) * dtheta = J^T * deltab, solve for dtheta

        Eigen::MatrixXd JT = J.transpose(); // J^T(nxm)
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(FKInputDim, FKInputDim); // I(nxn)
        double alpha = 0.01;

        dtheta = (JT * J + alpha * I).ldlt().solve(JT * db);
    } else
    {
        // pseudo inverse, same as solving the following:
        // dtheta = J_dagger * db, where J_dagger = J^T * (J * J^T)^(-1)

        Eigen::MatrixXd JT = J.transpose(); // J^T(nxm)
        dtheta = JT * (J * JT).inverse() * db;
    }

}

