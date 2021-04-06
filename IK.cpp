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
  // TODO: debug forwardKinematicsFunction()
  // Students should implement this.
  // The implementation of this function is very similar to function computeLocalAndGlobalTransforms in the FK class.
  // The recommended approach is to first implement FK::computeLocalAndGlobalTransforms.
  // Then, implement the same algorithm into this function. To do so,
  // you can use fk.getJointUpdateOrder(), fk.getJointRestTranslation(), and fk.getJointRotateOrder() functions.

  int numJoints = fk.getNumJoints();

  // compute local transforms of all joints
  vector<Mat3<real>> localR(numJoints);
  vector<Vec3<real>> localt(numJoints);
  for (int i = 0; i < numJoints; i++)
  {
      // rotation
      Mat3<real> RJoint, ROrient;
      RJoint = Euler2Rotation(&eulerAngles[i * 3], fk.getJointRotateOrder(i));
      real orient[3] = {fk.getJointOrient(i).data()[0], fk.getJointOrient(i).data()[1], fk.getJointOrient(i).data()[2]};
      ROrient = Euler2Rotation(orient, fk.getJointRotateOrder(i));
      // translation
      localR.push_back(ROrient * RJoint);
      Vec3<real> T(fk.getJointRestTranslation(i).data()[0], fk.getJointRestTranslation(i).data()[1], fk.getJointRestTranslation(i).data()[2]);
      localt.push_back(T);
  }
  // compute global transforms of all joints
  vector<Mat3<real>> globalR(numJoints);
  vector<Vec3<real>> globalt(numJoints);
  for (int i = 0; i < numJoints; i++)
  {
      Mat3<real> R1 = (fk.getJointParent(fk.getJointUpdateOrder(i)) == -1) ?
              Mat3<real>(1) : globalR[fk.getJointParent(fk.getJointUpdateOrder(i))];
      Vec3<real> t1 = (fk.getJointParent(fk.getJointUpdateOrder(i)) == -1) ?
              Vec3<real>(1) : globalt[fk.getJointParent(fk.getJointUpdateOrder(i))];
      Mat3<real> R2 = localR[fk.getJointUpdateOrder(i)];
      Vec3<real> t2 = localt[fk.getJointUpdateOrder(i)];

      Mat3<real> Rout;
      Vec3<real> tout;
      multiplyAffineTransform4ds(R1, t1, R2, t2, Rout, tout);
      globalR[fk.getJointUpdateOrder(i)] = Rout;
      globalt[fk.getJointUpdateOrder(i)] = tout;
  }

  // compute world position of all IK joints
  for (int i = 0; i < numIKJoints; i++)
  {
      Vec3<real> handlePos;
      handlePos = globalt[IKJointIDs[i]];
      for (int j = 0; j < 3; j++)
        handlePositions[i * 3 + j] = handlePos[j];
  }

  // Also useful is the multiplyAffineTransform4ds function in minivectorTemplate.h .
  // It would be in principle possible to unify this "forwardKinematicsFunction" and FK::computeLocalAndGlobalTransforms(),
  // so that code is only written once. We considered this; but it is actually not easily doable.
  // If you find a good approach, feel free to document it in the README file, for extra credit.
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
  // TODO: debug IK::train_adolc()
  // Students should implement this.
  // Here, you should setup adol_c:
  //   Define adol_c inputs and outputs.
  //   y = f(x), mapping R^m -> R^n (m = 3 * numJoints, n = 3 * numIKJoints)
  int n = 3 * fk->getNumJoints();   // input dimension is n
  int m = 3 * numIKJoints;          // output dimension is m

  // first, call trace_on to ask ADOL-C to begin recording how function f is implemented
  int adolc_tagID = 1; // This is an ID used in ADOL-C to represent each individual function.
  trace_on(adolc_tagID); // start tracking computation with ADOL-C

  vector<adouble> x(n); // define the input of the function f
  for(int i = 0; i < n; i++)
      x[i] <<= 0.0; // The <<= syntax tells ADOL-C that these are the input variables.

  vector<adouble> y(m); // define the output of the function f

  // The computation of f goes here:
  //   Use the "forwardKinematicsFunction" as the function that will be computed by adol_c.
  //   This will later make it possible for you to compute the gradient of this function in IK::doIK
  //   (in other words, compute the "Jacobian matrix" J).
  // See ADOLCExample.cpp .
  forwardKinematicsFunction(numIKJoints, IKJointIDs, *fk, x, y);

  vector<double> output(m);
  for(int i = 0; i < m; i++)
      y[i] >>= output[i]; // Use >>= to tell ADOL-C that y[i] are the output variables

  // Finally, call trace_off to stop recording the function f.
  trace_off(); // ADOL-C tracking finished

}

void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
  // You may find the following helpful:
  int numJoints = fk->getNumJoints(); // Note that is NOT the same as numIKJoints!

  // TODO: IK::doIK()
  // Students should implement this.
  // Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
  // Specifically, use ::function, and ::jacobian .
  // See ADOLCExample.cpp .
  //
  // Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
  // Note that at entry, "jointEulerAngles" contains the input Euler angles. 
  // Upon exit, jointEulerAngles should contain the new Euler angles.
}

