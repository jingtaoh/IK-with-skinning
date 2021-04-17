#include "skinning.h"
#include "vec3d.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

Skinning::Skinning(int numMeshVertices, const double * restMeshVertexPositions,
    const std::string & meshSkinningWeightsFilename)
{
  this->numMeshVertices = numMeshVertices;
  this->restMeshVertexPositions = restMeshVertexPositions;

  cout << "Loading skinning weights..." << endl;
  ifstream fin(meshSkinningWeightsFilename.c_str());
  assert(fin);
  int numWeightMatrixRows = 0, numWeightMatrixCols = 0;
  fin >> numWeightMatrixRows >> numWeightMatrixCols;
  assert(fin.fail() == false);
  assert(numWeightMatrixRows == numMeshVertices);
  int numJoints = numWeightMatrixCols;

  vector<vector<int>> weightMatrixColumnIndices(numWeightMatrixRows);
  vector<vector<double>> weightMatrixEntries(numWeightMatrixRows);
  fin >> ws;
  while(fin.eof() == false)
  {
    int rowID = 0, colID = 0;
    double w = 0.0;
    fin >> rowID >> colID >> w;
    weightMatrixColumnIndices[rowID].push_back(colID);
    weightMatrixEntries[rowID].push_back(w);
    assert(fin.fail() == false);
    fin >> ws;
  }
  fin.close();

  // Build skinning joints and weights.
  numJointsInfluencingEachVertex = 0;
  for (int i = 0; i < numMeshVertices; i++)
    numJointsInfluencingEachVertex = std::max(numJointsInfluencingEachVertex, (int)weightMatrixEntries[i].size());
  assert(numJointsInfluencingEachVertex >= 2);

  // Copy skinning weights from SparseMatrix into meshSkinningJoints and meshSkinningWeights.
  meshSkinningJoints.assign(numJointsInfluencingEachVertex * numMeshVertices, 0);
  meshSkinningWeights.assign(numJointsInfluencingEachVertex * numMeshVertices, 0.0);
  for (int vtxID = 0; vtxID < numMeshVertices; vtxID++)
  {
    vector<pair<double, int>> sortBuffer(numJointsInfluencingEachVertex);
    for (size_t j = 0; j < weightMatrixEntries[vtxID].size(); j++)
    {
      int frameID = weightMatrixColumnIndices[vtxID][j];
      double weight = weightMatrixEntries[vtxID][j];
      sortBuffer[j] = make_pair(weight, frameID);
    }
    sortBuffer.resize(weightMatrixEntries[vtxID].size());
    assert(sortBuffer.size() > 0);
    sort(sortBuffer.rbegin(), sortBuffer.rend()); // sort in descending order using reverse_iterators
    for(size_t i = 0; i < sortBuffer.size(); i++)
    {
      meshSkinningJoints[vtxID * numJointsInfluencingEachVertex + i] = sortBuffer[i].second;
      meshSkinningWeights[vtxID * numJointsInfluencingEachVertex + i] = sortBuffer[i].first;
    }

    // Note: When the number of joints used on this vertex is smaller than numJointsInfluencingEachVertex,
    // the remaining empty entries are initialized to zero due to vector::assign(XX, 0.0) .
  }
}

void Skinning::applySkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions, bool useLBS) const
{
    if (useLBS) // linear blend skinning
        applyLBS(jointSkinTransforms, newMeshVertexPositions);
   else // dual quaternion skinning
        applyDQS(jointSkinTransforms, newMeshVertexPositions);
}

void Skinning::applyLBS(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions) const
{
    // linear blend skinning
    for(int i = 0; i < numMeshVertices; i++)
    {
        RigidTransform4d transformMatrix(Mat4d::Zero);
        Vec4d restMeshVertexPos(restMeshVertexPositions[3 * i + 0], restMeshVertexPositions[3 * i + 1], restMeshVertexPositions[3 * i + 2], 1);
        Vec4d newMeshVertexPos(0.0);
        for (int j = 0; j < numJointsInfluencingEachVertex; j++)
        {
            newMeshVertexPos += meshSkinningWeights[i * numJointsInfluencingEachVertex + j]
                                * jointSkinTransforms[meshSkinningJoints[i * numJointsInfluencingEachVertex + j]]
                                * restMeshVertexPos;
        }
        for (int j = 0; j < 3; j++)
            newMeshVertexPositions[3 * i + j] = newMeshVertexPos[j];
    }
}

void Skinning::applyDQS(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions) const
{
    // double quaternion skinning

    // form a dual quaternion q = (q0, q1)
    for(int i = 0; i < numMeshVertices; i++)
    {
        Eigen::Quaterniond q0(0, 0, 0, 0);
        Eigen::Quaterniond q1(0, 0, 0, 0);

        for (int j = 0; j < numJointsInfluencingEachVertex; j++)
        {
            int idx = i * numJointsInfluencingEachVertex + j;
            int jointIdx = meshSkinningJoints[idx];

            // form a dual quaternion for each joint q_j = (q0_j, q1_j)

            // q0_j <- R
            Eigen::Matrix3d R;
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                    R(r, c) = jointSkinTransforms[jointIdx][r][c];
            Eigen::Quaterniond q0_j(R);
            // negate q0 if q0.s < 0
            if (q0_j.w() < 0) q0_j.w() *= -1.0;

            // q1_j <- 0.5 * T * q1
            Vec3d t = jointSkinTransforms[jointIdx].getTranslation();
            Eigen::Quaterniond T(0, t[0], t[1], t[2]);
            Eigen::Quaterniond q1_j = T * q0_j;
            q1_j.coeffs() *= 0.5;

            // normalize
            q0_j.normalized();
            q1_j.normalized();

            // sum over each joint's dual quaternion with weight: q = sum(w_j * q_j)
            q0.coeffs() += meshSkinningWeights[idx] * q0_j.coeffs();
            q1.coeffs() += meshSkinningWeights[idx] * q1_j.coeffs();
        }

        // normalize
        q0.normalized();
        q1.normalized();

        // R <- q0
        Eigen::Matrix3d R = q0.toRotationMatrix();

        // T <- 2 * q1 * q0^(-1)
        Eigen::Quaterniond T = q1 * q0.inverse();
        T.coeffs() *= 2.0;
        Eigen::Vector3d t = {T.x(), T.y(), T.z()};

        Eigen::Vector3d restMeshVertexPos = {restMeshVertexPositions[3 * i + 0], restMeshVertexPositions[3 * i + 1], restMeshVertexPositions[3 * i + 2]};
        Eigen::Vector3d newMeshVertexPos = R * restMeshVertexPos + t;

        for (int j = 0; j < 3; j++)
            newMeshVertexPositions[3 * i + j] = newMeshVertexPos[j];
    }
}

