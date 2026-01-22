#include "../headers/inout.hpp"
#include "../headers/mesh_n_model.hpp"
#include "../headers/concurrent_compute_flow.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include "omp.h"
#include <chrono>

void initialConditionsOMP(
  std::vector<double>& feature, const uint initFuncInd, 
  const model_data& params, const double startTime)
{
  const uint domainPartition = params.domainPartition;
  const uint offsetY = domainPartition + 1;
  const uint offsetZ = (domainPartition+1) * (domainPartition+1);

  const double xStep (params.xLen / domainPartition);
  const double yStep (params.yLen / domainPartition);
  const double zStep (params.zLen / domainPartition);
  
  ABC_Flow functionSet;
  functionContainer<ABC_Flow> fC(&functionSet);
  fC.setFunctions();
  #pragma omp parallel for collapse(2) schedule(static)
  for (size_t k = 0; k < domainPartition + 1; k++)      // Z-Axis
  {
    for (size_t j = 0; j < domainPartition + 1; j++)    // Y-Axis
    {
      for (size_t i = 0; i < domainPartition + 1; i++)  // X-Axis
      {
        uint index (k*offsetZ + j*offsetY + i);
        feature[index] = (functionSet.*fC.indexedFunc[initFuncInd])
          (i*xStep, j*yStep, k*zStep, startTime, params.Reyn) ;
      }
    }
  }
}

void conc_compute_cube_FDA1_3( const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p0)
{
  bool isFDA1 (params.fdaNumber == 1);
  int code = 0;  // exit code
  uint tick = 1; // number of time step
  const auto dimSize ( params.domainPartition );  // 1-dimension size
  const uint offsetY = dimSize + 1;               // j+1 component (i+1 component is just [index+1])
  const uint offsetZ = (dimSize+1) * (dimSize+1); // k+1 component

  const double tau (  params.duration / params.timePartition ); // time step
  const double hX ( params.xLen / dimSize ); // x-step
  const double hY ( params.yLen / dimSize ); // y-step
  const double hZ ( params.zLen / dimSize ); // z-step
  std::ofstream outputFile( params.PATH + params.PATH_log, std::ios::out);
  std::ofstream outputResidualFile( params.PATH + params.PATH_residual, std::ios::out);
  const std::string outputFuncFile =  params.PATH;
  
  // vector size
  const size_t vecSize = (dimSize + 1) * (dimSize + 1) * (dimSize + 1);
  std::vector<double> u1(vecSize);
  std::vector<double> v1(vecSize);
  std::vector<double> w1(vecSize);
  
  Eigen::VectorXd p(vecSize);
  for (size_t i = 0; i < vecSize; i++) p[i] = p0[i];
  //! initialization of equation entities (Ax = b)
  Eigen::SparseMatrix<double> A(vecSize, vecSize);
  Eigen::VectorXd B(vecSize);
  std::vector<Eigen::Triplet<double>> triplets; // entities for filling a sparse matrix

  int numProcs = omp_get_num_procs();
  omp_set_num_threads(numProcs);
  while (tick < 2)//params.timePartition + 1)
  {
    //! velocity exact
    std::vector<double> uExac(vecSize);
    initialConditionsOMP(uExac, 0, params, tick*tau);
    std::vector<double> vExac(vecSize);
    initialConditionsOMP(vExac, 1, params, tick*tau);
    std::vector<double> wExac(vecSize);
    initialConditionsOMP(wExac, 2, params, tick*tau);
    std::vector<double> pExac(vecSize);
    initialConditionsOMP(pExac, 3, params, tick*tau);

    #pragma omp parallel
    {
      //! velocity compute
      //---------------------------- inner knots --------------------------------
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 1; k < dimSize; k++)      // Z-Axis
      {
        for (size_t j = 1; j < dimSize; j++)    // Y-Axis
        {
          for (size_t i = 1; i < dimSize; i++)  // X-Axis
          {
            uint index (k*offsetZ + j*offsetY + i);

            u1[index] =  u[index] - tau * (
              ( u[index+1]* u[index+1] -  u[index-1]* u[index-1]) / (2.0*hX) + 
              ( u[index+offsetY]* v[index+offsetY] -  u[index-offsetY]* v[index-offsetY]) / (2.0*hY) + 
              ( u[index+offsetZ]* w[index+offsetZ] -  u[index-offsetZ]* w[index-offsetZ]) / (2.0*hZ) +
              (p[index+1] - p[index-1]) / (2.0*hX) -
              ( (u[index+1] - 2*u[index] + u[index-1]) / (hX*hX) +
                (u[index+offsetY] - 2*u[index] + u[index-offsetY]) / (hY*hY) +
                (u[index+offsetZ] - 2*u[index] + u[index-offsetZ]) / (hZ*hZ)
              ) / params.Reyn
            );
            v1[index] = v[index] - tau * (
              (v[index+offsetY]*v[index+offsetY] - v[index-offsetY]*v[index-offsetY]) / (2.0*hY) + 
              (u[index+1]*v[index+1] - u[index-1]*v[index-1]) / (2.0*hX) + 
              (v[index+offsetZ]*w[index+offsetZ] - v[index-offsetZ]*w[index-offsetZ]) / (2.0*hZ) +
              (p[index+offsetY] - p[index-offsetY]) / (2.0*hY) -
              ( (v[index+1] - 2*v[index] + v[index-1]) / (hX*hX) +
                (v[index+offsetY] - 2*v[index] + v[index-offsetY]) / (hY*hY) +
                (v[index+offsetZ] - 2*v[index] + v[index-offsetZ]) / (hZ*hZ)
              ) / params.Reyn
            );
            w1[index] = w[index] - tau * (
              (w[index+offsetZ]*w[index+offsetZ] - w[index-offsetZ]*w[index-offsetZ]) / (2.0*hZ) + 
              (u[index+1]*w[index+1] - u[index-1]*w[index-1]) / (2.0*hX) + 
              (v[index+offsetY]*w[index+offsetY] - v[index-offsetY]*w[index-offsetY]) / (2.0*hY) +
              (p[index+offsetZ] - p[index-offsetZ]) / (2.0*hZ) -
              ( (w[index+1] - 2*w[index] + w[index-1]) / (hX*hX) +
                (w[index+offsetY] - 2*w[index] + w[index-offsetY]) / (hY*hY) +
                (w[index+offsetZ] - 2*w[index] + w[index-offsetZ]) / (hZ*hZ)
              ) / params.Reyn
            );
          }
        }
      }
      //-------------------------------------------------------------------------
      #pragma omp barrier
      //---------------------------- border knots -------------------------------
      //! XY-plane Z = 0 / Z = MAX; Neiman's condition dw/dz = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t j = 1; j < dimSize; j++)    // Y-Axis
      {
        for (size_t i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index1 (j*offsetY + i);
          uint index2 (dimSize*offsetZ + j*offsetY + i);

          u1[index1] = uExac[index1];
          u1[index2] = uExac[index2];
          v1[index1] = vExac[index1];
          v1[index2] = vExac[index2];
          w1[index1] = wExac[index1];
          w1[index2] = wExac[index2];
        }
      }
      //! XZ-plane Y = 0 / Y = MAX; Neiman's condition dv/dy = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 1; k < dimSize; k++)    // Z-Axis
      {
        for (size_t i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index1 (k*offsetZ + i);
          uint index2 (k*offsetZ + dimSize*offsetY + i);

          w1[index1] = wExac[index1];
          w1[index2] = wExac[index2];
          u1[index1] = uExac[index1];
          u1[index2] = uExac[index2];
          v1[index1] = vExac[index1];
          v1[index2] = vExac[index2];
        }
      }
      //! YZ-plane X = 0 / X = MAX; Neiman's condition du/dx = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 1; k < dimSize; k++)    // Z-Axis
      {
        for (size_t j = 1; j < dimSize; j++)  // X-Axis
        {
          uint index1 (k*offsetZ + j*offsetY);
          uint index2 (k*offsetZ + j*offsetY + dimSize);

          w1[index1] = wExac[index1];
          w1[index2] = wExac[index2];
          v1[index1] = vExac[index1];
          v1[index2] = vExac[index2];
          u1[index1] = uExac[index1];
          u1[index2] = uExac[index2];
        }
      }
      //-------------------------------------------------------------------------

      //---------------------------- edge knots ---------------------------------
      //! Y = Z = (0 || MAX); dv/dy = 0; dw/dz = 0
      #pragma omp for schedule(static)
      for (size_t i = 1; i < dimSize; i++)
      {
        // Y = Z = 0
        uint index1 (i);
        // Y = MAX, Z = 0
        uint index2 (dimSize*offsetY + i);
        // Y = 0, Z = MAX
        uint index3 (dimSize*offsetZ + i);
        // Y = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        u1[index3] = uExac[index3];
        u1[index4] = uExac[index4];

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
        v1[index3] = vExac[index3];
        v1[index4] = vExac[index4];

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        w1[index3] = wExac[index3];
        w1[index4] = wExac[index4];

      }
      //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
      #pragma omp for schedule(static)
      for (size_t j = 1; j < dimSize; j++)
      {
        // X = Z = 0
        uint index1 (j*offsetY);
        // X = MAX, Z = 0
        uint index2 (j*offsetY + dimSize);
        // X = 0, Z = MAX
        uint index3 (dimSize*offsetZ + j*offsetY);
        // X = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + j*offsetY + dimSize);

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
        v1[index3] = vExac[index3];
        v1[index4] = vExac[index4];

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        u1[index3] = uExac[index3];
        u1[index4] = uExac[index4];

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        w1[index3] = wExac[index3];
        w1[index4] = wExac[index4];
      }
      //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
      #pragma omp for schedule(static)
      for (size_t k = 1; k < dimSize; k++)
      {
        // X = Y = 0
        uint index1 (k*offsetZ);
        // X = MAX, Y = 0
        uint index2 (k*offsetZ + dimSize);
        // X = 0, Y = MAX
        uint index3 (k*offsetZ + dimSize*offsetY);
        // X = MAX, Y = MAX
        uint index4 (k*offsetZ + dimSize*offsetY + dimSize);

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        w1[index3] = wExac[index3];
        w1[index4] = wExac[index4];

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        u1[index3] = uExac[index3];
        u1[index4] = uExac[index4];

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
        v1[index3] = vExac[index3];
        v1[index4] = vExac[index4];
      }
      //-------------------------------------------------------------------------

    } // end of OMP velocity concurrent section -------------------------------

    //------------------------- cube vertices ---------------------------------
    uint index (0);
    // X = Y = Z = 0
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    
    // velocity residual
    //-------------------------------------------------------------------------
    const double velResidual 
      = conc_velocity_residual_FDA1_3(params, u1, v1,w1,p0, uExac,vExac,wExac,pExac);
    outputResidualFile << velResidual << std::endl;

    //std::cout << velResidual << std::endl;
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // move data from upper time layer
    u.clear(); u = std::move(u1);
    v.clear(); v = std::move(v1);
    w.clear(); w = std::move(w1);
    u1.resize(vecSize); v1.resize(vecSize); w1.resize(vecSize);
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    #pragma omp parallel
    {
      //! pressure compute
      //---------------------------- inner knots --------------------------------
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)      // Z-Axis
      {
        for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
        {
          for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
          {
            uint index (k*offsetZ + j*offsetY + i); // central index
            // values at indexes
            const double indexVal_000 = - (1.0/(2.0*hX*hX) + 1.0/(2.0*hY*hY) + 1.0/(2.0*hZ*hZ)); // centre
            const double indexVal_R00 = 1.0/(4.0*hX*hX); // +offsetX
            const double indexVal_L00 = 1.0/(4.0*hX*hX); // -offsetX
            const double indexVal_0R0 = 1.0/(4.0*hY*hY); // +offsetY
            const double indexVal_0L0 = 1.0/(4.0*hY*hY); // -offsetY
            const double indexVal_00R = 1.0/(4.0*hZ*hZ); // +offsetZ
            const double indexVal_00L = 1.0/(4.0*hZ*hZ); // -offsetZ
            // matrix construct via EIGEN triplets 
            std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
              triplets_local.emplace_back(index,index,indexVal_000);
              triplets_local.emplace_back(index,index+2,indexVal_R00);
              triplets_local.emplace_back(index,index-2,indexVal_L00);
              triplets_local.emplace_back(index,index+2*offsetY,indexVal_0R0);
              triplets_local.emplace_back(index,index-2*offsetY,indexVal_0L0);
              triplets_local.emplace_back(index,index+2*offsetZ,indexVal_00R);
              triplets_local.emplace_back(index,index-2*offsetZ,indexVal_00L);
            #pragma omp critical
            triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
            // right side assignement
            B[index] = 
              -(u[index+2]*u[index+2] - 2*u[index]*u[index] + u[index-2]*u[index-2])
                /(4.0*hX*hX)
              -(v[index+2*offsetY]*v[index+2*offsetY] - 2*v[index]*v[index] + v[index-2*offsetY]*v[index-2*offsetY])
                /(4.0*hY*hY)
              -(w[index+2*offsetZ]*w[index+2*offsetZ] - 2*w[index]*w[index] + w[index-2*offsetZ]*w[index-2*offsetZ])
                /(4.0*hZ*hZ)
              -2.0 * (
                (u[index+offsetY+1]*v[index+offsetY+1] 
                  - u[index-offsetY+1]*v[index-offsetY+1]
                  - u[index+offsetY-1]*v[index+offsetY-1] 
                  + u[index-offsetY-1]*v[index-offsetY-1]) 
                    / (4.0*hX*hY)
                + (u[index+offsetZ+1]*w[index+offsetZ+1] 
                  - u[index-offsetZ+1]*w[index-offsetZ+1]
                  - u[index+offsetZ-1]*w[index+offsetZ-1] 
                  + u[index-offsetZ-1]*w[index-offsetZ-1]) 
                    / (4.0*hX*hZ)
                + (v[index+offsetZ+offsetY]*w[index+offsetZ+offsetY] 
                  - v[index-offsetZ+offsetY]*w[index-offsetZ+offsetY]
                  - v[index+offsetZ-offsetY]*w[index+offsetZ-offsetY] 
                  + v[index-offsetZ-offsetY]*w[index-offsetZ-offsetY])
                    / (4.0*hY*hZ)
              );
            //! for FDA1 additional term is required
            if (isFDA1) 
            {
              const double fda1AdditionalTerm ( (
                (u[index+2]-2*u[index+1]+2*u[index-1]-u[index-2]) / (2.0*hX*hX*hX) +
                (u[index+offsetY+1]-2*u[index+1]+u[index-offsetY+1]
                  -u[index+offsetY-1]+2*u[index-1]-u[index-offsetY-1]) 
                    / (2*hX*hY*hY) +
                (u[index+offsetZ+1]-2*u[index+1]+u[index-offsetZ+1]
                  -u[index+offsetZ-1]+2*u[index-1]-u[index-offsetZ-1]) 
                    / (2*hX*hZ*hZ)
              ) / params.Reyn
              + (
                (v[index+2*offsetY]-2*v[index+offsetY]+2*v[index-offsetY]-v[index-2*offsetY]) 
                  / (2.0*hY*hY*hY) +
                (v[index+offsetY+1]-2*v[index+offsetY]+v[index+offsetY-1]
                  -v[index-offsetY+1]+2*v[index-offsetY]-v[index-offsetY-1]) 
                    / (2*hY*hX*hX) +
                (v[index+offsetY+offsetZ]-2*v[index+offsetY]+v[index+offsetY-offsetZ]
                  -v[index-offsetY+offsetZ]+2*v[index-offsetY]-v[index-offsetY-offsetZ]) 
                    / (2*hY*hZ*hZ)
              ) / params.Reyn
              + (
                (w[index+2*offsetZ]-2*w[index+offsetZ]+2*w[index-offsetZ]-w[index-2*offsetZ]) 
                  / (2.0*hZ*hZ*hZ) +
                (w[index+offsetZ+1]-2*w[index+offsetZ]+w[index+offsetZ-1]
                  -w[index-offsetZ+1]+2*w[index-offsetZ]-w[index-offsetZ-1]) 
                    / (2*hZ*hX*hX) +
                (w[index+offsetZ+offsetY]-2*w[index+offsetZ]+w[index+offsetZ-offsetY]
                  -w[index-offsetZ+offsetY]+2*w[index-offsetZ]-w[index-offsetZ-offsetY]) 
                    / (2*hZ*hY*hY)
              ) / params.Reyn );
              
              B[index] += fda1AdditionalTerm;
            }
          }
        }
      }
      //-------------------------------------------------------------------------

      //---------------------------- border knots -------------------------------
      //! XY-plane Z = 0,1 / Z = MAX-1,MAX; Neiman's condition dp/dz = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
      {
        for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
        {
          uint index1 (j*offsetY + i);
          uint index2 (dimSize*offsetZ + j*offsetY + i);
          uint index3 (offsetZ + j*offsetY + i);
          uint index4 ((dimSize-1)*offsetZ + j*offsetY + i);
          // matrix construct via EIGEN triplets 
          std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
            // behind low border
            triplets_local.emplace_back(index3,index3, 1.0);
            // low border
            triplets_local.emplace_back(index1,index1, 1.0);
            // behind up border
            triplets_local.emplace_back(index4,index4, 1.0);
            // up border
            triplets_local.emplace_back(index2,index2, 1.0);
          #pragma omp critical
          triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
          B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
        }
      }
      //! XZ-plane Y = 0,1 / Y = MAX-1,MAX; Neiman's condition dp/dy = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
      {
        for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
        {
          uint index1 (k*offsetZ + i);
          uint index2 (k*offsetZ + dimSize*offsetY + i);
          uint index3 (k*offsetZ + offsetY + 1);
          uint index4 (k*offsetZ + (dimSize-1)*offsetY + i);
          // matrix construct via EIGEN triplets 
          std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
            // behind Y0 border
            triplets_local.emplace_back(index3,index3, 1.0);
            // low border
            triplets_local.emplace_back(index1,index1, 1.0);
            // behind up border
            triplets_local.emplace_back(index4,index4, 1.0);
            // up border
            triplets_local.emplace_back(index2,index2, 1.0);
          #pragma omp critical
          triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
          B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
        }
      }
      
      //! YZ-plane X = 0,1 / X = MAX-1,MAX; Neiman's condition dp/dx = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
      {
        for (size_t j = 2; j < dimSize-1; j++)  // X-Axis
        {
          uint index1 (k*offsetZ + j*offsetY);
          uint index2 (k*offsetZ + j*offsetY + dimSize);
          uint index3 (k*offsetZ + j*offsetY + 1);
          uint index4 (k*offsetZ + j*offsetY + dimSize-1);
          // matrix construct via EIGEN triplets (flow out)
          std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
            // behind X0 border
            triplets_local.emplace_back(index3,index3, 1.0);
            // low border
            triplets_local.emplace_back(index1,index1, 1.0);
            // behind up border
            triplets_local.emplace_back(index4,index4, 1.0);
            // up border
            triplets_local.emplace_back(index2,index2, 1.0);
          #pragma omp critical
          triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
          B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
        }
      }
      //-------------------------------------------------------------------------

      //---------------------------- edge knots ---------------------------------
      //! Y = Z = (0,1 || MAX,MAX-1); 
      #pragma omp for schedule(static)
      for (size_t i = 2; i < dimSize-1; i++)
      {
        // Y = Z = 0
        uint index1 (i);
        // Y = MAX, Z = 0
        uint index2 (dimSize*offsetY + i);
        // Y = 0, Z = MAX
        uint index3 (dimSize*offsetZ + i);
        // Y = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

        std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];

        // Y = Z = 1
        index1 = offsetZ + offsetY + i;
        // Y = MAX-1, Z = 1
        index2 = (dimSize-1)*offsetY + offsetZ + i;
        // Y = 1, Z = MAX-1
        index3 = (dimSize-1)*offsetZ + offsetY + i;
        // Y = MAX-1, Z = MAX-1
        index4 = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + i;
        
        triplets_local.clear(); // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
      }
      //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
      #pragma omp for schedule(static)
      for (size_t j = 2; j < dimSize-1; j++)
      {
        // X = Z = 0
        uint index1 (j*offsetY);
        // X = MAX, Z = 0
        uint index2 (j*offsetY + dimSize);
        // X = 0, Z = MAX
        uint index3 (dimSize*offsetZ + j*offsetY);
        // X = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + j*offsetY + dimSize);
        
        std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];

        // X = Z = 1
        index1 = offsetZ + j*offsetY + 1;
        // X = MAX-1, Z = 1
        index2 = offsetZ + j*offsetY + dimSize-1;
        // X = 1, Z = MAX-1
        index3 = (dimSize-1)*offsetZ + j*offsetY + 1;
        // X = MAX-1, Z = MAX-1
        index4 = (dimSize-1)*offsetZ + j*offsetY + dimSize-1;
        
        triplets_local.clear();
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
      }
      //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
      #pragma omp for schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)
      {
        // X = Y = 0
        uint index1 (k*offsetZ);
        // X = MAX, Y = 0
        uint index2 (k*offsetZ + dimSize);
        // X = 0, Y = MAX
        uint index3 (k*offsetZ + dimSize*offsetY);
        // X = MAX, Y = MAX
        uint index4 (k*offsetZ + dimSize*offsetY + dimSize);
        
        std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];

        // X = Y = 1
        index1 = k*offsetZ + offsetY + 1;
        // X = MAX-1, Y = 1
        index2 = k*offsetZ + offsetY + dimSize-1;
        // X = 1, Y = MAX-1
        index3 = k*offsetZ + (dimSize-1)*offsetY + 1;
        // X = MAX-1, Y = MAX-1
        index4 = k*offsetZ + (dimSize-1)*offsetY + dimSize-1;
        
        triplets_local.clear();
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
      }
      //-------------------------------------------------------------------------

    } // end of OMP pressure concurrent section -------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 1
    index = offsetZ + offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = MAX-1, Y = Z = 1
    index = offsetZ + offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = 1, Y = MAX-1
    index = offsetZ + (dimSize-1)*offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = 1, Z = MAX-1
    index = (dimSize-1)*offsetZ + offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = MAX-1, Z = 1
    index = offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = MAX-1, Y = 1
    index = (dimSize-1)*offsetZ + offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // Y = Z = MAX-1, X = 1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = Z = MAX-1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];

    // X = Y = Z = 0
    index = 0;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // call Eigen solver for Pressure (idk, does is support concurrency or not)
    //-------------------------------------------------------------------------
    A.setFromTriplets (triplets.begin(), triplets.end());

    // biconjugate gradient stabilized algorithm
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, 
      Eigen::DiagonalPreconditioner<double>> solver(A);
    if (solver.info() != Eigen::Success)
    {
      outputFile << "Can not build preconditioner" << std::endl;
      std::cout << "Can not build preconditioner" << std::endl;
      return;
    }
    Eigen::VectorXd pHat(vecSize);
    pHat = solver.solveWithGuess(B, p);
    if (solver.info() != Eigen::Success)
    {
      outputFile << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      std::cout << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      return;
    }
    //-------------------------------------------------------------------------
    // refresh Pressure values (transfer Eigen-vector to vector)
    //-------------------------------------------------------------------------
    p = std::move(pHat);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < vecSize; ++i) p0[i] = p[i];
    tick += 1;
    if ((tick % 10 == 0)) 
    {
      std::cout << "tick: " << tick << '\n';
      /*code = funcOutput(outputFuncFile, "v1", std::to_string(tick), ".txt", u, params, false);
      code = funcOutput(outputFuncFile, "v2", std::to_string(tick), ".txt", v, params, false);
      code = funcOutput(outputFuncFile, "v3", std::to_string(tick), ".txt", w, params, false);
      code = funcOutput(outputFuncFile, "p", std::to_string(tick), ".txt", p0, params, false);*/
    }
  }
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < vecSize; ++i) p0[i] = p[i];
  std::cout << "final tick: " << tick << '\n';
}

void conc_compute_cube_FDA2_4( const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p0)
{
  bool isFDA2 (params.fdaNumber == 1);
  int code = 0;  // exit code
  uint tick = 1; // number of time step
  const auto dimSize ( params.domainPartition );  // 1-dimension size
  const uint offsetY = dimSize + 1;               // j+1 component (i+1 component is just [index+1])
  const uint offsetZ = (dimSize+1) * (dimSize+1); // k+1 component

  const double tau (  params.duration / params.timePartition ); // time step
  const double hX ( params.xLen / dimSize ); // x-step
  const double hY ( params.yLen / dimSize ); // y-step
  const double hZ ( params.zLen / dimSize ); // z-step
  std::ofstream outputFile( params.PATH + params.PATH_log, std::ios::out);
  std::ofstream outputResidualFile( params.PATH + params.PATH_residual, std::ios::out);
  const std::string outputFuncFile =  params.PATH;
  
  // vector size
  const size_t vecSize = (dimSize + 1) * (dimSize + 1) * (dimSize + 1);
  std::vector<double> u1(vecSize);
  std::vector<double> v1(vecSize);
  std::vector<double> w1(vecSize);
  
  Eigen::VectorXd p(vecSize);
  for (size_t i = 0; i < vecSize; i++) p[i] = p0[i];
  //! initialization of equation entities (Ax = b)
  Eigen::SparseMatrix<double> A(vecSize, vecSize);
  Eigen::VectorXd B(vecSize);
  std::vector<Eigen::Triplet<double>> triplets; // entities for filling a sparse matrix

  int numProcs = omp_get_num_procs();
  omp_set_num_threads(numProcs);
  while (tick < 2)//params.timePartition + 1)
  {
    //! velocity exact
    std::vector<double> uExac(vecSize);
    initialConditionsOMP(uExac, 0, params, tick*tau);
    std::vector<double> vExac(vecSize);
    initialConditionsOMP(vExac, 1, params, tick*tau);
    std::vector<double> wExac(vecSize);
    initialConditionsOMP(wExac, 2, params, tick*tau);
    std::vector<double> pExac(vecSize);
    initialConditionsOMP(pExac, 3, params, tick*tau);

    #pragma omp parallel
    {
      //! velocity compute
      //---------------------------- inner knots --------------------------------
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 1; k < dimSize; k++)      // Z-Axis
      {
        for (size_t j = 1; j < dimSize; j++)    // Y-Axis
        {
          for (size_t i = 1; i < dimSize; i++)  // X-Axis
          {
            uint index (k*offsetZ + j*offsetY + i);

            u1[index] = u[index] - tau * ( 
            u[index]*(u[index+1]-u[index-1])/(2*hX) 
            + v[index]*(u[index+offsetY]-u[index-offsetY])/(2*hY)
            + w[index]*(u[index+offsetZ]-u[index-offsetZ])/(2*hZ) 
            + (p[index+1]-p[index-1])/(2*hX)
            - ( 
              (u[index+1] - 2*u[index] + u[index-1]) / (hX*hX) +
              (u[index+offsetY] - 2*u[index] + u[index-offsetY]) / (hY*hY) +
              (u[index+offsetZ] - 2*u[index] + u[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
          );
          v1[index] = v[index] - tau * ( 
            u[index]*(v[index+1]-v[index-1])/(2*hX) 
            + v[index]*(v[index+offsetY]-v[index-offsetY])/(2*hY)
            + w[index]*(v[index+offsetZ]-v[index-offsetZ])/(2*hZ) 
            + (p[index+offsetY]-p[index-offsetY])/(2*hY)
            - ( 
              (v[index+1] - 2*v[index] + v[index-1]) / (hX*hX) +
              (v[index+offsetY] - 2*v[index] + v[index-offsetY]) / (hY*hY) +
              (v[index+offsetZ] - 2*v[index] + v[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
           );
          w1[index] = w[index] - tau * ( 
            u[index]*(w[index+1]-w[index-1])/(2*hX) 
            + v[index]*(w[index+offsetY]-w[index-offsetY])/(2*hY)
            + w[index]*(w[index+offsetZ]-w[index-offsetZ])/(2*hZ) 
            + (p[index+offsetZ]-p[index-offsetZ])/(2*hZ)
            - ( 
              (w[index+1] - 2*w[index] + w[index-1]) / (hX*hX) +
              (w[index+offsetY] - 2*w[index] + w[index-offsetY]) / (hY*hY) +
              (w[index+offsetZ] - 2*w[index] + w[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
           );
          }
        }
      }
      //-------------------------------------------------------------------------
      #pragma omp barrier
      //---------------------------- border knots -------------------------------
      //! XY-plane Z = 0 / Z = MAX; Neiman's condition dw/dz = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t j = 1; j < dimSize; j++)    // Y-Axis
      {
        for (size_t i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index1 (j*offsetY + i);
          uint index2 (dimSize*offsetZ + j*offsetY + i);

          u1[index1] = uExac[index1];
          u1[index2] = uExac[index2];
          v1[index1] = vExac[index1];
          v1[index2] = vExac[index2];
          w1[index1] = wExac[index1];
          w1[index2] = wExac[index2];
        }
      }
      //! XZ-plane Y = 0 / Y = MAX; Neiman's condition dv/dy = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 1; k < dimSize; k++)    // Z-Axis
      {
        for (size_t i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index1 (k*offsetZ + i);
          uint index2 (k*offsetZ + dimSize*offsetY + i);

          w1[index1] = wExac[index1];
          w1[index2] = wExac[index2];
          u1[index1] = uExac[index1];
          u1[index2] = uExac[index2];
          v1[index1] = vExac[index1];
          v1[index2] = vExac[index2];
        }
      }
      //! YZ-plane X = 0 / X = MAX; Neiman's condition du/dx = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 1; k < dimSize; k++)    // Z-Axis
      {
        for (size_t j = 1; j < dimSize; j++)  // X-Axis
        {
          uint index1 (k*offsetZ + j*offsetY);
          uint index2 (k*offsetZ + j*offsetY + dimSize);

          w1[index1] = wExac[index1];
          w1[index2] = wExac[index2];
          v1[index1] = vExac[index1];
          v1[index2] = vExac[index2];
          u1[index1] = uExac[index1];
          u1[index2] = uExac[index2];
        }
      }
      //-------------------------------------------------------------------------

      //---------------------------- edge knots ---------------------------------
      //! Y = Z = (0 || MAX); dv/dy = 0; dw/dz = 0
      #pragma omp for schedule(static)
      for (size_t i = 1; i < dimSize; i++)
      {
        // Y = Z = 0
        uint index1 (i);
        // Y = MAX, Z = 0
        uint index2 (dimSize*offsetY + i);
        // Y = 0, Z = MAX
        uint index3 (dimSize*offsetZ + i);
        // Y = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        u1[index3] = uExac[index3];
        u1[index4] = uExac[index4];

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
        v1[index3] = vExac[index3];
        v1[index4] = vExac[index4];

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        w1[index3] = wExac[index3];
        w1[index4] = wExac[index4];

      }
      //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
      #pragma omp for schedule(static)
      for (size_t j = 1; j < dimSize; j++)
      {
        // X = Z = 0
        uint index1 (j*offsetY);
        // X = MAX, Z = 0
        uint index2 (j*offsetY + dimSize);
        // X = 0, Z = MAX
        uint index3 (dimSize*offsetZ + j*offsetY);
        // X = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + j*offsetY + dimSize);

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
        v1[index3] = vExac[index3];
        v1[index4] = vExac[index4];

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        u1[index3] = uExac[index3];
        u1[index4] = uExac[index4];

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        w1[index3] = wExac[index3];
        w1[index4] = wExac[index4];
      }
      //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
      #pragma omp for schedule(static)
      for (size_t k = 1; k < dimSize; k++)
      {
        // X = Y = 0
        uint index1 (k*offsetZ);
        // X = MAX, Y = 0
        uint index2 (k*offsetZ + dimSize);
        // X = 0, Y = MAX
        uint index3 (k*offsetZ + dimSize*offsetY);
        // X = MAX, Y = MAX
        uint index4 (k*offsetZ + dimSize*offsetY + dimSize);

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        w1[index3] = wExac[index3];
        w1[index4] = wExac[index4];

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        u1[index3] = uExac[index3];
        u1[index4] = uExac[index4];

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
        v1[index3] = vExac[index3];
        v1[index4] = vExac[index4];
      }
      //-------------------------------------------------------------------------

    } // end of OMP velocity concurrent section -------------------------------

    //------------------------- cube vertices ---------------------------------
    uint index (0);
    // X = Y = Z = 0
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    
    // velocity residual
    //-------------------------------------------------------------------------
    const double velResidual 
      = conc_velocity_residual_FDA2_4(params, u1, v1,w1,p0, uExac,vExac,wExac,pExac);
    outputResidualFile << velResidual << std::endl;

    //std::cout << velResidual << std::endl;
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // move data from upper time layer
    u.clear(); u = std::move(u1);
    v.clear(); v = std::move(v1);
    w.clear(); w = std::move(w1);
    u1.resize(vecSize); v1.resize(vecSize); w1.resize(vecSize);
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    #pragma omp parallel
    {
      //! pressure compute
      //---------------------------- inner knots --------------------------------
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)      // Z-Axis
      {
        for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
        {
          for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
          {
            uint index (k*offsetZ + j*offsetY + i); // central index
            // values at indexes
            const double indexVal_000 = - (1.0/(2.0*hX*hX) + 1.0/(2.0*hY*hY) + 1.0/(2.0*hZ*hZ)); // centre
            const double indexVal_R00 = 1.0/(4.0*hX*hX); // +offsetX
            const double indexVal_L00 = 1.0/(4.0*hX*hX); // -offsetX
            const double indexVal_0R0 = 1.0/(4.0*hY*hY); // +offsetY
            const double indexVal_0L0 = 1.0/(4.0*hY*hY); // -offsetY
            const double indexVal_00R = 1.0/(4.0*hZ*hZ); // +offsetZ
            const double indexVal_00L = 1.0/(4.0*hZ*hZ); // -offsetZ
            // matrix construct via EIGEN triplets 
            std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
              triplets_local.emplace_back(index,index,indexVal_000);
              triplets_local.emplace_back(index,index+2,indexVal_R00);
              triplets_local.emplace_back(index,index-2,indexVal_L00);
              triplets_local.emplace_back(index,index+2*offsetY,indexVal_0R0);
              triplets_local.emplace_back(index,index-2*offsetY,indexVal_0L0);
              triplets_local.emplace_back(index,index+2*offsetZ,indexVal_00R);
              triplets_local.emplace_back(index,index-2*offsetZ,indexVal_00L);
            #pragma omp critical
            triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
            // right side assignement
            B[index] = 0.0;
            const double FDA2_right_part ( -(
              (u[index+1]-u[index-1])*(u[index+1]-u[index-1])/(4.0*hX*hX) +
              (v[index+offsetY]-v[index-offsetY])*(v[index+offsetY]-v[index-offsetY])/(4.0*hY*hY) +
              (w[index+offsetZ]-w[index-offsetZ])*(w[index+offsetZ]-w[index-offsetZ])/(4.0*hZ*hZ) +
              2.0 * ( (v[index+1]-v[index-1])*(u[index+offsetY]-u[index-offsetY])/(hX*hY)
                + (w[index+1]-w[index-1])*(u[index+offsetZ]-u[index-offsetZ])/(hX*hZ)
                + (v[index+offsetZ]-v[index-offsetZ])*(w[index+offsetY]-w[index-offsetY])/(hZ*hY)
              )
            ) );
            if (isFDA2) {
              B[index] += FDA2_right_part;
              continue;
            }
            const double FDA4_right_part ( 
            -(
              (u[index+2]*u[index+1]-u[index+1]*u[index]-u[index]*u[index-1]+u[index-1]*u[index-2])
                /(4.0*hX*hX) +
              (v[index+2*offsetY]*v[index+offsetY]-v[index+offsetY]*v[index]
                -v[index]*v[index-offsetY]+v[index-offsetY]*v[index-2*offsetY]) 
                / (4.0*hY*hY) +
              (w[index+2*offsetZ]*w[index+offsetZ]-w[index+offsetZ]*w[index]
                -w[index]*w[index-offsetZ]+w[index-offsetZ]*w[index-2*offsetZ]) 
                / (4.0*hZ*hZ) +
              (v[index+1]*u[index+offsetY+1]-v[index+1]*u[index-offsetY+1]
                -v[index-1]*u[index+offsetY-1]+v[index-1]*u[index-offsetY-1]) 
                / (4.0*hX*hY) +
              (u[index+offsetZ]*w[index+offsetZ+1]-u[index+offsetZ]*w[index+offsetZ-1]
                -u[index-offsetZ]*w[index-offsetZ+1]+u[index-offsetZ]*w[index-offsetZ-1]) 
                / (4.0*hZ*hX) +
              (w[index+offsetY]*v[index+offsetY+offsetZ]-w[index+offsetY]*v[index+offsetY-offsetZ]
                -w[index-offsetY]*v[index-offsetY+offsetZ]+w[index-offsetY]*v[index-offsetY-offsetZ]) 
                / (4.0*hY*hZ) -
              ( (
              (u[index+2]-2*u[index+1]+2*u[index-1]-u[index-2]) / (2.0*hX*hX*hX) +
              (u[index+offsetY+1]-2*u[index+1]+u[index-offsetY+1]
                -u[index+offsetY-1]+2*u[index-1]-u[index-offsetY-1]) 
                  / (2*hX*hY*hY) +
              (u[index+offsetZ+1]-2*u[index+1]+u[index-offsetZ+1]
                -u[index+offsetZ-1]+2*u[index-1]-u[index-offsetZ-1]) 
                  / (2*hX*hZ*hZ)
                ) / params.Reyn
              + (
              (v[index+2*offsetY]-2*v[index+offsetY]+2*v[index-offsetY]-v[index-2*offsetY]) 
                / (2.0*hY*hY*hY) +
              (v[index+offsetY+1]-2*v[index+offsetY]+v[index+offsetY-1]
                -v[index-offsetY+1]+2*v[index-offsetY]-v[index-offsetY-1]) 
                  / (2*hY*hX*hX) +
              (v[index+offsetY+offsetZ]-2*v[index+offsetY]+v[index+offsetY-offsetZ]
                -v[index-offsetY+offsetZ]+2*v[index-offsetY]-v[index-offsetY-offsetZ]) 
                  / (2*hY*hZ*hZ)
                ) / params.Reyn
              + (
              (w[index+2*offsetZ]-2*w[index+offsetZ]+2*w[index-offsetZ]-w[index-2*offsetZ]) 
                / (2.0*hZ*hZ*hZ) +
              (w[index+offsetZ+1]-2*w[index+offsetZ]+w[index+offsetZ-1]
                -w[index-offsetZ+1]+2*w[index-offsetZ]-w[index-offsetZ-1]) 
                  / (2*hZ*hX*hX) +
              (w[index+offsetZ+offsetY]-2*w[index+offsetZ]+w[index+offsetZ-offsetY]
                -w[index-offsetZ+offsetY]+2*w[index-offsetZ]-w[index-offsetZ-offsetY]) 
                  / (2*hZ*hY*hY)
                ) / params.Reyn )
            ) );
            B[index] += FDA4_right_part;
          }
        }
      }
      //-------------------------------------------------------------------------

      //---------------------------- border knots -------------------------------
      //! XY-plane Z = 0,1 / Z = MAX-1,MAX; Neiman's condition dp/dz = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
      {
        for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
        {
          uint index1 (j*offsetY + i);
          uint index2 (dimSize*offsetZ + j*offsetY + i);
          uint index3 (offsetZ + j*offsetY + i);
          uint index4 ((dimSize-1)*offsetZ + j*offsetY + i);
          // matrix construct via EIGEN triplets 
          std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
            // behind low border
            triplets_local.emplace_back(index3,index3, 1.0);
            // low border
            triplets_local.emplace_back(index1,index1, 1.0);
            // behind up border
            triplets_local.emplace_back(index4,index4, 1.0);
            // up border
            triplets_local.emplace_back(index2,index2, 1.0);
          #pragma omp critical
          triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
          B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
        }
      }
      //! XZ-plane Y = 0,1 / Y = MAX-1,MAX; Neiman's condition dp/dy = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
      {
        for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
        {
          uint index1 (k*offsetZ + i);
          uint index2 (k*offsetZ + dimSize*offsetY + i);
          uint index3 (k*offsetZ + offsetY + 1);
          uint index4 (k*offsetZ + (dimSize-1)*offsetY + i);
          // matrix construct via EIGEN triplets 
          std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
            // behind Y0 border
            triplets_local.emplace_back(index3,index3, 1.0);
            // low border
            triplets_local.emplace_back(index1,index1, 1.0);
            // behind up border
            triplets_local.emplace_back(index4,index4, 1.0);
            // up border
            triplets_local.emplace_back(index2,index2, 1.0);
          #pragma omp critical
          triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
          B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
        }
      }
      
      //! YZ-plane X = 0,1 / X = MAX-1,MAX; Neiman's condition dp/dx = 0
      #pragma omp for collapse(2) schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
      {
        for (size_t j = 2; j < dimSize-1; j++)  // X-Axis
        {
          uint index1 (k*offsetZ + j*offsetY);
          uint index2 (k*offsetZ + j*offsetY + dimSize);
          uint index3 (k*offsetZ + j*offsetY + 1);
          uint index4 (k*offsetZ + j*offsetY + dimSize-1);
          // matrix construct via EIGEN triplets (flow out)
          std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
            // behind X0 border
            triplets_local.emplace_back(index3,index3, 1.0);
            // low border
            triplets_local.emplace_back(index1,index1, 1.0);
            // behind up border
            triplets_local.emplace_back(index4,index4, 1.0);
            // up border
            triplets_local.emplace_back(index2,index2, 1.0);
          #pragma omp critical
          triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
          B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
        }
      }
      //-------------------------------------------------------------------------

      //---------------------------- edge knots ---------------------------------
      //! Y = Z = (0,1 || MAX,MAX-1); 
      #pragma omp for schedule(static)
      for (size_t i = 2; i < dimSize-1; i++)
      {
        // Y = Z = 0
        uint index1 (i);
        // Y = MAX, Z = 0
        uint index2 (dimSize*offsetY + i);
        // Y = 0, Z = MAX
        uint index3 (dimSize*offsetZ + i);
        // Y = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

        std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];

        // Y = Z = 1
        index1 = offsetZ + offsetY + i;
        // Y = MAX-1, Z = 1
        index2 = (dimSize-1)*offsetY + offsetZ + i;
        // Y = 1, Z = MAX-1
        index3 = (dimSize-1)*offsetZ + offsetY + i;
        // Y = MAX-1, Z = MAX-1
        index4 = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + i;
        
        triplets_local.clear(); // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
      }
      //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
      #pragma omp for schedule(static)
      for (size_t j = 2; j < dimSize-1; j++)
      {
        // X = Z = 0
        uint index1 (j*offsetY);
        // X = MAX, Z = 0
        uint index2 (j*offsetY + dimSize);
        // X = 0, Z = MAX
        uint index3 (dimSize*offsetZ + j*offsetY);
        // X = MAX, Z = MAX
        uint index4 (dimSize*offsetZ + j*offsetY + dimSize);
        
        std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];

        // X = Z = 1
        index1 = offsetZ + j*offsetY + 1;
        // X = MAX-1, Z = 1
        index2 = offsetZ + j*offsetY + dimSize-1;
        // X = 1, Z = MAX-1
        index3 = (dimSize-1)*offsetZ + j*offsetY + 1;
        // X = MAX-1, Z = MAX-1
        index4 = (dimSize-1)*offsetZ + j*offsetY + dimSize-1;
        
        triplets_local.clear();
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
      }
      //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
      #pragma omp for schedule(static)
      for (size_t k = 2; k < dimSize-1; k++)
      {
        // X = Y = 0
        uint index1 (k*offsetZ);
        // X = MAX, Y = 0
        uint index2 (k*offsetZ + dimSize);
        // X = 0, Y = MAX
        uint index3 (k*offsetZ + dimSize*offsetY);
        // X = MAX, Y = MAX
        uint index4 (k*offsetZ + dimSize*offsetY + dimSize);
        
        std::vector<Eigen::Triplet<double>> triplets_local; // Private to each thread
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];

        // X = Y = 1
        index1 = k*offsetZ + offsetY + 1;
        // X = MAX-1, Y = 1
        index2 = k*offsetZ + offsetY + dimSize-1;
        // X = 1, Y = MAX-1
        index3 = k*offsetZ + (dimSize-1)*offsetY + 1;
        // X = MAX-1, Y = MAX-1
        index4 = k*offsetZ + (dimSize-1)*offsetY + dimSize-1;
        
        triplets_local.clear();
          triplets_local.emplace_back(index1,index1, 1.0);
          triplets_local.emplace_back(index2,index2, 1.0);
          triplets_local.emplace_back(index3,index3, 1.0);
          triplets_local.emplace_back(index4,index4, 1.0);
        #pragma omp critical
        triplets.insert(triplets.end(), triplets_local.begin(), triplets_local.end());
        B[index1] = pExac[index1]; B[index2] = pExac[index2]; B[index3] = pExac[index3]; B[index4] = pExac[index4];
      }
      //-------------------------------------------------------------------------

    } // end of OMP pressure concurrent section -------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 1
    index = offsetZ + offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = MAX-1, Y = Z = 1
    index = offsetZ + offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = 1, Y = MAX-1
    index = offsetZ + (dimSize-1)*offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = 1, Z = MAX-1
    index = (dimSize-1)*offsetZ + offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = MAX-1, Z = 1
    index = offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = MAX-1, Y = 1
    index = (dimSize-1)*offsetZ + offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // Y = Z = MAX-1, X = 1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + 1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = Z = MAX-1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];

    // X = Y = Z = 0
    index = 0;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    triplets.emplace_back(index,index, 1.0);
    B[index] = pExac[index];
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // call Eigen solver for Pressure (idk, does is support concurrency or not)
    //-------------------------------------------------------------------------
    A.setFromTriplets (triplets.begin(), triplets.end());

    // biconjugate gradient stabilized algorithm
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, 
      Eigen::DiagonalPreconditioner<double>> solver(A);
    if (solver.info() != Eigen::Success)
    {
      outputFile << "Can not build preconditioner" << std::endl;
      std::cout << "Can not build preconditioner" << std::endl;
      return;
    }
    Eigen::VectorXd pHat(vecSize);
    pHat = solver.solveWithGuess(B, p);
    if (solver.info() != Eigen::Success)
    {
      outputFile << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      std::cout << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      return;
    }
    //-------------------------------------------------------------------------
    // refresh Pressure values (transfer Eigen-vector to vector)
    //-------------------------------------------------------------------------
    p = std::move(pHat);
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < vecSize; ++i) p0[i] = p[i];
    tick += 1;
    if ((tick % 10 == 0)) 
    {
      std::cout << "tick: " << tick << '\n';
      /*code = funcOutput(outputFuncFile, "v1", std::to_string(tick), ".txt", u, params, false);
      code = funcOutput(outputFuncFile, "v2", std::to_string(tick), ".txt", v, params, false);
      code = funcOutput(outputFuncFile, "v3", std::to_string(tick), ".txt", w, params, false);
      code = funcOutput(outputFuncFile, "p", std::to_string(tick), ".txt", p0, params, false);*/
    }
  }
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < vecSize; ++i) p0[i] = p[i];
  std::cout << "final tick: " << tick << '\n';
}
//=================================================================================================



//======================================== RESIDUALS ==============================================
double conc_velocity_residual_FDA1_3(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac)
{
  const auto dimSize ( params.domainPartition );
  const uint offsetY = dimSize + 1;
  const uint offsetZ = (dimSize+1) * (dimSize+1);

  const double tau ( params.duration / params.timePartition ); // time step
  const double hX (params.xLen / dimSize );
  const double hY (params.yLen / dimSize );
  const double hZ (params.zLen / dimSize );

  double vectorResidual (0.0);

  int numProcs = omp_get_num_procs();
  omp_set_num_threads(numProcs);
  #pragma omp parallel
  {
  //---------------------------- inner knots --------------------------------
  #pragma omp for collapse(2) schedule(static)
  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index (k*offsetZ + j*offsetY + i);
        double resTerm (0.0); // residual term
        //! Insert precise values to the scheme, take the difference with the estimated values
        resTerm = uEst[index] - ( uExac[index] - tau * (
          (uExac[index+1]*uExac[index+1] - uExac[index-1]*uExac[index-1]) / (2.0*hX) + 
          (uExac[index+offsetY]*vExac[index+offsetY] - uExac[index-offsetY]*vExac[index-offsetY]) / (2.0*hY) + 
          (uExac[index+offsetZ]*wExac[index+offsetZ] - uExac[index-offsetZ]*wExac[index-offsetZ]) / (2.0*hZ) +
          (pExac[index+1] - pExac[index-1]) / (2.0*hX) -
          ( (uExac[index+1] - 2*uExac[index] + uExac[index-1]) / (hX*hX) +
            (uExac[index+offsetY] - 2*uExac[index] + uExac[index-offsetY]) / (hY*hY) +
            (uExac[index+offsetZ] - 2*uExac[index] + uExac[index-offsetZ]) / (hZ*hZ)
          ) / params.Reyn
        ) );
        vectorResidual += resTerm*resTerm;

        resTerm = vEst[index] - ( vExac[index] - tau * (
          (vExac[index+offsetY]*vExac[index+offsetY] - vExac[index-offsetY]*vExac[index-offsetY]) / (2.0*hY) + 
          (uExac[index+1]*vExac[index+1] - uExac[index-1]*vExac[index-1]) / (2.0*hX) + 
          (vExac[index+offsetZ]*wExac[index+offsetZ] - vExac[index-offsetZ]*wExac[index-offsetZ]) / (2.0*hZ) +
          (pExac[index+offsetY] - pExac[index-offsetY]) / (2.0*hY) -
          ( (vExac[index+1] - 2*vExac[index] + vExac[index-1]) / (hX*hX) +
            (vExac[index+offsetY] - 2*vExac[index] + vExac[index-offsetY]) / (hY*hY) +
            (vExac[index+offsetZ] - 2*vExac[index] + vExac[index-offsetZ]) / (hZ*hZ)
          ) / params.Reyn
        ) );
        vectorResidual += resTerm*resTerm;

        resTerm = wEst[index] - ( wExac[index] - tau * (
          (wExac[index+offsetZ]*wExac[index+offsetZ] - wExac[index-offsetZ]*wExac[index-offsetZ]) / (2.0*hZ) + 
          (uExac[index+1]*wExac[index+1] - uExac[index-1]*wExac[index-1]) / (2.0*hX) + 
          (vExac[index+offsetY]*wExac[index+offsetY] - vExac[index-offsetY]*wExac[index-offsetY]) / (2.0*hY) +
          (pExac[index+offsetZ] - pExac[index-offsetZ]) / (2.0*hZ) -
          ( (wExac[index+1] - 2*wExac[index] + wExac[index-1]) / (hX*hX) +
            (wExac[index+offsetY] - 2*wExac[index] + wExac[index-offsetY]) / (hY*hY) +
            (wExac[index+offsetZ] - 2*wExac[index] + wExac[index-offsetZ]) / (hZ*hZ)
          ) / params.Reyn
        ) );
        vectorResidual += resTerm*resTerm;
      }
    }
  }
  //-------------------------------------------------------------------------
  }
  return std::sqrt(vectorResidual);
}

double conc_velocity_residual_FDA2_4(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac)
{
  const auto dimSize ( params.domainPartition );
  const uint offsetY = dimSize + 1;
  const uint offsetZ = (dimSize+1) * (dimSize+1);

  const double tau ( params.duration / params.timePartition ); // time step
  const double hX (params.xLen / dimSize );
  const double hY (params.yLen / dimSize );
  const double hZ (params.zLen / dimSize );

  double vectorResidual (0.0);

  int numProcs = omp_get_num_procs();
  omp_set_num_threads(numProcs);
  #pragma omp parallel
  {
  //---------------------------- inner knots --------------------------------
  #pragma omp for collapse(2) schedule(static)
  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index (k*offsetZ + j*offsetY + i);
        double resTerm (0.0); // residual term
        //! Insert precise values to the scheme, take the difference with the estimated values
        resTerm = uEst[index] - ( uExac[index] - tau * (
            uExac[index]*(uExac[index+1]-uExac[index-1])/(2*hX) 
            + vExac[index]*(uExac[index+offsetY]-uExac[index-offsetY])/(2*hY)
            + wExac[index]*(uExac[index+offsetZ]-uExac[index-offsetZ])/(2*hZ) 
            + (pExac[index+1]-pExac[index-1])/(2*hX)
            - ( 
              (uExac[index+1] - 2*uExac[index] + uExac[index-1]) / (hX*hX) +
              (uExac[index+offsetY] - 2*uExac[index] + uExac[index-offsetY]) / (hY*hY) +
              (uExac[index+offsetZ] - 2*uExac[index] + uExac[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
        ) );
        vectorResidual += resTerm*resTerm;

        resTerm = vEst[index] - ( vExac[index] - tau * (
            uExac[index]*(vExac[index+1]-vExac[index-1])/(2*hX) 
            + vExac[index]*(vExac[index+offsetY]-vExac[index-offsetY])/(2*hY)
            + wExac[index]*(vExac[index+offsetZ]-vExac[index-offsetZ])/(2*hZ) 
            + (pExac[index+offsetY]-pExac[index-offsetY])/(2*hY)
            - ( 
              (vExac[index+1] - 2*vExac[index] + vExac[index-1]) / (hX*hX) +
              (vExac[index+offsetY] - 2*vExac[index] + vExac[index-offsetY]) / (hY*hY) +
              (vExac[index+offsetZ] - 2*vExac[index] + vExac[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
        ) );
        vectorResidual += resTerm*resTerm;

        resTerm = wEst[index] - ( wExac[index] - tau * (
            uExac[index]*(wExac[index+1]-wExac[index-1])/(2*hX) 
            + vExac[index]*(wExac[index+offsetY]-wExac[index-offsetY])/(2*hY)
            + wExac[index]*(wExac[index+offsetZ]-wExac[index-offsetZ])/(2*hZ) 
            + (pExac[index+offsetZ]-pExac[index-offsetZ])/(2*hZ)
            - ( 
              (wExac[index+1] - 2*wExac[index] + wExac[index-1]) / (hX*hX) +
              (wExac[index+offsetY] - 2*wExac[index] + wExac[index-offsetY]) / (hY*hY) +
              (wExac[index+offsetZ] - 2*wExac[index] + wExac[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
        ) );
        vectorResidual += resTerm*resTerm;
      }
    }
  }
  //-------------------------------------------------------------------------
  }
  return std::sqrt(vectorResidual);
}
//=================================================================================================