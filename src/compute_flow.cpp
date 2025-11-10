#include "D:\Projects\C++\ns_fda\headers\compute_flow.hpp"

//==================================== FLOW COMPUTATION ===========================================
void compute(
  model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p)
{
  uint tick = 0;
  const auto dimSize ( params.domainPartition );
  uint offsetY = dimSize + 1;
  uint offsetZ = (dimSize+1) * (dimSize+1);

  const double tau ( params.duration / params.timePartition ); // time step
  const double hX (params.xLen / dimSize );
  const double hY (params.yLen / dimSize );
  const double hZ (params.zLen / dimSize );

  const size_t vecSize = (dimSize + 1) * (dimSize + 1) * (dimSize + 1);
  std::vector<double> u1(vecSize);
  std::vector<double> v1(vecSize);
  std::vector<double> w1(vecSize);
  
  while (tick < params.timePartition + 1)
  {
    //! velocity compute
    //---------------------------- inner knots --------------------------------
    for (auto k = 1; k < dimSize; k++)      // Z-Axis
    {
      for (auto j = 1; j < dimSize; j++)    // Y-Axis
      {
        for (auto i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index (k*offsetZ + j*offsetY + i);

          u1[index] = u[index] - tau * (
            (u[index+1]*u[index+1] - u[index-1]*u[index-1]) / (2*hX) + 
            (u[index+offsetY]*v[index+offsetY] - u[index-offsetY]*v[index-offsetY]) / (2*hY) + 
            (u[index+offsetZ]*w[index+offsetZ] - u[index-offsetZ]*w[index-offsetZ]) / (2*hZ) +
            (p[index+1] - p[index-1]) / (2*hX) -
            ( (u[index+1] - 2*u[index] + u[index-1]) / (hX*hX) +
              (u[index+offsetY] - 2*u[index] + u[index-offsetY]) / (hY*hY) +
              (u[index+offsetZ] - 2*u[index] + u[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
          );
          v1[index] = v[index] - tau * (
            (v[index+offsetY]*v[index+offsetY] - v[index-offsetY]*v[index-offsetY]) / (2*hY) + 
            (u[index+1]*v[index+1] - u[index-1]*v[index-1]) / (2*hX) + 
            (v[index+offsetZ]*w[index+offsetZ] - v[index-offsetZ]*w[index-offsetZ]) / (2*hZ) +
            (p[index+offsetY] - p[index-offsetY]) / (2*hY) -
            ( (v[index+1] - 2*v[index] + v[index-1]) / (hX*hX) +
              (v[index+offsetY] - 2*v[index] + v[index-offsetY]) / (hY*hY) +
              (v[index+offsetZ] - 2*v[index] + v[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
          );
          w1[index] = w[index] - tau * (
            (w[index+offsetZ]*w[index+offsetZ] - w[index-offsetZ]*w[index-offsetZ]) / (2*hZ) + 
            (u[index+1]*w[index+1] - u[index-1]*w[index-1]) / (2*hX) + 
            (v[index+offsetY]*w[index+offsetY] - v[index-offsetY]*w[index-offsetY]) / (2*hY) +
            (p[index+offsetZ] - p[index-offsetZ]) / (2*hZ) -
            ( (w[index+1] - 2*w[index] + w[index-1]) / (hX*hX) +
              (w[index+offsetY] - 2*w[index] + w[index-offsetY]) / (hY*hY) +
              (w[index+offsetZ] - 2*w[index] + w[index-offsetZ]) / (hZ*hZ)
            ) / params.Reyn
          );
        }
      }
    }
    //-------------------------------------------------------------------------

    //---------------------------- border knots -------------------------------
    //! XY-plane Z = 0 / Z = MAX; Neiman's condition dw/dz = 0
    for (auto j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (auto i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (j*offsetY + i);
        uint index2 (dimSize*offsetZ + j*offsetY + i);
        u1[index1] = v1[index1] = 0.0;
        u1[index2] = v1[index2] = 0.0;
        w1[index1] = w1[index1 + offsetZ];
        w1[index2] = w1[index2 - offsetZ];
      }
    }
    //! XZ-plane Y = 0 / Y = MAX; Neiman's condition dv/dy = 0
    for (auto k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (auto i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (k*offsetZ + i);
        uint index2 (k*offsetZ + dimSize*offsetY + i);
        u1[index1] = w1[index1] = 0.0;
        u1[index2] = w1[index2] = 0.0;
        v1[index1] = v1[index1 + offsetY];
        v1[index2] = v1[index2 - offsetY];
      }
    }
    //! YZ-plane X = 0 / X = MAX; Neiman's condition du/dx = 0
    for (auto k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (auto j = 1; j < dimSize; j++)  // X-Axis
      {
        uint index1 (k*offsetZ + j*offsetY);
        uint index2 (k*offsetZ + j*offsetY + dimSize);
        v1[index1] = w1[index1] = 0.0;
        v1[index2] = w1[index2] = 0.0;
        u1[index1] = u1[index1 + 1];
        u1[index2] = u1[index2 - 1];
      }
    }
    //-------------------------------------------------------------------------

    //---------------------------- edge knots ---------------------------------
    //! Y = Z = (0 || MAX); 
    for (auto i = 1; i < dimSize; i++)
    {
      // Y = Z = 0
      uint index1 (i);
      // Y = MAX, Z = 0
      uint index2 (dimSize*offsetY + i);
      // Y = 0, Z = MAX
      uint index3 (dimSize*offsetZ + i);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

      u1[index1] = u1[index2] = u1[index3] = u1[index4] = 0.0;

      v1[index1] = v1[index1 + offsetY];
      v1[index2] = v1[index2 - offsetY];
      v1[index3] = v1[index3 + offsetY];
      v1[index4] = v1[index4 - offsetY];

      w1[index1] = w1[index1 + offsetZ];
      w1[index2] = w1[index2 + offsetZ];
      w1[index3] = w1[index3 - offsetZ];
      w1[index4] = w1[index4 - offsetZ];

    }
    //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
    for (auto j = 1; j < dimSize; j++)
    {
      // X = Z = 0
      uint index1 (j*offsetY);
      // X = MAX, Z = 0
      uint index2 (j*offsetY + dimSize);
      // X = 0, Z = MAX
      uint index3 (dimSize*offsetZ + j*offsetY);
      // X = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + j*offsetY + dimSize);

      v1[index1] = v1[index2] = v1[index3] = v1[index4] = 0.0;

      u1[index1] = u1[index1 + 1];
      u1[index2] = u1[index2 - 1];
      u1[index3] = u1[index3 + 1];
      u1[index4] = u1[index4 - 1];

      w1[index1] = w1[index1 + offsetZ];
      w1[index2] = w1[index2 + offsetZ];
      w1[index3] = w1[index3 - offsetZ];
      w1[index4] = w1[index4 - offsetZ];
    }
    //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
    for (auto k = 1; k < dimSize; k++)
    {
      // X = Y = 0
      uint index1 (k*offsetZ);
      // X = MAX, Y = 0
      uint index2 (k*offsetZ + dimSize);
      // X = 0, Y = MAX
      uint index3 (k*offsetZ + dimSize*offsetY);
      // X = MAX, Y = MAX
      uint index4 (k*offsetZ + dimSize*offsetY + dimSize);

      w1[index1] = w1[index2] = w1[index3] = w1[index4] = 0.0;

      u1[index1] = u1[index1 + 1];
      u1[index2] = u1[index2 - 1];
      u1[index3] = u1[index3 + 1];
      u1[index4] = u1[index4 - 1];

      v1[index1] = v1[index1 + offsetY];
      v1[index2] = v1[index2 + offsetY];
      v1[index3] = v1[index3 - offsetY];
      v1[index4] = v1[index4 - offsetY];
    }
    //-------------------------------------------------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 0
    uint index (0);
    u1[index] = u1[index + 1];
    v1[index] = v1[index + offsetY];
    w1[index] = w1[index + offsetZ];
    // X = MAX, Y = Z = 0
    index = dimSize;
    u1[index] = u1[index - 1];
    v1[index] = v1[index + offsetY];
    w1[index] = w1[index + offsetZ];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    u1[index] = u1[index + 1];
    v1[index] = v1[index - offsetY];
    w1[index] = w1[index + offsetZ];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    u1[index] = u1[index + 1];
    v1[index] = v1[index + offsetY];
    w1[index] = w1[index - offsetZ];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    u1[index] = u1[index - 1];
    v1[index] = v1[index - offsetY];
    w1[index] = w1[index + offsetZ];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    u1[index] = u1[index - 1];
    v1[index] = v1[index + offsetY];
    w1[index] = w1[index - offsetZ];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    u1[index] = u1[index + 1];
    v1[index] = v1[index - offsetY];
    w1[index] = w1[index - offsetZ];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    u1[index] = u1[index - 1];
    v1[index] = v1[index - offsetY];
    w1[index] = w1[index - offsetZ];
    //-------------------------------------------------------------------------
    // move data from upper time layer
    u.clear(); u = std::move(u1);
    v.clear(); v = std::move(v1);
    w.clear(); w = std::move(w1);
    u1.resize(vecSize); v1.resize(vecSize); w1.resize(vecSize);
    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //! pressure compute
    //---------------------------- inner knots --------------------------------
    for (auto k = 1; k < dimSize; k++)      // Z-Axis
    {
      for (auto j = 1; j < dimSize; j++)    // Y-Axis
      {
        for (auto i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index (k*offsetZ + j*offsetY + i);
          // matrix construct
          // call solver
        }
      }
    }
    //-------------------------------------------------------------------------

    //---------------------------- border knots -------------------------------
    //! XY-plane Z = 0 / Z = MAX; Neiman's condition dp/dz = 0
    for (auto j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (auto i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (j*offsetY + i);
        uint index2 (dimSize*offsetZ + j*offsetY + i);
        p[index1] = p[index1 + offsetZ];
        p[index2] = p[index2 - offsetZ];
      }
    }
    //! XZ-plane Y = 0 / Y = MAX; Neiman's condition dp/dy = 0
    for (auto k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (auto i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (k*offsetZ + i);
        uint index2 (k*offsetZ + dimSize*offsetY + i);
        p[index1] = p[index1 + offsetY];
        p[index2] = p[index2 - offsetY];
      }
    }
    //! YZ-plane X = 0 / X = MAX; Neiman's condition dp/dx = 0
    for (auto k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (auto j = 1; j < dimSize; j++)  // X-Axis
      {
        uint index1 (k*offsetZ + j*offsetY);
        uint index2 (k*offsetZ + j*offsetY + dimSize);
        p[index1] = p[index1 + 1];
        p[index2] = p[index2 - 1];
      }
    }
    //-------------------------------------------------------------------------

    //---------------------------- edge knots ---------------------------------
    //! Y = Z = (0 || MAX); 
    for (auto i = 1; i < dimSize; i++)
    {
      // Y = Z = 0
      uint index1 (i);
      // Y = MAX, Z = 0
      uint index2 (dimSize*offsetY + i);
      // Y = 0, Z = MAX
      uint index3 (dimSize*offsetZ + i);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

      p[index1] = p[index1 + offsetY];
      p[index2] = p[index2 - offsetY];
      p[index3] = p[index3 + offsetY];
      p[index4] = p[index4 - offsetY];

    }
    //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
    for (auto j = 1; j < dimSize; j++)
    {
      // X = Z = 0
      uint index1 (j*offsetY);
      // X = MAX, Z = 0
      uint index2 (j*offsetY + dimSize);
      // X = 0, Z = MAX
      uint index3 (dimSize*offsetZ + j*offsetY);
      // X = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + j*offsetY + dimSize);

      p[index1] = p[index1 + 1];
      p[index2] = p[index2 - 1];
      p[index3] = p[index3 + 1];
      p[index4] = p[index4 - 1];
    }
    //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
    for (auto k = 1; k < dimSize; k++)
    {
      // X = Y = 0
      uint index1 (k*offsetZ);
      // X = MAX, Y = 0
      uint index2 (k*offsetZ + dimSize);
      // X = 0, Y = MAX
      uint index3 (k*offsetZ + dimSize*offsetY);
      // X = MAX, Y = MAX
      uint index4 (k*offsetZ + dimSize*offsetY + dimSize);

      p[index1] = p[index1 + 1];
      p[index2] = p[index2 - 1];
      p[index3] = p[index3 + 1];
      p[index4] = p[index4 - 1];
    }
    //-------------------------------------------------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 0
    index = 0;
    p[index] = p[index + 1];
    // X = MAX, Y = Z = 0
    index = dimSize;
    p[index] = p[index - 1];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    p[index] = p[index + 1];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    p[index] = p[index + 1];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    p[index] = p[index - 1];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    p[index] = p[index - 1];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    p[index] = p[index + 1];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    p[index] = p[index - 1];
    //-------------------------------------------------------------------------
    tick += 1;
  }
  std::cout << "final tick: " << tick << '\n';
}
//=================================================================================================