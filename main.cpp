#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <vector>
#include <cmath>


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;

typedef std::pair<size_t, size_t> harm_pair;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Plot the 3D mesh
    viewer.data().set_mesh(V,F);
    viewer.core.align_camera_center(V,F);
  }
  else if (key == '2')
  {
    // Plot the mesh in 2D using the UV coordinates as vertex coordinates
    viewer.data().set_mesh(V_uv,F);
    viewer.core.align_camera_center(V_uv,F);
  }

  viewer.data().compute_normals();

  return false;
}

bool check_eigen_vector(Eigen::VectorXi &bnd, int key)
{
  for (size_t i = 0; i < bnd.size(); i++) {
    if(bnd(i) == key)
      return true;
  }
  return false;
}

bool find_2_common_vertices(std::vector<double> &src, std::vector<double> &dst, harm_pair &comm_2)
{
  std::vector<double> common;
  for (size_t i = 0; i < src.size(); i++) {
    for (size_t j = 0; j < dst.size(); j++) {
      if(src[i] == dst[j])
      {
        common.push_back(src[i]);
      }
    }
  }
  if(common.size() >= 2)
  {
    comm_2.first = common[0];
    comm_2.second = common[1];
    return true;
  }
  else
    false;
}

double cot_2_vectors(Eigen::Vector3d &vec_1, Eigen::Vector3d &vec_2)
{
  double dot = vec_1.dot(vec_2);   //between [x1, y1, z1] and [x2, y2, z2]
  // Eigen::Vector3d cross;
  // cross = vec_1.cross(vec_2);
  // double cross_norm = cross.norm();
  double v_1_norm = vec_1.norm();
  double v_2_norm = vec_2.norm();
  double angle = acos(dot/(v_1_norm*v_2_norm));

  return 1/tan(angle);

}

double tan_over_2_vectors(Eigen::Vector3d &vec_1, Eigen::Vector3d &vec_2)
{
  double dot = vec_1.dot(vec_2);   //between [x1, y1, z1] and [x2, y2, z2]
  // Eigen::Vector3d cross;
  // cross = vec_1.cross(vec_2);
  // double cross_norm = cross.norm();
  double v_1_norm = vec_1.norm();
  double v_2_norm = vec_2.norm();
  double angle = acos(dot/(v_1_norm*v_2_norm));

  return tan(angle/2);

}

void boundry_map_circle(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd,
  Eigen::MatrixXd& bnd_uv)
{
  // Get sorted list of boundary vertices
  std::vector<int> interior,map_ij;
  map_ij.resize(V.rows());

  std::vector<bool> isOnBnd(V.rows(),false);
  for (int i = 0; i < bnd.size(); i++)
  {
    isOnBnd[bnd[i]] = true;
    map_ij[bnd[i]] = i;
  }

  for (int i = 0; i < (int)isOnBnd.size(); i++)
  {
    if (!isOnBnd[i])
    {
      map_ij[i] = interior.size();
      interior.push_back(i);
    }
  }

  // Map boundary to unit circle
  std::vector<double> circle_len(bnd.size());
  circle_len[0] = 0.;

  for (int i = 1; i < bnd.size(); i++)
  {
    circle_len[i] = circle_len[i-1] + (V.row(bnd[i-1]) - V.row(bnd[i])).norm();
  }
  double total_circle_len = circle_len[circle_len.size()-1] + (V.row(bnd[0]) - V.row(bnd[bnd.size()-1])).norm();

  bnd_uv.resize(bnd.size(),2);
  for (int i = 0; i < bnd.size(); i++)
  {
    double len_frac = circle_len[i] * 2. * igl::PI / total_circle_len;
    bnd_uv.row(map_ij[bnd[i]]) << cos(len_frac), sin(len_frac);
  }

}

void parametrization_uniform(Eigen::MatrixXd &V, Eigen::VectorXi &bnd,
  Eigen::MatrixXd &bnd_uv, std::vector<std::vector<double>> &A,
  Eigen::MatrixXd &V_uv_uni)
{
  Eigen::MatrixXd W;
  W = Eigen::MatrixXd::Zero(V.rows(), V.rows());

  // boundry mask
  Eigen::MatrixXd bnd_mask ;
  bnd_mask = Eigen::MatrixXd::Zero(V.rows(), 3);

  for (size_t i = 0; i < bnd.size(); i++) {
    bnd_mask(bnd(i), 0) = 1;
    bnd_mask(bnd(i), 1) = bnd_uv(i,0);
    bnd_mask(bnd(i), 2) = bnd_uv(i,1);
  }

  Eigen::MatrixXd b;
  b = Eigen::MatrixXd::Zero(V.rows(), 2);

  for (size_t i = 0; i < V.rows(); i++) {
    if(bnd_mask(i, 0) == 1)
    {
      W(i,i) = 1;
      b(i, 0) = bnd_mask(i, 1);
      b(i, 1) = bnd_mask(i, 2);
    }
    else
    {
      std::vector<double> temp_adj = A[i];
      int temp_len = temp_adj.size();
      // std::cout << "Diagonal value is " << temp_len << '\n';
      for (size_t j = 0; j < temp_adj.size(); j++) {
        int k =  int(temp_adj[j]);
        W(i, k) = 1;
      }
      W(i,i) = -1*temp_len;
      // std::cout << "Inserted element is " << W.coeffRef(i,i) << '\n';
    }
  }
  V_uv_uni = W.colPivHouseholderQr().solve(b);

}

void parametrization_harmonic_2(Eigen::MatrixXd &V, Eigen::VectorXi &bnd,
  Eigen::MatrixXd &bnd_uv, std::vector<std::vector<double>> &A,
  Eigen::MatrixXd &V_uv_uni)
{
  Eigen::MatrixXd W;
  W = Eigen::MatrixXd::Zero(V.rows(), V.rows());

  // boundry mask
  Eigen::MatrixXd bnd_mask ;
  bnd_mask = Eigen::MatrixXd::Zero(V.rows(), 3);

  for (size_t i = 0; i < bnd.size(); i++) {
    bnd_mask(bnd(i), 0) = 1;
    bnd_mask(bnd(i), 1) = bnd_uv(i,0);
    bnd_mask(bnd(i), 2) = bnd_uv(i,1);
  }

  Eigen::MatrixXd b;
  b = Eigen::MatrixXd::Zero(V.rows(), 2);

  for (size_t i = 0; i < V.rows(); i++) {
    if(bnd_mask(i, 0) == 1)
    {
      W(i,i) = 1;
      b(i, 0) = bnd_mask(i, 1);
      b(i, 1) = bnd_mask(i, 2);
    }
    else
    {
      std::vector<double> temp_adj = A[i];
      int temp_len = temp_adj.size();
      // std::cout << "Diagonal value is " << temp_len << '\n';
      for (size_t j = 0; j < temp_adj.size(); j++) {
        int k =  int(temp_adj[j]);
        // W(i, k) = 1;
        // method 2 beginning
        harm_pair temp_pair;
        if(find_2_common_vertices(A[i], A[k], temp_pair))
        {
          // std::cout << "i and k are " << i << " " << k << '\n';
          // std::cout << "Pair values are " << temp_pair.first << " " << temp_pair.second << '\n';
          Eigen::Vector3d vec_1 = V.row(i) - V.row(temp_pair.first);
          Eigen::Vector3d vec_2 = V.row(k) - V.row(temp_pair.first);
          Eigen::Vector3d vec_3 = V.row(i) - V.row(temp_pair.second);
          Eigen::Vector3d vec_4 = V.row(k) - V.row(temp_pair.second);
          double cot_1 = cot_2_vectors(vec_1, vec_2);
          double cot_2 = cot_2_vectors(vec_3, vec_4);
          // std::cout << "cot values are " << cot_1 << " " << cot_2 << '\n';
          W(i, k) = (cot_1 + cot_2)/2;
          W(i,i) -= (cot_1 + cot_2)/2;
        }
      }
    }
  }
  V_uv_uni = W.colPivHouseholderQr().solve(b);
}

void parametrization_harmonic_3(Eigen::MatrixXd &V, Eigen::VectorXi &bnd,
  Eigen::MatrixXd &bnd_uv, std::vector<std::vector<double>> &A,
  Eigen::MatrixXd &V_uv_uni)
{
  Eigen::MatrixXd W;
  W = Eigen::MatrixXd::Zero(V.rows(), V.rows());

  // boundry mask
  Eigen::MatrixXd bnd_mask ;
  bnd_mask = Eigen::MatrixXd::Zero(V.rows(), 3);

  for (size_t i = 0; i < bnd.size(); i++) {
    bnd_mask(bnd(i), 0) = 1;
    bnd_mask(bnd(i), 1) = bnd_uv(i,0);
    bnd_mask(bnd(i), 2) = bnd_uv(i,1);
  }

  Eigen::MatrixXd b;
  b = Eigen::MatrixXd::Zero(V.rows(), 2);

  for (size_t i = 0; i < V.rows(); i++) {
    if(bnd_mask(i, 0) == 1)
    {
      W(i,i) = 1;
      b(i, 0) = bnd_mask(i, 1);
      b(i, 1) = bnd_mask(i, 2);
    }
    else
    {
      std::vector<double> temp_adj = A[i];
      int temp_len = temp_adj.size();
      // std::cout << "Diagonal value is " << temp_len << '\n';
      for (size_t j = 0; j < temp_adj.size(); j++) {
        int k =  int(temp_adj[j]);
        // W(i, k) = 1;
        // method 2 beginning
        harm_pair temp_pair;
        if(find_2_common_vertices(A[i], A[k], temp_pair))
        {
          // std::cout << "i and k are " << i << " " << k << '\n';
          // std::cout << "Pair values are " << temp_pair.first << " " << temp_pair.second << '\n';
          Eigen::Vector3d vec_1 = V.row(temp_pair.first) - V.row(k);
          Eigen::Vector3d vec_2 = V.row(i) - V.row(k);
          Eigen::Vector3d vec_3 = V.row(temp_pair.second) - V.row(k);

          double tan_1 = tan_over_2_vectors(vec_1, vec_2);
          double tan_2 = tan_over_2_vectors(vec_2, vec_3);
          // std::cout << "cot values are " << tan_1 << " " << tan_2 << '\n';

          W(i, k) = (tan_1 + tan_2)/(2*vec_2.norm());
          W(i,i) -= (tan_1 + tan_2)/(2*vec_2.norm());
        }
      }
    }
  }
  V_uv_uni = W.colPivHouseholderQr().solve(b);
}




int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../../faces/face-low.off", V, F);

  std::cout << "V size is " << V.rows() << '\n';
  std::cout << "F size is " << F.rows() << '\n';

  // Find the open boundary
  Eigen::VectorXi bnd;
  igl::boundary_loop(F,bnd);

  std::cout << "bnd size is " << bnd.size() << '\n';

  Eigen::MatrixXd P(bnd.size(), 3);
  for (size_t i = 0; i < bnd.size(); i++) {
    P.row(i) = V.row(bnd(i));
  }

  // Map the boundary to a circle, preserving edge proportions
  Eigen::MatrixXd bnd_uv;
  // igl::map_vertices_to_circle(V,bnd,bnd_uv);
  boundry_map_circle(V,bnd,bnd_uv);

  std::cout << "bnd_uv size is " << bnd_uv.rows() << '\n';

  // Eigen::MatrixXd P_uv(bnd_uv.size(), 3);
  // for (size_t i = 0; i < 10; i++) {
  //   P_uv.row(i) = V.row(bnd_uv(i));
  // }

  // Eigen::SparseMatrix<double> W(V.rows(), V.rows());
  // Eigen::MatrixXd W;
  // W = Eigen::MatrixXd::Zero(V.rows(), V.rows());
  //
  // // boundry mask
  // Eigen::MatrixXd bnd_mask ;
  // bnd_mask = Eigen::MatrixXd::Zero(V.rows(), 3);
  //
  // for (size_t i = 0; i < bnd.size(); i++) {
  //   bnd_mask(bnd(i), 0) = 1;
  //   bnd_mask(bnd(i), 1) = bnd_uv(i,0);
  //   bnd_mask(bnd(i), 2) = bnd_uv(i,1);
  // }

  // find the adjacency list
  std::vector<std::vector<double>> A;
  igl::adjacency_list(F,A, true);
  // Eigen::MatrixXd b;
  // b = Eigen::MatrixXd::Zero(V.rows(), 2);
  //
  // for (size_t i = 0; i < A[10].size(); i++) {
  //   std::cout << A[10][i] << " ";
  // }
  // std::cout << " end of A[0]" << '\n';
  // std::cout << "size of A " << A.size() << '\n';
  //
  // for (size_t i = 0; i < V.rows(); i++) {
  //   if(bnd_mask(i, 0) == 1)
  //   {
  //     W(i,i) = 1;
  //     b(i, 0) = bnd_mask(i, 1);
  //     b(i, 1) = bnd_mask(i, 2);
  //   }
  //   else
  //   {
  //     std::vector<double> temp_adj = A[i];
  //     int temp_len = temp_adj.size();
  //     // std::cout << "Diagonal value is " << temp_len << '\n';
  //     for (size_t j = 0; j < temp_adj.size(); j++) {
  //       int k =  int(temp_adj[j]);
  //       W(i, k) = 1;
  //     }
  //     W(i,i) = -1*temp_len;
  //     // std::cout << "Inserted element is " << W.coeffRef(i,i) << '\n';
  //   }
  // }

  // std::cout << "Non zero element number is " << W.nonZeros() << '\n';
  //
  // for (size_t i = 0; i < V.rows(); i++) {
  //   if(bnd_mask(i, 0) == 1)
  //   {
  //     W(i,i) = 1;
  //     b(i, 0) = bnd_mask(i, 1);
  //     b(i, 1) = bnd_mask(i, 2);
  //   }
  //   else
  //   {
  //     std::vector<double> temp_adj = A[i];
  //     int temp_len = temp_adj.size();
  //     // std::cout << "Diagonal value is " << temp_len << '\n';
  //     for (size_t j = 0; j < temp_adj.size(); j++) {
  //       int k =  int(temp_adj[j]);
  //       // W(i, k) = 1;
  //       // method 2 beginning
  //       harm_pair temp_pair;
  //       if(find_2_common_vertices(A[i], A[k], temp_pair))
  //       {
  //         std::cout << "i and k are " << i << " " << k << '\n';
  //         std::cout << "Pair values are " << temp_pair.first << " " << temp_pair.second << '\n';
  //         Eigen::Vector3d vec_1 = V.row(i) - V.row(temp_pair.first);
  //         Eigen::Vector3d vec_2 = V.row(k) - V.row(temp_pair.first);
  //         Eigen::Vector3d vec_3 = V.row(i) - V.row(temp_pair.second);
  //         Eigen::Vector3d vec_4 = V.row(k) - V.row(temp_pair.second);
  //         double cot_1 = cot_2_vectors(vec_1, vec_2);
  //         double cot_2 = cot_2_vectors(vec_3, vec_4);
  //         std::cout << "cot values are " << cot_1 << " " << cot_2 << '\n';
  //         W(i, k) = (cot_1 + cot_2)/2;
  //         W(i,i) -= (cot_1 + cot_2)/2;
  //       }
  //     }
  //   }
  // }

  // for (size_t i = 0; i < V.rows(); i++) {
  //   if(bnd_mask(i, 0) == 1)
  //   {
  //     W(i,i) = 1;
  //     b(i, 0) = bnd_mask(i, 1);
  //     b(i, 1) = bnd_mask(i, 2);
  //   }
  //   else
  //   {
  //     std::vector<double> temp_adj = A[i];
  //     int temp_len = temp_adj.size();
  //     // std::cout << "Diagonal value is " << temp_len << '\n';
  //     for (size_t j = 0; j < temp_adj.size(); j++) {
  //       int k =  int(temp_adj[j]);
  //       // W(i, k) = 1;
  //       // method 2 beginning
  //       harm_pair temp_pair;
  //       if(find_2_common_vertices(A[i], A[k], temp_pair))
  //       {
  //         std::cout << "i and k are " << i << " " << k << '\n';
  //         std::cout << "Pair values are " << temp_pair.first << " " << temp_pair.second << '\n';
  //         Eigen::Vector3d vec_1 = V.row(temp_pair.first) - V.row(k);
  //         Eigen::Vector3d vec_2 = V.row(i) - V.row(k);
  //         Eigen::Vector3d vec_3 = V.row(temp_pair.second) - V.row(k);
  //
  //         double tan_1 = tan_over_2_vectors(vec_1, vec_2);
  //         double tan_2 = tan_over_2_vectors(vec_2, vec_3);
  //         std::cout << "cot values are " << tan_1 << " " << tan_2 << '\n';
  //
  //         W(i, k) = (tan_1 + tan_2)/(2*vec_2.norm());
  //         W(i,i) -= (tan_1 + tan_2)/(2*vec_2.norm());
  //       }
  //     }
  //   }
  // }



  Eigen::VectorXd x_res, y_res;
  // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  // x_res =solver.compute(W).solve(b.col(0));
  // y_res =solver.compute(W).solve(b.col(1));
  // x_res = W.colPivHouseholderQr().solve(b.col(0));
  // y_res = W.colPivHouseholderQr().solve(b.col(1));
  //
  Eigen::MatrixXd V_uv_uni(V.rows(), 2);
  parametrization_harmonic_3(V, bnd, bnd_uv, A, V_uv_uni);
  // V_uv_uni.col(0) = x_res;
  // V_uv_uni.col(1) = y_res;
  // V_uv_uni =  W.colPivHouseholderQr().solve(b);

  // Harmonic parametrization for the internal vertices
  igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

  // Scale UV to make the texture more clear
  V_uv *= 5;
  V_uv_uni *= 5;
  bnd_uv *= 5;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_uv_uni, F);
  // viewer.data().set_uv(V_uv);
  viewer.data().add_points(bnd_uv, Eigen::RowVector3d(0,0,1));
  // viewer.callback_key_down = &key_down;

  // Disable wireframe
  // viewer.data().show_lines = false;

  // Draw checkerboard texture
  // viewer.data().show_texture = true;

  // Launch the viewer
  viewer.launch();
}
