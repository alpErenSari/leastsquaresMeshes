#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <vector>
#include <cmath>
#include <random>
#include "fiboheap.h"
#define MAX 100000


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;

typedef std::pair<double, int> dij_Pair;

struct dij_out
{
  double min_dist;
  std::vector<int> prev;
  std::vector<double> dist;
  std::vector<int> parent;
};

double l2_norm(Eigen::MatrixXd &V, long int i, long int j)
{
  return sqrt(pow(V(i,0) - V(j,0), 2) + pow(V(i,1) - V(j,1), 2) +
    pow(V(i,2) - V(j,2), 2));
}

dij_out dijkstra_fibo_heap(std::vector<std::vector<dij_Pair>> adj_vec, int src, int des)
{
      int V_m = adj_vec.size();
      // std::cout << "V_m is " << V_m << std::endl;
     std::vector<double> dist(V_m);     // The output array.  dist[i] will hold the shortest
     std::vector<int> prev(V_m);
     std::vector<bool> sptSet(V_m); // sptSet[i] will be true if vertex i is included in shortest

                      // distance from src to i
     std::vector<int> parent;
     FibHeap <dij_Pair> pq;

                     // path tree or shortest distance from src to i is finalized
     // Distance of source vertex from itself is always 0

     // Initialize all distances as INFINITE and stpSet[] as false
     for (int i = 0; i < V_m; i++)
     {
       if(i == src)
       {
         pq.push(std::make_pair(0, i));
         dist[i] = 0;
       }
       else
       {
         pq.push(std::make_pair(INT_MAX, i)), dist[i] = INT_MAX;
       }
       prev[i] = -1;
     }

     // Find shortest path for all vertices
     while (!pq.empty())
     {
       // std::cout << "count number " << count << std::endl;
       // Pick the minimum distance vertex from the set of vertices not
       // yet processed. u is always equal to src in the first iteration.
       int u = pq.top().second;
       pq.pop();
       // std::cout << "min vertex is " << u << std::endl;
       // if(u == des)
       //   break;

       std::vector<dij_Pair>::iterator i;
       // Update dist value of the adjacent vertices of the picked vertex.
       for (i = adj_vec[u].begin(); i != adj_vec[u].end(); i++)
       {
         // Update dist[v] only if is not in sptSet, there is an edge from
         // u to v, and total weight of path from src to  v through u is
         // smaller than current value of dist[v]
         int v = (*i).second;
         double weight = (*i).first;

         // If there is shorted path to v through u.
         if (dist[v] > dist[u] + weight)
         {
             // Updating distance of v
             dist[v] = dist[u] + weight;
             pq.push(std::make_pair(dist[v], v));
             prev[v] = u;
         }
       }
     }

     int u = des;
     if(prev[u]>=0 || u==src)
     {
       while(u>=0)
       {
         parent.push_back(u);
         u = prev[u];
       }
     }

     dij_out retval;
     retval.min_dist = dist[des];
     retval.prev = prev;
     retval.dist = dist;
     retval.parent = parent;

     return retval;
}

std::vector<int> furthest_sample(std::vector<std::vector<double>> &A, int n_samples, int m_v){
  // this the FPS part

  std::vector<std::vector<dij_Pair>> my_pair(A.size());
  for(size_t i=0; i<A.size(); i++)
  {
    for(size_t j=0; j<A[i].size(); j++)
    {
      long int k = A[i][j];
      double weight = l2_norm(V, i, k);
      // V_cost(i, k) = weight;
      my_pair[i].push_back(std::make_pair(weight, k));
    }
  }

  std::vector<bool> isVertexUsed(m_v);
  Eigen::MatrixXd fps_geo_dist(n_samples, m_v);
  for (size_t i = 0; i < m_v; i++) {
    isVertexUsed[i] = false;
  }
  isVertexUsed[0] = true;
  std::vector<int> fps(n_samples);
  fps[0] = 0;
  for (size_t i = 1; i < n_samples; i++) {
    // int len_try = m_v-i;
    double max = 0;
    int max_place = 0;
    dij_out my_dij = dijkstra_fibo_heap(my_pair, fps[i-1], 0);
    for (size_t j = 0; j < m_v; j++) {
      fps_geo_dist(i-1, j) = my_dij.dist[j];
    }
    // std::vector<double> distances(len_try);
    for (size_t j = 0; j < m_v; j++) {
      if(isVertexUsed[j])
        continue;

      double min_inner = INT_MAX;
      for (size_t k = 0; k < i; k++) {
        // dij_out my_dij = dijkstra(V_cost, fps[k], j);
        size_t vert_index = fps[k];
        double curr_dist = fps_geo_dist(k, j);
        if(curr_dist < min_inner)
          min_inner = curr_dist;
      }
      if(max < min_inner)
      {
        max = min_inner;
        max_place = j;
      }

    }
    fps[i] = max_place;
    isVertexUsed[max_place] = true;
    std::cout << "The iteration number is " << i << '\n';
    std::cout << "The selected vertex is " << max_place << '\n';
  }

  return fps;
}



int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  if(argc >= 2)
    igl::readOFF(argv[1], V, F);
  else
  {
    std::cerr << "You need to enter the mesh file path!" << '\n';
    std::exit(1);
  }

  int n_samples = 100;

  try
  {
    if(argc >= 3)
      n_samples = std::stoi(argv[2]);
  }
  catch(...)
  {
    std::cout << "Method input cannot be converted to int" <<
    std::endl;
  }


  std::cout << "V size is " << V.rows() << '\n';
  std::cout << "F size is " << F.rows() << '\n';

  // // Find the open boundary
  // Eigen::VectorXi bnd;
  // igl::boundary_loop(F,bnd);
  //
  // std::cout << "bnd size is " << bnd.size() << '\n';

  // std::cout << "bnd vertices are" << '\n';
  // for (size_t i = 0; i < bnd.size(); i++) {
  //   std::cout << bnd[i] << " " << '\n';
  // }
  // std::cout << '\n';

  // Eigen::MatrixXd P(bnd.size(), 3);
  // for (size_t i = 0; i < bnd.size(); i++) {
  //   P.row(i) = V.row(bnd(i));
  // }

  // Map the boundary to a circle, preserving edge proportions
  // Eigen::MatrixXd bnd_uv;
  // igl::map_vertices_to_circle(V,bnd,bnd_uv);
  // boundry_map_circle(V,bnd,bnd_uv);
  //
  // std::cout << "bnd_uv size is " << bnd_uv.rows() << '\n';

  std::random_device rd; // obtain a random number from hardware
  std::mt19937 eng(rd()); // seed the gen#include<Eigen/SparseCholesky>erator
  std::uniform_int_distribution<> distr(0, V.rows()); // define the range

// Eigen::MatrixXd L = Eigen::MatrixXd::Zero(V.rows(), V.rows());

// Eigen::SparseMatrix<double> L(V.rows(), V.rows());
Eigen::SparseMatrix<double> A_sp(V.rows() + n_samples, V.rows());
  // find the adjacency list
std::vector<std::vector<double>> A;
igl::adjacency_list(F,A, true);



// std::vector<Eigen::Triplet<double>> tripletList;
// tripletList.reserve(estimation_of_entries);
std::cout << "Creating L matrix" << '\n';
for (size_t i = 0; i < A.size(); i++) {
  A_sp.coeffRef(i,i) = 1.0;
  float A_size = A[i].size();
  for (size_t k = 0; k < A[i].size(); k++) {
    int j = int(A[i][k]);
    A_sp.coeffRef(i,j) = -1/A_size;
  }
}

std::vector<int> f_list(n_samples);

// for (size_t i = 0; i < n_samples; i++) {
//   int new_vertex = distr(eng);
//   bool is_vertex_element = true;
//   while (is_vertex_element) {
//     is_vertex_element = false;
//     for (size_t j = 0; j < i; j++) {
//       if(f_list[j]==new_vertex) {
//         is_vertex_element = true;
//         new_vertex = distr(eng);
//       }
//     }
//   }
//   f_list[i] = new_vertex;
// }

// compute the all permutations
std::vector<bool> v(V.rows());
std::fill(v.end() - n_samples, v.end(), true);

double min_rec_cost = MAX;
int count_main_loop = 0;
Eigen::MatrixXd V_lse(V.rows(), V.cols());

do {
  int count_p = 0;
  std::cout << "Count number is " << count_main_loop << '\n';
    for (int i = 0; i < V.rows(); ++i) {
        if (v[i]) {
            // std::cout << (i + 1) << " ";
            f_list[count_p] = i;
            count_p++;
        }
    }
    // this part is about computing the matrix
    for (size_t i = 0; i < n_samples; i++) {
      long long int j = f_list[i];
      A_sp.coeffRef(i + V.rows(),j) = 1.0;
    }

    A_sp.makeCompressed();

    Eigen::VectorXd x_s = Eigen::VectorXd::Zero(V.rows() + n_samples);
    Eigen::VectorXd y_s = Eigen::VectorXd::Zero(V.rows() + n_samples);
    Eigen::VectorXd z_s = Eigen::VectorXd::Zero(V.rows() + n_samples);
    for (size_t i = 0; i < n_samples; i++) {
      long long int j = f_list[i];
      x_s(i + V.rows()) = V(j, 0);
      y_s(i + V.rows()) = V(j, 1);
      z_s(i + V.rows()) = V(j, 2);
    }


    Eigen::VectorXd x_res, y_res, z_res;
    // Eigen::SparseMatrix<double> AT(A_sp.transpose());
    // std::cout << "A_sp size is " << A_sp.rows() << " " << A_sp.cols() << '\n';
    // std::cout << "AT size is " << AT.rows() << " " << AT.cols() << '\n';
    // Eigen::SparseMatrix<double> ATA(AT * A_sp);
    // Eigen::SparseMatrix<double> ATA_Inv(ATA.inverse());
    // Eigen::COLAMDOrdering<int>
    std::cout << "Computing A_sp at step " << count_main_loop << '\n';
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.compute(A_sp);
    if(solver.info() != Eigen::Success) {
        // decomposition failed
        std::cerr << "Decompostion failed for sparse matrix S_p" << '\n';
        return -1;
      }
    // std::cout << "Computing x_res" << '\n';
    x_res = solver.solve(x_s);
    // std::cout << "Computing y_res" << '\n';
    y_res = solver.solve(y_s);
    // std::cout << "Computing z_res" << '\n';
    z_res = solver.solve(z_s);




    // Compute cost of reconstruction
    Eigen::MatrixXd delta_V = V - V_lse;
    double rec_cost = 0;
    int m_v = delta_V.rows();
    for (size_t i = 0; i < delta_V.rows(); i++) {
      Eigen::VectorXd delta_v_i = delta_V.row(i);
      rec_cost += delta_v_i.norm();
    }
    rec_cost /= delta_V.rows();
    if(rec_cost < min_rec_cost)
    {
      min_rec_cost = rec_cost;
      V_lse.col(0) = x_res;
      V_lse.col(1) = y_res;
      V_lse.col(2) = z_res;
    }
    count_main_loop++;
    // std::cout << "Total reconstruction loss is " << rec_cost << '\n';
} while (std::next_permutation(v.begin(), v.end()));


// get ready for fps
// std::vector<int> f_list = furthest_sample(A, n_samples, V.rows());

// for (size_t i = 0; i < f_list.size(); i++) {
//   if(i == 0) {std::cout << "Random vertices are :" << '\n';}
//   else if(i == f_list.size()-1) {std::cout << '\n';}
//   else {std::cout << f_list[i] << " ";}
// }

// Eigen::SparseMatrix<double> F_sample(n_samples, V.rows());
// insert 1's to selected sample vertices' places


  // Eigen::VectorXd x_res, y_res;
  // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  // x_res =solver.compute(W).solve(b.col(0));
  // y_res =solver.compute(W).solve(b.col(1));
  // x_res = W.colPivHouseholderQr().solve(b.col(0));
  // y_res = W.colPivHouseholderQr().solve(b.col(1));
  //
  // Eigen::MatrixXd V_uv_uni(V.rows(), 2);
  // if(method_no==1)
  //   parametrization_uniform(V, bnd, bnd_uv, A, V_uv_uni);
  // else if(method_no==2)
  //   parametrization_harmonic_2(V, bnd, bnd_uv, A, V_uv_uni);
  // else if(method_no==3)
  //   parametrization_harmonic_3(V, bnd, bnd_uv, A, V_uv_uni);
  // else
  // {
  //   std::cerr << "Invalid method number!" << '\n';
  //   std::exit(1);
  // }
  // V_uv_uni.col(0) = x_res;
  // V_uv_uni.col(1) = y_res;
  // V_uv_uni =  W.colPivHouseholderQr().solve(b);

  // Harmonic parametrization for the internal vertices
  // igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

  // Scale UV to make the texture more clear
  // V_uv *= 5;
  // V_uv_uni *= 5;
  // bnd_uv *= 5;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_lse, F);
  // viewer.data().set_uv(V_uv);
  // viewer.data().add_points(bnd_uv, Eigen::RowVector3d(0,0,1));
  // viewer.callback_key_down = &key_down;

  // Disable wireframe
  // viewer.data().show_lines = false;

  // Draw checkerboard texture
  // viewer.data().show_texture = true;

  // Launch the viewer
  viewer.launch();
}
