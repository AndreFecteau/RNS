#ifndef GRID1D_H
#define GRID1D_H

template <typename scalar_type_in, typename size_type_in, typename global_solution_vector_type_in, typename matrix_type_in>
struct Grid1D {

 public:

  using scalar_type = scalar_type_in;
  using size_type = size_type_in;
  using global_solution_vector_type = global_solution_vector_type_in;
  using matrix_type = matrix_type_in;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Grid1D() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Grid1D(const Grid1D&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Grid1D(Grid1D&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Grid1D& operator=(const Grid1D&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Grid1D& operator=(Grid1D&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor setting up required inputs.
  Grid1D(scalar_type x_min_in, scalar_type x_max_in,
          global_solution_vector_type global_solution_vector_in) :
          x_min(x_min_in), x_max(x_max_in),
          global_solution_vector(global_solution_vector_in) {}

  scalar_type x_min;
  scalar_type x_max;
  global_solution_vector_type global_solution_vector;

  size_type number_of_cells() const {return global_solution_vector.size();}
  scalar_type domaine_length() const {return x_max - x_min;}
  scalar_type dx() const {return domaine_length() / static_cast<scalar_type>(number_of_cells());}
  scalar_type per_FL() const {return number_of_cells() / (x_max - x_min);}

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(x_min, x_max, global_solution_vector);
  }

};


#endif //#ifndef NON_DIMENSIONAL_NAVIER_STOKES_H
