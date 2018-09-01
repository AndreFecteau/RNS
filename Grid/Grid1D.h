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
  Grid1D( scalar_type x_min_in, scalar_type x_max_in,
          global_solution_vector_type global_solution_vector_in) :
          x_min(x_min_in), x_max(x_max_in),
          global_solution_vector(global_solution_vector_in) {}

  Grid1D( scalar_type x_min_in, scalar_type x_max_in,
          scalar_type number_of_cells_in) :
          x_min(x_min_in), x_max(x_max_in),
          global_solution_vector(global_solution_vector_type(number_of_cells_in,
          global_solution_vector_type::value_type::Zero())) {}

  scalar_type x_min;
  scalar_type x_max;
  global_solution_vector_type global_solution_vector;

  size_type number_of_cells() const {return global_solution_vector.size();}
  scalar_type domaine_length() const {return x_max - x_min;}
  scalar_type dx() const {return domaine_length() / static_cast<scalar_type>(number_of_cells());}
  size_type per_FL() const {return number_of_cells() / (x_max - x_min);}

  void add_space_in_back(scalar_type space);
  void add_space_in_front(scalar_type space);
  void refine(size_type per_FL);

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(x_min, x_max, global_solution_vector);
  }

};

///////////////////////////////////////////////////////////////////////////////
// add_space_in_back
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type_in, typename size_type_in, typename global_solution_vector_type_in, typename matrix_type_in>
void Grid1D<scalar_type_in, size_type_in, global_solution_vector_type_in, matrix_type_in>::
add_space_in_back(scalar_type space) {
  global_solution_vector.resize(number_of_cells() + space/dx());
  x_max += space;
  for(size_t i = number_of_cells() - space/dx(); i < number_of_cells(); ++i) {
    global_solution_vector[i] = global_solution_vector[i-1];
  }
  std::cout << "Added space in domaine, New domaine size: " << domaine_length() << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// add_space_in_front
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type_in, typename size_type_in, typename global_solution_vector_type_in, typename matrix_type_in>
void Grid1D<scalar_type_in, size_type_in, global_solution_vector_type_in, matrix_type_in>::
add_space_in_front(scalar_type space) {
  global_solution_vector_type temp_global_solution_vector;
  temp_global_solution_vector.resize(number_of_cells() + space/dx());
  for(size_t i = 0; i < space/dx(); ++i){
    temp_global_solution_vector[i] = global_solution_vector[0];
  }
  for(size_t i = space/dx(); i < temp_global_solution_vector.size(); ++i) {
    temp_global_solution_vector[i] = global_solution_vector[i - space/dx()];
  }
  x_max += space;
  global_solution_vector = temp_global_solution_vector;
  std::cout << "Added space in domaine, New domaine size: " << domaine_length() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// refine
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type_in, typename size_type_in, typename global_solution_vector_type_in, typename matrix_type_in>
void Grid1D<scalar_type_in, size_type_in, global_solution_vector_type_in, matrix_type_in>::
refine(size_type per_flame_length) {
  global_solution_vector_type temp_global_solution_vector;
  if(per_flame_length == per_FL()){
    temp_global_solution_vector = global_solution_vector;
  } else if(per_flame_length > per_FL()){
    if(per_flame_length % per_FL() != 0) {
      std::cout << "cant refine to asked amount" << std::endl;
    } else {
      temp_global_solution_vector.resize(domaine_length() * per_flame_length);
      for(size_type i = 0; i < number_of_cells(); ++i) {
        for(size_type j = 0; j < per_flame_length / per_FL(); ++j) {
          temp_global_solution_vector[i*per_flame_length / per_FL()+j] = global_solution_vector[i];
        }
      }
    }
  } else if(per_flame_length < per_FL()) {
    if(per_FL() % per_flame_length != 0) {
      std::cout << "cant refine to asked amount" << std::endl;
    } else {
      // global_solution_vector_type temp_global_solution_vector;
      temp_global_solution_vector.resize(domaine_length() * per_flame_length);
      for(size_type i = 0; i < temp_global_solution_vector.size(); ++i) {
        temp_global_solution_vector[i] = global_solution_vector[i * per_FL() / per_flame_length];
      }
      temp_global_solution_vector[0] = global_solution_vector[0];
      temp_global_solution_vector[temp_global_solution_vector.size()-1] = global_solution_vector[number_of_cells()-1];
    }
  } else {
      std::cout << "Something went wrong in refinement" << std::endl;
  }
  global_solution_vector = temp_global_solution_vector;
  std::cout << "Number_cells: " << number_of_cells() << std::endl;
}


#endif //#ifndef NON_DIMENSIONAL_NAVIER_STOKES_H
