#ifndef TEXT_OUTPUT_H
#define TEXT_OUTPUT_H

#include <ostream>

namespace Color {
  enum Code {
      FG_RED      = 31,
      FG_GREEN    = 32,
      FG_BLUE     = 34,
      FG_BLACK     = 30,
      FG_DEFAULT  = 39,
  };

class Modifier {
    Code code;
public:
    Modifier(Code pCode) : code(pCode) {}
    friend std::ostream&
    operator<<(std::ostream& os, const Modifier& mod) {
        return os << "\033[" << mod.code << "m";
    }
};
}

void loading_text(){
  Color::Modifier red(Color::FG_RED);
  Color::Modifier black(Color::FG_BLACK);
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);
  std::cout << black <<
"//////////////////////////////////////////////////////////////////////////////////////////////////\n"
"//////////////////////////////////////////////////////////////////////////////////////////////////\n"
 << red <<
"\n"
"                        |||||||||        ||          ||        ||||||||                           \n"
"                        ||       ||      ||||        ||      ||        ||                         \n"
"                        ||       ||      ||  ||      ||      ||                                   \n"
"                        |||||||||        ||    ||    ||        ||||||||                           \n"
"                        ||   ||          ||      ||  ||                ||                         \n"
"                        ||     ||        ||        ||||      ||        ||                         \n"
"                        ||       ||      ||          ||        ||||||||                           \n"
 << black <<
"Things are only impossible until they are not. - Jean-Luc Picard                     Andre Fecteau\n"
"//////////////////////////////////////////////////////////////////////////////////////////////////\n"
"//////////////////////////////////////////////////////////////////////////////////////////////////"
  << def << std::endl;
}

template <typename scalar_type>
void lambda_itteration_print(scalar_type lambda){
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);
  std::cout << green << "\rThe Eigen Value Current Guess: " << def << lambda << std::flush;
}

template <typename scalar_type>
void lambda_final_print(scalar_type lambda){
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier def(Color::FG_DEFAULT);
  Color::Modifier black(Color::FG_BLACK);
  std::cout << green << "\rThe Eigen Value From Shooting Method: " << def << lambda << std::endl;
  std::cout << black <<
"//////////////////////////////////////////////////////////////////////////////////////////////////"
  << def << std::endl;
}
#endif
