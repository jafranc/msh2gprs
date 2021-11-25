#include <vector>
#include <string>

namespace gprs_data
{

void invert_matrix_3x3(const std::vector<double> &A, std::vector<double> &A_inv);

double determinant_3x3(const std::vector<double> &m);

std::vector<double> transpose3x3(const std::vector<double> &mat);

std::string str_toupper(std::string s);


}  // end namespace gprs_data
