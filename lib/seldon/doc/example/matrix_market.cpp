#define SELDON_DEBUG_LEVEL_2

#include "SeldonFlag.hxx"
#include "Seldon.hxx"

using namespace Seldon;

#include "matrix_sparse/IOMatrixMarket.cxx"

int main(int argc, char** argv)
{

  if (argc != 2)
    {
      cout << "Enter the file name where the matrix is stored" << endl;
      cout << "Usage : ./a.out matrice.rua" << endl;
      abort();
    }
  
  TRY;

  Matrix<double, General, ColSparse> A;

  ReadHarwellBoeing(argv[1], A);

  A.Print();

  WriteHarwellBoeing(A, "result.rua");

  END;

  return 0;

}
