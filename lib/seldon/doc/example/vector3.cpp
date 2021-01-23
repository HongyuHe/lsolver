#define SELDON_DEBUG_LEVEL_4

#include "SeldonFlag.hxx"
#include "Seldon.hxx"

using namespace Seldon;

#include "vector/Vector3.cxx"

int main()
{

  TRY;

  Vector<Vector<int> > length(2);
  length(0).Reallocate(2);
  length(0)(0) = 4;
  length(0)(1) = 5;
  length(1).Reallocate(3);
  length(1)(0) = 7;
  length(1)(1) = 2;
  length(1)(2) = 3;
  Vector3<double> V(length);

  // Fills all inner vectors with 2.
  V.Fill(2.);
  // Access to the second inner vector, to fill it with 5.
  V(1).Fill(5.);
  // Access to an inner vector at lowest level, and fill it
  V(0, 1).Fill(3.0);
  V.Print();
  cout << "First vector of the second inner vector: " << V(1, 0) << endl;
  // Note that V(1)(0) would have returned the same element.
  
  Vector<double> inner_vector(4);
  inner_vector.Fill();
  // Appends a new inner vector to second vector.
  V.PushBack(1, inner_vector);
  V.Print();

  cout << "After setting to -10 the second vector of the last inner vector:"
       << endl;
  V(1, 1) = -10.;
  V.Print();

  END;

  return 0;

}
