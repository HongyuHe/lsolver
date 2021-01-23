// Copyright (C) 2001-2009 Vivien Mallet
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_ERRORS_CXX

#include "Errors.hxx"

namespace Seldon
{


  ///////////
  // ERROR //
  ///////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  Error::Error(string function, string comment):
    description_("ERROR!\nAn undefined error occurred"),
    function_(function), comment_(comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Alternative constructor.
  /*! Error associated with a description, a function and a comment.
    \param description short description of the error.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  Error::Error(string description, string function, string comment):
    description_("ERROR!\n" + description),
    function_(function), comment_(comment)
  {
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  /*!
    \note Empty.
  */
  Error::~Error()
  {
  }


  /***********
   * METHODS *
   ***********/


  //! Delivers information about the error.
  /*! Displays available information, i.e. the error description, the function
    and/or a comment.
  */
  string Error::What()
  {
    string message(description_);
    if (!function_.empty())
      message += " in " + function_;
    message += ".\n";
    if (!comment_.empty())
      message += "   " + comment_;
    return message;
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e. the error description, the function
    and/or a comment.
  */
  void Error::CoutWhat()
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
    cout << "Error on processor " << rank_proc << endl;
#endif

    cout << this->What() << endl;

#ifdef SELDON_WITH_MPI
    // waiting 5 seconds so that each processor should have the time
    // to display its own errors before aborting the processus
    //sleep(5);

    // other approach : imposing a Barrier
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }



  ///////////////
  // UNDEFINED //
  ///////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param[in] function function with which the error is associated.
    \param[in] comment comment associated with the error.
  */
  Undefined::Undefined(string function, string comment):
  Error("", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  Undefined::~Undefined()
  {
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  string Undefined::What()
  {
    string message;
    if (!this->function_.empty())
      message = this->function_;
    else
      message = "A function or a method";
    message += " is undefined.\nEither its implementation is missing,"
      + string(" or it does not make sense or it is impossible ")
      + "to implement it.\n";
    if (!this->comment_.empty())
      message += "   " + this->comment_;
    return message;
  }


  ///////////////////
  // WRONGARGUMENT //
  ///////////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongArgument::WrongArgument(string function, string comment):
  Error("Wrong argument given to ", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  WrongArgument::~WrongArgument()
  {
  }
  
  
  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  string WrongArgument::What()
  {
    string message(this->description_);
    if (!this->function_.empty())
      message += this->function_;
    message += ".\n";
    if (!this->comment_.empty())
      message += "   " + this->comment_;
    return message;
  }


  //////////////
  // NOMEMORY //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  NoMemory::NoMemory(string function, string comment):
    Error("Out of memory", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  NoMemory::~NoMemory()
  {
  }


  //////////////
  // WRONGDIM //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongDim::WrongDim(string function, string comment):
    Error("Wrong dimensions involved", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }

  
  //! Destructor (present to avoid duplicate symbols in .o files)
  WrongDim::~WrongDim()
  {
  }


  ////////////////
  // WRONGINDEX //
  ////////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongIndex::WrongIndex(string function, string comment):
    Error("Index out of range", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  WrongIndex::~WrongIndex()
  {
  }


  //////////////
  // WRONGROW //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongRow::WrongRow(string function, string comment):
    Error("Row index out of range", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  WrongRow::~WrongRow()
  {
  }


  //////////////
  // WRONGCOL //
  //////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  WrongCol::WrongCol(string function, string comment):
    Error("Column index out of range", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  WrongCol::~WrongCol()
  {
  }


  /////////////
  // IOERROR //
  /////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  IOError::IOError(string function, string comment):
    Error("Error while performing an I/O operation", function, comment)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  IOError::~IOError()
  {
  }

  
  //////////////////
  // other errors //
  //////////////////


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  SolverMaximumIterationError::SolverMaximumIterationError(string function, string comment):
    Error("Maximum number of iterations reached by the solver", function, comment)
  {
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  SolverMaximumIterationError::~SolverMaximumIterationError()
  {
  }


  //! Main constructor.
  /*! Error associated with both a function and a comment.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  SolverDivergenceError::SolverDivergenceError(string function, string comment):
    Error("Maximum number of iterations reached by the solver", function, comment)
  {
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  SolverDivergenceError::~SolverDivergenceError()
  {
  }
  

  /////////////////
  // CheckBounds //
  /////////////////


  void CheckBounds(long i, long length1_, string nom)
  {
    if (i < 0 || i >= length1_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
  }


  void CheckBounds(int i, int j, int length1_, int length2_, string nom)
  {
    if (i < 0 || i >= length1_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length2_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
  }


  void CheckBoundsSym(int i, int j, int length1_, int length2_, string nom)
  {
    CheckBounds(i, j, length1_, length2_, nom);
    
    if (i > j)
      throw WrongRow(nom+"::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
  }


  void CheckBoundsTriang(int i, int j, int length1_, int length2_,
			 bool uplo, string nom)
  {
    CheckBounds(i, j, length1_, length2_, nom);
    
    if (uplo)
      {
	if (i > j)
	  throw WrongRow(nom + "::Val(int, int)",
			 string("Attempted to access to element (")
			 + to_str(i) + ", " + to_str(j) + string(") but row")
			 + string(" index should not be strictly more")
			 + " than column index (upper triangular matrix).");
      }
    else
      {
	if (j > i)
	  throw WrongCol(nom + "::Val(int, int)",
			 string("Attempted to access to element (")
			 + to_str(i) + ", " + to_str(j) + string(") but")
			 + string(" column index should not be strictly more")
			 + " than row index (lower triangular matrix).");
      }
  }


  void CheckBounds(int i, int j, int k,
		   int length1_, int length2_, int length3_, string nom)
  {
    if (i < 0 || i >= length1_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length2_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length3_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length3_-1) + "], but is equal to "
		       + to_str(k) + ".");
  }


  void CheckBounds(int i, int j, int k, int l, int length1_, int length2_,
		   int length3_, int length4_, string nom)
  {
    if (i < 0 || i >= length1_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length2_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length3_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length3_-1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length4_)
      throw WrongIndex(nom + "::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length4_-1) + "], but is equal to "
		       + to_str(l) + ".");
  }


  void CheckBounds(int i, int j, int k, int l, int m,
		   int length0, int length1, int length2, int length3,
		   int length4, string nom)
  {
    if (i < 0 || i >= length0)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length0 - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length1)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length1 - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length2)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length2 - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length3)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length3 - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length4)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length4 - 1) + "], but is equal to "
		       + to_str(m) + ".");
  }
  

  void CheckBounds(int i, int j, int k, int l, int m, int n,
		   int length0, int length1, int length2, int length3,
		   int length4, int length5, string nom)
  {
    if (i < 0 || i >= length0)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length0 - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length1)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length1 - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length2)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length2 - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length3)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length3 - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length4)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length4 - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length5)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length5 - 1) + "], but is equal to "
		       + to_str(n) + ".");
  }


  void CheckBounds(int i, int j, int k, int l, int m, int n, int o,
		   int length0, int length1, int length2, int length3,
		   int length4, int length5, int length6, string nom)
  {
    if (i < 0 || i >= length0)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length0 - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length1)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length1 - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length2)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length2 - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length3)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length3 - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length4)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length4 - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length5)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length5 - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length6)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length6 - 1) + "], but is equal to "
		       + to_str(o) + ".");
  }


  void CheckBounds(int i, int j, int k, int l, int m, int n, int o, int p,
		   int length0, int length1, int length2, int length3,
		   int length4, int length5, int length6, int length7,
		   string nom)
  {
    if (i < 0 || i >= length0)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length0 - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length1)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length1 - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length2)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length2 - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length3)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length3 - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length4)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length4 - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length5)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length5 - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length6)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length6 - 1) + "], but is equal to "
		       + to_str(o) + ".");
    if (p < 0 || p >= length7)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #8 should be in [0, ")
		       + to_str(length7 - 1) + "], but is equal to "
		       + to_str(p) + ".");
  }


  void CheckBounds(int i, int j, int k, int l, int m, int n, int o, int p,
		   int q, int length0, int length1, int length2, int length3,
		   int length4, int length5, int length6, int length7,
		   int length8, string nom)
  {
    if (i < 0 || i >= length0)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length0 - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length1)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length1 - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length2)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length2 - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length3)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length3 - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length4)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length4 - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length5)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length5 - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length6)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length6 - 1) + "], but is equal to "
		       + to_str(o) + ".");
    if (p < 0 || p >= length7)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #8 should be in [0, ")
		       + to_str(length7 - 1) + "], but is equal to "
		       + to_str(p) + ".");
    if (q < 0 || q >= length8)
      throw WrongIndex(nom+"::operator()",
		       string("Index along dimension #9 should be in [0, ")
		       + to_str(length8 - 1) + "], but is equal to "
		       + to_str(q) + ".");
  }
  
  
  /////////////////
  // LAPACKERROR //
  /////////////////


  //! Main constructor.
  /*! Error associated with a diagnostic integer, a function and a comment.
    \param info Lapack diagnostic integer.
    \param function function with which the error is associated.
    \param comment comment associated with the error.
  */
  LapackError::LapackError(int info, string function, string comment):
  Error("Error returned by Lapack", function, comment), info_(info)
  {
#ifdef SELDON_WITH_ABORT
    this->CoutWhat();
    abort();
#endif
  }


  //! Destructor (present to avoid duplicate symbols in .o files)
  LapackError::~LapackError()
  {
  }


  //! Delivers information about the error.
  /*! Displays available information, i.e.
    the error description, the function and/or the comment.
  */
  string LapackError::What()
  {
    string message(this->description_);
    if (!this->function_.empty())
      message += " in " + this->function_;
    message += ".\n";
    if (!this->comment_.empty())
      message += "   " + this->comment_;
    message += "   Diagnostic integer (\"info\"): " + to_str(info_)
      + ".\n";
    return message;
  }


  ////////////////
  // LAPACKINFO //
  ////////////////


  LapackInfo::LapackInfo(int info): info_(info)
  {
  }


  LapackInfo::operator int ()
  {
    return info_;
  }


  int LapackInfo::GetInfo()
  {
    return info_;
  }


  int& LapackInfo::GetInfoRef()
  {
    return info_;
  }


#ifndef SELDON_WITH_COMPILED_LIBRARY
  LapackInfo lapack_info(0);
#endif


} // namespace Seldon.

#define SELDON_FILE_ERRORS_CXX
#endif
