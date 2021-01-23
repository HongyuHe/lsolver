// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_VECTOR_CXX

#include "Vector.hxx"

namespace Seldon
{

  //! Changes the length of the vector, and keeps previous values.
  /*!
    Reallocates the vector to size i. Previous values are kept.
    \param n new length of the vector.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::Resize(size_t n)
  {
    ResizeVector(n);
  }


  //! Changes the length of the vector, and keeps previous values.
  /*!
    Reallocates the vector to size i. Previous values are kept.
    \param n new length of the vector.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::ResizeVector(size_t n)
  {
    // function implemented in the aim that explicit specialization
    // of Resize can call ResizeVector
    if (n == this->m_)
      return;

    Vector<T, VectFull, Allocator> X_new(n);
    for (size_t i = 0; i < min(this->m_, n); i++)
      X_new(i) = this->data_[i];

    SetData(n, X_new.GetData());
    X_new.Nullify();
  }
  
  
  /*********
   * NORMS *
   *********/


  //! Returns the infinite norm.
  /*!
    \return The infinite norm.
  */
  template <class T, class Allocator>
  typename ClassComplexType<T>::Treal
  Vector<T, VectFull, Allocator>::GetNormInf() const
  {
    typename ClassComplexType<T>::Treal res(0);
    for (size_t i = 0; i < this->GetLength(); i++)
      res = max(res, ComplexAbs(this->data_[i]));
    
    return res;
  }


  //! Returns the index of the highest absolute value.
  /*!
    \return The index of the element that has the highest absolute value.
  */
  template <class T, class Allocator>
  long Vector<T, VectFull, Allocator>::GetNormInfIndex() const
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (this->GetLength() == 0)
      throw WrongDim("Vector<VectFull>::GetNormInfIndex()",
		     "Vector is null.");
#endif

    typename ClassComplexType<T>::Treal res(0), temp;
    long j = 0;
    for (size_t i = 0; i < this->GetLength(); i++)
      {
	temp = res;
	res = max(res, ComplexAbs(this->data_[i]));
	if (temp != res) j = i;
      }
    
    return j;
  }


  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/


  //! Writes the vector in a file.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::Write(string FileName, "
                    "bool with_size)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Write(ostream& FileStream, "
                    "bool with_size)",
                    "The stream is not ready.");
#endif
    int n = this->m_;
    if (with_size)
      FileStream.write(reinterpret_cast<char*>(&n),
		       sizeof(int));
    
    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Write(ostream& FileStream, "
                    "bool with_size)",
                    "Output operation failed.");
#endif

  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    // end of line to finish the file
    FileStream << endl;

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    if (this->GetLength() != 0)
      FileStream << (*this)(0);

    for (size_t i = 1; i < this->GetLength(); i++)
      FileStream << "\t" << (*this)(i);

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


#ifdef SELDON_WITH_HDF5
  //! Writes the vector in a HDF5 file.
  /*!
    All elements of the vector are stored in HDF5 format.
    \param FileName file name.
    \param group_name name of the group the vector must be stored in.
    \param dataset_name name of the dataset the vector must be stored in.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::WriteHDF5(string FileName,
                                                 string group_name,
                                                 string dataset_name) const
  {
    hid_t file_id = H5Fopen(FileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!file_id)
      throw IOError("Vector<VectFull>::WriteHDF5(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    hid_t dataset_id, dataspace_id, group_id;
    herr_t status;

    T x;
    SetComplexZero(x);
    hid_t datatype = GetH5Type(x);

    if (!H5Lexists(file_id, group_name.c_str(), H5P_DEFAULT))
      group_id = H5Gcreate(file_id, group_name.c_str(), 0);
    group_id = H5Gopen(file_id, group_name.c_str());

    if (!H5Lexists(group_id, dataset_name.c_str(), H5P_DEFAULT))
      {
        // Create the initial dataspace.
        hsize_t dim[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        dataspace_id = H5Screate_simple(1, dim, maxdims);

        // Define chunking parameters to allow extension of the dataspace.
        hid_t cparms = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {this->GetLength()};
        status = H5Pset_chunk(cparms, 1, chunk_dims);

        // Create the dataset.
        hid_t filetype = H5Tvlen_create(datatype);
        dataset_id = H5Dcreate(group_id, dataset_name.c_str(),
                               filetype, dataspace_id, cparms);
      }

    // Opens the dataset, and extend it to store a new vector.
    dataset_id = H5Dopen(group_id, dataset_name.c_str());
    dataspace_id = H5Dget_space(dataset_id);
    hsize_t dims_out[1];
    status  = H5Sget_simple_extent_dims(dataspace_id, dims_out, NULL);
    hsize_t new_dim[1]= {dims_out[0] + 1};
    status = H5Dextend(dataset_id, new_dim);

    // Selects the last memory part of the dataset, to store the vector.
    hsize_t offset[1] = {dims_out[0]};
    hsize_t dim2[1] = {1};
    dataspace_id = H5Dget_space(dataset_id);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET,
                                 offset, NULL,
                                 dim2, NULL);

    // Defines memory space.
    hid_t dataspace = H5Screate_simple(1, dim2, NULL);

    // Writes the data to the memory hyperslab selected.
    hid_t memtype = H5Tvlen_create(datatype);
    hvl_t wdata[1];
    wdata[0].len = this->GetLength();
    wdata[0].p = this->GetDataVoid();
    status = H5Dwrite(dataset_id, memtype, dataspace,
                      dataspace_id, H5P_DEFAULT, wdata);

    // Closes the dataset, group and file.
    status = H5Dclose(dataset_id);
    status = H5Gclose(group_id);
    status = H5Fclose(file_id);
  }
#endif


  //! Sets the vector from a file.
  /*!
    Sets the vector according to a binary file that stores the length of the
    vector (integer) and all elements.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not
    available in the file. In this case, the current size N of the vector is
    unchanged, and N elements are read in the file.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::Read(string FileName, bool with_size)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream, with_size);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets the vector according to a binary file stream that stores the length
    of the vector (integer) and all elements.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not
    available in the stream. In this case, the current size N of the vector is
    unchanged, and N elements are read in the stream.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Read(istream& FileStream, "
                    "bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int new_size;
        FileStream.read(reinterpret_cast<char*>(&new_size), sizeof(int));
        this->Reallocate(new_size);
      }

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    this->GetLength() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Read(istream& FileStream, "
                    "bool with_size)",
                    "Output operation failed.");
#endif

  }


  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::ReadText(istream& FileStream)
  {
    // Previous values of the vector are cleared.
    Clear();

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::ReadText(istream& FileStream)",
                    "The stream is not ready.");
#endif

    T entry;
    size_t number_element = 0;
    while (!FileStream.eof())
      {
	// Reads a new entry.
	FileStream >> entry;
	if (FileStream.fail())
	  break;
	else
	  {
	    number_element++;

	    // If needed, resizes the vector. Its size is already doubled so
	    // that the vector should be resized a limited number of times.
	    if (number_element > this->m_)
	      this->Resize(2 * number_element);

	    this->data_[number_element - 1] = entry;
	  }
      }

    // Resizes to the actual size.
    if (number_element > 0)
      this->Resize(number_element);
    else
      this->Clear();
  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Storage, Allocator>& V)
  {
    if (V.GetLength() >= 1)
      for (size_t i = 0; i < V.GetLength() - 1; i++)
        out << V(i) << '\t';
    
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_VECTOR_CXX
#endif
