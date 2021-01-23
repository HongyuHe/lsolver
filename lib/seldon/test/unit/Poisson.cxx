template<class MatrixSparse>
void ComputeHelmholtzMatrix(MatrixSparse& A, double dx, int Nx, int Ny, double w2)
{
  double coef = 0.25*w2*dx*dx;
  for (int i = 0; i < Nx-1; i++)
    for (int j = 0; j < Ny-1; j++)
      {
	int num_row = j*Nx + i;
	// interaction between node (i, j) and other nodes
        A.AddInteraction(num_row, num_row, 1.0-coef);
        A.AddInteraction(num_row, num_row+1, -0.5);
        A.AddInteraction(num_row, num_row+Nx, -0.5);
        	
	// interaction between node (i+1, j) and other nodes
        A.AddInteraction(num_row+1, num_row+1, 1.0-coef);
        A.AddInteraction(num_row+1, num_row, -0.5);
        A.AddInteraction(num_row+1, num_row+Nx+1, -0.5);

	// interaction between node (i+1, j+1) and other nodes
        A.AddInteraction(num_row+Nx+1, num_row+Nx+1, 1.0-coef);
        A.AddInteraction(num_row+Nx+1, num_row+1, -0.5);
        A.AddInteraction(num_row+Nx+1, num_row+Nx, -0.5);
        
	// interaction between node (i, j+1) and other nodes
        A.AddInteraction(num_row+Nx, num_row+Nx, 1.0-coef);
        A.AddInteraction(num_row+Nx, num_row, -0.5);
        A.AddInteraction(num_row+Nx, num_row+Nx+1, -0.5);
      }
}

void DistributePoissonDofs(int nb_proc_row, int nb_proc_col, int nb_points_glob, MPI_Comm comm,
			   IVect& GlobalRowNumbers, IVect& ListProcInterac,
			   Vector<IVect>& ListRowInterac, IVect& OverlapRowNumbers, IVect& OriginalProcNumbers,
			   int& i0_proc, int& i1_proc, int& j0_proc, int& j1_proc)
{
  int num_proc_row = rank_processor/nb_proc_col;
  int num_proc_col = rank_processor%nb_proc_col;
  //DISP(num_proc_row); DISP(num_proc_col);
  int nb_row = nb_points_glob/nb_proc_row;
  i0_proc = num_proc_row*nb_row;
  i1_proc = (num_proc_row+1)*nb_row+1;
  if (num_proc_row == nb_proc_row-1)
    i1_proc = nb_points_glob;

  j0_proc = num_proc_col*nb_row;
  j1_proc = (num_proc_col+1)*nb_row+1;
  if (num_proc_col == nb_proc_col-1)
    j1_proc = nb_points_glob;

  //DISP(i0_proc); DISP(i1_proc); DISP(j0_proc); DISP(j1_proc);
  
  int N = (i1_proc-i0_proc)*(j1_proc-j0_proc);
  //DISP(N); DISP(dx);
  


  // global row numbers
  GlobalRowNumbers.Reallocate(N);
  int Nx = i1_proc-i0_proc, Ny = j1_proc-j0_proc;
  for (int i = i0_proc; i < i1_proc; i++)
    for (int j = j0_proc; j < j1_proc; j++)
      GlobalRowNumbers((j-j0_proc)*Nx + i-i0_proc) = j*nb_points_glob + i;
  
  //DISP(GlobalRowNumbers);

  
  if (true)
    {
      // first choice : the arrays ListProcInterac, ListRowInterac, 
      // OverlapRowNumbers and OriginalProcNumbers
      // are recovered from the global row numbers
      
      // with this choice, you only need to give the global row numbers for each processor
      // however, more communications and memory will be needed to construct the arrays ListProcInterac, etc
      
      // this function also sorts GlobalRowNumber if not sorted, so be careful
      // if you are providing unsorted numbers on entry
      DistributedMatrix<double, General, ArrayRowSparse> Atest;
      Atest.Init(GlobalRowNumbers, OverlapRowNumbers, OriginalProcNumbers,
                 ListProcInterac, ListRowInterac, comm);
    }
  else
    {
      // second choice : the arrays ListProcInterac, ListRowInterac, 
      // OverlapRowNumbers and OriginalProcNumbers
      // are constructed manually with the knowledge of the geometry
      
      // Here, since we are using finite difference discretization
      // it is easy to know the numbers of overlapped rows
      
      // row numbers that are "overlapped" by another processor
      int Ns = 0;
      if (num_proc_row > 0)
        {
          Ns = Ny;
          if (num_proc_col > 0)
            Ns += Nx-1;
        }
      else if (num_proc_col > 0)
        Ns = Nx;
      
      OverlapRowNumbers.Reallocate(Ns); OriginalProcNumbers.Reallocate(Ns);
      OverlapRowNumbers.Fill(-1); OriginalProcNumbers.Fill(-1);
      if (num_proc_col > 0)
        {
          for (int i = 0; i < Nx; i++)
            {
              OverlapRowNumbers(i) = i;
              OriginalProcNumbers(i) = rank_processor - 1;
            }
          
          if (num_proc_row > 0)
            {
              OriginalProcNumbers(0) = rank_processor - 1 - nb_proc_col;
              for (int i = 1; i < Ny; i++)
                {
                  OverlapRowNumbers(Nx+i-1) = i*Nx;
                  OriginalProcNumbers(Nx+i-1) = rank_processor - nb_proc_col;
                }
            }
        }
      else if (num_proc_row > 0)
        {
          for (int i = 0; i < Ny; i++)
            {
              OverlapRowNumbers(i) = i*Nx;
              OriginalProcNumbers(i) = rank_processor - nb_proc_col;
            }      
        }
      
      //DISP(OverlapRowNumbers); DISP(OriginalProcNumbers);
      
      // dofs shared with other processors
      int nb_proc_interac = 0;
      if (num_proc_col > 0)
        {
          nb_proc_interac++;
          if (num_proc_row > 0)
            nb_proc_interac++;
          
          if (num_proc_row < nb_proc_row-1)
            nb_proc_interac++;      
        }
      
      if (num_proc_col < nb_proc_col-1)
        {
          nb_proc_interac++;
          if (num_proc_row > 0)
            nb_proc_interac++;
          
          if (num_proc_row < nb_proc_row-1)
            nb_proc_interac++;      
        }
      
      if (num_proc_row > 0)
        nb_proc_interac++;
      
      if (num_proc_row < nb_proc_row-1)
        nb_proc_interac++;
      
      ListProcInterac.Reallocate(nb_proc_interac);
      ListRowInterac.Reallocate(nb_proc_interac);
      
      nb_proc_interac = 0;
      if (num_proc_col > 0)
        {
          ListProcInterac(nb_proc_interac) = rank_processor-1;
          ListRowInterac(nb_proc_interac).Reallocate(Nx);
          for (int i = 0; i < Nx; i++)
            ListRowInterac(nb_proc_interac)(i) = i;
          
          nb_proc_interac++;
          
          if (num_proc_row > 0)
            {
              ListProcInterac(nb_proc_interac) = rank_processor-1-nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = 0;
              nb_proc_interac++;
            }
          
          if (num_proc_row < nb_proc_row-1)
            {
              ListProcInterac(nb_proc_interac) = rank_processor-1+nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = Nx-1;
              nb_proc_interac++;      
            }
        }
      
      if (num_proc_col < nb_proc_col-1)
        {
          ListProcInterac(nb_proc_interac) = rank_processor+1;
          ListRowInterac(nb_proc_interac).Reallocate(Nx);
          for (int i = 0; i < Nx; i++)
            ListRowInterac(nb_proc_interac)(i) = Nx*(Ny-1) + i;
          
          nb_proc_interac++;
          
          if (num_proc_row > 0)
            {
              ListProcInterac(nb_proc_interac) = rank_processor+1-nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = Nx*(Ny-1);
              nb_proc_interac++;
            }
          
          if (num_proc_row < nb_proc_row-1)
            {
              ListProcInterac(nb_proc_interac) = rank_processor+1+nb_proc_col;
              ListRowInterac(nb_proc_interac).Reallocate(1);
              ListRowInterac(nb_proc_interac)(0) = Nx*Ny-1;
              nb_proc_interac++;      
            }
        }
      
      if (num_proc_row > 0)
        {
          ListProcInterac(nb_proc_interac) = rank_processor-nb_proc_col;
          ListRowInterac(nb_proc_interac).Reallocate(Ny);
          for (int i = 0; i < Ny; i++)
            ListRowInterac(nb_proc_interac)(i) = Nx*i;
          
          nb_proc_interac++;
        }
      
      if (num_proc_row < nb_proc_row-1)
        {
          ListProcInterac(nb_proc_interac) = rank_processor+nb_proc_col;
          ListRowInterac(nb_proc_interac).Reallocate(Ny);
          for (int i = 0; i < Ny; i++)
            ListRowInterac(nb_proc_interac)(i) = Nx*i+Nx-1;
          
          nb_proc_interac++;
        }
      
      for (int j = 0; j < nb_proc_interac; j++)
        {
          //DISP(j);
          //DISP(ListProcInterac(j));
          //DISP(ListRowInterac(j));
        }
    }
}

