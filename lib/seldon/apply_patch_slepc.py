#!/usr/bin/env python

# - on rajouter matsolve dans include/slepc/private/stimpl.h dans _p_ST :
 #   PetscErrorCode (*matsolve)(Mat A, Vec b, Vec x);
#    PetscErrorCode (*matsolve_trans)(Mat A, Vec b, Vec x);
    
try:
    fid = open("include/slepc/private/stimpl.h")
    lignes = fid.readlines(); fid.close();
    test_struct = False;
    for i in range(len(lignes)):
        if (lignes[i].startswith('struct _p_ST')):
            test_struct = True;

        if (test_struct):
            # pour detecter si le patch a deja ete applique
            if (lignes[i].find("matsolve") >= 0):
                test_struct = False;
            
            if (lignes[i].startswith('};')):
                lignes.insert(i, "  PetscErrorCode (*matsolve_trans)(Mat A, Vec b, Vec x);\n");
                lignes.insert(i, "  PetscErrorCode (*matsolve)(Mat A, Vec b, Vec x);\n");
                break;
    
    fid = open("include/slepc/private/stimpl.h", 'w')
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify include/slepc/private/stimpl.h")


# - on modifie STMatSolve dans src/sys/classes/st/interface/stsles.c pour mettre
#    if (st->matsolve != NULL)
#    return (*st->matsolve)(st->A[0], b, x);

#  et dans STMatSolveTranspose
#    if (st->matsolve_trans != NULL)  
#    return (*st->matsolve_trans)(st->A[0], b, x);
try:
    fid = open("src/sys/classes/st/interface/stsles.c")
    lignes = fid.readlines(); fid.close();
    test_func = False; i = 0;
    while (i < len(lignes)):
        if (lignes[i].startswith("PetscErrorCode STMatSolve(ST st")):
            if (lignes[i+2].find("matsolve") == -1):
                lignes.insert(i+2, "    return (*st->matsolve)(st->A[0], b, x);\n");
                lignes.insert(i+2, "  if (st->matsolve != NULL)\n");

        if (lignes[i].startswith("PetscErrorCode STMatSolveTranspose(ST st")):
            if (lignes[i+2].find("matsolve") == -1):
                lignes.insert(i+2, "    return (*st->matsolve_trans)(st->A[0], b, x);\n");
                lignes.insert(i+2, "  if (st->matsolve_trans != NULL)\n");

        i += 1
    
    fid = open("src/sys/classes/st/interface/stsles.c", 'w')
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify src/sys/classes/st/interface/stsles.c")
                
# - on initialise matsolve dans STCreate (stfunc.c)
#    st->matsolve = NULL;
#    st->matsolve_trans = NULL;
try:
    fid = open("src/sys/classes/st/interface/stfunc.c")
    lignes = fid.readlines(); fid.close();
    i = 0;
    while (i < len(lignes)):
        if (lignes[i].startswith("PetscErrorCode STCreate")):
            j = i; test_solve = False; pos_solve = i+4;
            while (j < len(lignes)):
                if (lignes[j].startswith('}')):
                    break;

                if (lignes[j].find("st->data") >= 0):
                    pos_solve = j+1;
                
                if (lignes[j].find("matsolve") >= 0):
                    test_solve = True

                j += 1
            
            if (not test_solve):
                lignes.insert(pos_solve, "  st->matsolve_trans = NULL;\n");
                lignes.insert(pos_solve, "  st->matsolve = NULL;\n");

            break;

        i += 1
                
    fid = open("src/sys/classes/st/interface/stfunc.c", 'w')
    fid.writelines(lignes); fid.close();
except:
    print("Unable to modify src/sys/classes/st/interface/stfunc.c")
