namespace mymex {

    double* set_field_double(mxArray* my_struct, const char *name, size_t ndim, size_t*dims){
        
        auto var = mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL); 
        mxSetField(my_struct,0,name,var);
        return (double*) mxGetPr(var);    

    }

    int* set_field_int(mxArray* my_struct, const char *name, size_t ndim, size_t*dims){

        auto var = mxCreateNumericArray(ndim,dims,mxINT32_CLASS,mxREAL); 
        mxSetField(my_struct,0,name,var);
        return (int*) mxGetData(var);    

    }

    double** set_field_cell(mxArray * my_struct, const char *name, 
                            size_t ndim_cell, size_t *dims_cell, 
                            size_t *ndim, size_t **dims){
        
        // a. field
        auto cell = mxCreateCellArray(ndim_cell, dims_cell);
        mxSetField(my_struct,0,name,cell);

        // b. allocate memory
        auto out = new double*[dims_cell[1]*dims_cell[0]];
        for(int t = 0; t < dims_cell[1]; t++){
        for(int d = 0; d < dims_cell[0]; d++){

            int i_cell  = index::d2(t,d,dims_cell[0]);
            auto array  = mxCreateNumericArray(ndim[0],dims[0],mxDOUBLE_CLASS,mxREAL);
            mxSetCell(cell, i_cell, array);
            out[i_cell] = mxGetPr(array);

        } }

        return out;

    }

}