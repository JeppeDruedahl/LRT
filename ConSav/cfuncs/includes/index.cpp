namespace index {

int d2(int i_x1, int i_x2, int Nx2){

    int i_grid =   i_x1*Nx2
                 + i_x2;

    return i_grid;

}

int d3(int i_x1, int i_x2, int i_x3, int Nx2, int Nx3){

    int i_grid =   i_x1*Nx2*Nx3
                 + i_x2*Nx3
                 + i_x3;

    return i_grid;

}

} // index