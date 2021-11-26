program main
    use ex_vtkfortran_vts
    use ex_vtkfortran_vtr
    implicit none

    ! 3次元スカラ値（ノードデータ）をvtsファイルに出力する．
    call write_scalar_node_vts()
    print *, "output xml structured grid raw scalar vts"

    ! 3次元スカラ値（ノードデータ）をvtrファイルに出力する．
    call write_scalar_node_vtr()
    print *, "output xml rectilinear grid raw scalar vtr"

    ! 3次元ベクトル値（ノードデータ）をvtrファイルに出力する．
    call write_vector_node_vtr()
    print *, "output xml rectilinear grid raw vector vtr"

    ! 3次元スカラ値（セルデータ）をvtrファイルに出力する．
    call write_scalar_cell_vtr()
    print *, "output xml rectilinear grid raw scalar cell vtr"
end program main
