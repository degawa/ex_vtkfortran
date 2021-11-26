module ex_vtkfortran_vtr
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: write_scalar_node_vtr
    public :: write_vector_node_vtr
    public :: write_scalar_cell_vtr

contains
    !|RectilinearGrid Grid上のノードに定義されたスカラ値を出力する例．
    !
    ! <?xml version="1.0"?>
    ! <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">
    !   <RectilinearGrid WholeExtent="+1 +10 +1 +6 +1 +6">
    !     <Piece Extent="+1 +10 +1 +6 +1 +6">
    !       <Coordinates>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="X" format="appended" offset="0"/>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="appended" offset="84"/>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="Z" format="appended" offset="136"/>
    !       </Coordinates>
    !       <PointData>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="float64_scalar" format="appended" offset="188"/>
    !       </PointData>
    !     </Piece>
    !   </RectilinearGrid>
    !   <AppendedData encoding="raw">
    ! raw data...
    !   </AppendedData>
    ! </VTKFile>
    subroutine write_scalar_node_vtr()
        use :: vtk_fortran, only:vtk_file
        implicit none

        type(vtk_file) :: vtr !! VTK file
        integer(int32), parameter :: nx1 = 1 !! x方向の配列範囲の下限
        integer(int32), parameter :: nx2 = 10 !! x方向の配列範囲の上限
        integer(int32), parameter :: ny1 = 1 !! y方向の配列範囲の下限
        integer(int32), parameter :: ny2 = 6 !! y方向の配列範囲の上限
        integer(int32), parameter :: nz1 = 1 !! z方向の配列範囲の下限
        integer(int32), parameter :: nz2 = 6 !! z方向の配列範囲の上限
        real(real64) :: x(nx1:nx2) !! x座標値
        real(real64) :: y(ny1:ny2) !! y座標値
        real(real64) :: z(nz1:nz2) !! z座標値
        real(real64) :: v(nx1:nx2, ny1:ny2, nz1:nz2) !! スカラ値

        integer(int32) :: stat !! 手続の実行結果の状態

        block
            integer(int32) :: i, j, k

            x = [(dble(i), i=nx1, nx2)]
            y = [(dble(j), j=ny1, ny2)]
            z = [(dble(k), k=nz1, nz2)]

            do concurrent(k=nz1:nz2, j=ny1:ny2, i=nx1:nx2)
                v(i, j, k) = dble(i*j*k)
            end do
        end block

        ! ファイルを開き，VTKFileタグとStructuredGridタグをopen．
        ! <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">
        !   <RectilinearGrid WholeExtent="+1 +10 +1 +6 +1 +6">
        stat = vtr%initialize(format="raw", &
                              filename="xml_rectilinear_grid_raw_scalar.vtr", &
                              mesh_topology="RectilinearGrid", &
                              nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! Pieceタグ open．
        ! <Piece Extent="+1 +10 +1 +6 +1 +6">
        stat = vtr%xml_writer%write_piece(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! 各方向座標値の出力情報の書き出し．
        ! 実際のデータは </StructuredGrid>と</VTKFile>の間に，
        ! <AppendedData encoding="raw"></AppendedData>が設けられ，
        ! そこに書き出される．
        !
        ! <Coordinates>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="X" format="appended" offset="0"/>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="appended" offset="84"/>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="Z" format="appended" offset="136"/>
        ! </Coordinates>
        stat = vtr%xml_writer%write_geo(x=x, y=y, z=z)

        ! 配列の書き出し．
        ! 実際のデータはAppendedDataに追記される．
        !
        ! <PointData>
        stat = vtr%xml_writer%write_dataarray(location="node", action="open")
        !   <DataArray type="Float64" NumberOfComponents="1" Name="float64_scalar" format="appended" offset="188"/>
        stat = vtr%xml_writer%write_dataarray(x=v, data_name="float64_scalar", one_component=.true.)
        ! </PointData>
        stat = vtr%xml_writer%write_dataarray(location="node", action="close")

        ! Pieceタグ close
        ! </Piece>
        stat = vtr%xml_writer%write_piece()

        ! VTKFileタグとStructuredGridタグ close．
        ! vtrファイルを閉じる．
        !   </StructuredGrid>
        ! </VTKFile>
        stat = vtr%finalize()
    end subroutine write_scalar_node_vtr

    !|RectilinearGrid Grid上のノードに定義されたベクトル値を出力する例．
    !
    ! <?xml version="1.0"?>
    ! <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">
    !   <RectilinearGrid WholeExtent="+1 +10 +1 +6 +1 +6">
    !     <Piece Extent="+1 +10 +1 +6 +1 +6">
    !       <Coordinates>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="X" format="appended" offset="0"/>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="appended" offset="84"/>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="Z" format="appended" offset="136"/>
    !       </Coordinates>
    !       <PointData>
    !         <DataArray type="Float64" NumberOfComponents="3" Name="float64_vector" format="appended" offset="188"/>
    !       </PointData>
    !     </Piece>
    !   </RectilinearGrid>
    !   <AppendedData encoding="raw">
    ! raw data...
    !   </AppendedData>
    ! </VTKFile>
    subroutine write_vector_node_vtr()
        use :: vtk_fortran, only:vtk_file
        implicit none

        type(vtk_file) :: vtr !! VTK file
        integer(int32), parameter :: nx1 = 1 !! x方向の配列範囲の下限
        integer(int32), parameter :: nx2 = 10 !! x方向の配列範囲の上限
        integer(int32), parameter :: ny1 = 1 !! y方向の配列範囲の下限
        integer(int32), parameter :: ny2 = 6 !! y方向の配列範囲の上限
        integer(int32), parameter :: nz1 = 1 !! z方向の配列範囲の下限
        integer(int32), parameter :: nz2 = 6 !! z方向の配列範囲の上限
        real(real64) :: x(nx1:nx2) !! x座標値
        real(real64) :: y(ny1:ny2) !! y座標値
        real(real64) :: z(nz1:nz2) !! z座標値
        real(real64) :: u(nx1:nx2, ny1:ny2, nz1:nz2) !! ベクトル値のx成分
        real(real64) :: v(nx1:nx2, ny1:ny2, nz1:nz2) !! ベクトル値のy成分
        real(real64) :: w(nx1:nx2, ny1:ny2, nz1:nz2) !! ベクトル値のz成分

        integer(int32) :: stat !! 手続の実行結果の状態

        block
            integer(int32) :: i, j, k

            x = [(dble(i), i=nx1, nx2)]
            y = [(dble(j), j=ny1, ny2)]
            z = [(dble(k), k=nz1, nz2)]

            do concurrent(k=nz1:nz2, j=ny1:ny2, i=nx1:nx2)
                u(i, j, k) = 0.5d0*dble(i*j*k)
                v(i, j, k) = 1.0d0*dble(i*j*k)
                w(i, j, k) = 1.5d0*dble(i*j*k)
            end do
        end block

        ! ファイルを開き，VTKFileタグとStructuredGridタグをopen．
        ! <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">
        !   <RectilinearGrid WholeExtent="+1 +10 +1 +6 +1 +6">
        stat = vtr%initialize(format="raw", &
                              filename="xml_rectilinear_grid_raw_vector.vtr", &
                              mesh_topology="RectilinearGrid", &
                              nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! Pieceタグ open．
        ! <Piece Extent="+1 +10 +1 +6 +1 +6">
        stat = vtr%xml_writer%write_piece(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! 各方向座標値の出力情報の書き出し．
        ! 実際のデータは </StructuredGrid>と</VTKFile>の間に，
        ! <AppendedData encoding="raw"></AppendedData>が設けられ，
        ! そこに書き出される．
        !
        ! <Coordinates>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="X" format="appended" offset="0"/>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="appended" offset="84"/>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="Z" format="appended" offset="136"/>
        ! </Coordinates>
        stat = vtr%xml_writer%write_geo(x=x, y=y, z=z)

        ! 配列の書き出し．
        ! 実際のデータはAppendedDataに追記される．
        !
        ! <PointData>
        stat = vtr%xml_writer%write_dataarray(location="node", action="open")
        !   <DataArray type="Float64" NumberOfComponents="3" Name="float64_vector" format="appended" offset="188"/>
        stat = vtr%xml_writer%write_dataarray(x=u, y=v, z=w, data_name="float64_vector", is_tuples=.false.)
        ! </PointData>
        stat = vtr%xml_writer%write_dataarray(location="node", action="close")

        ! Pieceタグ close
        ! </Piece>
        stat = vtr%xml_writer%write_piece()

        ! VTKFileタグとStructuredGridタグ close．
        ! vtrファイルを閉じる．
        !   </StructuredGrid>
        ! </VTKFile>
        stat = vtr%finalize()
    end subroutine write_vector_node_vtr

    !|RectilinearGrid Grid上のセルに定義されたスカラ値を出力する例．
    !
    ! <?xml version="1.0"?>
    ! <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">
    !   <RectilinearGrid WholeExtent="+1 +16 +1 +16 +1 +16">
    !     <Piece Extent="+1 +16 +1 +16 +1 +16">
    !       <Coordinates>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="X" format="appended" offset="0"/>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="appended" offset="132"/>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="Z" format="appended" offset="264"/>
    !       </Coordinates>
    !       <CellData>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="float64_scalar_cell" format="appended" offset="396"/>
    !       </CellData>
    !     </Piece>
    !   </RectilinearGrid>
    !   <AppendedData encoding="raw">
    ! raw data...
    !   </AppendedData>
    ! </VTKFile>
    subroutine write_scalar_cell_vtr()
        use :: vtk_fortran, only:vtk_file
        implicit none

        type(vtk_file) :: vtr !! VTK file
        integer(int32), parameter :: nx1 = 1 !! x方向の配列範囲の下限
        integer(int32), parameter :: nx2 = 16 !! x方向の配列範囲の上限
        integer(int32), parameter :: ny1 = 1 !! y方向の配列範囲の下限
        integer(int32), parameter :: ny2 = 16 !! y方向の配列範囲の上限
        integer(int32), parameter :: nz1 = 1 !! z方向の配列範囲の下限
        integer(int32), parameter :: nz2 = 16 !! z方向の配列範囲の上限
        integer(int32), parameter :: ncell = (nx2 - nx1)*(ny2 - ny1)*(nz2 - nz1)
            !! セル数．ノードの数より各方向それぞれ1点少ない．
        real(real64) :: x(nx1:nx2) !! x座標値
        real(real64) :: y(ny1:ny2) !! y座標値
        real(real64) :: z(nz1:nz2) !! z座標値
        real(real64) :: v(1:ncell) !! セルのスカラ値

        integer(int32) :: stat !! 手続の実行結果の状態

        block
            integer(int32) :: i, j, k

            x = [(dble(i), i=nx1, nx2)]
            y = [(dble(j), j=ny1, ny2)]
            z = [(dble(k), k=nz1, nz2)]

            do concurrent(k=nz1:nz2 - 1, j=ny1:ny2 - 1, i=nx1:nx2 - 1)
                block
                    integer(int32) :: ijk
                    ijk = (k-1)*((nx2-1)*(ny2-1)) + (j-1)*(nx2-1) + i !&
                    v(ijk) = dble(i*j*k)
                end block
            end do
        end block

        stat = 0

        ! ファイルを開き，VTKFileタグとStructuredGridタグをopen．
        ! <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">
        !   <RectilinearGrid WholeExtent="+1 +16 +1 +16 +1 +16">
        stat = vtr%initialize(format="raw", &
                              filename="xml_rectilinear_grid_raw_scalar_cell.vtr", &
                              mesh_topology="RectilinearGrid", &
                              nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! Pieceタグ open．
        ! <Piece Extent="+1 +16 +1 +16 +1 +16">
        stat = vtr%xml_writer%write_piece(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! 各方向座標値の出力情報の書き出し．
        ! 実際のデータは </StructuredGrid>と</VTKFile>の間に，
        ! <AppendedData encoding="raw"></AppendedData>が設けられ，
        ! そこに書き出される．
        !
        ! <Coordinates>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="X" format="appended" offset="0"/>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="Y" format="appended" offset="132"/>
        !   <DataArray type="Float64" NumberOfComponents="1" Name="Z" format="appended" offset="264"/>
        ! </Coordinates>
        stat = vtr%xml_writer%write_geo(x=x, y=y, z=z)

        ! 配列の書き出し．
        ! 実際のデータはAppendedDataに追記される．
        ! <CellData>
        stat = vtr%xml_writer%write_dataarray(location="cell", action="open")
        !   <DataArray type="Float64" NumberOfComponents="1" Name="float64_scalar_cell" format="appended" offset="396"/>
        stat = vtr%xml_writer%write_dataarray(x=v, data_name="float64_scalar_cell")
        ! </CellData>!
        stat = vtr%xml_writer%write_dataarray(location="cell", action="close")

        ! Pieceタグ close
        ! </Piece>
        stat = vtr%xml_writer%write_piece()

        ! VTKFileタグとStructuredGridタグ close．
        ! vtrファイルを閉じる．
        !   </StructuredGrid>
        ! </VTKFile>
        stat = vtr%finalize()
    end subroutine write_scalar_cell_vtr
end module ex_vtkfortran_vtr
