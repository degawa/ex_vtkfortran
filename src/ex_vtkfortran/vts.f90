module ex_vtkfortran_vts
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: write_scalar_node_vts

contains
    !|Structured Grid上のノードに定義されたスカラ値を出力する例．
    !
    ! <?xml version="1.0"?>
    ! <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian">
    !   <StructuredGrid WholeExtent="+0 +9 +0 +5 +0 +5">
    !     <Piece Extent="+0 +9 +0 +5 +0 +5">
    !       <Points>
    !         <DataArray type="Float64" NumberOfComponents="3" Name="Points" format="appended" offset="0"/>
    !       </Points>
    !       <PointData>
    !         <DataArray type="Float64" NumberOfComponents="1" Name="float64_scalar" format="appended" offset="8644"/>
    !       </PointData>
    !     </Piece>
    !   </StructuredGrid>
    !   <AppendedData encoding="raw">
    ! raw data...
    !   </AppendedData>
    ! </VTKFile>
    subroutine write_scalar_node_vts()
        use :: vtk_fortran, only:vtk_file
        implicit none

        type(vtk_file) :: vts !! VTK file
        integer(int32), parameter :: nx1 = 0 !! x方向の配列範囲の下限
        integer(int32), parameter :: nx2 = 9 !! x方向の配列範囲の上限
        integer(int32), parameter :: ny1 = 0 !! y方向の配列範囲の下限
        integer(int32), parameter :: ny2 = 5 !! y方向の配列範囲の上限
        integer(int32), parameter :: nz1 = 0 !! z方向の配列範囲の下限
        integer(int32), parameter :: nz2 = 5 !! z方向の配列範囲の上限
        integer(int32), parameter :: nn = (nx2 - nx1 + 1)*(ny2 - ny1 + 1)*(nz2 - nz1 + 1) !! 全要素数
        real(real64) :: x(nx1:nx2, ny1:ny2, nz1:nz2) !! x座標値
        real(real64) :: y(nx1:nx2, ny1:ny2, nz1:nz2) !! y座標値
        real(real64) :: z(nx1:nx2, ny1:ny2, nz1:nz2) !! z座標値
        real(real64) :: v(nx1:nx2, ny1:ny2, nz1:nz2) !! スカラ値

        integer(int32) :: stat !! 手続の実行結果の状態

        block
            integer(int32) :: i, j, k
            real(real64) :: x_, y_, z_
            do concurrent(k=nz1:nz2, j=ny1:ny2, i=nx1:nx2)
                x_ = nx2*(dble(i)/dble(nx2))**2
                y_ = ny2*(dble(j)/dble(ny2))**2
                z_ = nz2*(dble(k)/dble(nz2))**2
                x(i, j, k) = x_
                y(i, j, k) = y_
                z(i, j, k) = z_
                v(i, j, k) = dble(i*j*k)
            end do
        end block

        stat = 0

        ! ファイルを開き，VTKFileタグとStructuredGridタグをopen．
        ! <VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian">
        !   <StructuredGrid WholeExtent="+0 +9 +0 +5 +0 +5">
        stat = vts%initialize(format='raw', &
                              filename='xml_structured_grid_raw_scalar.vts', &
                              mesh_topology='StructuredGrid', &
                              nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! Pieceタグ open．
        ! <Piece Extent="+0 +9 +0 +5 +0 +5">
        stat = vts%xml_writer%write_piece(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2)

        ! 各方向座標値の出力情報の書き出し．
        ! 実際のデータは </StructuredGrid>と</VTKFile>の間に，
        ! <AppendedData encoding="raw"></AppendedData>が設けられ，
        ! そこに書き出される．
        !
        ! <Points>
        !     <DataArray type="Float64" NumberOfComponents="3" Name="Points" format="appended" offset="0"/>
        ! </Points>
        stat = vts%xml_writer%write_geo(n=nn, x=x, y=y, z=z)

        ! 配列の書き出し．
        ! 実際のデータはAppendedDataに追記される．
        !
        ! <PointData>
        stat = vts%xml_writer%write_dataarray(location='node', action='open')
        ! <DataArray type="Float64" NumberOfComponents="1" Name="float64_scalar" format="appended" offset="8644"/>
        stat = vts%xml_writer%write_dataarray(data_name='float64_scalar', x=v, one_component=.true.)
        ! </PointData>
        stat = vts%xml_writer%write_dataarray(location='node', action='close')

        ! Pieceタグ close
        ! </Piece>
        stat = vts%xml_writer%write_piece()

        ! VTKFileタグとStructuredGridタグ close．
        ! vtsファイルを閉じる．
        !   </StructuredGrid>
        ! </VTKFile>
        stat = vts%finalize()
    end subroutine write_scalar_node_vts
end module ex_vtkfortran_vts
