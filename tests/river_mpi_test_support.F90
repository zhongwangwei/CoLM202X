MODULE MOD_Precision
   IMPLICIT NONE
   integer, parameter :: r8 = selected_real_kind(12)
END MODULE MOD_Precision

MODULE MOD_DataType
   USE MOD_Precision
   IMPLICIT NONE
   type :: pointer_int32_1d
      integer, allocatable :: val(:)
   END type pointer_int32_1d
   type :: pointer_int32_2d
      integer, allocatable :: val(:,:)
   END type pointer_int32_2d
   type :: pointer_real8_1d
      real(r8), allocatable :: val(:)
   END type pointer_real8_1d
END MODULE MOD_DataType

MODULE MOD_SPMD_Task
   IMPLICIT NONE
   include 'mpif.h'
   logical :: p_is_worker = .true.
   integer :: p_comm_worker = MPI_COMM_NULL
   integer :: p_iam_worker = -1
   integer :: p_np_worker = 0
   integer :: p_np_io = 0
   integer :: p_stat(MPI_STATUS_SIZE)
   integer :: p_err = MPI_SUCCESS
CONTAINS
   SUBROUTINE CoLM_stop(message)
      character(len=*), intent(in), optional :: message
      IF (present(message)) write(*,'(A)') trim(message)
      CALL MPI_Abort(MPI_COMM_WORLD, 1, p_err)
   END SUBROUTINE CoLM_stop
END MODULE MOD_SPMD_Task

MODULE MOD_Utils
   IMPLICIT NONE
   INTERFACE insert_into_sorted_list1
      MODULE procedure insert_into_sorted_list1_int32
   END INTERFACE
   INTERFACE find_in_sorted_list1
      MODULE procedure find_in_sorted_list1_int32
   END INTERFACE
   INTERFACE quicksort
      MODULE procedure quicksort_int32
   END INTERFACE
CONTAINS
   SUBROUTINE insert_into_sorted_list1_int32(x, n, list, iloc, is_new_out)
      integer, intent(in) :: x
      integer, intent(inout) :: n
      integer, intent(inout) :: list(:)
      integer, intent(out) :: iloc
      logical, intent(out), optional :: is_new_out
      integer :: ileft, iright
      logical :: is_new

      IF (n == 0) THEN
         iloc = 1
         is_new = .true.
      ELSEIF (x <= list(1)) THEN
         iloc = 1
         is_new = x /= list(1)
      ELSEIF (x > list(n)) THEN
         iloc = n + 1
         is_new = .true.
      ELSEIF (x == list(n)) THEN
         iloc = n
         is_new = .false.
      ELSE
         ileft = 1
         iright = n
         DO
            IF (iright - ileft <= 1) THEN
               iloc = iright
               is_new = .true.
               EXIT
            ENDIF
            iloc = (ileft + iright) / 2
            IF (x > list(iloc)) THEN
               ileft = iloc
            ELSEIF (x < list(iloc)) THEN
               iright = iloc
            ELSE
               is_new = .false.
               EXIT
            ENDIF
         ENDDO
      ENDIF
      IF (is_new) THEN
         IF (iloc <= n) list(iloc+1:n+1) = list(iloc:n)
         list(iloc) = x
         n = n + 1
      ENDIF
      IF (present(is_new_out)) is_new_out = is_new
   END SUBROUTINE insert_into_sorted_list1_int32

   integer FUNCTION find_in_sorted_list1_int32(x, n, list) RESULT(iloc)
      integer, intent(in) :: x, n
      integer, intent(in) :: list(n)
      integer :: ileft, iright

      iloc = 0
      IF (n <= 0 .or. x < list(1) .or. x > list(n)) RETURN
      ileft = 1
      iright = n
      DO WHILE (ileft <= iright)
         iloc = (ileft + iright) / 2
         IF (x < list(iloc)) THEN
            iright = iloc - 1
         ELSEIF (x > list(iloc)) THEN
            ileft = iloc + 1
         ELSE
            RETURN
         ENDIF
      ENDDO
      iloc = 0
   END FUNCTION find_in_sorted_list1_int32

   SUBROUTINE insert_into_sorted_list2(x, y, n, xlist, ylist, iloc, is_new_out)
      integer, intent(in) :: x, y
      integer, intent(inout) :: n
      integer, intent(inout) :: xlist(:), ylist(:)
      integer, intent(out) :: iloc
      logical, intent(out), optional :: is_new_out
      iloc = n + 1
      n = n + 1
      xlist(iloc) = x
      ylist(iloc) = y
      IF (present(is_new_out)) is_new_out = .true.
   END SUBROUTINE insert_into_sorted_list2

   integer FUNCTION find_in_sorted_list2(x, y, n, xlist, ylist) RESULT(iloc)
      integer, intent(in) :: x, y, n
      integer, intent(in) :: xlist(n), ylist(n)
      integer :: i
      iloc = 0
      DO i = 1, n
         IF (xlist(i) == x .and. ylist(i) == y) THEN
            iloc = i
            RETURN
         ENDIF
      ENDDO
   END FUNCTION find_in_sorted_list2

   recursive SUBROUTINE quicksort_int32(n, values, order)
      integer, intent(in) :: n
      integer, intent(inout) :: values(n), order(n)
      integer :: left, right, pivot, marker, value_tmp, order_tmp

      IF (n <= 1) RETURN
      pivot = values(n / 2)
      left = 0
      right = n + 1
      DO WHILE (left < right)
         right = right - 1
         DO WHILE (values(right) > pivot)
            right = right - 1
         ENDDO
         left = left + 1
         DO WHILE (values(left) < pivot)
            left = left + 1
         ENDDO
         IF (left < right) THEN
            value_tmp = values(left)
            values(left) = values(right)
            values(right) = value_tmp
            order_tmp = order(left)
            order(left) = order(right)
            order(right) = order_tmp
         ENDIF
      ENDDO
      marker = right
      CALL quicksort_int32(marker, values(1:marker), order(1:marker))
      CALL quicksort_int32(n-marker, values(marker+1:n), order(marker+1:n))
   END SUBROUTINE quicksort_int32
END MODULE MOD_Utils

MODULE MOD_Pixelset
   IMPLICIT NONE
   type :: subset_type
      integer, allocatable :: substt(:), subend(:)
   END type subset_type
   type :: pixelset_type
      integer :: unused = 0
   END type pixelset_type
END MODULE MOD_Pixelset

MODULE MOD_Grid
   IMPLICIT NONE
   type :: grid_type
      integer :: nlon = 0
   END type grid_type
   type :: grid_list_type
      integer :: ng = 0
      integer, allocatable :: ilon(:), ilat(:)
   END type grid_list_type
END MODULE MOD_Grid

MODULE MOD_SpatialMapping
   USE MOD_Precision
   USE MOD_DataType
   USE MOD_Grid
   USE MOD_Pixelset
   IMPLICIT NONE
   type :: spatial_mapping_type
      type(grid_list_type), allocatable :: glist(:)
      integer :: npset = 0
      integer, allocatable :: npart(:)
      type(pointer_int32_2d), allocatable :: address(:)
      type(pointer_real8_1d), allocatable :: areapart(:)
   CONTAINS
      procedure :: build_arealweighted
   END type spatial_mapping_type
CONTAINS
   SUBROUTINE build_arealweighted(this, grid, pixelset)
      class(spatial_mapping_type), intent(inout) :: this
      type(grid_type), intent(in) :: grid
      type(pixelset_type), intent(in) :: pixelset
   END SUBROUTINE build_arealweighted
END MODULE MOD_SpatialMapping
